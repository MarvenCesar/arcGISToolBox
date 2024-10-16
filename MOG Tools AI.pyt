import random
import arcpy
import pandas as pd
import numpy as np
import os
import math
import traceback
import joblib

class Toolbox(object):
    def __init__(self):
        self.label = "Aircraft MOG Optimization"
        self.alias = "AircraftMOG"
        self.tools = [CalculateMaximumOnGround]

class CalculateMaximumOnGround(object):
    def __init__(self):
        self.label = "Calculate Maximum On Ground"
        self.description = "Calculate the maximum number of aircraft that can park on a specific apron using AI for optimized parking."

    def getParameterInfo(self):
        params = [
            arcpy.Parameter(
                displayName="Input Aircraft Table",
                name="in_table",
                datatype="DETable",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Airfield Layer",
                name="airfield_layer",
                datatype="GPFeatureLayer",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Airfield Name",
                name="airfield_name",
                datatype="GPString",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Aircraft Clearance (ft)",
                name="aircraft_clearance",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input"),  # Clearance around each aircraft
            arcpy.Parameter(
                displayName="Select Aircraft (MDS)",
                name="selected_aircraft",
                datatype="GPString",
                parameterType="Required",
                direction="Input",
                multiValue=True),  # Allow multiple aircraft to be selected
            arcpy.Parameter(
                displayName="Aircraft Quantities (Comma Separated)",
                name="aircraft_quantities",
                datatype="GPString",
                parameterType="Required",
                direction="Input"),  # Comma-separated values for the quantities
            arcpy.Parameter(
                displayName="Apply LCN Compatibility?",
                name="apply_lcn",
                datatype="GPBoolean",
                parameterType="Optional",
                direction="Input")
        ]
        params[1].filter.list = ["Polygon"]  # Only allow polygon-type airfield layers
        return params

    def updateParameters(self, parameters):
        # Populate Airfield Dropdown
        if parameters[1].altered and not parameters[2].altered:
            airfield_layer = parameters[1].valueAsText
            if airfield_layer:
                airfield_names = [row[0] for row in arcpy.da.SearchCursor(airfield_layer, "AFLD_NAME")]
                parameters[2].filter.list = sorted(set(airfield_names))
        
        # Populate Aircraft Dropdown based on the Input Aircraft Table
        if parameters[0].altered and not parameters[4].altered:
            aircraft_table = parameters[0].valueAsText
            if aircraft_table:
                try:
                    # Fetch unique aircraft names (MDS) from the table
                    aircraft_names = set()
                    with arcpy.da.SearchCursor(aircraft_table, ["MDS"]) as cursor:
                        for row in cursor:
                            aircraft_names.add(row[0])
                    parameters[4].filter.list = sorted(aircraft_names)  # Populate dropdown with unique MDS
                except Exception as e:
                    arcpy.AddError(f"Error loading aircraft names: {e}")

        return

    def execute(self, parameters, messages):
        in_table = parameters[0].valueAsText
        airfield_layer = parameters[1].valueAsText
        airfield_name = parameters[2].valueAsText
        aircraft_clearance = float(parameters[3].valueAsText)  # Clearance around each aircraft
        selected_aircraft = parameters[4].values  # List of selected aircraft MDS
        aircraft_quantities = parameters[5].valueAsText.split(',')  # Quantities entered as comma-separated values
        apply_lcn = parameters[6].value  # Whether to apply LCN compatibility constraint

        # Ensure that the number of aircraft matches the number of quantities
        if len(selected_aircraft) != len(aircraft_quantities):
            arcpy.AddError(f"Number of selected aircraft does not match the number of quantities provided.")
            return

        # Convert aircraft quantities to integers
        try:
            aircraft_quantities = [int(q.strip()) for q in aircraft_quantities]
        except ValueError:
            arcpy.AddError(f"Invalid input for aircraft quantities. Ensure all quantities are valid integers.")
            return

        try:
            # Load the AI model
            model_path = r"aircraft_parking_model.pkl"
            aircraft_parking_model = joblib.load(model_path)
            arcpy.AddMessage("Loaded AI model for optimized aircraft parking.")

            # Get airfield data based on the selected name
            where_clause = f"AFLD_NAME = '{airfield_name}'"
            with arcpy.da.SearchCursor(airfield_layer, ["SHAPE@", "LENGTH", "WIDTH", "LCN"], where_clause) as cursor:
                for row in cursor:
                    airfield_shape, apron_length, apron_width, apron_lcn = row
                    break
                else:
                    arcpy.AddError(f"Airfield '{airfield_name}' not found.")
                    return

            arcpy.AddMessage(f"Airfield dimensions: Length={apron_length} ft, Width={apron_width} ft")
            arcpy.AddMessage(f"Airfield LCN: {apron_lcn}")

            # Prepare input features for the AI model
            input_features = []
            for i, mds in enumerate(selected_aircraft):
                with arcpy.da.SearchCursor(in_table, ["MDS", "LENGTH", "WING_SPAN", "ACFT_LCN"]) as cursor:
                    for row in cursor:
                        aircraft_mds, length, wingspan, aircraft_lcn = row
                        if aircraft_mds == mds:
                            length = float(length) if length is not None else 0
                            wingspan = float(wingspan) if wingspan is not None else 0
                            aircraft_lcn = float(aircraft_lcn) if aircraft_lcn is not None else 0
                            if length == 0 or wingspan == 0 or aircraft_lcn == 0:
                                arcpy.AddError(f"Invalid data for {aircraft_mds}: LENGTH={length}, WING_SPAN={wingspan}, ACFT_LCN={aircraft_lcn}")
                                return
                            # Add data as features for prediction
                            input_features.append([apron_length, apron_width, length, wingspan, aircraft_clearance, aircraft_quantities[i], apply_lcn])

            # Convert input features to numpy array for model prediction
            input_features = np.array(input_features)
            
            # Predict the optimal parking arrangement using the AI model
            predicted_parking = aircraft_parking_model.predict(input_features)
            arcpy.AddMessage(f"Predicted parking arrangement: {predicted_parking}")

            # Use the AI model's output to place aircraft
            # Here, we assume the model returns the optimal positions or arrangement that we can integrate with ArcGIS.
            for i, (mds, quantity) in enumerate(zip(selected_aircraft, aircraft_quantities)):
                arcpy.AddMessage(f"Place {quantity} of {mds} as per AI optimization.")
                # The placement logic can be implemented here based on the prediction output `predicted_parking[i]`.

            arcpy.AddMessage(f"Aircraft parking optimization completed using AI model.")

        except Exception as e:
            arcpy.AddError(f"An error occurred: {str(e)}")
            arcpy.AddError(traceback.format_exc())
