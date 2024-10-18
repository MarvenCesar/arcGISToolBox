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
        self.canRunInBackground = False
        self.optimization_attempts = 20
        self.aircraft_data = None  # Add this line

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
                displayname="Interior Taxiway Width (ft)",
                name="interior_taxiway_width",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayname="Peripheral Taxiway Width (ft)",
                name="peripheral_taxiway_width",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Aircraft Wingtip Clearance (ft)",
                name="wingtip_clearance",
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
                direction="Input"),
            arcpy.Parameter(
                displayName="Optimize Quantities",
                name="optimize_quantities",
                datatype="GPBoolean",
                parameterType="Optional",
                direction="Input"),
            arcpy.Parameter(
                displayName="Minimum Aircraft Quantities (Comma Separated)",
                name="min_quantities",
                datatype="GPString",
                parameterType="Optional",
                direction="Input",
                enabled=False),
            arcpy.Parameter(
                displayName="Number of Optimization Attempts",
                name="optimization_attempts",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
                enabled=False),
            arcpy.Parameter(
                displayName="Output Feature Class",
                name="output_fc",
                datatype="DEFeatureClass",
                parameterType="Required",
                direction="Output")
        ]
        params[1].filter.list = ["Polygon"]  # Only allow polygon-type airfield layers
        params[9].value = False  # Default value for optimization
        params[11].value = 20 # Default value for optimization attempts
        return params

    def updateParameters(self, parameters):
        # Populate Airfield Dropdown
        if parameters[1].altered and not parameters[2].altered:
            airfield_layer = parameters[1].valueAsText
            if airfield_layer:
                airfield_names = [row[0] for row in arcpy.da.SearchCursor(airfield_layer, "AFLD_NAME")]
                parameters[2].filter.list = sorted(set(airfield_names))
        
        # Populate Aircraft Dropdown based on the Input Aircraft Table
        if parameters[0].altered and not parameters[6].altered:
            aircraft_table = parameters[0].valueAsText
            if aircraft_table:
                try:
                    # Fetch unique aircraft names (MDS) from the table
                    aircraft_names = set(row[0] for row in arcpy.da.SearchCursor(aircraft_table, ["MDS"]))
                    parameters[6].filter.list = sorted(aircraft_names)  # Populate dropdown with unique MDS
                except Exception as e:
                    arcpy.AddError(f"Error loading aircraft names: {e}")

        optimize_quantities_param = parameters[9].value
        min_aircraft_quantities_param = parameters[10]
        optimization_attempts_param = parameters[11]

        optimize_quantities_param.value = False if optimize_quantities_param.value is None else optimize_quantities_param.value
        min_aircraft_quantities_param.enabled = optimize_quantities_param.value
        optimization_attempts_param.enabled = optimize_quantities_param.value

    def execute(self, parameters, messages):
        # Extract parameter values
        aircraft_table = parameters[0].valueAsText
        airfield_layer = parameters[1].valueAsText
        airfield_layer_object = parameters[1].value # added to work for the airfield layer generation
        airfield_name = parameters[2].valueAsText
        interior_taxi_width = float(parameters[3].valueAsText)
        peripheral_taxi_width = float(parameters[4].valueAsText)
        wingtip_clearance = float(parameters[5].valueAsText)
        selected_aircraft = parameters[6].valueAsText.split(';')
        aircraft_quantities = parameters[7].valueAsText
        apply_lcn = parameters[8].value
        optimize_quantities = parameters[9].value
        min_quantities = parameters[10].valueAsText

        # Process min quantities
        if parameters[10].value:
            min_quantities = parameters[10].valueAsText
        else:
            min_quantities = ','.join(['1'] * len(selected_aircraft))  # Default to 1 for each aircraft type

        # Get airfield data
        apron_length, apron_width, apron_lcn = self.get_airfield_data(airfield_layer, airfield_name)

        # Get aircraft data
        aircraft_data = self.get_aircraft_data(aircraft_table, selected_aircraft)

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
                with arcpy.da.SearchCursor(aircraft_table, ["MDS", "LENGTH", "WING_SPAN", "ACFT_LCN"]) as cursor:
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
                            input_features.append([apron_length, apron_width, length, wingspan, interior_taxi_width, peripheral_taxi_width, wingtip_clearance, aircraft_lcn, aircraft_quantities[i], apply_lcn])

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
