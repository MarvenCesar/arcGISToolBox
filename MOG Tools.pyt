# -*- coding: utf-8 -*-
import arcpy
import pandas as pd
import os
import math
import traceback


class Toolbox(object):
    def __init__(self):
        self.label = "Aircraft MOG Optimization"
        self.alias = "AircraftMOG"
        self.tools = [ImportAircraftData, CalculateAircraftFootprint, CalculateMaximumOnGround]

# Class for importing aircraft data
class ImportAircraftData(object):
    def __init__(self):
        self.label = "Import Aircraft Data"
        self.description = "Import aircraft specifications from a CSV file"

    def getParameterInfo(self):
        params = [
            arcpy.Parameter(
                displayName="Input CSV File", 
                name="in_csv", 
                datatype="DEFile", 
                parameterType="Required", 
                direction="Input"),
            arcpy.Parameter(
                displayName="Output Aircraft Table", 
                name="out_table", 
                datatype="DETable", 
                parameterType="Required", 
                direction="Output")
        ]
        return params

    def execute(self, parameters, messages):
        in_csv = parameters[0].valueAsText
        out_table = parameters[1].valueAsText
        try:
            if not os.path.exists(in_csv):
                arcpy.AddError(f"Input CSV file does not exist: {in_csv}")
                return

            df = pd.read_csv(in_csv)
            df = df.dropna(subset=['MDS'])  # Drop rows with missing MDS

            numeric_columns = ['LENGTH', 'WING_SPAN', 'HEIGHT', 'WING_HEIGHT', 'TURNING_RADIUS', 'MIN_RWY_LENGTH', 'ACFT_LCN']
            for col in numeric_columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')

            temp_csv = os.path.join(arcpy.env.scratchFolder, "temp.csv")
            df.to_csv(temp_csv, index=False)
            
            arcpy.conversion.TableToTable(temp_csv, os.path.dirname(out_table), os.path.basename(out_table))
            os.remove(temp_csv)
            
            arcpy.AddMessage(f"Successfully imported {len(df)} aircraft records.")
            arcpy.AddMessage(f"Table imported. Please open {out_table} to view data")
        
        except Exception as e:
            arcpy.AddError(f"An error occurred during import: {str(e)}")
            arcpy.AddError(arcpy.GetMessages())

class CalculateAircraftFootprint(object):
    def __init__(self):
        self.label = "Create Aircraft Symbol Layer"
        self.description = "Create points for selected aircraft on a specific airfield parking apron in a grid pattern"

    def getParameterInfo(self):
        params = [
            arcpy.Parameter(
                displayName="Input Aircraft Table",
                name="in_table",
                datatype="DETable",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Input Airfield Layer",
                name="airfield_layer",
                datatype="DEFeatureClass",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Output Feature Class",
                name="out_fc",
                datatype="DEFeatureClass",
                parameterType="Required",
                direction="Output"),
            arcpy.Parameter(
                displayName="Aircraft Names (MDS)",
                name="aircraft_names",
                datatype="GPString",
                parameterType="Required",
                direction="Input",
                multiValue=True),
            arcpy.Parameter(
                displayName="Quantity for Each Aircraft",
                name="quantity_per_aircraft",
                datatype="GPLong",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Airfield Name (AFLD_NAME)",
                name="afld_name",
                datatype="GPString",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Airfield ObjectID",
                name="object_id",
                datatype="GPLong",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Buffer Distance (in feet)",
                name="buffer_distance",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Maximum Aircraft per Row",
                name="max_per_row",
                datatype="GPLong",
                parameterType="Required",
                direction="Input")
        ]
        return params

    def execute(self, parameters, messages):
        in_table = parameters[0].valueAsText
        airfield_layer = parameters[1].valueAsText
        out_fc = parameters[2].valueAsText
        selected_aircraft = parameters[3].values
        quantity_per_aircraft = int(parameters[4].valueAsText)
        afld_name = parameters[5].valueAsText
        object_id = int(parameters[6].valueAsText)
        buffer_distance = float(parameters[7].valueAsText)
        max_per_row = int(parameters[8].valueAsText)

        try:
            # Validate and create output feature class
            workspace = os.path.dirname(out_fc)
            valid_name = arcpy.ValidateTableName(os.path.basename(out_fc), workspace)
            out_fc = os.path.join(workspace, valid_name)
            sr = arcpy.Describe(airfield_layer).spatialReference
            arcpy.CreateFeatureclass_management(workspace, valid_name, "POINT", spatial_reference=sr)

            # Add fields for aircraft properties
            arcpy.AddField_management(out_fc, "MDS", "TEXT")
            arcpy.AddField_management(out_fc, "LENGTH", "DOUBLE")
            arcpy.AddField_management(out_fc, "WINGSPAN", "DOUBLE")
            arcpy.AddField_management(out_fc, "Aircraft_Footprint", "DOUBLE")

            # Get the airfield data
            airfield_where_clause = f"AFLD_NAME = '{afld_name}' AND OBJECTID = {object_id}"
            with arcpy.da.SearchCursor(airfield_layer, ["SHAPE@", "LENGTH", "WIDTH", "LATITUDE", "LONGITUDE", "LCN"], airfield_where_clause) as cursor:
                for row in cursor:
                    airfield_shape, apron_length, apron_width, start_lat, start_lon, apron_lcn = row
                    # Convert to float and handle potential None values
                    apron_length = float(apron_length) if apron_length is not None else 0
                    apron_width = float(apron_width) if apron_width is not None else 0
                    start_lat = float(start_lat) if start_lat is not None else 0
                    start_lon = float(start_lon) if start_lon is not None else 0
                    apron_lcn = float(apron_lcn) if apron_lcn is not None else 0
                    break
                else:
                    arcpy.AddError(f"Airfield '{afld_name}' with ObjectID {object_id} not found.")
                    return

            arcpy.AddMessage(f"Airfield data: Length={apron_length}, Width={apron_width}, Lat={start_lat}, Lon={start_lon}, LCN={apron_lcn}")

            # Process each selected aircraft
            with arcpy.da.SearchCursor(in_table, ["MDS", "LENGTH", "WING_SPAN", "ACFT_LCN"]) as search_cursor:
                with arcpy.da.InsertCursor(out_fc, ["SHAPE@", "MDS", "LENGTH", "WINGSPAN", "Aircraft_Footprint"]) as insert_cursor:
                    for row in search_cursor:
                        mds, length, wingspan, aircraft_lcn = row
                        if mds in selected_aircraft:
                            # Convert to float and handle potential None values
                            length = float(length) if length is not None else 0
                            wingspan = float(wingspan) if wingspan is not None else 0
                            aircraft_lcn = float(aircraft_lcn) if aircraft_lcn is not None else 0

                            if length == 0 or wingspan == 0 or aircraft_lcn == 0:
                                arcpy.AddWarning(f"Skipping {mds}: Invalid LENGTH, WING_SPAN, or ACFT_LCN.")
                                continue

                            arcpy.AddMessage(f"Processing aircraft {mds}: Length={length}, Wingspan={wingspan}, LCN={aircraft_lcn}")

                            # Check LCN compatibility
                            if aircraft_lcn > apron_lcn:
                                arcpy.AddWarning(f"Aircraft {mds} LCN ({aircraft_lcn}) exceeds apron LCN ({apron_lcn}). Placement may not be suitable.")

                            aircraft_footprint = length * wingspan
                            arcpy.AddMessage(f"Placing {quantity_per_aircraft} {mds} aircraft (Footprint: {aircraft_footprint})")

                            points_placed = 0
                            row_index = 0
                            col_index = 0

                            while points_placed < quantity_per_aircraft:
                                # Calculate position based on grid
                                x = start_lon + (col_index * (wingspan + buffer_distance) / 111111)  # Approximate conversion
                                y = start_lat + (row_index * (length + buffer_distance) / 111111)  # Approximate conversion

                                # Check if we're still within a reasonable range of the apron
                                if (abs(x - start_lon) > (apron_width / 111111) * 2) or (abs(y - start_lat) > (apron_length / 111111) * 2):
                                    arcpy.AddWarning(f"Exceeded reasonable placement area. Placed {points_placed} out of {quantity_per_aircraft} {mds} aircraft.")
                                    break

                                point = arcpy.Point(x, y)
                                insert_cursor.insertRow([point, mds, length, wingspan, aircraft_footprint])
                                points_placed += 1
                                arcpy.AddMessage(f"Placed {mds} at ({x}, {y}). Point {points_placed} of {quantity_per_aircraft}")

                                col_index += 1
                                if col_index >= max_per_row:
                                    col_index = 0
                                    row_index += 1

                            if points_placed < quantity_per_aircraft:
                                arcpy.AddWarning(f"Could only place {points_placed} out of {quantity_per_aircraft} {mds} aircraft.")
                            else:
                                arcpy.AddMessage(f"Successfully placed all {quantity_per_aircraft} {mds} aircraft.")

                            # Add summary information
                            arcpy.AddMessage(f"Placement summary for {mds}:")
                            arcpy.AddMessage(f"  - Aircraft dimensions: Length={length}, Wingspan={wingspan}")
                            arcpy.AddMessage(f"  - Buffer distance: {buffer_distance}")
                            arcpy.AddMessage(f"  - Max per row: {max_per_row}")
                            arcpy.AddMessage(f"  - Apron dimensions: Length={apron_length}, Width={apron_width}")
                            arcpy.AddMessage(f"  - Total aircraft placed: {points_placed}")
                            arcpy.AddMessage(f"  - Rows used: {row_index + 1}")
                            arcpy.AddMessage(f"  - Columns in last row: {col_index}")

            arcpy.AddMessage(f"Aircraft points with footprints created in {out_fc}")

        except Exception as e:
            arcpy.AddError(f"An error occurred while creating the points: {str(e)}")
            arcpy.AddError(arcpy.GetMessages())
            arcpy.AddError(traceback.format_exc())

# Class for calculating Maximum On Ground (MOG)
class CalculateMaximumOnGround(object):
    def __init__(self):
        self.label = "Calculate Maximum On Ground"
        self.description = "Calculate the maximum number of aircraft that can park on a specific apron"

    def getParameterInfo(self):
        params = [
            arcpy.Parameter(
                displayName="Input Aircraft Table", 
                name="in_table", 
                datatype="DETable", 
                parameterType="Required", 
                direction="Input"),
            arcpy.Parameter(
                displayName="Airfield Geodatabase", 
                name="gdb", 
                datatype="DEWorkspace", 
                parameterType="Required", 
                direction="Input"),
            arcpy.Parameter(
                displayName="Movement Surface Feature Class", 
                name="surface_fc", 
                datatype="GPFeatureLayer", 
                parameterType="Required", 
                direction="Input"),
            arcpy.Parameter(
                displayName="Movement Surface OBJECTID", 
                name="surface_id", 
                datatype="GPLong", 
                parameterType="Required", 
                direction="Input"),
            arcpy.Parameter(
                displayName="Aircraft MDS", 
                name="mds", 
                datatype="GPString", 
                parameterType="Required", 
                direction="Input"),
            arcpy.Parameter(
                displayName="Wingtip Clearance (ft)", 
                name="wingtip_clearance", 
                datatype="GPDouble", 
                parameterType="Required", 
                direction="Input")
        ]
        params[2].filter.list = ["Polygon"]
        return params

    def execute(self, parameters, messages):
        in_table = parameters[0].valueAsText
        gdb = parameters[1].valueAsText
        surface_fc = parameters[2].valueAsText
        surface_id = parameters[3].value
        mds = parameters[4].valueAsText
        wingtip_clearance = float(parameters[5].valueAsText)

        try:
            # Get airfield data
            airfield_table = os.path.join(gdb, surface_fc)
            required_fields = ["OBJECTID", "LENGTH", "WIDTH", "LCN"]
            
            existing_fields = [f.name for f in arcpy.ListFields(airfield_table)]
            missing_fields = [f for f in required_fields if f not in existing_fields]
            if missing_fields:
                arcpy.AddError(f"Missing required fields in airfield table: {', '.join(missing_fields)}")
                return

            with arcpy.da.SearchCursor(airfield_table, required_fields, f"OBJECTID = {surface_id}") as cursor:
                for row in cursor:
                    apron_length, apron_width, lcn = row[1], row[2], row[3]
                    break
                else:
                    arcpy.AddError(f"Movement Surface with OBJECTID {surface_id} not found.")
                    return

            # Get aircraft data
            with arcpy.da.SearchCursor(in_table, ["MDS", "LENGTH", "WING_SPAN", "ACFT_LCN"]) as cursor:
                for row in cursor:
                    if row[0] == mds:
                        aircraft_length, aircraft_wingspan, aircraft_lcn = row[1], row[2], row[3]
                        break
                else:
                    arcpy.AddError(f"Aircraft {mds} not found in the table.")
                    return

            # Check LCN compatibility
            if aircraft_lcn > lcn:
                arcpy.AddWarning(f"Aircraft LCN ({aircraft_lcn}) exceeds apron LCN ({lcn}). Parking may not be suitable.")

            # Calculate MOG
            interior_taxi_width = 30 + aircraft_wingspan + 30
            peripheral_taxi_width = 50 + aircraft_wingspan / 2 + 37.5
            mog = self.calculate_mog(apron_length, apron_width, aircraft_length, aircraft_wingspan,
                                     interior_taxi_width, peripheral_taxi_width, wingtip_clearance)
            
            arcpy.AddMessage(f"Maximum On Ground (MOG) for {mds} on selected apron: {mog}")
        
        except Exception as e:
            arcpy.AddError(f"An error occurred while calculating MOG: {str(e)}")
            arcpy.AddError(arcpy.GetMessages())

    def calculate_mog(self, apron_length, apron_width, aircraft_length, aircraft_wingspan,
                      interior_taxi_width, peripheral_taxi_width, wingtip_clearance):
        available_length = apron_length
        rows = math.floor((available_length + interior_taxi_width) / (aircraft_length + interior_taxi_width))
        available_width = max(0, apron_width - peripheral_taxi_width)
        cols = math.floor((available_width + wingtip_clearance) / (aircraft_wingspan + wingtip_clearance))
        parking_available_1 = rows * cols

        available_length = apron_width
        rows = math.floor((available_length + interior_taxi_width) / (aircraft_length + interior_taxi_width))
        available_width = max(0, apron_length - peripheral_taxi_width)
        cols = math.floor((available_width + wingtip_clearance) / (aircraft_wingspan + wingtip_clearance))
        parking_available_2 = rows * cols

        return max(parking_available_1, parking_available_2)