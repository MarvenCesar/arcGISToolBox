# -*- coding: utf-8 -*-
import arcpy
import pandas as pd
import os
import math

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

# Class for calculating aircraft footprints
class CalculateAircraftFootprint(object):
    def __init__(self):
        self.label = "Create Aircraft Symbol Layer"
        self.description = "Create points for selected aircraft at a specific location with corresponding properties, including aircraft footprint"

    def getParameterInfo(self):
        params = [
            arcpy.Parameter(
                displayName="Input Aircraft Table", 
                name="in_table", 
                datatype="DETable", 
                parameterType="Required", 
                direction="Input"),
            arcpy.Parameter(
                displayName="Output Feature Class", 
                name="out_fc", 
                datatype="DEFeatureClass", 
                parameterType="Required", 
                direction="Output"),
            arcpy.Parameter(
                displayName="Longitude for Points",  # X-coordinate, comes first
                name="longitude", 
                datatype="GPDouble", 
                parameterType="Required", 
                direction="Input"),
            arcpy.Parameter(
                displayName="Latitude for Points",  # Y-coordinate, comes second
                name="latitude", 
                datatype="GPDouble", 
                parameterType="Required", 
                direction="Input"),
            arcpy.Parameter(
                displayName="Aircraft Names (MDS)",  # List of aircraft names to filter
                name="aircraft_names", 
                datatype="GPString", 
                parameterType="Required", 
                direction="Input",
                multiValue=True),  # Allow multiple aircraft names
            arcpy.Parameter(
                displayName="Quantity for Each Aircraft",  # Quantity of each aircraft to display
                name="quantity_per_aircraft", 
                datatype="GPLong", 
                parameterType="Required", 
                direction="Input")
        ]
        return params

    def execute(self, parameters, messages):
        in_table = parameters[0].valueAsText
        out_fc = parameters[1].valueAsText
        base_longitude = float(parameters[2].valueAsText)  # Longitude for initial point placement (X)
        base_latitude = float(parameters[3].valueAsText)  # Latitude for initial point placement (Y)
        selected_aircraft = parameters[4].values  # List of selected aircraft names (MDS)
        quantity_per_aircraft = int(parameters[5].valueAsText)  # Quantity of each aircraft to display

        try:
            # Validate the output feature class name
            workspace = os.path.dirname(out_fc)
            valid_name = arcpy.ValidateTableName(os.path.basename(out_fc), workspace)
            out_fc = os.path.join(workspace, valid_name)

            # Create the feature class as a point layer
            sr = arcpy.SpatialReference(4326)  # WGS 1984 (or any other appropriate CRS)
            arcpy.CreateFeatureclass_management(workspace, valid_name, "POINT", spatial_reference=sr)

            # Add fields for the aircraft properties and footprint
            arcpy.AddField_management(out_fc, "MDS", "TEXT")
            arcpy.AddField_management(out_fc, "LENGTH", "DOUBLE")
            arcpy.AddField_management(out_fc, "WINGSPAN", "DOUBLE")
            arcpy.AddField_management(out_fc, "Aircraft_Footprint", "DOUBLE")  # New field for footprint

            # Filter the selected aircraft and calculate footprint
            with arcpy.da.SearchCursor(in_table, ["MDS", "LENGTH", "WING_SPAN"]) as search_cursor:
                with arcpy.da.InsertCursor(out_fc, ["SHAPE@", "MDS", "LENGTH", "WINGSPAN", "Aircraft_Footprint"]) as insert_cursor:
                    i = 0  # Counter for placing the points slightly apart

                    for row in search_cursor:
                        mds, length, wingspan = row

                        # Only process if MDS is in the list of selected aircraft
                        if mds in selected_aircraft:
                            if length is None or wingspan is None:
                                arcpy.AddWarning(f"Skipping {mds}: Missing LENGTH or WING_SPAN.")
                                continue  # Skip records with missing values

                            # Calculate the footprint (area of the rectangle)
                            aircraft_footprint = length * wingspan

                            # For each aircraft, insert the specified number of points
                            for _ in range(quantity_per_aircraft):
                                # Offset each point slightly from the base latitude/longitude to avoid overlap
                                lat_offset = base_latitude + (i * 0.0001)  # Offset by a small value for each aircraft
                                lon_offset = base_longitude + (i * 0.0001)  # Offset similarly in longitude

                                # Create a point geometry at the specific location
                                point = arcpy.Point(lon_offset, lat_offset)

                                # Insert the point along with its properties and footprint
                                insert_cursor.insertRow([point, mds, length, wingspan, aircraft_footprint])

                                # Debugging message
                                arcpy.AddMessage(f"Created point for {mds} at ({lon_offset}, {lat_offset}) with LENGTH: {length}, WING_SPAN: {wingspan}, and Footprint: {aircraft_footprint}")
                                
                                i += 1  # Increment offset for the next point

            arcpy.AddMessage(f"Aircraft points with footprints created in {out_fc}")

        except Exception as e:
            arcpy.AddError(f"An error occurred while creating the points: {str(e)}")
            arcpy.AddError(arcpy.GetMessages())







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

