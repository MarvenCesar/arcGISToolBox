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
        self.tools = ['''ImportAircraftData''',
                      CalculateAircraftFootprint,
                      CalculateMaximumOnGround
                      ]

class CalculateAircraftFootprint(object):
    def __init__(self):
        self.label = "Create Aircraft Symbol Layer (Aircraft-Shaped Polygons)"
        self.description = "Create polygon footprints resembling aircraft shapes for a selected aircraft at specified airfield location."

    def getParameterInfo(self):
        params = [
            arcpy.Parameter(
                displayName="Input Aircraft Table",
                name="in_table",
                datatype="GPTableView",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Input Airfield Layer",
                name="airfield_layer",
                datatype="GPLayer",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Airfield ObjectID",
                name="object_id",
                datatype="GPLong",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Aircraft Name (MDS)",
                name="aircraft_name",
                datatype="GPString",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Quantity of Aircraft",
                name="quantity_of_aircraft",
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
                displayName="Output Feature Class",
                name="out_fc",
                datatype="DEFeatureClass",
                parameterType="Required",
                direction="Output")
        ]

        params[2].parameterDependencies = [params[1].name]
        
        return params
    
    def updateParameters(self, parameters):
        airfield_layer = parameters[1]
        object_id_param = parameters[2]

        if airfield_layer.altered and not airfield_layer.hasBeenValidated:
            if airfield_layer.valueAsText:
                try:
                    # Get the list of OBJECTIDs from the airfield layer
                    with arcpy.da.SearchCursor(airfield_layer.valueAsText, ["OBJECTID"]) as cursor:
                        object_ids = [row[0] for row in cursor]
                    
                    object_id_param.filter.list = object_ids
                    if not object_id_param.value:
                        object_id_param.value = object_ids[0] if object_ids else None
                except Exception as e:
                    arcpy.AddError(f"An error occurred while retrieving OBJECTIDs: {str(e)}")
        return

    def create_aircraft_shape(self, x_start, y_start, length, wingspan, angle):
        # Define proportions
        fuselage_width = length * 0.1
        nose_length = length * 0.2
        tail_length = length * 0.15
        wing_sweep = length * 0.1
        tail_sweep = length * 0.05

        corners = [
            arcpy.Point(x_start, y_start + length/2),  # Nose tip
            arcpy.Point(x_start - fuselage_width/2, y_start + length/2 - nose_length),  # Nose left
            arcpy.Point(x_start - wingspan/2, y_start + wing_sweep),  # Left wingtip front
            arcpy.Point(x_start - wingspan/2, y_start),  # Left wingtip middle
            arcpy.Point(x_start - wingspan/2, y_start - wing_sweep),  # Left wingtip rear
            arcpy.Point(x_start - fuselage_width/2, y_start - length/2 + tail_length),  # Fuselage left before tail
            arcpy.Point(x_start - wingspan/4, y_start - length/2),  # Left tail tip
            arcpy.Point(x_start, y_start - length/2 - tail_sweep),  # Tail bottom tip
            arcpy.Point(x_start + wingspan/4, y_start - length/2),  # Right tail tip
            arcpy.Point(x_start + fuselage_width/2, y_start - length/2 + tail_length),  # Fuselage right before tail
            arcpy.Point(x_start + wingspan/2, y_start - wing_sweep),  # Right wingtip rear
            arcpy.Point(x_start + wingspan/2, y_start),  # Right wingtip middle
            arcpy.Point(x_start + wingspan/2, y_start + wing_sweep),  # Right wingtip front
            arcpy.Point(x_start + fuselage_width/2, y_start + length/2 - nose_length),  # Nose right
            arcpy.Point(x_start, y_start + length/2)  # Back to nose tip
        ]

        # Rotate the aircraft shape based on the angle
        angle_rad = math.radians(angle)
        cos_angle = math.cos(angle_rad)
        sin_angle = math.sin(angle_rad)

        # rotated_corners = []
        # for corner in corners:
        #     x_shifted = corner.X - x_start
        #     y_shifted = corner.Y - y_start
        #     x_rotated = x_shifted * cos_angle - y_shifted * sin_angle
        #     y_rotated = x_shifted * sin_angle + y_shifted * cos_angle
        #     rotated_corners.append(arcpy.Point(x_start + x_rotated, y_start + y_rotated))

        return corners

    def execute(self, parameters, messages):
        # Retrieve parameters
        in_table = parameters[0].valueAsText
        airfield_layer = parameters[1].valueAsText
        object_id = int(parameters[2].valueAsText)
        selected_aircraft = parameters[3].valueAsText
        quantity_of_aircraft = int(parameters[4].valueAsText)
        buffer_distance = float(parameters[5].valueAsText)
        out_fc = parameters[6].valueAsText

        try:
            # Validate and create output feature class
            workspace = os.path.dirname(out_fc)
            valid_name = arcpy.ValidateTableName(os.path.basename(out_fc), workspace)
            out_fc = os.path.join(workspace, valid_name)
            sr = arcpy.Describe(airfield_layer).spatialReference
            arcpy.CreateFeatureclass_management(workspace, valid_name, "POLYGON", spatial_reference=sr)

            # Add fields for aircraft properties
            arcpy.AddField_management(out_fc, "MDS", "TEXT")
            arcpy.AddField_management(out_fc, "LENGTH", "DOUBLE")
            arcpy.AddField_management(out_fc, "WINGSPAN", "DOUBLE")
            arcpy.AddField_management(out_fc, "FOOTPRINT", "DOUBLE")

            # Check if column with shape angles exists in the airfield layer

            if "SHAPE_ANGLE" not in [columns.name for columns in arcpy.ListFields(airfield_layer)]:
                arcpy.AddField_management(airfield_layer, "SHAPE_ANGLE", "DOUBLE")
                arcpy.AddMessage("Field 'SHAPE_ANGLE' added to the input table.")
            else:
                arcpy.AddMessage("Field 'SHAPE_ANGLE' already exists in the input table.")

                            
            # Calculate the angle of the selected object of the airfield
            arcpy.CalculatePolygonMainAngle_cartography(airfield_layer, "SHAPE_ANGLE", "GEOGRAPHIC")

            # Get the airfield data (location and size)
            with arcpy.da.SearchCursor(airfield_layer, ["LENGTH", "WIDTH", "SHAPE_ANGLE", "LATITUDE", "LONGITUDE", "LCN"], f"OBJECTID = {object_id}") as search_cursor:
                for row in search_cursor:
                    apron_length, apron_width, apron_angle, start_lat, start_lon, apron_lcn = row
                    apron_length = float(apron_length) if apron_length is not None else 0
                    apron_width = float(apron_width) if apron_width is not None else 0
                    apron_angle = float(apron_angle) if apron_angle is not None else 0
                    start_lat = float(start_lat) if start_lat is not None else 0
                    start_lon = float(start_lon) if start_lon is not None else 0
                    apron_lcn = float(apron_lcn) if apron_lcn is not None else 0
                    break

            arcpy.AddMessage(f"Airfield data: Length = {apron_length}, Width = {apron_width}, Angle = {apron_angle}, Lat = {start_lat}, Lon = {start_lon}, LCN = {apron_lcn}")

            # Process the selected aircraft
            with arcpy.da.SearchCursor(in_table, ["MDS", "LENGTH", "WING_SPAN", "ACFT_LCN"]) as search_cursor:
                for row in search_cursor:
                    mds, length, wingspan, aircraft_lcn = row
                    if mds == selected_aircraft:
                        length = float(length) if length is not None else 0
                        wingspan = float(wingspan) if wingspan is not None else 0
                        aircraft_lcn = float(aircraft_lcn) if aircraft_lcn is not None else 0

                        if length == 0 or wingspan == 0 or aircraft_lcn == 0:
                            arcpy.AddError(f"Invalid data for {mds}: LENGTH = {length}, WING_SPAN = {wingspan}, ACFT_LCN = {aircraft_lcn}")
                            return

                        arcpy.AddMessage(f"Processing aircraft {mds}: Length = {length}, Wingspan = {wingspan}, LCN = {aircraft_lcn}")

                        if aircraft_lcn > apron_lcn:
                            arcpy.AddWarning(f"Aircraft {mds} LCN ({aircraft_lcn}) exceeds apron LCN ({apron_lcn}). Placement may not be suitable.")

                        length_in_degrees = length / 364000  # Approximate conversion from feet to degrees latitude
                        wingspan_in_degrees = wingspan / 364000

                        max_per_row = int(apron_width / (wingspan + buffer_distance))
                        max_rows = int(apron_length / (length + buffer_distance))

                        arcpy.AddMessage(f"Calculated maximum per row: {max_per_row}, and maximum rows: {max_rows}")

                        # Calculate the total width and height of the aircraft layout
                        total_width = max_per_row * (wingspan_in_degrees + buffer_distance / 364000)
                        total_height = max_rows * (length_in_degrees + buffer_distance / 364000)

                        arcpy.AddMessage(f"Width taken by aircrafts: {total_width}, and height: {total_height}")

                        # Calculate the offset to center the layout
                        x_offset = total_width / 2
                        #y_offset = total_height / 2

                        with arcpy.da.InsertCursor(out_fc, ["SHAPE@", "MDS", "LENGTH", "WINGSPAN", "FOOTPRINT"]) as insert_cursor:
                            points_placed = 0
                            row_index = 0
                            col_index = 0

                            while points_placed < quantity_of_aircraft:
                                x_start = start_lon - x_offset + (col_index * (wingspan_in_degrees + buffer_distance / 364000))
                                y_start = start_lat + (row_index * (length_in_degrees + buffer_distance / 364000))

                                # Create the aircraft shape
                                aircraft_shapes = self.create_aircraft_shape(x_start, y_start, length_in_degrees, wingspan_in_degrees, apron_angle)

                                # Create the polygon
                                polygon = arcpy.Polygon(arcpy.Array(aircraft_shapes), sr)
                                insert_cursor.insertRow([polygon, mds, length, wingspan, length * wingspan])
                                points_placed += 1

                                arcpy.AddMessage(f"Placed {mds} at ({x_start}, {y_start}) - {points_placed}/{quantity_of_aircraft}")

                                # Move to the next column or row
                                col_index += 1
                                if col_index >= max_per_row:
                                    col_index = 0
                                    row_index += 1

                        break
                else:
                    arcpy.AddError(f"Aircraft {selected_aircraft} not found in the input table.")
                    return

            arcpy.AddMessage(f"Aircraft polygons created in {out_fc}")

        except Exception as e:
            arcpy.AddError(f"An error occurred while creating the polygons: {str(e)}")
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