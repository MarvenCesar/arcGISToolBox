import arcpy
import pandas as pd
import os
import math
import traceback
import random

class Toolbox(object):
    def __init__(self):
        self.label = "Aircraft MOG Optimization"
        self.alias = "AircraftMOG"
        self.tools = [CalculateMaximumOnGround] # ImportAircraftData, CalculateAircraftFootprint, 

class CalculateMaximumOnGround(object):
    def __init__(self):
        self.label = "Calculate Maximum On Ground"
        self.description = "Calculates the maximum number of aircraft that can be parked on an airfield, with optional manual input of airfield dimensions and taxiway widths."
        self.canRunInBackground = False

    def getParameterInfo(self):
        params = [
            # Input Aircraft Table
            arcpy.Parameter(
                displayName="Input Aircraft Table",
                name="aircraft_table",
                datatype="GPTableView",
                parameterType="Required",
                direction="Input"),
            # Airfield Layer
            arcpy.Parameter(
                displayName="Airfield Layer",
                name="airfield_layer",
                datatype="GPFeatureLayer",
                parameterType="Required",
                direction="Input"),
            # Airfield Name
            # arcpy.Parameter(
            #     displayName="Airfield Name",
            #     name="airfield_name",
            #     datatype="GPString",
            #     parameterType="Required",
            #     direction="Input"),
            # Select Aircraft (MDS)
            arcpy.Parameter(
                displayName="Select Aircraft (MDS)",
                name="selected_aircraft",
                datatype="GPString",
                parameterType="Required",
                direction="Input"),
            # Enable Advanced Settings
            arcpy.Parameter(
                displayName="Enable Advanced Settings",
                name="enable_advanced_settings",
                datatype="GPBoolean",
                parameterType="Optional",
                direction="Input"),
            # Manual Airfield Length (ft)
            arcpy.Parameter(
                displayName="Manual Airfield Length (ft)",
                name="manual_length",
                datatype="GPDouble",
                parameterType="Optional",
                direction="Input"),
            # Manual Airfield Width (ft)
            arcpy.Parameter(
                displayName="Manual Airfield Width (ft)",
                name="manual_width",
                datatype="GPDouble",
                parameterType="Optional",
                direction="Input"),
            # Manual Interior Taxiway Width (ft)
            arcpy.Parameter(
                displayName="Manual Interior Taxiway Width (ft)",
                name="manual_interior_taxi_width",
                datatype="GPDouble",
                parameterType="Optional",
                direction="Input"),
            # Manual Peripheral Taxiway Width (ft)
            arcpy.Parameter(
                displayName="Manual Peripheral Taxiway Width (ft)",
                name="manual_peripheral_taxi_width",
                datatype="GPDouble",
                parameterType="Optional",
                direction="Input"),
            # Output Aircraft Positions Feature Class
            arcpy.Parameter(
                displayName="Output Aircraft Positions",
                name="out_aircraft_fc",
                datatype="DEFeatureClass",
                parameterType="Required",
                direction="Output"),
            # Output Constraint Polygons Feature Class
            arcpy.Parameter(
                displayName="Output Constraint Polygons",
                name="out_constraint_fc",
                datatype="DEFeatureClass",
                parameterType="Required",
                direction="Output"),
        ]

        # Disable advanced settings parameters by default
        for param in params[4:8]:
            param.enabled = False

        return params

    def updateParameters(self, parameters):
        # Populate the Aircraft Name (MDS) dropdown
        if parameters[0].altered or not parameters[2].altered:
            aircraft_table = parameters[0].valueAsText
            if aircraft_table:
                try:
                    aircraft_names = [row[0] for row in arcpy.da.SearchCursor(aircraft_table, "MDS") if row[0]] # Check if row[0] is valid and not <Null>. If <null>, filter.list breaks
                    parameters[2].filter.list = aircraft_names
                except Exception as e:
                    arcpy.AddError(f"Error loading aircraft names: {e}")

        enable_advanced_settings = parameters[3].value
        for param in parameters[4:8]:
            param.enabled = enable_advanced_settings

        return

    def calculate_taxiway_widths(self, aircraft_length, aircraft_wingspan):
        # Initialize taxiway widths
        interior_taxi_width = 0
        peripheral_taxi_width = 0

        # Wingspan is assumed to be in feet
        wingspan_ft = aircraft_wingspan

        # Determine taxiway widths based on wingspan
        if wingspan_ft >= 110:  # Aircraft with wingspan >= 110 ft
            peripheral_taxi_width = 50  # Wingtip clearance for moving aircraft on peripheral
            interior_taxi_width = 30 + 30  # Wingtip clearance on each side for moving aircraft between parked
        else:  # Aircraft with wingspan < 110 ft
            peripheral_taxi_width = 30  # Wingtip clearance for moving aircraft on peripheral
            interior_taxi_width = 20 + 20  # Wingtip clearance on each side for moving aircraft between parked

        return interior_taxi_width, peripheral_taxi_width

    def calculate_parking_available(self, apron_length, apron_width, aircraft_length, aircraft_wingspan):
        # Calculate the taxiway dimensions based on the aircraft specifications
        if aircraft_wingspan >= 110:  # Condition for larger aircraft
            interior_taxi_width = 30 + aircraft_wingspan + 30  # Larger aircraft calculation
        else:
            interior_taxi_width = 20 + aircraft_wingspan + 20  # Smaller aircraft calculation

        peripheral_taxi_width = 50 + (aircraft_wingspan / 2) + 37.5  # Peripheral Taxi Width
        wingtip_between_parked = 25  # Space between parked aircraft

        # I. Standard apron configuration
        # 1) Determine number of rows
        available_length = apron_length
        num_rows = max(1, math.floor((available_length + interior_taxi_width) / (aircraft_length + interior_taxi_width)))
        arcpy.AddMessage(f"Standard Configuration - Available Length: {available_length} ft, Aircraft Length: {aircraft_length} ft, Interior Taxi Width: {interior_taxi_width} ft")
        arcpy.AddMessage(f"Calculated Rows: {num_rows}")

        # 2) Determine number of columns
        available_width = max(0, apron_width - peripheral_taxi_width)
        num_cols = max(1, math.floor((available_width + wingtip_between_parked) / (aircraft_wingspan + wingtip_between_parked)))
        arcpy.AddMessage(f"Standard Configuration - Available Width: {available_width} ft, Aircraft Wingspan: {aircraft_wingspan} ft, Wingtip Between Parked: {wingtip_between_parked} ft")
        arcpy.AddMessage(f"Calculated Columns: {num_cols}")

        # 3) Parking available in standard configuration
        parking_available_I = num_rows * num_cols

        # II. Rotated apron configuration
        # 1) Determine number of rows
        available_length_rotated = apron_width
        num_rows_rotated = max(1, math.floor((available_length_rotated + interior_taxi_width) / (aircraft_length + interior_taxi_width)))
        arcpy.AddMessage(f"Rotated Configuration - Available Length: {available_length_rotated} ft, Aircraft Length: {aircraft_length} ft, Interior Taxi Width: {interior_taxi_width} ft")
        arcpy.AddMessage(f"Calculated Rotated Rows: {num_rows_rotated}")

        # 2) Determine number of columns
        available_width_rotated = max(0, apron_length - peripheral_taxi_width)
        num_cols_rotated = max(1, math.floor((available_width_rotated + wingtip_between_parked) / (aircraft_wingspan + wingtip_between_parked)))
        arcpy.AddMessage(f"Rotated Configuration - Available Width: {available_width_rotated} ft, Aircraft Wingspan: {aircraft_wingspan} ft, Wingtip Between Parked: {wingtip_between_parked} ft")
        arcpy.AddMessage(f"Calculated Rotated Columns: {num_cols_rotated}")

        # 3) Parking available in rotated configuration
        parking_available_II = num_rows_rotated * num_cols_rotated

        # III. Final Parking Available
        final_parking_available = max(parking_available_I, parking_available_II)
        arcpy.AddMessage(f"Final Parking Available (Maximum of Standard and Rotated): {final_parking_available}")

        return final_parking_available, num_rows, num_cols, num_rows_rotated, num_cols_rotated, parking_available_I, parking_available_II

    def calculate_airfield_orientation(self, airfield_shape, sr):
        # Calculate the orientation of the airfield based on the longest edge
        # max_length = 0
        # orientation_angle = 0
        # for row in arcpy.da.SearchCursor(airfield_shape, ["SHAPE@"], spatial_reference=sr):
        #     for part in row:
        #         for i in range(len(part) - 1):
        #             if part[i] and part[i + 1]:
        #                 p1 = part[i]
        #                 p2 = part[i + 1]
        #                 dx = p2.X - p1.X
        #                 dy = p2.Y - p1.Y
        #                 length = math.hypot(dx, dy)
        #                 if length > max_length:
        #                     max_length = length
        #                     orientation_angle = math.degrees(math.atan2(dy, dx))

        # Commented out the above code and replaced with CalculatePolygonMainAngle tool

        if "SHAPE_ANGLE" not in [columns.name for columns in arcpy.ListFields(airfield_shape)]:
            arcpy.management.AddField(airfield_shape, "Shape_Angle", "FLOAT")
            arcpy.AddMessage("Field 'SHAPE_ANGLE' added to the input table.")
        else:
            arcpy.AddMessage("Field 'SHAPE_ANGLE' already exists in the input table.")

        arcpy.cartography.CalculatePolygonMainAngle(airfield_shape, "Shape_Angle", "ARITHMETIC")

        with arcpy.da.SearchCursor(airfield_shape, ["Shape_Angle"]) as cursor:
            for row in cursor:
                orientation_angle = row[0]

        return orientation_angle

    def execute(self, parameters, messages):
        # Extract parameter values
        aircraft_table = parameters[0].valueAsText
        airfield_layer = parameters[1].valueAsText
        selected_aircraft = parameters[2].valueAsText
        enable_advanced_settings = parameters[3].value
        manual_length = parameters[4].value
        manual_width = parameters[5].value
        manual_interior_taxi_width = parameters[6].value
        manual_peripheral_taxi_width = parameters[7].value
        out_aircraft_fc = parameters[8].valueAsText  # Output Aircraft Feature Class
        out_constraint_fc = parameters[9].valueAsText  # Output Constraint Feature Class

        try:
            # Get the spatial reference of the airfield layer
            airfield_sr = arcpy.Describe(airfield_layer).spatialReference

            # Ensure the airfield layer is in a projected coordinate system
            if airfield_sr.type != 'Projected':
                arcpy.AddError("Airfield layer must be in a projected coordinate system with linear units (e.g., meters or feet).")
                return

            sr = airfield_sr

            # Get airfield geometry
            with arcpy.da.SearchCursor(airfield_layer, ["SHAPE@", "LENGTH", "WIDTH"], spatial_reference=sr) as cursor:
                for row in cursor:
                    airfield_shape, default_length, default_width = row
                    break
                else:
                    arcpy.AddError(f"Airfield '{airfield_layer}' not found.")
                    return

            # Determine orientation angle
            angle = self.calculate_airfield_orientation(airfield_layer, sr)
            arcpy.AddMessage(f"Calculated airfield orientation: {angle} degrees")

            # Get airfield dimensions
            if enable_advanced_settings and manual_length is not None and manual_width is not None:
                apron_length = float(manual_length)
                apron_width = float(manual_width)
                arcpy.AddMessage(f"Using manual airfield dimensions: Length={apron_length} ft, Width={apron_width} ft")
            else:
                apron_length = float(default_length)
                apron_width = float(default_width)
                arcpy.AddMessage(f"Retrieved airfield dimensions: Length={apron_length} ft, Width={apron_width} ft")

            # Get aircraft dimensions
            aircraft_data = self.get_aircraft_data(aircraft_table, selected_aircraft)
            if not aircraft_data:
                arcpy.AddError(f"Aircraft '{selected_aircraft}' not found.")
                return
            mds, aircraft_length, aircraft_wingspan, _ = aircraft_data[0]
            arcpy.AddMessage(f"Selected Aircraft: {mds}, Length={aircraft_length} ft, Wingspan={aircraft_wingspan} ft")

            # Use manual taxiway widths if provided
            if enable_advanced_settings and manual_interior_taxi_width is not None and manual_peripheral_taxi_width is not None:
                interior_taxi_width = float(manual_interior_taxi_width)
                peripheral_taxi_width = float(manual_peripheral_taxi_width)
                arcpy.AddMessage(f"Using manual taxiway widths: Interior Taxi Width={interior_taxi_width} ft, Peripheral Taxi Width={peripheral_taxi_width} ft")
            else:
                # Calculate taxiway widths based on aircraft dimensions
                interior_taxi_width, peripheral_taxi_width = self.calculate_taxiway_widths(aircraft_length, aircraft_wingspan)
                arcpy.AddMessage(f"Calculated taxiway widths: Interior Taxi Width={interior_taxi_width} ft, Peripheral Taxi Width={peripheral_taxi_width} ft")

            # Wingtip clearance between parked aircraft
            wingtip_between_parked = 25  # Adjust as necessary
            arcpy.AddMessage(f"Using wingtip clearance between parked aircraft: {wingtip_between_parked} ft")

            # Calculate parking availability
            (parking_available, num_rows_standard, num_cols_standard, num_rows_rotated, num_cols_rotated,
             parking_available_I, parking_available_II) = self.calculate_parking_available(
                apron_length, apron_width, aircraft_length, aircraft_wingspan)

            arcpy.AddMessage(f"Final Parking Available for {mds}: {parking_available}")

            # Decide whether to use standard or rotated configuration based on which has more parking available
            if parking_available == parking_available_I:
                # Standard configuration
                num_rows_to_place = num_rows_standard
                num_cols_to_place = num_cols_standard
                dx = aircraft_wingspan + wingtip_between_parked
                dy = aircraft_length + interior_taxi_width
                arcpy.AddMessage("Using standard configuration for aircraft placement.")
            else:
                # Rotated configuration
                num_rows_to_place = num_rows_rotated
                num_cols_to_place = num_cols_rotated
                dx = aircraft_wingspan + wingtip_between_parked
                dy = aircraft_length + interior_taxi_width
                # Swap dx and dy for rotated configuration
                dx, dy = dy, dx
                angle += 90  # Adjust angle by 90 degrees for rotation
                arcpy.AddMessage("Using rotated configuration for aircraft placement.")

            # Convert dimensions from feet to spatial reference units
            unit_name = sr.linearUnitName.lower()
            if 'foot' in unit_name:
                unit_factor = 1.0  # Units are already in feet
            elif 'meter' in unit_name:
                unit_factor = 0.3048  # Convert feet to meters
            else:
                arcpy.AddError(f"Unsupported linear unit in spatial reference: {sr.linearUnitName}")
                return

            dx_units = dx * unit_factor
            dy_units = dy * unit_factor

            # Create separate feature classes for aircraft and constraints
            # Aircraft Feature Class already defined as out_aircraft_fc
            # Constraint Feature Class defined as out_constraint_fc

            # Ensure feature classes do not already exist. If they do, delete them.
            for fc in [out_aircraft_fc, out_constraint_fc]:
                if arcpy.Exists(fc):
                    arcpy.Delete_management(fc)
                    arcpy.AddMessage(f"Existing feature class '{fc}' deleted.")

            # Create Aircraft Feature Class
            arcpy.CreateFeatureclass_management(os.path.dirname(out_aircraft_fc), os.path.basename(out_aircraft_fc), "POLYGON", spatial_reference=sr)
            # Add necessary fields to Aircraft Feature Class
            arcpy.AddField_management(out_aircraft_fc, "MDS", "TEXT")
            arcpy.AddField_management(out_aircraft_fc, "AircraftID", "LONG")
            arcpy.AddField_management(out_aircraft_fc, "Rotation", "DOUBLE")
            arcpy.AddField_management(out_aircraft_fc, "Type", "TEXT")  # New field to identify type

            # Create Constraint Feature Class
            arcpy.CreateFeatureclass_management(os.path.dirname(out_constraint_fc), os.path.basename(out_constraint_fc), "POLYGON", spatial_reference=sr)
            # Add necessary fields to Constraint Feature Class
            arcpy.AddField_management(out_constraint_fc, "AircraftID", "LONG")  # To associate with aircraft
            arcpy.AddField_management(out_constraint_fc, "Type", "TEXT")  # New field to identify type

            # Determine the starting point (use airfield centroid)
            start_point = airfield_shape.centroid
            arcpy.AddMessage(f"Using airfield centroid as starting point: ({start_point.X}, {start_point.Y})")

            # Generate grid coordinates centered at (0, 0)
            grid_points = []
            for i in range(num_rows_to_place):
                for j in range(num_cols_to_place):
                    x = j * dx_units - ((num_cols_to_place - 1) * dx_units) / 2
                    y = i * dy_units - ((num_rows_to_place - 1) * dy_units) / 2
                    grid_points.append((x, y))

            # Rotate and translate grid points
            angle_rad = math.radians(angle)
            cos_angle = math.cos(angle_rad)
            sin_angle = math.sin(angle_rad)

            aircraft_id = 1
            aircraft_records = []  # To store aircraft data for constraint creation

            with arcpy.da.InsertCursor(out_aircraft_fc, ["SHAPE@", "MDS", "AircraftID", "Rotation", "Type"]) as aircraft_cursor:
                for x_local, y_local in grid_points:
                    x_rotated = x_local * cos_angle - y_local * sin_angle
                    y_rotated = x_local * sin_angle + y_local * cos_angle
                    x_global = start_point.X + x_rotated
                    y_global = start_point.Y + y_rotated

                    # Create aircraft rectangle corners
                    half_length = (aircraft_length * unit_factor) / 2
                    half_wingspan = (aircraft_wingspan * unit_factor) / 2
                    corners = [
                        arcpy.Point(x_global - half_wingspan, y_global + half_length),
                        arcpy.Point(x_global + half_wingspan, y_global + half_length),
                        arcpy.Point(x_global + half_wingspan, y_global - half_length),
                        arcpy.Point(x_global - half_wingspan, y_global - half_length),
                        arcpy.Point(x_global - half_wingspan, y_global + half_length)
                    ]

                    # Create aircraft polygon with rotation
                    polygon_array = arcpy.Array()
                    for corner in corners:
                        # Rotate corner around center point
                        dx = corner.X - x_global
                        dy = corner.Y - y_global
                        rotated_x = x_global + (dx * cos_angle - dy * sin_angle)
                        rotated_y = y_global + (dx * sin_angle + dy * cos_angle)
                        polygon_array.add(arcpy.Point(rotated_x, rotated_y))

                    aircraft_polygon = arcpy.Polygon(polygon_array, sr)

                    # Check if aircraft polygon is within airfield
                    if airfield_shape.contains(aircraft_polygon.centroid):
                        # Insert aircraft polygon
                        aircraft_cursor.insertRow([aircraft_polygon, mds, aircraft_id, angle % 360, "Aircraft"])
                        arcpy.AddMessage(f"Placed {mds} at ({x_global}, {y_global}) - ID: {aircraft_id}/{parking_available}")
                        # Store aircraft data for constraint creation
                        aircraft_records.append({
                            "AircraftID": aircraft_id,
                            "Centroid": aircraft_polygon.centroid,
                            "Length": aircraft_length * unit_factor,
                            "Wingspan": aircraft_wingspan * unit_factor,
                            "Rotation": angle % 360
                        })
                        aircraft_id += 1
                    else:
                        arcpy.AddWarning(f"Aircraft at position ({x_global}, {y_global}) is outside the airfield boundary.")

            arcpy.AddMessage(f"Aircraft positions created in {out_aircraft_fc}, total aircraft placed: {aircraft_id - 1}")

            # After aircraft polygons are created, create constraint polygons
            self.create_constraint_polygons(out_constraint_fc, sr, aircraft_records)

        except Exception as e:
            arcpy.AddError(f"An error occurred during the MOG calculation: {str(e)}")
            arcpy.AddError(traceback.format_exc())

    def create_constraint_polygons(self, out_constraint_fc, sr, aircraft_records):
        """
        Creates constraint polygons based on the aircraft polygons.
        Each constraint polygon is a larger rectangle around the aircraft.

        :param out_constraint_fc: Path to the Constraint Feature Class
        :param sr: Spatial Reference
        :param aircraft_records: List of dictionaries containing aircraft data
        """
        try:
            with arcpy.da.InsertCursor(out_constraint_fc, ["SHAPE@", "AircraftID", "Type"]) as constraint_cursor:
                for record in aircraft_records:
                    aircraft_id = record["AircraftID"]
                    centroid = record["Centroid"]
                    length = record["Length"]
                    wingspan = record["Wingspan"]
                    rotation = record["Rotation"]  # Get the rotation angle

                    # Define the size of the constraint polygon
                    clearance_factor = 1.2  # 20% larger
                    constraint_length = length * clearance_factor
                    constraint_wingspan = wingspan * clearance_factor

                    # Create constraint rectangle corners (unrotated)
                    half_length = constraint_length / 2
                    half_wingspan = constraint_wingspan / 2
                    corners = [
                        arcpy.Point(centroid.X - half_wingspan, centroid.Y + half_length),
                        arcpy.Point(centroid.X + half_wingspan, centroid.Y + half_length),
                        arcpy.Point(centroid.X + half_wingspan, centroid.Y - half_length),
                        arcpy.Point(centroid.X - half_wingspan, centroid.Y - half_length),
                        arcpy.Point(centroid.X - half_wingspan, centroid.Y + half_length)
                    ]

                    # Rotate corners around the centroid
                    angle_rad = math.radians(rotation)
                    cos_angle = math.cos(angle_rad)
                    sin_angle = math.sin(angle_rad)
                    rotated_corners = []
                    for corner in corners:
                        dx = corner.X - centroid.X
                        dy = corner.Y - centroid.Y
                        rotated_x = centroid.X + (dx * cos_angle - dy * sin_angle)
                        rotated_y = centroid.Y + (dx * sin_angle + dy * cos_angle)
                        rotated_corners.append(arcpy.Point(rotated_x, rotated_y))

                    # Create constraint polygon with rotation
                    constraint_polygon = arcpy.Polygon(arcpy.Array(rotated_corners), sr)

                    # Insert constraint polygon
                    constraint_cursor.insertRow([constraint_polygon, aircraft_id, "Constraint"])
                    arcpy.AddMessage(f"Constraint polygon created for Aircraft ID: {aircraft_id}")

            arcpy.AddMessage(f"Constraint polygons created in {out_constraint_fc}")

        except Exception as e:
            arcpy.AddError(f"An error occurred while creating constraint polygons: {str(e)}")
            arcpy.AddError(traceback.format_exc())

    def get_aircraft_data(self, aircraft_table, selected_aircraft):
        aircraft_data = []
        with arcpy.da.SearchCursor(aircraft_table, ["MDS", "LENGTH", "WING_SPAN", "ACFT_LCN"]) as cursor:
            for row in cursor:
                if row[0] == selected_aircraft:
                    aircraft_data.append(row)
                    break
        return aircraft_data

    def get_airfield_dimension(self, airfield_layer, dimension_field):
        with arcpy.da.SearchCursor(airfield_layer, [dimension_field]) as cursor:
            for row in cursor:
                return float(row[0])  # Return the requested dimension as float
        arcpy.AddError(f"Airfield '{airfield_layer}' not found or dimension '{dimension_field}' not available.")
        return None  # Return None if not found