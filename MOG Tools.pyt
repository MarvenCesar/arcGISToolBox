import arcpy
import os
import math
import traceback

class Toolbox(object):
    def __init__(self):
        self.label = "Aircraft MOG Optimization"
        self.alias = "AircraftMOG"
        self.tools = [CalculateMaximumOnGround]

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
            # Use Genetic Algorithm
            arcpy.Parameter(
                displayName="Use Genetic Algorithm",
                name="use_genetic_algorithm",
                datatype="GPBoolean",
                parameterType="Optional",
                direction="Input"),
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

        # Determine taxiway widths based on wingspan
        if aircraft_wingspan >= 110:  # Condition for larger aircraft
            interior_taxi_width = 30 + aircraft_wingspan + 30  # Larger aircraft calculation
        else:
            interior_taxi_width = 20 + aircraft_wingspan + 20  # Smaller aircraft calculation

        peripheral_taxi_width = 50 + (aircraft_wingspan / 2) + 37.5  # Peripheral Taxi Width

        return interior_taxi_width, peripheral_taxi_width

    def calculate_parking_available(self, apron_length, apron_width, aircraft_length, aircraft_wingspan, interior_taxi_width, peripheral_taxi_width, wingtip_between_parked):

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
        # Use the CalculatePolygonMainAngle tool to determine the airfield orientation
        if "Shape_Angle" not in [columns.name for columns in arcpy.ListFields(airfield_shape)]:
            arcpy.management.AddField(airfield_shape, "Shape_Angle", "FLOAT")
            arcpy.AddMessage("Field 'Shape_Angle' added to the input table.")
        else:
            arcpy.AddMessage("Field 'Shape_Angle' already exists in the input table.")

        arcpy.cartography.CalculatePolygonMainAngle(airfield_shape, "Shape_Angle", "ARITHMETIC")

        with arcpy.da.SearchCursor(airfield_shape, ["Shape_Angle"]) as cursor:
            for row in cursor:
                orientation_angle = row[0]

        return orientation_angle

    def convert_and_verify_dimensions(self, airfield_layer, airfield_name, target_sr_code=32617):
        try:
            # Define the target spatial reference
            target_sr = arcpy.SpatialReference(target_sr_code)

            # Get the current spatial reference of the airfield layer
            current_sr = arcpy.Describe(airfield_layer).spatialReference
            arcpy.AddMessage(f"Current spatial reference: {current_sr.name}, Units: {current_sr.linearUnitName}")

            # Define the where_clause to filter by airfield name
            where_clause = f"AFLD_NAME = '{airfield_name}'"

            # Retrieve original dimensions using the where_clause
            original_length = None
            original_width = None
            with arcpy.da.SearchCursor(airfield_layer, ["LENGTH", "WIDTH"], where_clause) as cursor:
                for row in cursor:
                    original_length, original_width = row
                    arcpy.AddMessage(f"Retrieved dimensions: Length={original_length}, Width={original_width}")
                    break
                else:
                    arcpy.AddError(f"Airfield '{airfield_name}' not found or dimensions not available.")
                    return

            # Check if conversion is necessary
            if current_sr.factoryCode != target_sr_code:
                arcpy.AddMessage(f"Converting airfield layer from {current_sr.name} to {target_sr.name}...")

                # Define the output path for the converted layer
                converted_layer = arcpy.env.scratchGDB + "/converted_airfield"

                # Perform the conversion
                arcpy.management.Project(airfield_layer, converted_layer, target_sr)

                # Verify the conversion
                converted_sr = arcpy.Describe(converted_layer).spatialReference
                arcpy.AddMessage(f"Converted spatial reference: {converted_sr.name}, Units: {converted_sr.linearUnitName}")

                # Check if units need conversion
                if current_sr.linearUnitName != target_sr.linearUnitName:
                    arcpy.AddMessage("Units differ between current and target spatial reference. Applying conversion factor...")

                    # Example conversion factor from feet to meters
                    conversion_factor = 0.3048 if 'foot' in current_sr.linearUnitName.lower() else 1.0

                    # Adjust dimensions based on conversion factor
                    with arcpy.da.UpdateCursor(converted_layer, ["LENGTH", "WIDTH"]) as cursor:
                        for row in cursor:
                            if row[0] is not None:
                                row[0] = original_length * conversion_factor  # Convert length
                            if row[1] is not None:
                                row[1] = original_width * conversion_factor  # Convert width
                            cursor.updateRow(row)

                arcpy.AddMessage("Conversion and verification successful. Using converted airfield layer.")
                return converted_layer
            else:
                # If no conversion is needed, return the original layer
                return airfield_layer

        except Exception as e:
            arcpy.AddError(f"An error occurred: {str(e)}")
            arcpy.AddError(traceback.format_exc())
            return airfield_layer

    def genetic_algorithm(self, apron_length, apron_width, aircraft_length, aircraft_wingspan, interior_taxi_width, peripheral_taxi_width, wingtip_between_parked, population_size=50, generations=100, mutation_rate=0.01):
        import random

        def create_chromosome():
            num_rows = max(1, math.floor((apron_length + interior_taxi_width) / (aircraft_length + interior_taxi_width)))
            num_cols = max(1, math.floor((apron_width + wingtip_between_parked) / (aircraft_wingspan + wingtip_between_parked)))
            return [(random.randint(0, num_rows - 1), random.randint(0, num_cols - 1)) for _ in range(num_rows * num_cols)]

        def fitness(chromosome):
            occupied_positions = set()
            for row, col in chromosome:
                if (row, col) in occupied_positions:
                    continue
                occupied_positions.add((row, col))
            return len(occupied_positions)

        def crossover(parent1, parent2):
            crossover_point = random.randint(0, len(parent1) - 1)
            child1 = parent1[:crossover_point] + parent2[crossover_point:]
            child2 = parent2[:crossover_point] + parent1[crossover_point:]
            return child1, child2

        def mutate(chromosome):
            if random.random() < mutation_rate:
                index = random.randint(0, len(chromosome) - 1)
                num_rows = max(1, math.floor((apron_length + interior_taxi_width) / (aircraft_length + interior_taxi_width)))
                num_cols = max(1, math.floor((apron_width + wingtip_between_parked) / (aircraft_wingspan + wingtip_between_parked)))
                chromosome[index] = (random.randint(0, num_rows - 1), random.randint(0, num_cols - 1))
            return chromosome

        population = [create_chromosome() for _ in range(population_size)]

        for generation in range(generations):
            population = sorted(population, key=fitness, reverse=True)
            new_population = population[:population_size // 2]

            while len(new_population) < population_size:
                parent1, parent2 = random.sample(population[:population_size // 2], 2)
                child1, child2 = crossover(parent1, parent2)
                new_population.append(mutate(child1))
                new_population.append(mutate(child2))

            population = new_population

        best_solution = max(population, key=fitness)
        return best_solution, fitness(best_solution)

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
        use_genetic_algorithm = parameters[10].value  # Use Genetic Algorithm parameter

        try:
            # Define the target spatial reference
            target_sr_code = 26917
            target_sr = arcpy.SpatialReference(target_sr_code)

            # Get the current spatial reference of the airfield layer
            current_sr = arcpy.Describe(airfield_layer).spatialReference
            arcpy.AddMessage(f"Current spatial reference: {current_sr.name}, Units: {current_sr.linearUnitName}")

            # Retrieve original dimensions using the where_clause
            original_length = None
            original_width = None
            with arcpy.da.SearchCursor(airfield_layer, ["LENGTH", "WIDTH"]) as cursor:
                for row in cursor:
                    original_length, original_width = row
                    arcpy.AddMessage(f"Original dimensions before conversion: Length={original_length}, Width={original_width}")
                    break
                else:
                    arcpy.AddError(f"Airfield '{airfield_layer}' not found.")
                    return

            # Check if conversion is necessary
            if current_sr.factoryCode != target_sr_code:
                arcpy.AddMessage(f"Converting airfield layer from {current_sr.name} to {target_sr.name}...")

                # Define the output path for the converted layer
                converted_layer = arcpy.env.scratchGDB + "/converted_airfield"

                # Perform the conversion
                arcpy.Project_management(airfield_layer, converted_layer, target_sr)

                # Verify the conversion
                converted_sr = arcpy.Describe(converted_layer).spatialReference
                arcpy.AddMessage(f"Converted spatial reference: {converted_sr.name}, Units: {converted_sr.linearUnitName}")

                # Check if units need conversion
                if current_sr.linearUnitName != target_sr.linearUnitName:
                    arcpy.AddMessage("Units differ between current and target spatial reference. Applying conversion factor...")

                    # Example conversion factor from feet to meters
                    conversion_factor = 0.3048 if 'foot' in current_sr.linearUnitName.lower() else 1.0

                    # Adjust dimensions based on conversion factor
                    with arcpy.da.UpdateCursor(converted_layer, ["LENGTH", "WIDTH"]) as cursor:
                        for row in cursor:
                            if row[0] is not None:
                                row[0] = original_length * conversion_factor  # Convert length
                            if row[1] is not None:
                                row[1] = original_width * conversion_factor  # Convert width
                            cursor.updateRow(row)

                arcpy.AddMessage("Conversion and verification successful. Using converted airfield layer.")
                airfield_layer = converted_layer

            # Ensure the airfield layer is in a projected coordinate system
            final_sr = arcpy.Describe(airfield_layer).spatialReference
            if final_sr.type != 'Projected':
                arcpy.AddError("Airfield layer must be in a projected coordinate system with linear units (e.g., meters or feet).")
                return

            sr = final_sr

            # Get airfield geometry
            with arcpy.da.SearchCursor(airfield_layer, ["SHAPE@", "LENGTH", "WIDTH"], spatial_reference=sr) as cursor:
                for row in cursor:
                    airfield_shape, default_length, default_width = row
                    arcpy.AddMessage(f"Retrieved airfield dimensions: Length={default_length}, Width={default_width}")
                    break
                else:
                    arcpy.AddError(f"Airfield '{airfield_layer}' not found.")
                    return

            # Determine orientation angle using the new method
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

            if use_genetic_algorithm:
                arcpy.AddMessage("Running genetic algorithm for MOG calculation...")
                best_solution, best_fitness = self.genetic_algorithm(
                    apron_length, apron_width, aircraft_length, aircraft_wingspan,
                    interior_taxi_width, peripheral_taxi_width, wingtip_between_parked
                )
                arcpy.AddMessage(f"Genetic Algorithm Result: Best Fitness = {best_fitness}")
                return

            # Calculate parking availability
            (parking_available, num_rows_standard, num_cols_standard, num_rows_rotated, num_cols_rotated,
             parking_available_I, parking_available_II) = self.calculate_parking_available(
                apron_length, apron_width, aircraft_length, aircraft_wingspan, interior_taxi_width, peripheral_taxi_width, wingtip_between_parked)

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
                    if airfield_shape.contains(aircraft_polygon):
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
            self.create_constraint_polygons(out_constraint_fc, sr, aircraft_records, dx_units, dy_units)

        except Exception as e:
            arcpy.AddError(f"An error occurred during the MOG calculation: {str(e)}")
            arcpy.AddError(traceback.format_exc())

    def create_constraint_polygons(self, out_constraint_fc, sr, aircraft_records, constraint_wingspan, constraint_length):
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
