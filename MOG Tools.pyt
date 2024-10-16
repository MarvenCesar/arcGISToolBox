# -*- coding: utf-8 -*-
import random
import arcpy
import pandas as pd
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
        self.description = "Calculates the maximum number of aircraft that can be parked on an airfield."
        self.canRunInBackground = False
        self.optimization_attempts = 20
        self.aircraft_data = None  # Add this line

    def getParameterInfo(self):
        params = [
            arcpy.Parameter( # 0
                displayName="Input Aircraft Table",
                name="in_table",
                datatype="GPTableView",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter( # 1
                displayName="Airfield Layer",
                name="airfield_layer",
                datatype="GPFeatureLayer",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter( # 2
                displayName="Origin Point",
                name="origin_point",
                datatype="GPString",
                parameterType="Optional",
                direction="Input"),
            arcpy.Parameter( # 3
                displayName="Interior Taxiway Width (ft)",
                name="interior_taxi_width",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter( # 4
                displayName="Peripheral Taxiway Width (ft)",
                name="peripheral_taxi_width",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter( # 5
                displayName="Aircraft Wingtip Clearance (ft)",
                name="wingtip_clearance",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter( # 6
                displayName="Select Aircraft (MDS)",
                name="selected_aircraft",
                datatype="GPString",
                parameterType="Required",
                direction="Input",
                multiValue=True),
            arcpy.Parameter( # 7
                displayName="Aircraft Quantities (Comma Separated)",
                name="aircraft_quantities",
                datatype="GPString",
                parameterType="Optional",
                direction="Input"),
            arcpy.Parameter( # 8
                displayName="Apply LCN Compatibility?",
                name="apply_lcn",
                datatype="GPBoolean",
                parameterType="Optional",
                direction="Input"),
            arcpy.Parameter( # 9
                displayName="Optimize Quantities",
                name="optimize_quantities",
                datatype="GPBoolean",
                parameterType="Optional",
                direction="Input"),
            arcpy.Parameter( # 10
                displayName="Minimum Aircraft Quantities (Comma Separated)",
                name="min_quantities",
                datatype="GPString",
                parameterType="Optional",
                direction="Input",
                enabled=False),
            arcpy.Parameter( # 11
                displayName="Number of Optimization Attempts",
                name="optimization_attempts",
                datatype="GPLong",
                parameterType="Optional",
                direction="Input",
                enabled=False),
            arcpy.Parameter( # 12
                displayName="Output Feature Class",
                name="out_fc",
                datatype="DEFeatureClass",
                parameterType="Required",
                direction="Output")
        ]
        params[1].filter.list = ["Polygon"]
        params[2].parameterDependencies = [params[1].name]
        params[6].parameterDependencies = [params[0].name]
        params[9].value = False  # Default to not optimizing
        params[11].value = 20  # Use class attribute for default
        return params

    def updateParameters(self, parameters):

        if parameters[0].altered:
            aircraft_table = parameters[0].valueAsText
            if aircraft_table:
                aircraft_names = [row[0] for row in arcpy.da.SearchCursor(aircraft_table, "MDS") if row[0]] # Check if row[0] is valid and not <Null>. If <null>, filter.list breaks
                parameters[6].filter.list = aircraft_names
        
        if parameters[1].altered:
            airfield_layer_object = parameters[1].value
            airfield_layer_object_name = parameters[1].valueAsText
            if airfield_layer_object:
                # Check if the original airfield_layer_object has more than 4 vertices
                vertex_count = 0
                with arcpy.da.SearchCursor(airfield_layer_object, ["SHAPE@"]) as cursor:
                    for row in cursor:
                        for part in row[0]:
                            vertex_count += len(part)

                if vertex_count > 4:
                    # Create a copy of the airfield_layer_object as a new feature class
                    workspace = os.path.dirname(arcpy.Describe(airfield_layer_object).catalogPath)
                    valid_name = arcpy.ValidateTableName(airfield_layer_object_name + "_generalized", workspace)
                    airfield_layer_object_generalized = os.path.join(workspace, valid_name)
                    sr = arcpy.Describe(airfield_layer_object).spatialReference
                    arcpy.management.CreateFeatureclass(workspace, valid_name, "POLYGON", spatial_reference=sr)

                    # Copy all data from the original to the new feature class
                    arcpy.management.CopyFeatures(airfield_layer_object, airfield_layer_object_generalized)

                    # Generalize the copied feature class
                    arcpy.edit.Generalize(airfield_layer_object_generalized, "10 Feet")
                else:
                    airfield_layer_object_generalized = airfield_layer_object

                # Proceed with gathering the points from the generalized feature class
                points_list = []
                for polygon in arcpy.da.SearchCursor(airfield_layer_object_generalized, ["SHAPE@"]):
                    for part in polygon[0].getPart():
                        i_point = 0
                        for point in part:
                            if point:
                                i_point += 1
                                points_list.append(f"{i_point}: X: {point.X},\tY: {point.Y}")
                            else:
                                break # Break is needed to ensure no inner rings


                parameters[2].filter.list = list(points_list)
    
        optimize_quantities_param = parameters[9]
        min_aircraft_quantities_param = parameters[10]
        optimization_attempts_param = parameters[11]

        optimize_quantities_param.value = False if optimize_quantities_param.value is None else optimize_quantities_param.value
        min_aircraft_quantities_param.enabled = optimize_quantities_param.value
        optimization_attempts_param.enabled = optimize_quantities_param.value

    def execute(self, parameters, messages):
        # Extract parameter values
        aircraft_table = parameters[0].valueAsText
        airfield_layer_object = parameters[1].value # added to work for the airfield layer generation

        airfield_layer_object_origin = parameters[2].valueAsText

        interior_taxi_width = float(parameters[3].value)
        peripheral_taxi_width = float(parameters[4].value)
        wingtip_clearance = float(parameters[5].value)
        selected_aircraft = parameters[6].valueAsText
        aircraft_quantities = parameters[7].valueAsText
        apply_lcn = parameters[8].value
        optimize_quantities = parameters[9].value
        min_quantities = parameters[10].valueAsText

        out_fc = parameters[12].valueAsText

        # Process min_quantities
        if parameters[10].value:  # Assuming min_quantities is the 11th parameter
            min_quantities = parameters[11].valueAsText
        else:
            min_quantities = ','.join(['1'] * len(selected_aircraft))  # Default to 1 for each aircraft type

        # Get airfield data
        apron_length, apron_width, apron_lcn = self.get_airfield_data(airfield_layer_object)

        # Get aircraft data
        aircraft_data = self.get_aircraft_data(aircraft_table, selected_aircraft)

        # Calculate MOG
        result = self.calculate_mog(aircraft_data, apron_length, apron_width, apron_lcn,
                                    interior_taxi_width, peripheral_taxi_width, 
                                    wingtip_clearance, apply_lcn, aircraft_quantities, 
                                    min_quantities, optimize_quantities)

        # Display results
        self.display_results(result, apron_length, apron_width, apron_lcn,
                             interior_taxi_width, peripheral_taxi_width, wingtip_clearance)
        
        self.create_aircraft_layer(airfield_layer_object, aircraft_data, result, out_fc, wingtip_clearance, interior_taxi_width)

    def get_airfield_data(self, airfield_layer_object):
        
        with arcpy.da.SearchCursor(airfield_layer_object, ["LENGTH", "WIDTH", "LCN", "OID@", "SHAPE@"]) as cursor:
            for row in cursor:
                # Print the current polygon or polyline's ID
                arcpy.AddMessage(f"Feature {row[3]}:")

                # Step through each part of the feature
                for part in row[4]:

                    # # Check if the part has more than 4 vertices
                    # if len(part) > 4:

                    #     arcpy.AddMessage("More than 4 points detected. Generalizing the airfield layer object...")
                    #     # Copy the original feature class
                    #     airfield_layer_object_generalized = airfield_layer_object
                    #     arcpy.management.CopyFeatures(airfield_layer_object, airfield_layer_object_generalized)

                    #     # Generalize the copied feature class
                    #     arcpy.edit.Generalize(airfield_layer_object_generalized, "10 Feet")

                    #     with arcpy.da.SearchCursor(airfield_layer_object_generalized, ["SHAPE@"]) as gen_cursor:
                    #         for gen_row in gen_cursor:

                    #             partnum = 0

                    #             # Print the generalized part's points
                    #             for gen_part in gen_row[0]:
                    #                 for gen_pnt in gen_part:
                    #                     if gen_pnt:
                    #                         arcpy.AddMessage("Generalized Point {}: {}, {}".format(gen_part.index(gen_pnt), gen_pnt.X, gen_pnt.Y))
                    #                     else:
                    #                         arcpy.AddMessage("Interior Ring:")

                    #         partnum += 1
                    # else:
                    partnum = 0
                    pointnum = 0
                    # Print the part number
                    arcpy.AddMessage(f"Part {partnum}:")
                    # Print the part's points
                    for pnt in part:
                        if pnt:
                            arcpy.AddMessage(f"Point {pointnum}: {pnt.X}, {pnt.Y}")
                        else:
                            arcpy.AddMessage("Interior Ring:")
                        pointnum += 1

                    partnum += 1

                return float(row[0]), float(row[1]), float(row[2])
            
        arcpy.AddError(f"Object '{airfield_layer_object}' not found.")
        return None, None, None

    def get_aircraft_data(self, aircraft_table, selected_aircraft):
        aircraft_data = []
        with arcpy.da.SearchCursor(aircraft_table, ["MDS", "LENGTH", "WING_SPAN", "ACFT_LCN"]) as cursor:
            for row in cursor:
                if row[0] and row[0] in selected_aircraft:
                    aircraft_data.append(row)
        return aircraft_data

    def calculate_mog(self, aircraft_data, apron_length, apron_width, apron_lcn,
                      interior_taxi_width, peripheral_taxi_width, wingtip_clearance,
                      apply_lcn=False, aircraft_quantities=None, min_quantities=None, optimize_quantities=False):
        self.aircraft_data = aircraft_data  # Store aircraft data for later use
        usable_length = apron_length - (2 * peripheral_taxi_width)
        usable_width = apron_width - (2 * peripheral_taxi_width)
        usable_area = usable_length * usable_width

        if aircraft_quantities:
            requested_quantities = [int(q.strip()) for q in aircraft_quantities.split(',')]
            total_requested_area = sum(qty * length * wingspan 
                                       for (mds, length, wingspan, _), qty in zip(aircraft_data, requested_quantities))
            
            if total_requested_area > usable_area:
                arcpy.AddWarning("WARNING: The requested aircraft quantities exceed the usable apron area.")
                arcpy.AddWarning(f"Total area required by aircraft: {total_requested_area:.2f} sq ft")
                arcpy.AddWarning(f"Usable apron area: {usable_area:.2f} sq ft")
                arcpy.AddWarning("The tool will attempt to place as many aircraft as possible within the available space.")

        if optimize_quantities:
            return self.optimize_placement(aircraft_data, apron_length, apron_width, apron_lcn,
                                           interior_taxi_width, peripheral_taxi_width, wingtip_clearance,
                                           apply_lcn, min_quantities)
        else:
            results = []
            for i, (mds, length, wingspan, aircraft_lcn) in enumerate(aircraft_data):
                if apply_lcn and aircraft_lcn > apron_lcn:
                    arcpy.AddWarning(f"Aircraft {mds} LCN ({aircraft_lcn}) exceeds apron LCN ({apron_lcn}). Skipping this aircraft.")
                    continue

                max_count = self.calculate_max_aircraft(length, wingspan, usable_length, usable_width,
                                                        interior_taxi_width, wingtip_clearance)

                if aircraft_quantities:
                    requested_count = requested_quantities[i]
                    if requested_count > max_count:
                        arcpy.AddWarning(f"Requested quantity for {mds} ({requested_count}) exceeds maximum possible ({max_count}). Placing {max_count} aircraft.")
                        count = max_count
                    else:
                        count = requested_count
                else:
                    count = max_count

                num_rows, per_row, orientation = self.calculate_layout(count, length, wingspan, 
                                                                       usable_length, usable_width,
                                                                       interior_taxi_width, wingtip_clearance)
                results.append((mds, count, num_rows, per_row, orientation, aircraft_lcn))

            return results

    def calculate_max_aircraft(self, length, wingspan, usable_length, usable_width, 
                               interior_taxi_width, wingtip_clearance):
        max_normal = self.calculate_num_rows(usable_width, wingspan, wingtip_clearance, interior_taxi_width) * \
                     self.calculate_aircraft_per_row(usable_length, length)
        max_rotated = self.calculate_num_rows(usable_width, length, wingtip_clearance, interior_taxi_width) * \
                      self.calculate_aircraft_per_row(usable_length, wingspan)
        return max(max_normal, max_rotated)

    def calculate_layout(self, count, length, wingspan, usable_length, usable_width,
                        interior_taxi_width, wingtip_clearance):
        normal_rows = self.calculate_num_rows(usable_width, wingspan, wingtip_clearance, interior_taxi_width)
        normal_per_row = self.calculate_aircraft_per_row(usable_length, length)
        normal_total = normal_rows * normal_per_row

        rotated_rows = self.calculate_num_rows(usable_width, length, wingtip_clearance, interior_taxi_width)
        rotated_per_row = self.calculate_aircraft_per_row(usable_length, wingspan)
        rotated_total = rotated_rows * rotated_per_row

        arcpy.AddMessage(f"Normal Layout: Rows = {normal_rows}, Per Row = {normal_per_row}, Total = {normal_total}")
        arcpy.AddMessage(f"Rotated Layout: Rows = {rotated_rows}, Per Row = {rotated_per_row}, Total = {rotated_total}")

        if rotated_total > normal_total:
            return rotated_rows, rotated_per_row, "rotated"
        else:
            return normal_rows, normal_per_row, "normal"

    def calculate_num_rows(self, usable_width, wingspan, wingtip_clearance, interior_taxi_width):
        space_per_row = wingspan + wingtip_clearance
        total_row_space = usable_width - interior_taxi_width  # Reserve space for one less taxiway than rows
        num_rows = math.floor(total_row_space / space_per_row)
        return max(1, num_rows)  # Ensure at least one row

    def calculate_aircraft_per_row(self, usable_length, aircraft_length):
        return math.floor(usable_length / aircraft_length)

    def display_results(self, results, apron_length, apron_width, apron_lcn,
                        interior_taxi_width, peripheral_taxi_width, wingtip_clearance):
        arcpy.AddMessage("\nSelected Aircraft Dimensions:")
        for mds, length, wingspan, aircraft_lcn in self.aircraft_data:
            arcpy.AddMessage(f"{mds}: Length = {length:.2f} ft, Wingspan = {wingspan:.2f} ft")

        arcpy.AddMessage("\nMaximum On Ground (MOG) Results:")
        total_aircraft = 0
        total_aircraft_area = 0
        for mds, count, rows, per_row, orientation, aircraft_lcn in results:
            arcpy.AddMessage(f"{mds}: {count} aircraft ({rows} rows, {per_row} per row, {orientation} orientation, LCN: {aircraft_lcn})")
            total_aircraft += count
            # Find the corresponding aircraft dimensions
            aircraft_dims = next((ac for ac in self.aircraft_data if ac[0] == mds), None)
            if aircraft_dims:
                _, length, wingspan, _ = aircraft_dims
                total_aircraft_area += count * length * wingspan

        arcpy.AddMessage(f"\nTotal aircraft: {total_aircraft}")
        arcpy.AddMessage(f"Apron dimensions: {apron_length:.2f} x {apron_width:.2f} ft")
        arcpy.AddMessage(f"Apron LCN: {apron_lcn}")
        arcpy.AddMessage(f"Interior taxiway width: {interior_taxi_width:.2f} ft")
        arcpy.AddMessage(f"Peripheral taxiway width: {peripheral_taxi_width:.2f} ft")
        arcpy.AddMessage(f"Wingtip clearance: {wingtip_clearance:.2f} ft")

        usable_length = apron_length - (2 * peripheral_taxi_width)
        usable_width = apron_width - (2 * peripheral_taxi_width)
        usable_area = usable_length * usable_width
        arcpy.AddMessage(f"Usable apron area: {usable_length:.2f} x {usable_width:.2f} ft")

        # Calculate additional metrics
        total_apron_area = apron_length * apron_width
        taxiway_area = total_apron_area - usable_area
        total_used_area = total_aircraft_area + taxiway_area
        remaining_area = max(0, total_apron_area - total_used_area)
        remaining_space_percentage = max(0, (remaining_area / total_apron_area) * 100)
        space_utilization_efficiency = min(100, (total_used_area / total_apron_area) * 100)

        arcpy.AddMessage(f"\nTotal airfield area: {total_apron_area:.2f} sq ft")
        arcpy.AddMessage(f"Area occupied by aircraft: {total_aircraft_area:.2f} sq ft")
        arcpy.AddMessage(f"Area used for taxiways: {taxiway_area:.2f} sq ft")
        arcpy.AddMessage(f"Total used area: {total_used_area:.2f} sq ft")
        arcpy.AddMessage(f"Remaining area: {remaining_area:.2f} sq ft")
        arcpy.AddMessage(f"Remaining space: {remaining_space_percentage:.2f}% of total airfield area")
        arcpy.AddMessage(f"Space utilization efficiency: {space_utilization_efficiency:.2f}%")

        if total_used_area > total_apron_area:
            arcpy.AddWarning("WARNING: The total used area exceeds the total airfield area.")
            arcpy.AddWarning("This may indicate that some aircraft are overlapping or extending beyond the apron boundaries.")

    def optimize_placement(self, aircraft_data, apron_length, apron_width, apron_lcn,
                           interior_taxi_width, peripheral_taxi_width, wingtip_clearance,
                           apply_lcn, min_quantities):
        population_size = 50
        generations = 20
        mutation_rate = 0.1

        usable_length = apron_length - (2 * peripheral_taxi_width)
        usable_width = apron_width - (2 * peripheral_taxi_width)

        # Process min_quantities
        if min_quantities:
            min_counts = [int(q) for q in min_quantities.split(',')]
        else:
            min_counts = [1] * len(aircraft_data)  # Default to 1 if not specified

        # Initialize population
        population = []
        for _ in range(population_size):
            individual = [random.randint(min_count, 20) for min_count in min_counts]
            population.append(individual)

        best_solution = None
        best_score = float('-inf')

        for generation in range(generations):
            arcpy.AddMessage(f"\nOptimization attempt {generation + 1}/{generations}")
            
            # Evaluate fitness
            fitness_scores = []
            for individual in population:
                score = self.evaluate_fitness(individual, aircraft_data, usable_length, usable_width, apron_lcn,
                                              interior_taxi_width, wingtip_clearance,
                                              apply_lcn, min_counts)
                fitness_scores.append(score)
                if score > best_score:
                    best_score = score
                    best_solution = individual

            # Selection
            new_population = []
            for _ in range(population_size):
                parent1 = self.tournament_selection(population, fitness_scores)
                parent2 = self.tournament_selection(population, fitness_scores)
                child = self.crossover(parent1, parent2)
                child = self.mutate(child, mutation_rate, min_counts)
                new_population.append(child)

            population = new_population

        # Convert best solution to results format
        results = []
        for i, (mds, length, wingspan, aircraft_lcn) in enumerate(aircraft_data):
            count = best_solution[i]
            if count > 0:
                num_rows, per_row, orientation = self.calculate_layout(count, length, wingspan, 
                                                                       usable_length, usable_width,
                                                                       interior_taxi_width, wingtip_clearance)
                results.append((mds, count, num_rows, per_row, orientation, aircraft_lcn))

        return results

    def evaluate_fitness(self, individual, aircraft_data, usable_length, usable_width, apron_lcn,
                         interior_taxi_width, wingtip_clearance,
                         apply_lcn, min_counts):
        total_area = 0
        total_aircraft = 0

        for i, (mds, length, wingspan, aircraft_lcn) in enumerate(aircraft_data):
            count = individual[i]
            if count < min_counts[i]:
                return float('-inf')  # Invalid solution if below minimum count

            if apply_lcn and aircraft_lcn > apron_lcn:
                return float('-inf')  # Invalid solution

            area = count * length * wingspan
            total_area += area
            total_aircraft += count

        if total_area > usable_length * usable_width:
            return float('-inf')  # Invalid solution

        return total_aircraft  # Fitness is the total number of aircraft

    def tournament_selection(self, population, fitness_scores, tournament_size=3):
        selected = random.sample(range(len(population)), tournament_size)
        winner = max(selected, key=lambda i: fitness_scores[i])
        return population[winner]

    def crossover(self, parent1, parent2):
        crossover_point = random.randint(1, len(parent1) - 1)
        child = parent1[:crossover_point] + parent2[crossover_point:]
        return child

    def mutate(self, individual, mutation_rate, min_counts):
        return [max(min_count, gene + random.randint(-1, 1)) if random.random() < mutation_rate else gene 
                for gene, min_count in zip(individual, min_counts)]

    def create_aircraft_layer(self, airfield_layer_object, aircraft_data, mog_result, out_fc, wingtip_clearance, interior_taxi_width):

        try:
            # Validate and create output feature class
            workspace = os.path.dirname(out_fc)
            valid_name = arcpy.ValidateTableName(os.path.basename(out_fc), workspace)
            out_fc = os.path.join(workspace, valid_name)
            sr = arcpy.Describe(airfield_layer_object).spatialReference
            arcpy.management.CreateFeatureclass(workspace, valid_name, "POLYGON", spatial_reference = sr)

            # Add fields for aircraft properties
            arcpy.management.AddField(out_fc, "MDS", "TEXT")
            arcpy.management.AddField(out_fc, "LENGTH", "DOUBLE")
            arcpy.management.AddField(out_fc, "WINGSPAN", "DOUBLE")
            arcpy.management.AddField(out_fc, "FOOTPRINT", "DOUBLE")

            # Check if column with shape angles exists in the airfield layer

            if "SHAPE_ANGLE" not in [columns.name for columns in arcpy.ListFields(airfield_layer_object)]:
                arcpy.management.AddField(airfield_layer_object, "SHAPE_ANGLE", "DOUBLE")
                arcpy.AddMessage("Field 'SHAPE_ANGLE' added to the input table.")
            else:
                arcpy.AddMessage("Field 'SHAPE_ANGLE' already exists in the input table.")

                            
            # Calculate the angle of the selected object of the airfield
            arcpy.cartography.CalculatePolygonMainAngle(airfield_layer_object, "SHAPE_ANGLE", "GEOGRAPHIC")

            # Get the airfield data (location and size)
            with arcpy.da.SearchCursor(airfield_layer_object, ["LENGTH", "WIDTH", "SHAPE_ANGLE", "LATITUDE", "LONGITUDE", "LCN"]) as search_cursor:
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
            arcpy.AddMessage(f"Creatig objects for the following aircrafts: {aircraft_data}")

            x_start = start_lon
            y_start = start_lat

            # Process the selected aircraft
            for row in aircraft_data:
                mds, length, wingspan, aircraft_lcn = row

                length = float(length) if length is not None else 0
                wingspan = float(wingspan) if wingspan is not None else 0

                length_in_degrees = length / 364000  # Approximate conversion from feet to degrees latitude
                wingspan_in_degrees = wingspan / 364000

                result_tuple = next((item for item in mog_result if item[0] == mds), None)
                if result_tuple is None:
                    arcpy.AddMessage(f"No result found for MDS: {mds}")
                    continue
                
                quantity = result_tuple[1]
                max_rows = result_tuple[2]
                max_per_row = result_tuple[3]

                arcpy.AddMessage(f"A total of {quantity} aircraft will be placed over {max_per_row} per row with {max_rows} in each row.")

                with arcpy.da.InsertCursor(out_fc, ["SHAPE@", "MDS", "LENGTH", "WINGSPAN", "FOOTPRINT"]) as insert_cursor:
                    points_placed = 0
                    row_index = 0
                    col_index = 0

                    while points_placed < quantity:
                        x = x_start + (col_index * (wingspan_in_degrees + wingtip_clearance / 364000))
                        y = y_start + (row_index * (length_in_degrees + interior_taxi_width / 364000))

                        # Create the aircraft shape
                        aircraft_shapes = self.create_aircraft_shape(x, y, length_in_degrees, wingspan_in_degrees, apron_angle)

                        # Create the polygon
                        polygon = arcpy.Polygon(arcpy.Array(aircraft_shapes), sr)
                        insert_cursor.insertRow([polygon, mds, length, wingspan, length * wingspan])
                        points_placed += 1

                        arcpy.AddMessage(f"Placed {mds} at ({x_start}, {y_start}) - {points_placed}/{quantity}")

                        # Move to the next column or row
                        col_index += 1
                        if col_index >= max_per_row:
                            col_index = 0
                            row_index += 1
                            
                y_start += interior_taxi_width / 364000

            arcpy.AddMessage(f"Aircraft polygons created in {out_fc}")

        except Exception as e:
            arcpy.AddError(f"An error occurred while creating the polygons: {str(e)}")
            arcpy.AddError(arcpy.GetMessages())
            arcpy.AddError(traceback.format_exc())

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