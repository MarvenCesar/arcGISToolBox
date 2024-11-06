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
        self.population_size = 50
        self.mutation_rate = 0.01
        self.generations = 100

    def getParameterInfo(self):
        params = [
            arcpy.Parameter(
                displayName="Input Aircraft Table",
                name="aircraft_table",
                datatype="GPTableView",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Airfield Layer",
                name="airfield_layer",
                datatype="GPFeatureLayer",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Select Aircraft (MDS)",
                name="selected_aircraft",
                datatype="GPString",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Enable Advanced Settings",
                name="enable_advanced_settings",
                datatype="GPBoolean",
                parameterType="Optional",
                direction="Input"),
            arcpy.Parameter(
                displayName="Manual Airfield Length (ft)",
                name="manual_length",
                datatype="GPDouble",
                parameterType="Optional",
                direction="Input"),
            arcpy.Parameter(
                displayName="Manual Airfield Width (ft)",
                name="manual_width",
                datatype="GPDouble",
                parameterType="Optional",
                direction="Input"),
            arcpy.Parameter(
                displayName="Manual Interior Taxiway Width (ft)",
                name="manual_interior_taxi_width",
                datatype="GPDouble",
                parameterType="Optional",
                direction="Input"),
            arcpy.Parameter(
                displayName="Manual Peripheral Taxiway Width (ft)",
                name="manual_peripheral_taxi_width",
                datatype="GPDouble",
                parameterType="Optional",
                direction="Input"),
            arcpy.Parameter(
                displayName="Output Aircraft Positions",
                name="out_aircraft_fc",
                datatype="DEFeatureClass",
                parameterType="Required",
                direction="Output"),
            arcpy.Parameter(
                displayName="Output Constraint Polygons",
                name="out_constraint_fc",
                datatype="DEFeatureClass",
                parameterType="Required",
                direction="Output"),
        ]

        for param in params[4:8]:
            param.enabled = False

        return params

    def updateParameters(self, parameters):
        if parameters[0].altered or not parameters[2].altered:
            aircraft_table = parameters[0].valueAsText
            if aircraft_table:
                try:
                    aircraft_names = [row[0] for row in arcpy.da.SearchCursor(aircraft_table, "MDS") if row[0]]
                    parameters[2].filter.list = aircraft_names
                except Exception as e:
                    arcpy.AddError(f"Error loading aircraft names: {e}")

        enable_advanced_settings = parameters[3].value
        for param in parameters[4:8]:
            param.enabled = enable_advanced_settings

        return

    def initialize_population(self, apron_length, apron_width, aircraft_length, aircraft_wingspan, max_per_row):
        if apron_length <= 0 or apron_width <= 0 or aircraft_length <= 0 or aircraft_wingspan <= 0 or max_per_row <= 0:
            raise ValueError("All dimensions and max_per_row must be greater than zero.")
        
        population = []
        for _ in range(self.population_size):
            solution = []
            for _ in range(max_per_row):
                x = random.uniform(0, apron_width - aircraft_wingspan)
                y = random.uniform(0, apron_length - aircraft_length)
                solution.append((x, y))
            population.append(solution)
        return population

    def evaluate_fitness(self, solution, apron_length, apron_width, aircraft_length, aircraft_wingspan, buffer_distance):
        fitness = 0
        for (x, y) in solution:
            if 0 <= x <= apron_width - aircraft_wingspan and 0 <= y <= apron_length - aircraft_length:
                fitness += 1
        return fitness

    def select_parents(self, population, fitnesses):
        selected = random.choices(population, weights=fitnesses, k=2)
        return selected

    def crossover(self, parent1, parent2):
        crossover_point = random.randint(0, len(parent1) - 1)
        child1 = parent1[:crossover_point] + parent2[crossover_point:]
        child2 = parent2[:crossover_point] + parent1[crossover_point:]
        return child1, child2

    def mutate(self, solution, apron_length, apron_width, aircraft_length, aircraft_wingspan):
        for i in range(len(solution)):
            if random.random() < self.mutation_rate:
                x = random.uniform(0, apron_width - aircraft_wingspan)
                y = random.uniform(0, apron_length - aircraft_length)
                solution[i] = (x, y)
        return solution

    def run_genetic_algorithm(self, apron_length, apron_width, aircraft_length, aircraft_wingspan, buffer_distance, max_per_row):
        population = self.initialize_population(apron_length, apron_width, aircraft_length, aircraft_wingspan, max_per_row)
        for generation in range(self.generations):
            fitnesses = [self.evaluate_fitness(solution, apron_length, apron_width, aircraft_length, aircraft_wingspan, buffer_distance) for solution in population]
            new_population = []
            for _ in range(self.population_size // 2):
                parent1, parent2 = self.select_parents(population, fitnesses)
                child1, child2 = self.crossover(parent1, parent2)
                child1 = self.mutate(child1, apron_length, apron_width, aircraft_length, aircraft_wingspan)
                child2 = self.mutate(child2, apron_length, apron_width, aircraft_length, aircraft_wingspan)
                new_population.extend([child1, child2])
            population = new_population
        best_solution = max(population, key=lambda sol: self.evaluate_fitness(sol, apron_length, apron_width, aircraft_length, aircraft_wingspan, buffer_distance))
        return best_solution

    def execute(self, parameters, messages):
        aircraft_table = parameters[0].valueAsText
        airfield_layer = parameters[1].valueAsText
        selected_aircraft = parameters[2].valueAsText
        enable_advanced_settings = parameters[3].value
        manual_length = parameters[4].value
        manual_width = parameters[5].value
        manual_interior_taxi_width = parameters[6].value
        manual_peripheral_taxi_width = parameters[7].value
        out_aircraft_fc = parameters[8].valueAsText
        out_constraint_fc = parameters[9].valueAsText

        try:
            airfield_sr = arcpy.Describe(airfield_layer).spatialReference

            if airfield_sr.type != 'Projected':
                arcpy.AddError("Airfield layer must be in a projected coordinate system with linear units (e.g., meters or feet).")
                return

            sr = airfield_sr

            # Determine unit factor based on the spatial reference units
            if sr.linearUnitName == "Meter":
                unit_factor = 0.3048  # Convert feet to meters
            elif sr.linearUnitName == "Foot":
                unit_factor = 1.0  # No conversion needed
            else:
                arcpy.AddError("Unsupported spatial reference units.")
                return

            with arcpy.da.SearchCursor(airfield_layer, ["SHAPE@", "LENGTH", "WIDTH"], spatial_reference=sr) as cursor:
                for row in cursor:
                    airfield_shape, default_length, default_width = row
                    break
                else:
                    arcpy.AddError(f"Airfield '{airfield_layer}' not found.")
                    return

            angle = self.calculate_airfield_orientation(airfield_layer, sr)
            arcpy.AddMessage(f"Calculated airfield orientation: {angle} degrees")

            if enable_advanced_settings and manual_length is not None and manual_width is not None:
                apron_length = float(manual_length)
                apron_width = float(manual_width)
                arcpy.AddMessage(f"Using manual airfield dimensions: Length={apron_length} ft, Width={apron_width} ft")
            else:
                apron_length = float(default_length)
                apron_width = float(default_width)
                arcpy.AddMessage(f"Retrieved airfield dimensions: Length={apron_length} ft, Width={apron_width} ft")

            aircraft_data = self.get_aircraft_data(aircraft_table, selected_aircraft)
            if not aircraft_data:
                arcpy.AddError(f"Aircraft '{selected_aircraft}' not found.")
                return
            mds, aircraft_length, aircraft_wingspan, _ = aircraft_data[0]
            arcpy.AddMessage(f"Selected Aircraft: {mds}, Length={aircraft_length} ft, Wingspan={aircraft_wingspan} ft")

            if enable_advanced_settings and manual_interior_taxi_width is not None and manual_peripheral_taxi_width is not None:
                interior_taxi_width = float(manual_interior_taxi_width)
                peripheral_taxi_width = float(manual_peripheral_taxi_width)
                arcpy.AddMessage(f"Using manual taxiway widths: Interior Taxi Width={interior_taxi_width} ft, Peripheral Taxi Width={peripheral_taxi_width} ft")
            else:
                interior_taxi_width, peripheral_taxi_width = self.calculate_taxiway_widths(aircraft_length, aircraft_wingspan)
                arcpy.AddMessage(f"Calculated taxiway widths: Interior Taxi Width={interior_taxi_width} ft, Peripheral Taxi Width={peripheral_taxi_width} ft")

            wingtip_between_parked = 25
            arcpy.AddMessage(f"Using wingtip clearance between parked aircraft: {wingtip_between_parked} ft")

            max_per_row = math.floor(apron_width / (aircraft_wingspan + wingtip_between_parked))
            best_solution = self.run_genetic_algorithm(apron_length, apron_width, aircraft_length, aircraft_wingspan, wingtip_between_parked, max_per_row)
            arcpy.AddMessage(f"Best solution found: {best_solution}")

            for fc in [out_aircraft_fc, out_constraint_fc]:
                if arcpy.Exists(fc):
                    arcpy.Delete_management(fc)
                    arcpy.AddMessage(f"Existing feature class '{fc}' deleted.")

            arcpy.CreateFeatureclass_management(os.path.dirname(out_aircraft_fc), os.path.basename(out_aircraft_fc), "POLYGON", spatial_reference=sr)
            arcpy.AddField_management(out_aircraft_fc, "MDS", "TEXT")
            arcpy.AddField_management(out_aircraft_fc, "AircraftID", "LONG")
            arcpy.AddField_management(out_aircraft_fc, "Rotation", "DOUBLE")
            arcpy.AddField_management(out_aircraft_fc, "Type", "TEXT")

            arcpy.CreateFeatureclass_management(os.path.dirname(out_constraint_fc), os.path.basename(out_constraint_fc), "POLYGON", spatial_reference=sr)
            arcpy.AddField_management(out_constraint_fc, "AircraftID", "LONG")
            arcpy.AddField_management(out_constraint_fc, "Type", "TEXT")

            start_point = airfield_shape.centroid
            arcpy.AddMessage(f"Using airfield centroid as starting point: ({start_point.X}, {start_point.Y})")

            angle_rad = math.radians(angle)
            cos_angle = math.cos(angle_rad)
            sin_angle = math.sin(angle_rad)

            aircraft_id = 1
            aircraft_records = []

            with arcpy.da.InsertCursor(out_aircraft_fc, ["SHAPE@", "MDS", "AircraftID", "Rotation", "Type"]) as aircraft_cursor:
                for (x, y) in best_solution:
                    x_rotated = x * cos_angle - y * sin_angle
                    y_rotated = x * sin_angle + y * cos_angle
                    x_global = start_point.X + x_rotated
                    y_global = start_point.Y + y_rotated

                    half_length = (aircraft_length * unit_factor) / 2
                    half_wingspan = (aircraft_wingspan * unit_factor) / 2
                    corners = [
                        arcpy.Point(x_global - half_wingspan, y_global + half_length),
                        arcpy.Point(x_global + half_wingspan, y_global + half_length),
                        arcpy.Point(x_global + half_wingspan, y_global - half_length),
                        arcpy.Point(x_global - half_wingspan, y_global - half_length),
                        arcpy.Point(x_global - half_wingspan, y_global + half_length)
                    ]

                    polygon_array = arcpy.Array()
                    for corner in corners:
                        dx = corner.X - x_global
                        dy = corner.Y - y_global
                        rotated_x = x_global + (dx * cos_angle - dy * sin_angle)
                        rotated_y = y_global + (dx * sin_angle + dy * cos_angle)
                        polygon_array.add(arcpy.Point(rotated_x, rotated_y))

                    aircraft_polygon = arcpy.Polygon(polygon_array, sr)

                    if airfield_shape.contains(aircraft_polygon.centroid):
                        aircraft_cursor.insertRow([aircraft_polygon, mds, aircraft_id, angle % 360, "Aircraft"])
                        arcpy.AddMessage(f"Placed {mds} at ({x_global}, {y_global}) - ID: {aircraft_id}")
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

            self.create_constraint_polygons(out_constraint_fc, sr, aircraft_records)

        except Exception as e:
            arcpy.AddError(f"An error occurred during the MOG calculation: {str(e)}")
            arcpy.AddError(traceback.format_exc())

    def create_constraint_polygons(self, out_constraint_fc, sr, aircraft_records):
        try:
            with arcpy.da.InsertCursor(out_constraint_fc, ["SHAPE@", "AircraftID", "Type"]) as constraint_cursor:
                for record in aircraft_records:
                    aircraft_id = record["AircraftID"]
                    centroid = record["Centroid"]
                    length = record["Length"]
                    wingspan = record["Wingspan"]
                    rotation = record["Rotation"]

                    clearance_factor = 1.2
                    constraint_length = length * clearance_factor
                    constraint_wingspan = wingspan * clearance_factor

                    half_length = constraint_length / 2
                    half_wingspan = constraint_wingspan / 2
                    corners = [
                        arcpy.Point(centroid.X - half_wingspan, centroid.Y + half_length),
                        arcpy.Point(centroid.X + half_wingspan, centroid.Y + half_length),
                        arcpy.Point(centroid.X + half_wingspan, centroid.Y - half_length),
                        arcpy.Point(centroid.X - half_wingspan, centroid.Y - half_length),
                        arcpy.Point(centroid.X - half_wingspan, centroid.Y + half_length)
                    ]

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

                    constraint_polygon = arcpy.Polygon(arcpy.Array(rotated_corners), sr)

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

    def calculate_taxiway_widths(self, aircraft_length, aircraft_wingspan):
        interior_taxi_width = aircraft_length * 1.5
        peripheral_taxi_width = aircraft_wingspan * 1.5
        return interior_taxi_width, peripheral_taxi_width

    def calculate_airfield_orientation(self, airfield_layer, sr):
        with arcpy.da.SearchCursor(airfield_layer, ["SHAPE@"], spatial_reference=sr) as cursor:
            for row in cursor:
                airfield_shape = row[0]
                break
        angle = airfield_shape.extent.angle
        return angle

    def calculate_parking_available(self, apron_length, apron_width, aircraft_length, aircraft_wingspan):
        max_per_row = math.floor(apron_width / (aircraft_wingspan + 25))
        max_rows = math.floor(apron_length / (aircraft_length + 25))
        return max_per_row, max_rows