import arcpy
import pandas as pd
import os
import math
import traceback
import random

class GeneticAlgorithm(object):
    def __init__(self, population_size, mutation_rate, generations):
        self.population_size = population_size
        self.mutation_rate = mutation_rate
        self.generations = generations

    def initialize_population(self, apron_length, apron_width, aircraft_length, aircraft_wingspan, max_per_row):
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
        if all(f == 0 for f in fitnesses):
            raise ValueError("All fitness values are zero. Cannot proceed with selection.")
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

    def run(self, apron_length, apron_width, aircraft_length, aircraft_wingspan, buffer_distance, max_per_row):
        population = self.initialize_population(apron_length, apron_width, aircraft_length, aircraft_wingspan, max_per_row)
        for generation in range(self.generations):
            fitnesses = [self.evaluate_fitness(solution, apron_length, apron_width, aircraft_length, aircraft_wingspan, buffer_distance) for solution in population]
            if all(f == 0 for f in fitnesses):
                raise ValueError("All fitness values are zero. Cannot proceed with selection.")
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

class Toolbox(object):
    def __init__(self):
        self.label = "Aircraft MOG Optimization"
        self.alias = "AircraftMOG"
        self.tools = [ImportAircraftData, CalculateAircraftFootprint, CalculateMaximumOnGround]


# Class for importing aircraft data and saving it to GDB if required
class ImportAircraftData(object):
    def __init__(self):
        self.label = "Import Aircraft Data"
        self.description = "Import aircraft specifications from a CSV file and optionally save as a GDB table"

    def getParameterInfo(self):
        params = [
            arcpy.Parameter(
                displayName="Input CSV File", 
                name="in_csv", 
                datatype="DEFile", 
                parameterType="Required", 
                direction="Input"),
            arcpy.Parameter(
                displayName="Output Table (CSV or GDB)", 
                name="out_table", 
                datatype="DETable", 
                parameterType="Required", 
                direction="Output"),
            arcpy.Parameter(
                displayName="Save as Geodatabase Table", 
                name="save_as_gdb", 
                datatype="GPBoolean", 
                parameterType="Optional", 
                direction="Input")
        ]
        return params

    def execute(self, parameters, messages):
        in_csv = parameters[0].valueAsText
        out_table = parameters[1].valueAsText
        save_as_gdb = parameters[2].value   # Boolean for saving as GDB
        
        try:
            # Read CSV data
            if not os.path.exists(in_csv):
                arcpy.AddError(f"Input CSV file does not exist: {in_csv}")
                return

            df = pd.read_csv(in_csv)
            df = df.dropna(subset=['MDS'])  # Drop rows with missing MDS

            numeric_columns = ['LENGTH', 'WING_SPAN', 'HEIGHT', 'WING_HEIGHT', 'TURNING_RADIUS', 'MIN_RWY_LENGTH', 'ACFT_LCN']
            for col in numeric_columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')

            # Convert the CSV data to a temporary CSV file
            temp_csv = os.path.join(arcpy.env.scratchFolder, "temp.csv")
            df.to_csv(temp_csv, index=False)

            # Check if saving as GDB
            if save_as_gdb:
                # Save the table as a GDB table
                gdb_path = os.path.dirname(out_table)
                table_name = os.path.basename(out_table)
                arcpy.conversion.TableToTable(temp_csv, gdb_path, table_name)
                arcpy.AddMessage(f"Successfully saved as GDB table: {out_table}")
            else:
                # Save as a standard table (CSV)
                arcpy.conversion.TableToTable(temp_csv, os.path.dirname(out_table), os.path.basename(out_table))
                arcpy.AddMessage(f"Successfully imported aircraft records as CSV table.")
            
            # Clean up
            os.remove(temp_csv)

        except Exception as e:
            arcpy.AddError(f"An error occurred during import: {str(e)}")
            arcpy.AddError(arcpy.GetMessages())



class CalculateAircraftFootprint(object):
    def __init__(self):
        self.label = "Create Aircraft Symbol Layer (Aircraft-Shaped Polygons)"
        self.description = "Create polygon footprints resembling aircraft shapes for a selected aircraft at a specified airfield location."

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
                displayName="Clearance Distance (in feet)",
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

    def updateParameters(self, parameters):
        # Ensure parameters list is long enough before accessing
        if len(parameters) < 7:
            arcpy.AddError("Not enough parameters provided.")
            return

        # Populate the Aircraft Name dropdown based on the Input Aircraft Table
        if parameters[0].altered and not parameters[3].altered:
            aircraft_table = parameters[0].valueAsText
            if aircraft_table:
                try:
                    aircraft_names = set()
                    with arcpy.da.SearchCursor(aircraft_table, ["MDS"]) as cursor:
                        for row in cursor:
                            aircraft_names.add(row[0])
                    parameters[3].filter.list = sorted(aircraft_names)
                except Exception as e:
                    arcpy.AddError(f"Error loading aircraft names: {e}")
        
        # Populate the Airfield Name dropdown based on the Airfield Layer
        if parameters[1].altered and not parameters[5].altered:
            airfield_layer = parameters[1].valueAsText
            if airfield_layer:
                try:
                    airfield_names = set()
                    with arcpy.da.SearchCursor(airfield_layer, ["AFLD_NAME"]) as cursor:
                        for row in cursor:
                            airfield_names.add(row[0])
                    parameters[5].filter.list = sorted(airfield_names)
                except Exception as e:
                    arcpy.AddError(f"Error loading airfield names: {e}")
        
        # Populate the Airfield ObjectID dropdown based on the selected airfield name
        if parameters[5].altered and not parameters[6].altered:
            airfield_layer = parameters[1].valueAsText
            selected_airfield_name = parameters[5].valueAsText
            if airfield_layer and selected_airfield_name:
                try:
                    object_ids = set()
                    where_clause = f"AFLD_NAME = '{selected_airfield_name}'"
                    with arcpy.da.SearchCursor(airfield_layer, ["OBJECTID"], where_clause) as cursor:
                        for row in cursor:
                            object_ids.add(str(row[0]))  # Convert ObjectID to string for dropdown
                    parameters[6].filter.list = sorted(object_ids)
                except Exception as e:
                    arcpy.AddError(f"Error loading ObjectIDs: {e}")
        
        return

    def create_aircraft_shape(self, x_start, y_start, length, wingspan):
        # Define proportions for the aircraft shape
        fuselage_width = length * 0.1
        nose_length = length * 0.2
        tail_length = length * 0.15
        wing_sweep = length * 0.1
        tail_sweep = length * 0.05

        corners = [
            arcpy.Point(x_start, y_start + length / 2),  # Nose tip
            arcpy.Point(x_start - fuselage_width / 2, y_start + length / 2 - nose_length),  # Nose left
            arcpy.Point(x_start - wingspan / 2, y_start + wing_sweep),  # Left wingtip front
            arcpy.Point(x_start - wingspan / 2, y_start),  # Left wingtip middle
            arcpy.Point(x_start - wingspan / 2, y_start - wing_sweep),  # Left wingtip rear
            arcpy.Point(x_start - fuselage_width / 2, y_start - length / 2 + tail_length),  # Fuselage left before tail
            arcpy.Point(x_start - wingspan / 4, y_start - length / 2),  # Left tail tip
            arcpy.Point(x_start, y_start - length / 2 - tail_sweep),  # Tail bottom tip
            arcpy.Point(x_start + wingspan / 4, y_start - length / 2),  # Right tail tip
            arcpy.Point(x_start + fuselage_width / 2, y_start - length / 2 + tail_length),  # Fuselage right before tail
            arcpy.Point(x_start + wingspan / 2, y_start - wing_sweep),  # Right wingtip rear
            arcpy.Point(x_start + wingspan / 2, y_start),  # Right wingtip middle
            arcpy.Point(x_start + wingspan / 2, y_start + wing_sweep),  # Right wingtip front
            arcpy.Point(x_start + fuselage_width / 2, y_start + length / 2 - nose_length),  # Nose right
            arcpy.Point(x_start, y_start + length / 2)  # Back to nose tip
        ]
        return corners

    def execute(self, parameters, messages):
        # Retrieve parameters
        in_table = parameters[0].valueAsText
        airfield_layer = parameters[1].valueAsText
        out_fc = parameters[2].valueAsText
        selected_aircraft = parameters[3].valueAsText
        quantity_of_aircraft = int(parameters[4].valueAsText)
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
            arcpy.CreateFeatureclass_management(workspace, valid_name, "POLYGON", spatial_reference=sr)

            # Add fields for aircraft properties
            arcpy.AddField_management(out_fc, "MDS", "TEXT")
            arcpy.AddField_management(out_fc, "LENGTH", "DOUBLE")
            arcpy.AddField_management(out_fc, "WINGSPAN", "DOUBLE")
            arcpy.AddField_management(out_fc, "Aircraft_Footprint", "DOUBLE")

            # Get the airfield data (location and size)
            airfield_where_clause = f"AFLD_NAME = '{afld_name}' AND OBJECTID = {object_id}"
            with arcpy.da.SearchCursor(airfield_layer, ["SHAPE@", "LENGTH", "WIDTH", "LATITUDE", "LONGITUDE", "LCN"], airfield_where_clause) as cursor:
                for row in cursor:
                    airfield_shape, apron_length, apron_width, start_lat, start_lon, apron_lcn = row
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

            # Process the selected aircraft
            with arcpy.da.SearchCursor(in_table, ["MDS", "LENGTH", "WING_SPAN", "ACFT_LCN"]) as search_cursor:
                for row in search_cursor:
                    mds, length, wingspan, aircraft_lcn = row
                    if mds == selected_aircraft:
                        length = float(length) if length is not None else 0
                        wingspan = float(wingspan) if wingspan is not None else 0
                        aircraft_lcn = float(aircraft_lcn) if aircraft_lcn is not None else 0

                        if length == 0 or wingspan == 0 or aircraft_lcn == 0:
                            arcpy.AddError(f"Invalid data for {mds}: LENGTH={length}, WING_SPAN={wingspan}, ACFT_LCN={aircraft_lcn}")
                            return

                        arcpy.AddMessage(f"Processing aircraft {mds}: Length={length}, Wingspan={wingspan}, LCN={aircraft_lcn}")

                        if aircraft_lcn > apron_lcn:
                            arcpy.AddWarning(f"Aircraft {mds} LCN ({aircraft_lcn}) exceeds apron LCN ({apron_lcn}). Placement may not be suitable.")

                        length_in_degrees = length / 364000  # Approximate conversion from feet to degrees latitude
                        wingspan_in_degrees = wingspan / 364000

                        with arcpy.da.InsertCursor(out_fc, ["SHAPE@", "MDS", "LENGTH", "WINGSPAN", "Aircraft_Footprint"]) as insert_cursor:
                            points_placed = 0
                            row_index = 0
                            col_index = 0

                            while points_placed < quantity_of_aircraft:
                                x_start = start_lon + (col_index * (wingspan_in_degrees + buffer_distance / 364000))
                                y_start = start_lat + (row_index * (length_in_degrees + buffer_distance / 364000))

                                # Create the aircraft shape
                                corners = self.create_aircraft_shape(x_start, y_start, length_in_degrees, wingspan_in_degrees)

                                # Create the polygon
                                polygon = arcpy.Polygon(arcpy.Array(corners), sr)
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
            


class CalculateMaximumOnGround(object):
    def __init__(self):
        self.label = "Calculate Maximum On Ground"
        self.description = "Calculates the maximum number of aircraft that can be parked on an airfield using a genetic algorithm."
        self.canRunInBackground = False

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
                displayName="Airfield Name",
                name="airfield_name",
                datatype="GPString",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Select Aircraft (MDS)",
                name="selected_aircraft",
                datatype="GPString",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Population Size",
                name="population_size",
                datatype="GPLong",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Mutation Rate",
                name="mutation_rate",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Generations",
                name="generations",
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
                direction="Input"),
            arcpy.Parameter(
                displayName="Output Feature Class",
                name="out_fc",
                datatype="DEFeatureClass",
                parameterType="Required",
                direction="Output")
        ]
        return params

    def updateParameters(self, parameters):
        # Populate the Airfield Name drop-down
        if parameters[1].altered and not parameters[2].altered:
            airfield_layer = parameters[1].valueAsText
            if airfield_layer:
                airfield_names = set()
                with arcpy.da.SearchCursor(airfield_layer, ["AFLD_NAME"]) as cursor:
                    for row in cursor:
                        airfield_names.add(row[0])
                parameters[2].filter.list = sorted(airfield_names)

        # Populate the Aircraft Name drop-down
        if parameters[0].altered and not parameters[3].altered:
            aircraft_table = parameters[0].valueAsText
            if aircraft_table:
                aircraft_names = set()
                with arcpy.da.SearchCursor(aircraft_table, ["MDS"]) as cursor:
                    for row in cursor:
                        aircraft_names.add(row[0])
                parameters[3].filter.list = sorted(aircraft_names)
        return

    def execute(self, parameters, messages):
        aircraft_table = parameters[0].valueAsText
        airfield_layer = parameters[1].valueAsText
        airfield_name = parameters[2].valueAsText
        selected_aircraft = parameters[3].valueAsText
        population_size = int(parameters[4].valueAsText)
        mutation_rate = float(parameters[5].valueAsText)
        generations = int(parameters[6].valueAsText)
        buffer_distance = float(parameters[7].valueAsText)
        max_per_row = int(parameters[8].valueAsText)
        out_fc = parameters[9].valueAsText

        try:
            # Get airfield dimensions
            apron_length = self.get_airfield_dimension(airfield_layer, airfield_name, "LENGTH")
            apron_width = self.get_airfield_dimension(airfield_layer, airfield_name, "WIDTH")

            # Get aircraft dimensions
            aircraft_data = self.get_aircraft_data(aircraft_table, selected_aircraft)
            if not aircraft_data:
                arcpy.AddError(f"Aircraft {selected_aircraft} not found in the input table.")
                return

            aircraft_length, aircraft_wingspan = aircraft_data[0][1], aircraft_data[0][2]

            # Run genetic algorithm
            ga = GeneticAlgorithm(population_size, mutation_rate, generations)
            best_solution = ga.run(apron_length, apron_width, aircraft_length, aircraft_wingspan, buffer_distance, max_per_row)

            # Create output feature class
            workspace = os.path.dirname(out_fc)
            valid_name = arcpy.ValidateTableName(os.path.basename(out_fc), workspace)
            out_fc = os.path.join(workspace, valid_name)
            sr = arcpy.Describe(airfield_layer).spatialReference
            arcpy.CreateFeatureclass_management(workspace, valid_name, "POINT", spatial_reference=sr)

            # Add fields for aircraft properties
            arcpy.AddField_management(out_fc, "MDS", "TEXT")
            arcpy.AddField_management(out_fc, "LENGTH", "DOUBLE")
            arcpy.AddField_management(out_fc, "WINGSPAN", "DOUBLE")

            # Insert aircraft positions into the feature class
            with arcpy.da.InsertCursor(out_fc, ["SHAPE@", "MDS", "LENGTH", "WINGSPAN"]) as cursor:
                for (x, y) in best_solution:
                    point = arcpy.Point(x, y)
                    cursor.insertRow([point, selected_aircraft, aircraft_length, aircraft_wingspan])

            arcpy.AddMessage(f"Best solution visualized in {out_fc}")

        except Exception as e:
            arcpy.AddError(f"An error occurred: {str(e)}")
            arcpy.AddError(arcpy.GetMessages())
            arcpy.AddError(traceback.format_exc())

    def get_aircraft_data(self, aircraft_table, selected_aircraft):
        with arcpy.da.SearchCursor(aircraft_table, ["MDS", "LENGTH", "WING_SPAN"]) as cursor:
            for row in cursor:
                if row[0] == selected_aircraft:
                    return [row]
        arcpy.AddError(f"Aircraft '{selected_aircraft}' not found in the input table.")
        return None

    def get_airfield_dimension(self, airfield_layer, airfield_name, dimension_field):
        where_clause = f"AFLD_NAME = '{airfield_name}'"
        with arcpy.da.SearchCursor(airfield_layer, [dimension_field], where_clause) as cursor:
            for row in cursor:
                return row[0]
        arcpy.AddError(f"Airfield '{airfield_name}' not found or dimension '{dimension_field}' not available.")
        return None