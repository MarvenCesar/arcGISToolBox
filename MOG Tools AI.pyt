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
                direction="Input")

        return solution

    def run(self, apron_length, apron_width, aircraft_length, aircraft_wingspan, buffer_distance, max_per_row):
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
        self.description = "Calculates the maximum number of aircraft that can be parked on an airfield."
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
        ]

        # Set default values for advanced settings parameters
        for param in params[5:]:
            param.enabled = False  # Disable all advanced settings parameters by default

        return params

    def updateParameters(self, parameters):
        # Auto-populate airfield name based on selected airfield layer
        airfield_layer = parameters[1].valueAsText
        if airfield_layer:
            airfield_name_field = "AFLD_NAME"
            if airfield_name_field in [f.name for f in arcpy.ListFields(airfield_layer)]:
                with arcpy.da.SearchCursor(airfield_layer, [airfield_name_field]) as cursor:
                    for row in cursor:
                        parameters[2].value = row[0]  # Set the airfield name
            else:
                arcpy.AddWarning(f"The field '{airfield_name_field}' does not exist in the selected airfield layer.")

        # Populate the dropdown for Select Aircraft (MDS)
        aircraft_table = parameters[0].valueAsText
        if aircraft_table:
            mds_field = "MDS"  # Replace with the actual field name for MDS in your aircraft table
            mds_values = set()  # Use a set to avoid duplicates
            with arcpy.da.SearchCursor(aircraft_table, [mds_field]) as cursor:
                for row in cursor:
                    mds_values.add(row[0])  # Add MDS values to the set
            parameters[3].filter.list = list(mds_values)  # Set the dropdown list

        # Show/hide advanced settings based on checkbox
        enable_advanced_settings = parameters[4].value
        for param in parameters[5:]:
            param.enabled = enable_advanced_settings

        return

    def calculate_taxiway_widths(self, aircraft_length, aircraft_wingspan):
        # Initialize taxiway widths
        interior_taxi_width = 0
        peripheral_taxi_width = 0

        # Convert wingspan to feet for comparison
        wingspan_ft = aircraft_wingspan  # Assuming wingspan is already in feet

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
            interior_taxi_width = 30 + aircraft_wingspan + 30  # Example calculation for larger aircraft
        else:
            interior_taxi_width = 20 + aircraft_wingspan + 20  # Example calculation for smaller aircraft

        peripheral_taxi_width = 50 + (aircraft_wingspan / 2) + 37.5  # Peripheral Taxi Width
        wingtip_between_parked = 25  # Space between parked aircraft

        # I. Standard apron (Width = apron_width, Length = apron_length)
        # 1) Determine # rows
        available_length = apron_length
        num_rows = math.floor((available_length + interior_taxi_width) / (aircraft_length + interior_taxi_width))
        
        # Print the number of rows
        arcpy.AddMessage(f"Standard Configuration - Available Length: {available_length}, Interior Taxi Width: {interior_taxi_width}, Aircraft Length: {aircraft_length}")
        arcpy.AddMessage(f"Calculated Rows: {num_rows}")

        # 2) Determine # Cols
        available_width = max(0, apron_width - peripheral_taxi_width)
        num_cols = math.floor((available_width + wingtip_between_parked) / (aircraft_wingspan + wingtip_between_parked))
        
        # Print the number of columns
        arcpy.AddMessage(f"Standard Configuration - Available Width: {available_width}, Wingtip Between Parked: {wingtip_between_parked}, Aircraft Wingspan: {aircraft_wingspan}")
        arcpy.AddMessage(f"Calculated Columns: {num_cols}")

        # 3) Parking Available I
        parking_available_I = num_rows * num_cols

        # II) Rotated configuration (Width=Length, Length=Width)
        # 1) Determine # Rows for rotated configuration
        available_length_rotated = apron_width
        num_rows_rotated = math.floor((available_length_rotated + interior_taxi_width) / (aircraft_length + interior_taxi_width))

        # Run genetic algorithm
        ga = GeneticAlgorithm(population_size, mutation_rate, generations)
        best_solution = ga.run(apron_length, apron_width, aircraft_length, aircraft_wingspan, buffer_distance, max_per_row)

        # Output the best solution
        arcpy.AddMessage(f"Best solution: {best_solution}")
        # Print the number of rows for rotated configuration
        arcpy.AddMessage(f"Rotated Configuration - Available Length: {available_length_rotated}, Interior Taxi Width: {interior_taxi_width}, Aircraft Length: {aircraft_length}")
        arcpy.AddMessage(f"Calculated Rotated Rows: {num_rows_rotated}")

        # 2) Determine # Cols for rotated configuration
        available_width_rotated = max(0, apron_length - peripheral_taxi_width)
        num_cols_rotated = math.floor((available_width_rotated + wingtip_between_parked) / (aircraft_wingspan + wingtip_between_parked))
        
        # Print the number of columns for rotated configuration
        arcpy.AddMessage(f"Rotated Configuration - Available Width: {available_width_rotated}, Wingtip Between Parked: {wingtip_between_parked}, Aircraft Wingspan: {aircraft_wingspan}")
        arcpy.AddMessage(f"Calculated Rotated Columns: {num_cols_rotated}")

        # 3) Parking Available II
        parking_available_II = num_rows_rotated * num_cols_rotated

        # III) Final Parking Available
        final_parking_available = max(parking_available_I, parking_available_II)

        # Print final parking available
        arcpy.AddMessage(f"Final Parking Available: {final_parking_available}")

        return final_parking_available

    def execute(self, parameters, messages):
        # Extract parameter values
        aircraft_table = parameters[0].valueAsText
        airfield_layer = parameters[1].valueAsText
        airfield_name = parameters[2].valueAsText
        selected_aircraft = parameters[3].valueAsText.split(';')
        enable_advanced_settings = parameters[4].value

        # Check for manual airfield dimensions
        manual_length = parameters[5].value
        manual_width = parameters[6].value
        manual_interior_taxi_width = parameters[7].value
        manual_peripheral_taxi_width = parameters[8].value

        # Use manual inputs if provided, otherwise use defaults
        if manual_length is not None:
            apron_length = float(manual_length)
        else:
            # Retrieve apron length from airfield data
            apron_length = self.get_airfield_dimension(airfield_layer, airfield_name, "LENGTH")

        if manual_width is not None:
            apron_width = float(manual_width)
        else:
            # Retrieve apron width from airfield data
            apron_width = self.get_airfield_dimension(airfield_layer, airfield_name, "WIDTH")

        # Enable advanced settings if the checkbox is checked
        if enable_advanced_settings:
            if manual_interior_taxi_width is not None:
                interior_taxi_width = float(manual_interior_taxi_width)
            else:
                interior_taxi_width = None  # Will be calculated based on aircraft dimensions

            if manual_peripheral_taxi_width is not None:
                peripheral_taxi_width = float(manual_peripheral_taxi_width)
            else:
                peripheral_taxi_width = None  # Will be calculated based on aircraft dimensions
        else:
            interior_taxi_width = None  # Will be calculated based on aircraft dimensions
            peripheral_taxi_width = None  # Will be calculated based on aircraft dimensions

        # Get aircraft dimensions (assuming you have a method to retrieve these)
        aircraft_data = self.get_aircraft_data(aircraft_table, selected_aircraft)
        for mds, length, wingspan, _ in aircraft_data:
            # Calculate taxiway widths based on aircraft dimensions if not manually provided
            if interior_taxi_width is None or peripheral_taxi_width is None:
                interior_taxi_width, peripheral_taxi_width = self.calculate_taxiway_widths(length, wingspan)

            # Assuming a fixed wingtip clearance for parked aircraft
            wingtip_between_parked = 25  # Adjust as necessary based on your requirements

            # Calculate parking available for each aircraft
            # parking_available = self.calculate_parking_available(apron_length, apron_width, interior_taxi_width, peripheral_taxi_width, length, wingspan, wingtip_between_parked)
            parking_available = self.calculate_parking_available(apron_length, apron_width, length, wingspan)
            arcpy.AddMessage(f"Parking available for {mds}: {parking_available}")

        # Continue with the rest of the execution logic...

    def get_aircraft_data(self, aircraft_table, selected_aircraft):
        aircraft_data = []
        with arcpy.da.SearchCursor(aircraft_table, ["MDS", "LENGTH", "WING_SPAN", "ACFT_LCN"]) as cursor:
            for row in cursor:
                if row[0] in selected_aircraft:
                    aircraft_data.append(row)
        return aircraft_data

    def get_airfield_dimension(self, airfield_layer, airfield_name, dimension_field):
        where_clause = f"AFLD_NAME = '{airfield_name}'"
        with arcpy.da.SearchCursor(airfield_layer, [dimension_field], where_clause) as cursor:
            for row in cursor:
                return row[0]  # Return the requested dimension
        arcpy.AddError(f"Airfield '{airfield_name}' not found or dimension '{dimension_field}' not available.")
        return None  # Return None if not found