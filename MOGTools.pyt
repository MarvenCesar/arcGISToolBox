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
        save_as_gdb = parameters[2].value  # Boolean for saving as GDB
        
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
                displayName="Interior Taxiway Width (ft)",
                name="interior_taxi_width",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Peripheral Taxiway Width (ft)",
                name="peripheral_taxi_width",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Aircraft Wingtip Clearance (ft)",
                name="wingtip_clearance",
                datatype="GPDouble",
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Select Aircraft (MDS)",
                name="selected_aircraft",
                datatype="GPString",
                parameterType="Required",
                direction="Input",
                multiValue=True),
            arcpy.Parameter(
                displayName="Aircraft Quantities (Comma Separated)",
                name="aircraft_quantities",
                datatype="GPString",
                parameterType="Optional",
                direction="Input"),
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
        ]
        params[1].filter.list = ["Polygon"]
        params[9].value = False  # Default to not optimizing
        params[11].value = self.optimization_attempts  # Use class attribute for default
        return params

    def updateParameters(self, parameters):
        if parameters[1].altered and not parameters[2].altered:
            airfield_layer = parameters[1].valueAsText
            if airfield_layer:
                airfield_names = [row[0] for row in arcpy.da.SearchCursor(airfield_layer, "AFLD_NAME")]
                parameters[2].filter.list = sorted(set(airfield_names))
        
        if parameters[0].altered and not parameters[6].altered:
            aircraft_table = parameters[0].valueAsText
            if aircraft_table:
                try:
                    aircraft_names = set(row[0] for row in arcpy.da.SearchCursor(aircraft_table, ["MDS"]))
                    parameters[6].filter.list = sorted(aircraft_names)
                except Exception as e:
                    arcpy.AddError(f"Error loading aircraft names: {str(e)}")
        
        optimize_quantities_param = parameters[9]
        min_aircraft_quantities_param = parameters[10]
        optimization_attempts_param = parameters[11]

        optimize_quantities_param.value = False if optimize_quantities_param.value is None else optimize_quantities_param.value
        min_aircraft_quantities_param.enabled = optimize_quantities_param.value
        optimization_attempts_param.enabled = optimize_quantities_param.value

    def execute(self, parameters, messages):
        # Extract parameter values
        aircraft_table = parameters[0].valueAsText
        airfield_layer = parameters[1].valueAsText
        airfield_name = parameters[2].valueAsText
        interior_taxi_width = float(parameters[3].value)
        peripheral_taxi_width = float(parameters[4].value)
        wingtip_clearance = float(parameters[5].value)
        selected_aircraft = parameters[6].valueAsText.split(';')
        aircraft_quantities = parameters[7].valueAsText
        apply_lcn = parameters[8].value
        optimize_quantities = parameters[9].value
        min_quantities = parameters[10].valueAsText
        self.optimization_attempts = int(parameters[11].value) if optimize_quantities else 1

        # Get airfield data
        apron_length, apron_width = self.get_airfield_data(airfield_layer, airfield_name)

        # Get aircraft data
        self.aircraft_data = self.get_aircraft_data(aircraft_table, selected_aircraft, aircraft_quantities, min_quantities)

        # Calculate MOG
        result = self.calculate_mog(self.aircraft_data, apron_length, apron_width, 
                                    interior_taxi_width, peripheral_taxi_width, 
                                    wingtip_clearance, optimize_quantities)

        # Display results
        self.display_results(result, apron_length, apron_width, 
                             interior_taxi_width, peripheral_taxi_width, wingtip_clearance)

    def get_airfield_data(self, airfield_layer, airfield_name):
        with arcpy.da.SearchCursor(airfield_layer, ["LENGTH", "WIDTH"], f"AFLD_NAME = '{airfield_name}'") as cursor:
            for row in cursor:
                return float(row[0]), float(row[1])
        arcpy.AddError(f"Airfield '{airfield_name}' not found.")
        return None, None

    def get_aircraft_data(self, aircraft_table, selected_aircraft, quantities, min_quantities):
        aircraft_data = []
        mds_field, length_field, wingspan_field = "MDS", "LENGTH", "WING_SPAN"

        with arcpy.da.SearchCursor(aircraft_table, [mds_field, length_field, wingspan_field]) as cursor:
            for row in cursor:
                if row[0] in selected_aircraft:
                    mds, length, wingspan = row
                    qty = int(quantities.split(',')[selected_aircraft.index(mds)]) if quantities else 0
                    min_qty = int(min_quantities.split(',')[selected_aircraft.index(mds)]) if min_quantities else 0
                    aircraft_data.append((mds, length, wingspan, qty, min_qty))
        
        arcpy.AddMessage("Selected Aircraft Dimensions:")
        for ac in aircraft_data:
            arcpy.AddMessage(f"{ac[0]}: Length = {ac[1]:.2f} ft, Wingspan = {ac[2]:.2f} ft")
        
        return aircraft_data

    def calculate_mog(self, aircraft_data, apron_length, apron_width, 
                      interior_taxi_width, peripheral_taxi_width, 
                      wingtip_clearance, optimize_quantities):
        usable_length = apron_length - 2 * peripheral_taxi_width
        usable_width = apron_width - 2 * peripheral_taxi_width

        arcpy.AddMessage(f"\nTotal airfield dimensions: {apron_length:.2f} x {apron_width:.2f} ft")
        arcpy.AddMessage(f"Usable parking area: {usable_length:.2f} x {usable_width:.2f} ft")
        arcpy.AddMessage(f"Interior taxiway width: {interior_taxi_width:.2f} ft")
        arcpy.AddMessage(f"Peripheral taxiway width: {peripheral_taxi_width:.2f} ft")
        arcpy.AddMessage(f"Wingtip clearance: {wingtip_clearance:.2f} ft")

        if optimize_quantities:
            arcpy.AddMessage("\nUsing optimized placement algorithm")
            return self.optimize_placement(aircraft_data, usable_length, usable_width, 
                                           interior_taxi_width, wingtip_clearance)
        else:
            arcpy.AddMessage("\nUsing fixed placement algorithm")
            return self.fixed_placement(aircraft_data, usable_length, usable_width, 
                                        interior_taxi_width, wingtip_clearance)

    def optimize_placement(self, aircraft_data, usable_length, usable_width, 
                           interior_taxi_width, wingtip_clearance):
        arcpy.AddMessage(f"\nAvailable parking area: {usable_length:.2f} x {usable_width:.2f} ft")
        best_solution = None
        best_score = float('-inf')

        for attempt in range(self.optimization_attempts):
            arcpy.AddMessage(f"\nOptimization attempt {attempt + 1}/{self.optimization_attempts}")
            
            solution = []
            remaining_space = [(0, 0, usable_length, usable_width)]
            aircraft_counts = {ac[0]: 0 for ac in aircraft_data}

            while True:
                placed_aircraft = self.place_next_aircraft(aircraft_data, remaining_space, 
                                                           interior_taxi_width, wingtip_clearance, 
                                                           aircraft_counts)
                if not placed_aircraft:
                    break
                solution.append(placed_aircraft)
                aircraft_counts[placed_aircraft[0]] += 1

            score = self.calculate_score(aircraft_counts, aircraft_data)
            if score > best_score:
                best_solution = aircraft_counts
                best_score = score

            arcpy.AddMessage("After placement:")
            for mds, count in aircraft_counts.items():
                min_qty = next(ac[4] for ac in aircraft_data if ac[0] == mds)
                arcpy.AddMessage(f"Placed {count} {mds} aircraft (minimum required: {min_qty})")

        arcpy.AddMessage(f"\nBest optimization result: {best_solution}")
        return list(best_solution.items())

    def place_next_aircraft(self, aircraft_data, remaining_space, 
                            interior_taxi_width, wingtip_clearance, current_counts):
        # Prioritize aircraft types that haven't met their minimum quantities
        unmet_aircraft = [ac for ac in aircraft_data if current_counts[ac[0]] < ac[4]]
        if unmet_aircraft:
            aircraft_to_place = random.choice(unmet_aircraft)
        else:
            aircraft_to_place = random.choice(aircraft_data)

        mds, length, wingspan, _, _ = aircraft_to_place
        for i, (x, y, w, h) in enumerate(remaining_space):
            if w >= length + wingtip_clearance and h >= wingspan + wingtip_clearance:
                # Place aircraft
                new_spaces = [
                    (x + length + wingtip_clearance, y, w - length - wingtip_clearance, h),
                    (x, y + wingspan + wingtip_clearance, length + wingtip_clearance, h - wingspan - wingtip_clearance)
                ]
                remaining_space[i:i+1] = [s for s in new_spaces if s[2] > 0 and s[3] > 0]
                remaining_space.sort(key=lambda s: s[2] * s[3], reverse=True)
                return aircraft_to_place

        return None

    def calculate_score(self, aircraft_counts, aircraft_data):
        total_aircraft = sum(aircraft_counts.values())
        min_requirements_met = all(aircraft_counts[ac[0]] >= ac[4] for ac in aircraft_data)
        return total_aircraft + (1000 if min_requirements_met else 0)

    def fixed_placement(self, aircraft_data, usable_length, usable_width, 
                        interior_taxi_width, wingtip_clearance):
        arcpy.AddMessage(f"\nAvailable parking area: {usable_length:.2f} x {usable_width:.2f} ft")
        aircraft_counts = {}
        remaining_space = [(0, 0, usable_length, usable_width)]

        # Sort aircraft by area (length * wingspan) in descending order
        sorted_aircraft = sorted(aircraft_data, key=lambda ac: ac[1] * ac[2], reverse=True)

        for aircraft in sorted_aircraft:
            mds, length, wingspan, qty, _ = aircraft
            placed = 0

            while placed < qty and remaining_space:
                best_space = max(remaining_space, key=lambda s: s[2] * s[3])
                x, y, w, h = best_space

                # Try both orientations
                if w >= length + wingtip_clearance and h >= wingspan + wingtip_clearance:
                    orientation = "normal"
                elif w >= wingspan + wingtip_clearance and h >= length + wingtip_clearance:
                    orientation = "rotated"
                    length, wingspan = wingspan, length
                else:
                    break  # Can't fit any more of this aircraft type

                # Place aircraft
                placed += 1
                aircraft_counts[mds] = aircraft_counts.get(mds, 0) + 1

                # Update remaining space
                remaining_space.remove(best_space)
                new_spaces = [
                    (x + length + wingtip_clearance, y, w - length - wingtip_clearance, h),
                    (x, y + wingspan + wingtip_clearance, length + wingtip_clearance, h - wingspan - wingtip_clearance)
                ]
                remaining_space.extend([s for s in new_spaces if s[2] > 0 and s[3] > 0])

            arcpy.AddMessage(f"Placed {placed} {mds} aircraft (out of {qty} requested)")
            arcpy.AddMessage(f"  Orientation: {orientation}")

        return list(aircraft_counts.items())

    def display_results(self, result, apron_length, apron_width, 
                        interior_taxi_width, peripheral_taxi_width, wingtip_clearance):
        arcpy.AddMessage("\nMaximum On Ground (MOG) Results:")
        total_aircraft = 0
        total_area = 0
        for mds, count in result:
            arcpy.AddMessage(f"{mds}: {count}")
            total_aircraft += count
            aircraft_data = next((ac for ac in self.aircraft_data if ac[0] == mds), None)
            if aircraft_data:
                total_area += count * aircraft_data[1] * aircraft_data[2]

        arcpy.AddMessage(f"\nTotal aircraft placed: {total_aircraft}")
        
        estimated_rows = math.ceil(total_aircraft / 3)  # Assuming 3 columns, adjust as needed
        arcpy.AddMessage(f"Estimated rows: {estimated_rows}")

        arcpy.AddMessage(f"Airfield dimensions: {apron_length} x {apron_width} ft")
        total_area_sq_ft = apron_length * apron_width
        arcpy.AddMessage(f"Total airfield area: {total_area_sq_ft:.2f} sq ft")

        arcpy.AddMessage(f"Area occupied by aircraft: {total_area:.2f} sq ft")
        
        taxiway_area = (2 * peripheral_taxi_width * apron_width) + (interior_taxi_width * apron_length)
        arcpy.AddMessage(f"Area used for taxiways: {taxiway_area:.2f} sq ft")

        total_used_area = total_area + taxiway_area
        arcpy.AddMessage(f"Total used area: {total_used_area:.2f} sq ft")

        remaining_area = total_area_sq_ft - total_used_area
        arcpy.AddMessage(f"Remaining area: {remaining_area:.2f} sq ft")

        remaining_percentage = (remaining_area / total_area_sq_ft) * 100
        arcpy.AddMessage(f"Remaining space: {remaining_percentage:.2f}% of total airfield area")

        utilization_efficiency = ((total_used_area / total_area_sq_ft) * 100)
        arcpy.AddMessage(f"Space utilization efficiency: {utilization_efficiency:.2f}%")
