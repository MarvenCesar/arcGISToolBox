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

            if save_as_gdb:
                gdb_path = os.path.dirname(out_table)
                table_name = os.path.basename(out_table)
                arcpy.conversion.TableToTable(temp_csv, gdb_path, table_name)
                arcpy.AddMessage(f"Successfully saved as GDB table: {out_table}")
            else:
                arcpy.conversion.TableToTable(temp_csv, os.path.dirname(out_table), os.path.basename(out_table))
                arcpy.AddMessage(f"Successfully imported aircraft records as CSV table.")
            
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

    def create_aircraft_shape(self, x_start, y_start, length, wingspan):
        fuselage_width = length * 0.1
        nose_length = length * 0.2
        tail_length = length * 0.15
        wing_sweep = length * 0.1
        tail_sweep = length * 0.05

        corners = [
            arcpy.Point(x_start, y_start + length / 2),  
            arcpy.Point(x_start - fuselage_width / 2, y_start + length / 2 - nose_length),
            arcpy.Point(x_start - wingspan / 2, y_start + wing_sweep), 
            arcpy.Point(x_start - wingspan / 2, y_start),  
            arcpy.Point(x_start - wingspan / 2, y_start - wing_sweep),  
            arcpy.Point(x_start - fuselage_width / 2, y_start - length / 2 + tail_length),  
            arcpy.Point(x_start - wingspan / 4, y_start - length / 2),  
            arcpy.Point(x_start, y_start - length / 2 - tail_sweep),  
            arcpy.Point(x_start + wingspan / 4, y_start - length / 2),  
            arcpy.Point(x_start + fuselage_width / 2, y_start - length / 2 + tail_length),  
            arcpy.Point(x_start + wingspan / 2, y_start - wing_sweep),  
            arcpy.Point(x_start + wingspan / 2, y_start),  
            arcpy.Point(x_start + wingspan / 2, y_start + wing_sweep),  
            arcpy.Point(x_start + fuselage_width / 2, y_start + length / 2 - nose_length),  
            arcpy.Point(x_start, y_start + length / 2)  
        ]
        return corners

    def execute(self, parameters, messages):
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
            workspace = os.path.dirname(out_fc)
            valid_name = arcpy.ValidateTableName(os.path.basename(out_fc), workspace)
            out_fc = os.path.join(workspace, valid_name)
            sr = arcpy.Describe(airfield_layer).spatialReference
            arcpy.CreateFeatureclass_management(workspace, valid_name, "POLYGON", spatial_reference=sr)

            arcpy.AddField_management(out_fc, "MDS", "TEXT")
            arcpy.AddField_management(out_fc, "LENGTH", "DOUBLE")
            arcpy.AddField_management(out_fc, "WINGSPAN", "DOUBLE")
            arcpy.AddField_management(out_fc, "Aircraft_Footprint", "DOUBLE")

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

                        length_in_degrees = length / 364000
                        wingspan_in_degrees = wingspan / 364000

                        with arcpy.da.InsertCursor(out_fc, ["SHAPE@", "MDS", "LENGTH", "WINGSPAN", "Aircraft_Footprint"]) as insert_cursor:
                            points_placed = 0
                            row_index = 0
                            col_index = 0

                            while points_placed < quantity_of_aircraft:
                                x_start = start_lon + (col_index * (wingspan_in_degrees + buffer_distance / 364000))
                                y_start = start_lat + (row_index * (length_in_degrees + buffer_distance / 364000))

                                corners = self.create_aircraft_shape(x_start, y_start, length_in_degrees, wingspan_in_degrees)
                                polygon = arcpy.Polygon(arcpy.Array(corners), sr)
                                insert_cursor.insertRow([polygon, mds, length, wingspan, length * wingspan])
                                points_placed += 1

                                arcpy.AddMessage(f"Placed {mds} at ({x_start}, {y_start}) - {points_placed}/{quantity_of_aircraft}")

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
        self.description = "Calculate the maximum number of aircraft that can park on a specific apron with specified clearances and optional constraints (LCN compatibility) and taxi lane area."

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
                displayName="Aircraft Clearance (ft)",
                name="aircraft_clearance",
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
                parameterType="Required",
                direction="Input"),
            arcpy.Parameter(
                displayName="Apply LCN Compatibility?",
                name="apply_lcn",
                datatype="GPBoolean",
                parameterType="Optional",
                direction="Input"),
            arcpy.Parameter(
                displayName="Taxi Lane Width (ft)",
                name="taxi_lane_width",
                datatype="GPDouble",
                parameterType="Optional",
                direction="Input")
        ]
        params[1].filter.list = ["Polygon"]
        return params

    def updateParameters(self, parameters):
        if parameters[1].altered and not parameters[2].altered:
            airfield_layer = parameters[1].valueAsText
            if airfield_layer:
                airfield_names = [row[0] for row in arcpy.da.SearchCursor(airfield_layer, "AFLD_NAME")]
                parameters[2].filter.list = sorted(set(airfield_names))
        
        if parameters[0].altered and not parameters[4].altered:
            aircraft_table = parameters[0].valueAsText
            if aircraft_table:
                try:
                    aircraft_names = set()
                    with arcpy.da.SearchCursor(aircraft_table, ["MDS"]) as cursor:
                        for row in cursor:
                            aircraft_names.add(row[0])
                    parameters[4].filter.list = sorted(aircraft_names)
                except Exception as e:
                    arcpy.AddError(f"Error loading aircraft names: {e}")
        return

    def execute(self, parameters, messages):
        in_table = parameters[0].valueAsText
        airfield_layer = parameters[1].valueAsText
        airfield_name = parameters[2].valueAsText
        aircraft_clearance = float(parameters[3].valueAsText)
        selected_aircraft = parameters[4].values
        aircraft_quantities = parameters[5].valueAsText.split(',')
        apply_lcn = parameters[6].value
        taxi_lane_width = float(parameters[7].valueAsText) if parameters[7].valueAsText else 0

        if len(selected_aircraft) != len(aircraft_quantities):
            arcpy.AddError(f"Number of selected aircraft does not match the number of quantities provided.")
            return

        try:
            aircraft_quantities = [int(q.strip()) for q in aircraft_quantities]
        except ValueError:
            arcpy.AddError(f"Invalid input for aircraft quantities. Ensure all quantities are valid integers.")
            return

        try:
            where_clause = f"AFLD_NAME = '{airfield_name}'"
            with arcpy.da.SearchCursor(airfield_layer, ["SHAPE@", "LENGTH", "WIDTH", "LCN"], where_clause) as cursor:
                for row in cursor:
                    airfield_shape, apron_length, apron_width, apron_lcn = row
                    break
                else:
                    arcpy.AddError(f"Airfield '{airfield_name}' not found.")
                    return

            arcpy.AddMessage(f"Airfield dimensions: Length={apron_length} ft, Width={apron_width} ft")
            arcpy.AddMessage(f"Airfield LCN: {apron_lcn}")

            taxi_lane_area = taxi_lane_width * apron_length
            available_area = (float(apron_length) * float(apron_width)) - taxi_lane_area
            arcpy.AddMessage(f"Adjusted available area (excluding taxi lanes): {available_area} sq ft")

            aircraft_data = []
            
            with arcpy.da.SearchCursor(in_table, ["MDS", "LENGTH", "WING_SPAN", "ACFT_LCN"]) as cursor:
                for row in cursor:
                    mds, length, wingspan, aircraft_lcn = row
                    if mds in selected_aircraft:
                        footprint = float(length) * float(wingspan)
                        index = selected_aircraft.index(mds)
                        quantity = aircraft_quantities[index]

                        if apply_lcn and float(aircraft_lcn) > float(apron_lcn):
                            arcpy.AddWarning(f"Aircraft {mds} LCN ({aircraft_lcn}) exceeds apron LCN ({apron_lcn}). Skipping")
                            continue

                        aircraft_data.append((mds, length, wingspan, footprint, quantity))

            aircraft_data.sort(key=lambda x: x[3], reverse=True)

            total_area_used = 0
            results = []
            for mds, length, wingspan, footprint, quantity in aircraft_data:
                clearance_area = (length + aircraft_clearance) * (wingspan + aircraft_clearance)
                max_fit = math.floor((available_area - total_area_used) / clearance_area)
                actual_fit = min(max_fit, quantity)
                total_area_used += actual_fit * clearance_area
                results.append((mds, actual_fit))
                arcpy.AddMessage(f"Placed {actual_fit} {mds} aircraft")

            arcpy.AddMessage("\nMaximum On Ground (MOG) Results:")
            for mds, count in results:
                arcpy.AddMessage(f"{mds}: {count}")

            arcpy.AddMessage(f"\nTotal area used: {total_area_used:.2f} sq ft")
            arcpy.AddMessage(f"Remaining area: {available_area - total_area_used:.2f} sq ft")

        except Exception as e:
            arcpy.AddError(f"An error occurred while calculating MOG: {str(e)}")
            arcpy.AddError(arcpy.GetMessages())

        return
