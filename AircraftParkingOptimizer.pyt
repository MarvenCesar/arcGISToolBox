import arcpy
import pandas as pd
import logging
import heapq

logging.basicConfig(filename='AircraftParkingOptimizer.log', level=logging.INFO)

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the .pyt file)."""
        self.label = "Aircraft Parking Optimizer"
        self.alias = "MOG Optimization"

    def getParameterInfo(self):
        """Define parameter definitions"""
        param0 = arcpy.Parameter(
            displayName="Input Aircraft Data",
            name="aircraft_data",
            datatype="DETable",
            parameterType="Required",
            direction="Input")
        
        param1 = arcpy.Parameter(
            displayName="Output Optimized Data",
            name="output_optimized_data",
            datatype="DETable",
            parameterType="Required",
            direction="Output")
        
        params = [param0, param1]
        return params
    
    def optimize_aircraft_parking(aircraft_data, parking_spaces, constraints):
        # Sort using a heap to reduce time complexity
        aircraft_data = heapq.nlargest(len(aircraft_data), aircraft_data, key=lambda x: x.size)
        parking_spaces = heapq.nlargest(parking_spaces, key=lambda x: x.available_aera)

        # Initialize result list
        parking_plan = []

        # Greedy placement of aircraft
        for aircraft in aircraft_data:
            for parking_space in parking_spaces:
                if can_place_aircraft(aircraft, parking_space, constraints):
                    # Assign aircraft to parking space
                    parking_plan.append((aircraft.id, space.id))
                    # Mark parking space as occupied
                    occupy_space(space, aircraft)
                    break
        return parking_plan
    
    def meets_constraints(aircraft, parking_space, constraints):
        # Check if aircraft meets constraints
        return all(constraint(aircraft, parking_space) for constraint in constraints)
    
    def update_constraints(aircraft, parking_space):
        # Update constraints based on aircraft and parking space
        for other_space in parking_space.adjacent_spaces:
            other_space.available_area -= aircraft.size

    def can_place_aircraft(aircraft, parking_space, constraints):
        # Check if aircraft can be placed in 
        # parking space
        return (aircraft.space <= parking_space.available_area) and meets_constraints(aircraft, parking_space, constraints)
    
    def occupy_parking_space(parking_space, aircraft):
        # Update parking space area
        parking_space.available_area -= aircraft.size
        # Update adjacent parking spaces
        update_constraints(aircraft, parking_space)

    def execute(self, parameters, messages):
        """The source code of the tool."""
        # Initiate logging
        logging.info("Starting Aircraft Parking Optimizer tool...")

        input_file = parameters[0].valueAsText  # Path to input file from ArcGIS
        output_file = parameters[1].valueAsText  # Path to output file
        
        # Load the data using pandas
        aircraft_data = pd.read_csv(input_file)
        
        # Print progress update to user
        messages.addMessage("Loading aircraft data...")
        
        # Run the optimization logic
        optimized_data = optimize_aircraft_parking(aircraft_data)

        # Print progress update to user
        messages.addMessage("Running optimization algorithm...")
        
        # Save the optimized data to output location
        optimized_data.to_csv(output_file, index=False)

        # Print progress update to user
        messages.addMessage("Saving optimized data...")

        # Log completion message
        logging.info(f"Optimization completed! Results saved to {output_file}.")
