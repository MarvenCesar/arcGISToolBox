import arcpy
import pandas as pd
from pseudocode import optimize_aircraft_parking  # Importing from your AI logic

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

    def execute(self, parameters, messages):
        """The source code of the tool."""
        input_file = parameters[0].valueAsText  # Path to input file from ArcGIS
        output_file = parameters[1].valueAsText  # Path to output file
        
        # Load the data using pandas
        aircraft_data = pd.read_csv(input_file)
        
        # Run the optimization logic
        optimized_data = optimize_aircraft_parking(aircraft_data)
        
        # Save the optimized data to output location
        optimized_data.to_csv(output_file, index=False)
        
        messages.addMessage(f"Optimization completed! Results saved to {output_file}.")
