This Repo will be used for the development of our Python toolbox for ArcGIS. 

The MOG Tools.pyt file performs the calculation for the maximum-on-ground using a simple algorithm to see if the aircraft can be placed on the airfield if there is any space left and then it places polygons on the airfield to represent the aircraft and a bigger one surrounding it to represent the constraints. 

The MOG Tools AI.pyt performs the calculation using a genetic algorithm. It randomly places the aircraft on the airfield, chooses the 2 best solutions, merges them, adds some randomness to it and repeats this process for 100 times. Our group did not have time to fine tune the algorithm to be able to visualize the results on the map but it shows the text output with amount of aircraft placed, locations, and space available. 

For both files, we first select the ACFT_Characteristics(ACFT Characteristics) table from the database, which contains the aircraft data, select the section of the airfield we want to place the aircraft and the aircraft we want to place on the airfield from the drop down menu. 

For the toolboxes to work, we also need to import the database files included in the Sponsor Documents from Teams into a new project in ArcGIS by dragging and dropping the files. We can then import the toolboxes by right-clicking in the toolbox section on the right side of the screen and selecting to add a new toolbox. From there we can then navigate to the location where the file is saved