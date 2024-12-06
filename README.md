# Maximum On Ground (MOG) Calculation Toolbox

## Overview
The **Maximum On Ground (MOG) Calculation Toolbox** is a custom ArcGIS Python toolbox designed to optimize the parking layout of Air Mobility Command (AMC) aircraft on military airfields. The toolbox leverages ArcGIS's geoprocessing capabilities, powered by ArcPy, to calculate the maximum number of aircraft that can be parked on a given airfield while considering real-world constraints such as airfield dimensions, taxiway widths, and aircraft specifications. 

This tool is tailored for airfield operators, planners, and decision-makers, enabling them to make data-driven decisions about parking configurations that maximize efficiency and operational effectiveness.

---

## Features
### Core Features
1. **Standard and Rotated Parking Configurations**
   - Automatically evaluates both standard and rotated parking layouts to determine the most efficient configuration.

2. **Advanced Settings**
   - Allows users to manually specify airfield dimensions, taxiway widths, and other parameters for custom scenarios.

3. **Use Advanced AI**
   - Offers an optional **Genetic Algorithm Optimization** feature to provide AI-enhanced solutions for aircraft parking.

4. **Data Validation**
   - Built-in error handling to ensure input data integrity and guide users when errors are detected.

5. **Customizable Outputs**
   - Generates feature classes showing parking layouts and constraint areas.

6. **User-Friendly Integration**
   - Packaged as a Python toolbox (.pyt), making it easy to use within the ArcGIS Pro interface.

---

## Requirements

### Software
- **ArcGIS Pro** (Required)
  - Must be installed with an active license.
- **Python 3.x**
  - Included with ArcGIS Pro; ensure the ArcGIS Pro Python environment is active.

### Dependencies
- **ArcPy**
  - Part of the ArcGIS Pro Python environment; no separate installation is needed.
- **Standard Python Libraries**
  - Includes `os`, `math`, `traceback`, and `random`.

---

## Installation

### Adding the Toolbox to ArcGIS Pro
1. **Clone the Repository**
   ```bash
   git clone https://github.com/MarvenCesar/arcGISToolBox.git
   ```

2. **Add the Toolbox to ArcGIS Pro**
   - Open ArcGIS Pro.
   - In the **Catalog Pane**, right-click **Toolboxes** > **Add Toolbox**.
   - Navigate to the downloaded `.pyt` file and add it to your project.

3. **Configure Inputs**
   - Ensure your input data, such as airfield layers and aircraft tables, are ready and accessible.

---

## Usage

### Running the Tool in ArcGIS Pro
1. **Locate the Toolbox**
   - In the **Catalog Pane**, navigate to the toolbox and open the tool named **Calculate Maximum On Ground**.

2. **Input Parameters**
   - **Input Aircraft Table:** Select the table containing aircraft specifications.
   - **Airfield Layer:** Choose the airfield polygon layer.
   - **Select Aircraft (MDS):** Select a specific aircraft model designation (MDS).
   - **Enable Advanced Settings:** (Optional) Check this box to provide manual airfield dimensions and taxiway widths.
   - **Use Advanced AI Feature:** (Optional) Check this box to enable AI-based optimization using a genetic algorithm.
   - **Output Locations:** Specify paths for output feature classes.

3. **Execute the Tool**
   - Click **Run**. The tool will process the inputs and generate results based on your configuration.

### Outputs
- **Aircraft Positions Feature Class:**
  - Displays the calculated positions of parked aircraft on the airfield.
- **Constraint Polygons Feature Class:**
  - Highlights restricted areas based on aircraft dimensions and taxiway requirements.

---

## Configuration

### Parameters
- **Input Aircraft Table:** A table containing aircraft specifications such as length and wingspan.
- **Airfield Layer:** A feature class representing the airfield geometry.
- **Select Aircraft (MDS):** Dropdown to select a specific aircraft type.
- **Enable Advanced Settings:** Allows manual input of airfield dimensions and taxiway widths.
- **Use Advanced AI Feature:** Enables the genetic algorithm for AI-optimized parking configurations.

### Spatial Reference
- Ensure the airfield layer uses a **projected coordinate system** with linear units (e.g., feet or meters).

---

## Advanced AI Feature
The **Use Advanced AI Feature** checkbox activates a genetic algorithm to optimize parking configurations. This feature iteratively evaluates potential solutions to find the most efficient layout for aircraft parking, particularly useful for complex airfield layouts and constraints.

---

## Testing

### Unit Tests
- Test individual methods such as `calculate_taxiway_widths` and `calculate_parking_available`.

### Integration Tests
- Validate the interaction between ArcPy and Python toolbox logic.

### Functional Tests
- Perform end-to-end tests by running the tool with various input datasets to ensure outputs meet expectations.

---

## Contributing

### Steps to Contribute
1. Fork this repository.
2. Create a branch for your feature or bug fix.
3. Test your changes in ArcGIS Pro.
4. Submit a pull request with a detailed explanation of your changes.

We welcome contributions that improve functionality, usability, or documentation.

---

## License
Copyright [2024] [Marven Cesar, Alain Padron, Izzat Omar, Bruno De Nadai Mundim, Kyle Moore]

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
---

## Contact
### For questions or assistance:
- Marven Cesar (Interface developer) mcesar2022@fau.edu  

- Bruno De Nadai Mundim (AI-specialist) bdenadaimund2022@fau.edu  

- Alain Padron (Interface developer) apadron2020@fau.edu  

- Izzat Omar (Interface developer) iomar2021@fau.edu  

- Kyle Moore (AI-specialist) moorekyle2022@fau.edu  
- **GitHub Issues:** [MOG Calculation Toolbox Issues](https://github.com/your-repo/MOG-Calculation-Toolbox/issues)

---

## Acknowledgments
Special thanks to:
- The ArcGIS team for providing robust geoprocessing tools.
- Contributors and testers for refining this toolbox.
- Our sponsors, Department of Defense Air Force for their ongoing support and sponsorship.

This toolbox enhances airfield efficiency through advanced GIS and AI techniques, ensuring optimal space utilization and operational planning.
