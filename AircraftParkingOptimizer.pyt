import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import numpy as np
import joblib

# Function to process the aircraft characteristics and apply AI model logic
def optimize_aircraft_parking(data, model):
    # Prompt the user to enter the quantity of each aircraft
    quantities = []
    for i, row in data.iterrows():
        aircraft_type = f"Aircraft {i+1} (Wing Span: {row['WING_SPAN']}, Length: {row['LENGTH']})"
        quantity = int(input(f"Enter the quantity of {aircraft_type}: "))
        quantities.append(quantity)

    # Add the quantities to the data frame
    data['Quantity'] = quantities

    # Use the trained model to predict optimized parking
    X = data[['WING_SPAN', 'LENGTH']].values  # Input features
    predictions = model.predict(X)

    # Adjust predictions based on quantity
    data['OptimizedParking'] = predictions * data['Quantity']

    return data

def train_model(data):
    # Define input features and target variable
    X = data[['WING_SPAN', 'LENGTH']].values  # Input features
    y = data['MIN_TWY_WIDTH'].values # Target variable

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Initialize the linear regression model
    model = LinearRegression()

    # Train the model
    model.fit(X_train, y_train)

    # Evalueate the model
    y_pred = model.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)

    print(f"Model Training Completed. MSE: {mse}, R^2: {r2}")

    # Save the trained model to a file
    joblib.dump(model, "aircraft_parking_model.pkl")

    return model

if __name__ == "__main__":
    # Load the aircraft characteristics from the CSV file
    input_file = "ACFT_Characteristics(ACFT Characteristics).csv"
    aircraft_data = pd.read_csv(input_file)

    # Train the model with the dataset
    model = train_model(aircraft_data)

    # Load the trained model
    model = joblib.load("aircraft_parking_model.pkl")

    # Perform optimization on the aircraft data using the AI logic
    optimized_data = optimize_aircraft_parking(aircraft_data, model)

    # Save the optimized data to a new CSV file
    output_file = "Optimized_Aircraft_Parking.csv"
    optimized_data.to_csv(output_file, index=False)

    # Print completion message
    print(f"Optimization completed! Results saved to {output_file}.")