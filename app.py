from flask import Flask, render_template, request, jsonify
from calculations import get_water_geometries_near_point, calculate_proximity_to_water

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/calculate', methods=['POST'])
def calculate():
    try:
        # Get data from JSON instead of form
        data = request.get_json()
        lat = float(data['latitude'])
        lon = float(data['longitude'])

        # Fetch water geometries
        water_gdf = get_water_geometries_near_point(lat, lon)
        
        if water_gdf is None:
            result = "No water bodies found within the specified radius."
        else:
            # Calculate proximity
            proximity = calculate_proximity_to_water(lat, lon, water_gdf)
            if proximity == 9999:
                result = "Unable to calculate proximity due to an error."
            else:
                result = f"The nearest water body is approximately {proximity:.2f} meters away."

        return jsonify(success=True, result=result)

    except Exception as e:
        # Debugging information
        return jsonify(success=False, error=f"An error occurred: {e}")

if __name__ == '__main__':
    app.run(debug=True)

