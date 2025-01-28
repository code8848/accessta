import osmnx as ox
import geopandas as gpd
from shapely.geometry import Point

def find_best_epsg(lat, lon):
    """
    Determine the best UTM CRS for a given latitude and longitude.

    Parameters:
    - lat (float): Latitude of the location.
    - lon (float): Longitude of the location.

    Returns:
    - best_crs (CRS): The best UTM CRS for the given location.
    """
    try:
        # Create a GeoDataFrame for the given point
        point = gpd.GeoDataFrame(
            {'geometry': [Point(lon, lat)]},
            crs="EPSG:4326"  # WGS84
        )
        # Estimate the best UTM CRS
        best_crs = point.estimate_utm_crs()
        print(f"Best UTM CRS determined: {best_crs.to_string()}")
        return best_crs

    except Exception as e:
        print(f"An error occurred while determining the best EPSG: {e}")
        return None


def get_water_geometries_near_point(lat, lon, dist=5000):
    """
    Fetch water body geometries near a specific point using OSMnx.

    Parameters:
    - lat (float): Latitude of the center point.
    - lon (float): Longitude of the center point.
    - dist (int): Distance (in meters) to search for water bodies.

    Returns:
    - water_gdf (GeoDataFrame): GeoDataFrame of water geometries.
    """
    try:
        # Define tags for water features
        tags = {"natural": "water"}
        
        # Fetch water features using OSMnx
        water = ox.features.features_from_point((lat, lon), tags, dist=dist)

        if water.empty:
            print("No water bodies found within the specified radius.")
            return None

        # Set CRS to EPSG:4326 (WGS84)
        water = water.set_crs("EPSG:4326", allow_override=True)
        print("Water geometries fetched successfully.")
        return water

    except Exception as e:
        print(f"An error occurred while fetching water geometries: {e}")
        return None


def calculate_proximity_to_water(lat, lon, water_gdf):
    """
    Calculate the minimum distance between a user point and water geometries.

    Parameters:
    - lat (float): Latitude of the user's location.
    - lon (float): Longitude of the user's location.
    - water_gdf (GeoDataFrame): GeoDataFrame of water geometries.

    Returns:
    - distance (float): Minimum distance in meters to the nearest water body.
    """
    if water_gdf is None:
        print("No water geometries provided.")
        return 9999

    try:
        # Create a GeoDataFrame for the user's location
        user_point = gpd.GeoDataFrame(
            {'geometry': [Point(lon, lat)]},
            crs="EPSG:4326"  # WGS84
        )

        # Find the best UTM CRS
        best_crs = find_best_epsg(lat, lon)

        if not best_crs:
            print("Failed to determine the best EPSG. Using default.")
            return 9999

        # Reproject both the user point and water geometries to the best CRS
        user_point = user_point.to_crs(best_crs)
        water_gdf = water_gdf.to_crs(best_crs)

        # Combine all water geometries into a single geometry
        water_geometry = water_gdf.geometry.union_all()

        # Calculate the distance (in meters)
        distance = user_point.distance(water_geometry).min()
        return distance

    except Exception as e:
        print(f"An error occurred while calculating proximity: {e}")
        return 9999
