import pandas as pd
import requests
import time
import csv
import math

# Coordinates of Manteigas, Serra da Estrela, Portugal
latitude = 40.404139
longitude = -7.538167

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance between two points on Earth.
    Returns distance in kilometers.
    """
    R = 6378  # Earth radius in kilometers
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    d_phi = math.radians(lat2 - lat1)
    d_lambda = math.radians(lon2 - lon1)

    a = math.sin(d_phi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(d_lambda / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c

df = pd.read_csv("results/blast/ncbi_taxonomy_lookup.csv")
species_list = df["Species"].dropna().unique()

output_file = "results/gbif_occurrences.csv"

with open(output_file, "w", newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Species", "Lat", "Lon", "Date", "Country", "ID", "Distance"])

def fetch_occurrences(species_name, limit=5):
    """
    Query the GBIF API for occurrence records of a given species.

        Parameters:
            species_name (str): Scientific name of the species.
            limit (int): Maximum number of records to retrieve.

        Returns:
            list: A list of occurrence records (dictionaries).
    """
    url = "https://api.gbif.org/v1/occurrence/search"
    params = {
        "scientificName": species_name,
        "hasCoordinate": "true",
        "limit": limit
    }
    response = requests.get(url, params=params)
    data = response.json()
    return data.get("results", [])

# Loop over each species and fetch occurrence records
for species in species_list:
    print(f"\nSpecies: {species}")
    occurrences = fetch_occurrences(species)

    if not occurrences:
        print("  No occurrences found.")
        continue


    rows = []
    for occ in occurrences:
        lat = occ.get("decimalLatitude")
        lon = occ.get("decimalLongitude")
        date = occ.get("eventDate")
        country = occ.get("country")
        occ_id = occ.get("occurrenceID")

        if lat is not None and lon is not None:
            distance = haversine(latitude, longitude, lat, lon)
        else:
            distance = None

        print(
            f"Lat: {lat}, "
                f"Lon: {lon}, "
                f"Date: {date}, "
                f"Country: {country}, "
                f"ID: {occ_id}, "
                f"Distance: {distance:.2f} km" if distance is not None else "Distance: N/A"
            )

        rows.append([species, lat, lon, date, country, occ_id, distance])

    with open(output_file, mode='a', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(rows)
