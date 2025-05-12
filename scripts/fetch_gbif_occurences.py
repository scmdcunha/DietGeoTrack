#!/usr/bin/env python3
"""
Script to fetch the nearest GBIF occurrence with coordinates for a list of species,
calculate the distance to Serra da Estrela, and save the results to a CSV.
"""

from numpy.lib.index_tricks import diff
import pandas as pd
import requests
import csv
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from time import sleep

# Coordinates of Manteigas, Serra da Estrela, Portugal
latitude = 40.404139
longitude = -7.538167

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance between two points on Earth using the Haversine formula.

    Parameters:
        lat1, lon1 (float): Latitude and longitude of the first point.
        lat2, lon2 (float): Latitude and longitude of the second point.

    Returns:
        float: Distance between the two points in kilometers.
    """
    R = 6378
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    d_phi = math.radians(lat2 - lat1)
    d_lambda = math.radians(lon2 - lon1)

    # Haversine formula
    a = math.sin(d_phi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(d_lambda / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c

output_file = "results/blast/gbif_occurrences_nearest.csv"

# Write header to the output CSV file
with open(output_file, "w", newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Species", "Lat", "Lon", "Date", "Country", "ID", "Distance"])

# Function to fetch speciesKey from GBIF using the species name
def get_species_key(species_name):
    """
    Retrieve the GBIF speciesKey for a given species name using the GBIF species match API.

    Parameters:
        species_name (str): The scientific name of the species.

    Returns:
        int or None: The speciesKey if found, otherwise None.
    """
    url = "https://api.gbif.org/v1/species/match"
    params = {"name": species_name}
    try:
        response = requests.get(url, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()
        return data.get("usageKey")
    except Exception as e:
        print(f"  Error fetching speciesKey for {species_name}: {e}")
        return None

# Function to fetch the nearest occurrence of a species from GBIF
def fetch_nearest_occurrence(species_name):
    """
    Retrieve the nearest occurrence (with coordinates) of a species from GBIF and compute its distance
    to Manteigas (Serra da Estrela).

    Parameters:
        species_name (str): The scientific name of the species.

    Returns:
        list or None: List with species name, latitude, longitude, date, country, occurrenceID, and distance (km),
                      or None if no occurrence is found or an error occurs.
    """
    species_key = get_species_key(species_name)
    if not species_key:
        return None

    url = "https://api.gbif.org/v1/occurrence/search"
    params = {
        "speciesKey": species_key,
        "hasCoordinate": "true",
    }
    try:
        response = requests.get(url, params=params, timeout=15)
        response.raise_for_status()
        data = response.json()
        occurrences = data.get("results", [])

        nearest = None
        min_distance = float("inf")

        # Loop through occurrences and calculate distances
        for occ in occurrences:
            lat = occ.get("decimalLatitude")
            lon = occ.get("decimalLongitude")
            if lat is not None and lon is not None:
                distance = haversine(latitude, longitude, lat, lon)
                if distance < min_distance:
                    min_distance = distance
                    nearest = [
                        species_name,
                        lat,
                        lon,
                        occ.get("eventDate"),
                        occ.get("country"),
                        occ.get("occurrenceID"),
                        round(distance, 2)
                    ]

        return nearest

    except Exception as e:
        print(f"  Error fetching occurrences for {species_name}: {e}")
        return None

def main():
    """
    Main execution function:
    - Reads species list from CSV.
    - Fetches nearest GBIF occurrence for each species using parallel threads.
    - Saves the results to a new CSV file.
    """
    input_file = "results/blast/ncbi_taxonomy_lookup.csv"
    output_file = "results/blast/gbif_occurrences_nearest.csv"

    df_species = pd.read_csv(input_file)
    species_list = df_species["Species"].dropna().unique()

    with open(output_file, "w", newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Species", "Lat", "Lon", "Date", "Country", "ID", "Distance"])

    results = []
    with ThreadPoolExecutor(max_workers=10) as executor:
        future_to_species = {
            executor.submit(fetch_nearest_occurrence, species): species
            for species in species_list
        }

        for i, future in enumerate(as_completed(future_to_species)):
            species = future_to_species[future]
            try:
                result = future.result()
                if result:
                    results.append(result)
                    print(f"[{i+1}/{len(species_list)}] {species} → {result[-1]} km")
            except Exception as e:
                print(f"  Failed for {species}: {e}")
            sleep(0.1)  # Small delay to avoid rate limits

    with open(output_file, "a", newline='') as file:
        writer = csv.writer(file)
        writer.writerows(results)

    print(f"\nDone. Saved {len(results)} nearest occurrences to {output_file}")

if __name__ == "__main__":
    main()
