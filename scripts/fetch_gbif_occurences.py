#!/usr/bin/env python3
"""
Script to fetch the nearest GBIF occurrence with coordinates for a list of species,
calculate the distance to Serra da Estrela, and save the results to a CSV.
"""

from os import supports_effective_ids
from numpy.lib.index_tricks import diff
import pandas as pd
import requests
import csv
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from time import sleep
import re

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


def clean_species_name(species_name):
    """
    Cleans the species name by removing any extra identifiers or codes.

    Parameters:
        species_name (str): The species name to clean.

    Returns:
        str: Cleaned species name.
    """
    # Remove 'sp.', 'cf.', 'aff.', 'nr.' (case insensitive)
    species_name = re.sub(r'\b(sp|cf|aff|nr)\.?\b', '', species_name, flags=re.IGNORECASE)
    species_name = re.sub(r'\.+', '', species_name)
    species_name = re.sub(r'\s+', ' ', species_name).strip()

    # Reduce to a maximum of two terms (Genus + Species)
    parts = species_name.split()
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    else:
        return None


# Function to fetch speciesKey from GBIF using the species name
def get_species_key(species_name):
    """
    Tries to retrieve the GBIF speciesKey using the full cleaned name,
    and falls back to a simpler version if the match fails.

    Parameters:
        species_name (str): The original species name from input data.

    Returns:
        int or None: The GBIF speciesKey, or None if not found.
    """
    def fetch_key(name):
        """Helper function to request speciesKey from GBIF API."""
        try:
            response = requests.get("https://api.gbif.org/v1/species/match", params={"name": name}, timeout=10)
            response.raise_for_status()
            data = response.json()
            if data.get("matchType") != "NONE" and data.get("usageKey"):
                return data.get("usageKey")
        except Exception as e:
            print(f"  Error fetching speciesKey for '{name}': {e}")
        return None

    cleaned_name = clean_species_name(species_name)
    key = fetch_key(cleaned_name)

    if not key:
        # Fallback: Try only the first two words (Genus + Species)
        fallback_name = ' '.join(cleaned_name.strip().split()[:2])
        if fallback_name != cleaned_name:
            print(f"  '{cleaned_name}' → trying fallback '{fallback_name}'")
            key = fetch_key(fallback_name)

    if not key:
        print(f"speciesKey not found for '{species_name}'")

    return key

# Function to fetch the nearest occurrence of a species from GBIF
def fetch_nearest_occurrence(species_name):
    """
    Fetches the nearest occurrence with coordinates for a species from GBIF,
    and computes the distance to Serra da Estrela.

    Parameters:
        species_name (str): The full species name to search for.

    Returns:
        list or None: A list containing [Species, Lat, Lon, Date, Country, ID, Distance],
                      or None if no valid occurrence is found.
    """
    species_key = get_species_key(species_name)
    if species_key:
        params = {
            "speciesKey": species_key,
            "hasCoordinate": "true",
            "limit": 300
        }
    else:
        print(f"No speciesKey for '{species_name}', trying fallback with scientificName.")
        params = {
            "scientificName": species_name.strip(),
            "hasCoordinate": "true",
        }
    url = "https://api.gbif.org/v1/occurrence/search"

    try:
        response = requests.get(url, params=params, timeout=15)
        response.raise_for_status()
        data = response.json()
        occurrences = data.get("results", [])
        print(f"{species_name}: {len(occurrences)} occurrences found.")

        if not occurrences:
            print(f"No occurrences with coordinates for {species_name}")
            return None

        nearest = None
        min_distance = float("inf")

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

        if nearest:
            return nearest
        else:
            print(f"No valid coordinates found in occurrences for {species_name}")
            return None

    except Exception as e:
        print(f"Error fetching occurrences for {species_name}: {e}")
        return None

def main():
    """
    Main execution function:
    - Reads species list from CSV.
    - Fetches nearest GBIF occurrence for each species using parallel threads.
    - Saves the results to a CSV file.
    - Logs species with no valid occurrences to a separate CSV file.
    """
    input_file = "results/blast/ncbi_taxonomy_lookup.csv"
    output_file = "results/blast/gbif_occurrences_nearest.csv"
    failed_file = "results/blast/failed_species.csv"

    df_species = pd.read_csv(input_file)
    species_list = [
        name for s in df_species["Species"].dropna().unique()
        if (name := clean_species_name(s)) is not None
    ]

    # Initialize output files with headers
    with open(output_file, "w", newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Species", "Lat", "Lon", "Date", "Country", "ID", "Distance"])

    with open(failed_file, "w", newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Species"])

    results = []
    failed_species = []

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
                else:
                    failed_species.append([species])
            except Exception as e:
                print(f"  Failed for {species}: {e}")
                failed_species.append([species])
            sleep(0.1)  # Small delay to avoid rate limits

    # Write successful occurrences
    with open(output_file, "a", newline='') as file:
        writer = csv.writer(file)
        writer.writerows(results)

    # Write failed species
    with open(failed_file, "a", newline='') as file:
        writer = csv.writer(file)
        writer.writerows(failed_species)

    print(f"\nDone. Saved {len(results)} nearest occurrences to {output_file}")
    print(f"Logged {len(failed_species)} species with no valid occurrence to {failed_file}")

if __name__ == "__main__":
    main()
