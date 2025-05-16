#!/usr/bin/env python3
"""
Script to fetch the nearest GBIF occurrence with coordinates for a list of species,
calculate the distance to Serra da Estrela, and save the results to a CSV.
"""

import pandas as pd
import requests
import math
import csv
from time import sleep
from concurrent.futures import ThreadPoolExecutor, as_completed
import json


# Fixed coordinates for Serra da Estrela (reference location)
LATITUDE = 40.404139
LONGITUDE = -7.538167

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculates the great-circle distance (in kilometers) between two points
    on the Earth's surface using the Haversine formula.
    """
    r = 6378  # Earth's radius in kilometers
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    d_phi = math.radians(lat2 - lat1)
    d_lambda = math.radians(lon2 - lon1)

    a = (math.sin(d_phi / 2) ** 2
         + math.cos(phi1) * math.cos(phi2) * math.sin(d_lambda / 2) ** 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return r * c


def gbif_request_with_retry(url, params=None, max_retries=5, backoff_factor=1):
    for attempt in range(max_retries):
        try:
            response = requests.get(url, params=params, timeout=20)
            if response.status_code == 503:
                wait = backoff_factor * (2 ** attempt)
                print(f"503 error. Waiting {wait}s before retrying...")
                sleep(wait)
                continue

            response.raise_for_status()

            # Try to parse JSON to catch decode errors early and retry
            try:
                _ = response.json()
            except json.JSONDecodeError as e:
                wait = backoff_factor * (2 ** attempt)
                print(f"JSON decode error: {e}. Retrying in {wait}s...")
                sleep(wait)
                continue

            return response

        except requests.RequestException as e:
            wait = backoff_factor * (2 ** attempt)
            print(f"Request failed: {e}. Retrying in {wait}s...")
            sleep(wait)

    print(f"Failed to fetch data from GBIF after {max_retries} attempts. URL: {url}")
    return None

def get_species_key(species_name):
    """
    Retrieves the GBIF usageKey for a given scientific name via the species/match endpoint.
    Returns the usageKey and the type of match (exact, fuzzy, or none).
    """
    url = "https://api.gbif.org/v1/species/match"
    params = {"name": species_name}
    response = gbif_request_with_retry(url, params)
    if response is None:
        return None, "none"
    data = response.json()
    usage_key = data.get("usageKey")
    match_type = data.get("matchType", "").lower()
    if usage_key and match_type in ("exact", "fuzzy"):
        return usage_key, match_type
    return None, "none"

def fetch_nearest_occurrence(species_name):
    """
    Fetches the nearest GBIF occurrence of a species to Serra da Estrela.
    Returns a dictionary with metadata of the nearest occurrence (if any).
    """
    species_key, match_type = get_species_key(species_name)
    if match_type == "none":
        return None, species_name, False  # To check manually later

    if species_key:
        params = {
            "speciesKey": species_key,
            "hasCoordinate": "true",
            "limit": 300
        }
    else:
        params = {
            "scientificName": species_name,
            "hasCoordinate": "true",
            "limit": 300
        }

    url = "https://api.gbif.org/v1/occurrence/search"
    response = gbif_request_with_retry(url, params)
    if response is None:
        return None, species_name, False

    try:
        data = response.json()
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON for species {species_name}: {e}")
        return None, species_name, False

    occurrences = data.get("results", [])
    if not occurrences:
        return None, species_name, False

    nearest = None
    min_distance = float("inf")

    for occ in occurrences:
        lat = occ.get("decimalLatitude")
        lon = occ.get("decimalLongitude")
        if lat is not None and lon is not None:
            distance = haversine(LATITUDE, LONGITUDE, lat, lon)
            if distance < min_distance:
                min_distance = distance
                nearest = {
                    "species_name": species_name,
                    "lat": lat,
                    "lon": lon,
                    "date": occ.get("eventDate", ""),
                    "country": occ.get("country", ""),
                    "occurrenceID": occ.get("occurrenceID", ""),
                    "distance_km": round(distance, 2),
                    "match_type": match_type
                }
    if nearest:
        return nearest, None, False
    else:
        # No occurrences with coordinates but occurrences exist
        return None, species_name, True


def process_species_row(row, idx, total):
    """
    Processes a single row of the DataFrame to find the nearest GBIF occurrence.
    Returns a tuple: (result, unmatched_name, no_coordinates_flag)
    """
    species = str(row.get("Species", "")).strip()
    genus = str(row.get("Genus", "")).strip()

    if not genus or not species:
        print(f"Skipping incomplete row {idx}")
        return (None, None, False)

    if species.lower().startswith(genus.lower()):
        species = species[len(genus):].strip()

    species_name = f"{genus} {species}"
    print(f"({idx + 1}/{total}) Processing species: {species_name}")

    occurrence, failed_name, is_no_coordinates = fetch_nearest_occurrence(species_name)
    if occurrence:
        return (occurrence, None, False)
    elif failed_name:
        species_key, _ = get_species_key(species_name)
        if species_key:
            return (None, failed_name, True)  # No coordinates
        else:
            return (None, failed_name, False)  # Not found in GBIF

    return (None, None, False)

def main():
    input_csv = "results/blast/ncbi_taxonomy_lookup.csv"
    output_csv = "results/blast/gbif_nearest_occurrences.csv"
    manual_check_csv = "results/blast/gbif_names_to_check_manually.csv"

    df = pd.read_csv(input_csv)
    results = []
    to_check_manually = []
    no_coordinates = []

    total = len(df)

    with ThreadPoolExecutor(max_workers=5) as executor:
        futures = {
            executor.submit(process_species_row, row, idx, total): idx
            for idx, row in df.iterrows()
        }

        for future in as_completed(futures):
            occurrence, failed_name, is_no_coordinates = future.result()
            if occurrence:
                results.append(occurrence)
            elif failed_name:
                if is_no_coordinates:
                    no_coordinates.append({"species_name": failed_name})
                else:
                    to_check_manually.append({"unmatched_name": failed_name})

    results.sort(key=lambda x: x["distance_km"])
    # Saves occurrences with valid coordenates
    headers = ["species_name", "lat", "lon", "date", "country", "occurrenceID", "distance_km", "match_type"]
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    if to_check_manually:
        with open(manual_check_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["unmatched_name"])
            writer.writeheader()
            for item in to_check_manually:
                writer.writerow(item)

    if no_coordinates:
        with open("results/blast/gbif_species_with_no_coordinates.csv", "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["species_name"])
            writer.writeheader()
            for item in no_coordinates:
                writer.writerow(item)

    print(f"\nSaved {len(results)} valid nearest occurrences to {output_csv}")
    print(f"{len(to_check_manually)} unmatched names saved to {manual_check_csv}")
    print(f"{len(no_coordinates)} species had occurrences but no coordinates. Saved to gbif_species_with_no_coordinates.csv")

if __name__ == "__main__":
    main()
