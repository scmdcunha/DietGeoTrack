#!/usr/bin/env python3
"""
Script to find the closest GBIF occurrence(s) for each species listed in a CSV file.

For each scientific species name in the input CSV, the script queries the GBIF API
to retrieve occurrence records within a specified radius from a reference geographic
point (latitude and longitude). It caches API responses to avoid repeated queries,
calculates the geographic distance using GDAL between occurrences and the reference
point, and outputs the closest occurrence for each species.

Results are saved into a CSV file, and species without any occurrence found are saved
in a separate CSV for further manual review.
"""

import argparse
import pandas as pd
import requests
import time
import os
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from osgeo import ogr, osr
from math import radians, cos
from dateutil import parser as dateparser
from tqdm import tqdm

def calculate_distance_gdal(lat1, lon1, lat2, lon2):
    """
    Calculates the geodesic distance (in kilometers) between two points using GDAL.

    Parameters:
        lat1, lon1: Latitude and longitude of the first point.
        lat2, lon2: Latitude and longitude of the second point.

    Returns:
        Distance in kilometers.
    """
    point1 = ogr.Geometry(ogr.wkbPoint)
    point1.AddPoint(lon1, lat1)
    point2 = ogr.Geometry(ogr.wkbPoint)
    point2.AddPoint(lon2, lat2)

    sr = osr.SpatialReference()
    sr.ImportFromEPSG(4326)  # WGS84 coordinate system

    point1.AssignSpatialReference(sr)
    point2.AssignSpatialReference(sr)

    sr_m = osr.SpatialReference()
    sr_m.ImportFromEPSG(3857)  # Web Mercator projection (in meters)
    transform = osr.CoordinateTransformation(sr, sr_m)

    point1.Transform(transform)
    point2.Transform(transform)

    return point1.Distance(point2) / 1000  # Convert meters to kilometers

def bbox_from_radius(lat, lon, radius_km):
    """
    Generates a bounding box around a central point using an approximate radius.

    Parameters:
        lat, lon: Latitude and longitude of the central point.
        radius_km: Radius in kilometers.

    Returns:
        A tuple: (min_lat, max_lat, min_lon, max_lon)
    """
    delta_lat = radius_km / 111.0  # 1 degree latitude ≈ 111 km
    delta_lon = radius_km / (111.0 * cos(radians(lat)))
    return (lat - delta_lat, lat + delta_lat, lon - delta_lon, lon + delta_lon)

def get_closest_gbif(species, ref_lat, ref_lon, radius_km, cache_dir, top_n, min_year=None):
    """
    Searches the GBIF API for occurrences of a given species within a specified radius,
    and returns the closest ones to the reference point.

    Parameters:
        species: Scientific name of the species.
        ref_lat, ref_lon: Reference coordinates (e.g. sampling site).
        radius_km: Search radius in kilometers.
        cache_dir: Path to directory where JSON response will be cached.
        top_n: Number of closest occurrences to return.
        min_year: Minimum year of occurrence to consider (optional).

    Returns:
        A list of dictionaries with metadata of the closest occurrences, or None if not found.
        Also returns the species name if no occurrence was found.
    """
    cache_path = Path(cache_dir) / f"{species.replace(' ', '_')}.json"
    if cache_path.exists():
        with open(cache_path) as f:
            occurrences = json.load(f)
    else:
        min_lat, max_lat, min_lon, max_lon = bbox_from_radius(ref_lat, ref_lon, radius_km)
        base_url = "https://api.gbif.org/v1/occurrence/search"
        limit, offset = 300, 0
        occurrences = []

        # Paginate through GBIF API results
        while True:
            params = {
                "scientificName": species,
                "hasCoordinate": "true",
                "decimalLatitude": f"{min_lat},{max_lat}",
                "decimalLongitude": f"{min_lon},{max_lon}",
                "limit": limit,
                "offset": offset
            }

            tqdm.write(f"[DEBUG] Requesting {species} (offset={offset})")
            r = requests.get(base_url, params=params)
            if r.status_code != 200:
                tqdm.write(f"[Error] GBIF API call failed for {species} (status {r.status_code})")
                return None, species

            results = r.json().get("results", [])
            occurrences.extend(results)

            if offset + limit >= r.json().get("count", 0):
                break
            offset += limit
            time.sleep(0.5)  # Avoid overloading the API

        with open(cache_path, "w") as f:
            json.dump(occurrences, f)

    if not occurrences:
        return None, species

    if min_year is not None:
        filtered_occurrences = []
        for occ in occurrences:
            event_date = occ.get("eventDate", "")
            if event_date:
                try:
                    parsed_date = dateparser.parse(event_date)
                    if parsed_date.year >= min_year:
                        filtered_occurrences.append(occ)
                except (ValueError, OverflowError, TypeError):
                    pass
        occurrences = filtered_occurrences

    if not occurrences:
        return None, species

    # Compute distances and filter by radius
    occurrences_with_distance = []
    for occ in occurrences:
        lat = occ.get("decimalLatitude")
        lon = occ.get("decimalLongitude")
        if lat is None or lon is None:
            continue
        dist = calculate_distance_gdal(ref_lat, ref_lon, lat, lon)
        if dist <= radius_km:
            occ["distance_km"] = round(dist, 2)
            occurrences_with_distance.append(occ)

    if not occurrences_with_distance:
        return None, species

    sorted_occs = sorted(occurrences_with_distance, key=lambda x: x["distance_km"])
    closest_occurrences = sorted_occs[:top_n]

    formatted = []
    for occ in closest_occurrences:
        formatted.append({
            "species": species,
            "gbifKey": occ["key"],
            "lat": occ["decimalLatitude"],
            "lon": occ["decimalLongitude"],
            "distance_km": occ["distance_km"],
            "eventDate": occ.get("eventDate", ""),
            "country": occ.get("country", ""),
            "locality": occ.get("locality", ""),
            "datasetKey": occ.get("datasetKey", "")
        })
    return formatted, None

def main(args):
    """
    Main function that processes a list of species from a CSV file and finds the
    closest GBIF occurrence for each one.

    Results are saved to two CSVs: one with valid occurrences and one with missing data.
    """
    df = pd.read_csv(args.input_csv, sep=";")
    if args.column not in df.columns:
        raise ValueError(f"Column '{args.column}' not found in the input CSV.")

    species_list = df[args.column].dropna().drop_duplicates().tolist()
    print(f"Unique species to process: {len(species_list)}")

    os.makedirs(args.cache_dir, exist_ok=True)
    results, no_hits = [], []

    # If an output file already exists, load the existing results
    already_processed = set()
    if Path(args.output_csv).exists():
        existing_df = pd.read_csv(args.output_csv)
        already_processed = set(existing_df["species"].dropna().unique())
        results = existing_df.to_dict(orient="records")
        print(f"Resuming: {len(already_processed)} species already processed.")

    if Path(args.no_occurrences_csv).exists():
        no_hits_df = pd.read_csv(args.no_occurrences_csv)
        already_processed.update(no_hits_df["species"].dropna().unique())
        no_hits = no_hits_df.to_dict(orient="records")

    # Filter only species not already processed
    to_process = [sp for sp in species_list if sp not in already_processed]
    print(f"Species left to process: {len(to_process)}")

    # Parallel processing using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {
            executor.submit(get_closest_gbif, sp, args.ref_lat, args.ref_lon, args.radius, args.cache_dir, args.top_n, args.min_year): sp
            for sp in to_process
        }
        for i, future in enumerate(tqdm(as_completed(futures), total=len(futures), desc="Processing species", unit="species"), start=1):
            sp = futures[future]
            try:
                res, no_hit = future.result()
                if res:
                    results.extend(res)
                    tqdm.write(f"[{i}/{len(to_process)}] {sp}")
                elif no_hit:
                    no_hits.append({"species": no_hit})
                    tqdm.write(f"[{i}/{len(to_process)}] No occurrence: {sp}")
            except Exception as e:
                tqdm.write(f"[{i}/{len(to_process)}] Error with {sp}: {e}")
                no_hits.append({"species": sp})

            # Save progress each 10 species
            if i % 10 == 0 or i == len(to_process):
                pd.DataFrame(results).to_csv(args.output_csv, index=False)
                pd.DataFrame(no_hits) \
                    .drop_duplicates(subset=["species"]) \
                    .to_csv(args.no_occurrences_csv, index=False)
                tqdm.write(f"Progress saved after {i} species.")

    # Save everything again (including previous results)
    pd.DataFrame(results).to_csv(args.output_csv, index=False)
    print(f"\nResults saved to: {args.output_csv}")

    if no_hits:
        pd.DataFrame(no_hits).drop_duplicates(subset=["species"]).to_csv(args.no_occurrences_csv, index=False)
        print(f"Species without occurrences saved to: {args.no_occurrences_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find the closest GBIF occurrence(s) for each species.")
    parser.add_argument("--input_csv", required=True, help="Input CSV with species names.")
    parser.add_argument("--column", required=True, help="Column name with scientific names.")
    parser.add_argument("--output_csv", required=True, help="CSV file to save results.")
    parser.add_argument("--no_occurrences_csv", default="no_occurrences.csv", help="CSV for species with no occurrences.")
    parser.add_argument("--ref_lat", type=float, default=40.3397, help="Reference latitude.")
    parser.add_argument("--ref_lon", type=float, default=-7.6120, help="Reference longitude.")
    parser.add_argument("--radius", type=float, default=20, help="Search radius in km.")
    parser.add_argument("--top_n", type=int, default=1, choices=range(1, 11),
                        help="Number of closest occurrences to return per species (1-10).")
    parser.add_argument("--min_year", type=int, default=None, help="Minimum year for occurrences.")
    parser.add_argument("--cache_dir", default="cache", help="Cache directory.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads.")
    args = parser.parse_args()
    main(args)
