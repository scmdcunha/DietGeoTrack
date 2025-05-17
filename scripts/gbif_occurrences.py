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

def get_closest_gbif(species, ref_lat, ref_lon, radius_km, cache_dir):
    """
    Searches the GBIF API for occurrences of a given species within a specified radius,
    and returns the closest one to the reference point.

    Parameters:
        species: Scientific name of the species.
        ref_lat, ref_lon: Reference coordinates (e.g. sampling site).
        radius_km: Search radius in kilometers.
        cache_dir: Path to directory where JSON response will be cached.

    Returns:
        A dictionary with metadata of the closest occurrence, or None if not found.
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

            r = requests.get(base_url, params=params)
            if r.status_code != 200:
                print(f"[Error] GBIF API call failed for {species} (status {r.status_code})")
                return None, species

            results = r.json().get("results", [])
            occurrences += results

            if offset + limit >= r.json().get("count", 0):
                break

            offset += limit
            time.sleep(0.2)  # Avoid overloading the API

        with open(cache_path, "w") as f:
            json.dump(occurrences, f)

    # Identify closest occurrence to reference point
    best, best_dist = None, radius_km * 1.1  # Allow some margin over max radius
    for occ in occurrences:
        lat, lon = occ.get("decimalLatitude"), occ.get("decimalLongitude")
        if lat is None or lon is None:
            continue
        dist = calculate_distance_gdal(ref_lat, ref_lon, lat, lon)
        if dist < best_dist:
            best_dist = dist
            best = occ

    if best:
        return {
            "species": species,
            "gbifKey": best["key"],
            "lat": best["decimalLatitude"],
            "lon": best["decimalLongitude"],
            "distance_km": round(best_dist, 2),
            "eventDate": best.get("eventDate", ""),
            "country": best.get("country", ""),
            "locality": best.get("locality", ""),
            "datasetKey": best.get("datasetKey", "")
        }, None
    else:
        return None, species

def main(args):
    """
    Main function that processes a list of species from a CSV file and finds the
    closest GBIF occurrence for each one.

    Results are saved to two CSVs: one with valid occurrences and one with missing data.
    """
    df = pd.read_csv(args.input_csv)
    if args.column not in df.columns:
        raise ValueError(f"Column '{args.column}' not found in the input CSV.")

    species_list = df[args.column].dropna().drop_duplicates().tolist()
    print(f"Unique species to process: {len(species_list)}")

    os.makedirs(args.cache_dir, exist_ok=True)
    results, no_hits = [], []

    # Parallel processing using ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {
            executor.submit(get_closest_gbif, sp, args.ref_lat, args.ref_lon, args.radius, args.cache_dir): sp
            for sp in species_list
        }
        for i, future in enumerate(as_completed(futures), 1):
            sp = futures[future]
            try:
                res, no_hit = future.result()
                if res:
                    results.append(res)
                    print(f"[{i}/{len(species_list)}] {sp}")
                elif no_hit:
                    no_hits.append({"species": no_hit})
                    print(f"[{i}/{len(species_list)}] No occurrence: {sp}")
            except Exception as e:
                print(f"[{i}/{len(species_list)}] Error with {sp}: {e}")
                no_hits.append({"species": sp})

    pd.DataFrame(results).to_csv(args.output_csv, index=False)
    print(f"\nResults saved to: {args.output_csv}")

    if no_hits:
        pd.DataFrame(no_hits).to_csv(args.no_occurrences_csv, index=False)
        print(f"Species without occurrences saved to: {args.no_occurrences_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find the closest GBIF occurrence for each species.")
    parser.add_argument("--input_csv", required=True, help="Path to input CSV file with scientific names.")
    parser.add_argument("--column", required=True, help="Column name containing scientific names.")
    parser.add_argument("--output_csv", required=True, help="Path to save the results CSV.")
    parser.add_argument("--no_occurrences_csv", default="no_occurrences.csv", help="Path to save species with no matches.")
    parser.add_argument("--ref_lat", type=float, default=40.3397, help="Reference latitude (e.g., sample site).")
    parser.add_argument("--ref_lon", type=float, default=-7.6120, help="Reference longitude.")
    parser.add_argument("--radius", type=float, default=20, help="Maximum search radius in kilometers.")
    parser.add_argument("--cache_dir", default="cache", help="Directory to store cached GBIF responses.")
    parser.add_argument("--threads", type=int, default=5, help="Number of threads to use.")
    args = parser.parse_args()
    main(args)
