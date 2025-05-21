import pandas as pd
import folium
from folium.plugins import MarkerCluster
import argparse
from pathlib import Path

def plot_map(csv_path, output_html, ref_lat, ref_lon, zoom):
    df = pd.read_csv(csv_path)

    # Creates the map
    m = folium.Map(location=[ref_lat, ref_lon], zoom_start=zoom)

    # Point of reference
    folium.Marker(
        [ref_lat, ref_lon],
        popup="Ponto de Referência (Serra da Estrela)",
        icon=folium.Icon(color='red', icon='star')
    ).add_to(m)

    # Cluster of coordinates
    marker_cluster = MarkerCluster().add_to(m)

    # Markers of species
    for _, row in df.iterrows():
        species = row.get("species", "Espécie desconhecida")
        lat = row.get("lat")
        lon = row.get("lon")
        locality = row.get("locality", "Localidade desconhecida")
        distance = row.get("distance_km", "N/A")
        event_date = row.get("eventDate", "Data desconhecida")

        # Avoid invalid coordinates
        if pd.isna(lat) or pd.isna(lon):
            continue

        popup_text = f"<b>{species}</b><br><i>{locality}</i><br>{distance:.2f} km<br>{event_date}"
        folium.Marker(
            [lat, lon],
            popup=popup_text
        ).add_to(marker_cluster)

    # Create directory if it doesn't exist
    Path(output_html).parent.mkdir(parents=True, exist_ok=True)

    # Save
    m.save(output_html)
    print(f"Map saved at: {output_html}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generates a map with the closest occurrences from GBIF.")
    parser.add_argument("--input_csv", required=True, help="Path to the CSV with occurrences.")
    parser.add_argument("--output_html", default="occurrences_map.html", help="Path to save the HTML file.")
    parser.add_argument("--ref_lat", type=float, default=40.3397, help="Latitude of the reference point.")
    parser.add_argument("--ref_lon", type=float, default=-7.6120, help="Longitude of the reference point.")
    parser.add_argument("--zoom", type=int, default=8, help="Initial zoom level of the map.")

    args = parser.parse_args()
    plot_map(args.input_csv, args.output_html, args.ref_lat, args.ref_lon, args.zoom)
