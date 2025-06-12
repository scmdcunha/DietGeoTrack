#!/usr/bin/env python3
"""
This script calculates a validation score for taxonomic identifications by integrating:
- BLAST identity percentages,
- Geographic distances from GBIF occurrence data,
- Temporal information (dates of occurrence records).

It processes three input files:
- A BLAST output file (TSV format) with sequence alignment results.
- A taxonomy file (CSV format) linking query and accession IDs to taxonomic info.
- A GBIF occurrences file (CSV format) with species occurrences, geographic distances, and dates.

The output consists of:
- A scored results file with a combined validation score.
- A file containing rows with missing distance or date information.

The weights for each score component (identity, distance, date) can be customized via CLI arguments.

Usage:
python calculate_score.py \
    --blast blast_results.tsv \
    --taxonomy taxonomy.csv \
    --gbif gbif_occurrences.csv \
    --output final_scores.csv \
    --missing missing_data.tsv
"""


import pandas as pd
import argparse
from datetime import datetime
from dateutil import parser as dateparser
import numpy as np
import math

def load_data(blast_file, taxonomy_file, gbif_file):
    """
    Load and format the input data files.

    Parameters:
    blast_file (str): Path to the BLAST TSV file.
    taxonomy_file (str): Path to the taxonomy CSV file.
    gbif_file (str): Path to the GBIF occurrences CSV file.

    Returns:
    tuple: DataFrames for blast, taxonomy, and gbif data.
    """
    blast_df = pd.read_csv(blast_file, sep='\t', header=None)
    blast_df.columns = [
        "Query_ID", "Accession", "Percent Identity", "Alignment Length", "Mismatches",
        "Gap Opens", "Query Start", "Query End", "Subject Start", "Subject End", "E-value", "Bit Score"
    ]

    taxonomy_df = pd.read_csv(taxonomy_file, sep=';')
    taxonomy_df.rename(columns={
        'Query ID': 'Query_ID',
        'Accession ID': 'Accession'
    }, inplace=True)

    gbif_df = pd.read_csv(gbif_file, sep=',', quotechar='"')
    gbif_df.rename(columns={
        'species': 'Species',
        'distance_km': 'Distance_km',
        'eventDate': 'EventDate'
    }, inplace=True)

    return blast_df, taxonomy_df, gbif_df

def preprocess_data(blast_df, taxonomy_df, gbif_df):
    """
    Merges BLAST, taxonomy, and GBIF occurrence data into a unified DataFrame.

    Parameters:
    - blast_df (pd.DataFrame): DataFrame containing BLAST results.
    - taxonomy_df (pd.DataFrame): DataFrame with taxonomy annotations.
    - gbif_df (pd.DataFrame): DataFrame with GBIF occurrences, distances, and dates.

    Returns:
    - pd.DataFrame: Merged DataFrame with all relevant fields.
    """

    merged = pd.merge(blast_df, taxonomy_df, how='left', on=['Query_ID', 'Accession'])
    merged = pd.merge(merged, gbif_df, how='left', on='Species')
    return merged

def half_life_to_lambda(half_life):
    if half_life <= 0:
        raise ValueError("Half-life must be positive")
    return math.log(2) / half_life

def calculate_score(df, w_identity=1 / 3, w_distance=1 / 3, w_date=1 / 3, half_life_distance=20, half_life_date=1825):
    """
    Calculate a combined validation score based on percent identity,
    geographic distance, and date of occurrence.

    Parameters:
    df (DataFrame): Merged DataFrame with necessary columns.
    w_identity (float): Weight for identity score component.
    w_distance (float): Weight for distance score component.
    w_date (float): Weight for date score component.

    Returns:
    tuple: (DataFrame with scores for complete data rows, DataFrame with incomplete data rows)
    """
    df['Percent Identity'] = pd.to_numeric(df['Percent Identity'], errors='coerce')
    df['Distance_km'] = pd.to_numeric(df['Distance_km'], errors='coerce')

    def parse_date_safe(date_str):
        try:
            dt = dateparser.parse(date_str)
            if dt is None:
                return pd.NaT
            if dt.year < 1800 or dt > datetime.now():
                return pd.NaT
            return dt
        except Exception:
            return pd.NaT

    df['EventDate'] = df['EventDate'].apply(parse_date_safe)
    today = datetime.now()
    df['DaysSince'] = df['EventDate'].apply(lambda x: (today - x).days if pd.notnull(x) else np.nan)
    complete_rows = df[df['Distance_km'].notnull() & df['DaysSince'].notnull()].copy()
    complete_rows['DaysSince'] = complete_rows['DaysSince'].astype(float)
    incomplete_rows = df[~(df['Distance_km'].notnull() & df['DaysSince'].notnull())].copy()

    # Calculate lambda values from half-life parameters
    lambda_distance = half_life_to_lambda(half_life_distance)
    lambda_date = half_life_to_lambda(half_life_date)

    # Identity Score: 0–1
    complete_rows['identity_score'] = complete_rows['Percent Identity'].apply(
        lambda x: x / 100 if x > 1 else x
    )

    # Distance Score with exponential decay
    complete_rows['distance_score'] = np.exp(-lambda_distance * complete_rows['Distance_km'])
    complete_rows['distance_score'] = complete_rows['distance_score'].clip(lower=0.001)

    # Date Score with exponential decay
    complete_rows['date_score'] = np.exp(-lambda_date * complete_rows['DaysSince'])
    complete_rows['date_score'] = complete_rows['date_score'].clip(lower=0.001)

    # Final weighted score
    complete_rows['Score'] = (
        complete_rows['identity_score'] * w_identity
        + complete_rows['distance_score'] * w_distance
        + complete_rows['date_score'] * w_date
    )

    complete_rows = complete_rows.sort_values(by='Score', ascending=False)
    complete_rows['Score'] = complete_rows['Score'].round(3)

    return complete_rows, incomplete_rows


def save_output(scored_df, missing_df, scored_file, missing_file):
    """
    Save the scored data and incomplete data to output files.

    Parameters:
    scored_df (DataFrame): DataFrame with computed scores.
    missing_df (DataFrame): DataFrame with missing distance or date data.
    scored_file (str): Output file path for scored data.
    missing_file (str): Output file path for missing data.
    """
    scored_df = scored_df.rename(columns={'Query_ID': 'Query ID', 'Accession': 'Accession ID'})

    cols = ['Query ID', 'Accession ID', 'Percent Identity',
            'Order', 'Family', 'Genus', 'Species',
            'Distance_km', 'EventDate', 'Score']
    scored_df[cols].to_csv(scored_file, sep='\t', index=False)
    print(f"Final scored output saved to: {scored_file}")

    if not missing_df.empty:
        missing_df = missing_df.rename(columns={'Query_ID': 'Query ID', 'Accession': 'Accession ID'})
        cols_missing = ['Query ID', 'Accession ID', 'Percent Identity',
                        'Order', 'Family', 'Genus', 'Species',
                        'Distance_km', 'EventDate']
        missing_df[cols_missing].to_csv(missing_file, sep='\t', index=False)
        print(f"Rows with missing distance or date data saved to: {missing_file}")
    else:
        print("No rows with missing distance or date data.")

def main():
    parser = argparse.ArgumentParser(description="Calculate a validation score from identity, geographic distance and date")
    parser.add_argument("--blast", required=True, help="BLAST TSV input file")
    parser.add_argument("--taxonomy", required=True, help="Taxonomy CSV input file")
    parser.add_argument("--gbif", required=True, help="GBIF occurrences CSV input file")
    parser.add_argument("--output", default="final_scores.csv", help="Output file for scored results")
    parser.add_argument("--missing", default="missing_data.csv", help="Output file for missing data rows")
    parser.add_argument("--w_identity", type=float, default=1 / 3, help="Weight for identity score (default 1/3)")
    parser.add_argument("--w_distance", type=float, default=1 / 3, help="Weight for distance score (default 1/3)")
    parser.add_argument("--w_date", type=float, default=1 / 3, help="Weight for date score (default 1/3)")
    parser.add_argument(
        "--half_life_distance",
        type=float,
        default=20,
        help="Distance at which score decays to half (km). Default 20 km."
    )
    parser.add_argument(
        "--half_life_date",
        type=float,
        default=1825,
        help="Time in days at which score decays to half. Default 1825 days (5 years)."
    )

    args = parser.parse_args()

    blast_df, taxonomy_df, gbif_df = load_data(args.blast, args.taxonomy, args.gbif)
    merged_df = preprocess_data(blast_df, taxonomy_df, gbif_df)
    scored_df, missing_df = calculate_score(
        merged_df,
        w_identity=args.w_identity,
        w_distance=args.w_distance,
        w_date=args.w_date,
        half_life_distance=args.half_life_distance,
        half_life_date=args.half_life_date
    )

    save_output(scored_df, missing_df, args.output, args.missing)

if __name__ == "__main__":
    main()
