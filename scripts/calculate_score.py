import pandas as pd
import argparse
from datetime import datetime
import numpy as np

def load_data(blast_file, taxonomy_file, gbif_file):
    blast_df = pd.read_csv(blast_file, sep='\t', header=None)
    blast_df.columns = [
        'Query_ID', 'Accession', 'Percent Identity', 'Alignment Length', 'Mismatches',
        'Gap Opens', 'Q. Start', 'Q. End', 'S. Start', 'S. End', 'E-value', 'Bit Score'
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
    # Merge BLAST + Taxonomy based on Query ID and Accession ID
    merged = pd.merge(blast_df, taxonomy_df, how='left', left_on=['Query_ID', 'Accession'], right_on=['Query_ID', 'Accession'])

    # Merge GBIF data based on Species
    merged = pd.merge(merged, gbif_df, how='left', on='Species')

    return merged

def calculate_score(df):
    # Normalization of values
    df['Percent Identity'] = df['Percent Identity'].astype(float)
    df['Distance'] = pd.to_numeric(df['Distance_km'], errors='coerce')

    # Date: Convert and calculate recency
    today = datetime.now()
    df['EventDate'] = pd.to_datetime(df['EventDate'], errors='coerce')
    df['DaysSince'] = (today - df['EventDate']).dt.days
    df['DaysSince'] = df['DaysSince'].fillna(df['DaysSince'].max())

    # Normalize to 0-1
    df['identity_score'] = df['Percent Identity'] / 100
    df['distance_score'] = 1 -df['Distance_km'] / df['Distance_km'].max()
    df['date_score'] = 1 -df['DaysSince'] / df['DaysSince'].max()

    # Fill NaN values with 0
    df['distance_score'] = df['distance_score'].fillna(0)
    df['date_score'] = df['date_score'].fillna(0)

    # Calculate final score
    df['Score'] = df['identity_score'] * 0.5 + df['distance_score'] * 0.3 + df['date_score'] * 0.2

    return df

def save_output(df, output_file):
    df.rename(columns={'Query_ID': 'Query ID', 'Accession': 'Accession ID'}, inplace=True)
    columns = ['Query ID', 'Accession ID', 'Percent Identity',
               'Order', 'Family', 'Genus', 'Species',
               'Distance_km', 'EventDate', 'Score']
    df[columns].to_csv(output_file, sep='\t', index=False)
    print(f"Final output saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Calculates a score for species validation based on identity, distance, and date.")
    parser.add_argument("--blast", required=True, help="TSV file with BLAST results.")
    parser.add_argument("--taxonomy", required=True, help="CSV file with taxonomy.")
    parser.add_argument("--gbif", required=True, help="CSV file with GBIF occurrences.")
    parser.add_argument("--output", default="final_scores.csv", help="Name of the output file.")

    args = parser.parse_args()

    blast_df, taxonomy_df, gbif_df = load_data(args.blast, args.taxonomy, args.gbif)
    merged_df = preprocess_data(blast_df, taxonomy_df, gbif_df)
    scored_df = calculate_score(merged_df)
    save_output(scored_df, args.output)

if __name__ == "__main__":
    main()
