import pandas as pd
import argparse
from datetime import datetime
from dateutil import parser as dateparser

def load_data(blast_file, taxonomy_file, gbif_file):
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
    merged = pd.merge(blast_df, taxonomy_df, how='left', on=['Query_ID', 'Accession'])

    merged = pd.merge(merged, gbif_df, how='left', on='Species')
    return merged

def calculate_score(df, w_identity=1/3, w_distance=1/3, w_date=1/3):
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
    df['DaysSince'] = df['EventDate'].apply(lambda x: (today - x).days if pd.notnull(x) else pd.NA)

    complete_rows = df[df['Distance_km'].notnull() & df['DaysSince'].notnull()].copy()
    incomplete_rows = df[~(df['Distance_km'].notnull() & df['DaysSince'].notnull())].copy()

    complete_rows['identity_score'] = complete_rows['Percent Identity'] / 100
    complete_rows['distance_score'] = 1 - (complete_rows['Distance_km'] / complete_rows['Distance_km'].max())
    complete_rows['date_score'] = 1 - (complete_rows['DaysSince'] / complete_rows['DaysSince'].max())

    complete_rows['distance_score'] = complete_rows['distance_score'].fillna(0)
    complete_rows['date_score'] = complete_rows['date_score'].fillna(0)

    complete_rows['Score'] = (
        complete_rows['identity_score'] * w_identity +
        complete_rows['distance_score'] * w_distance +
        complete_rows['date_score'] * w_date
    )

    complete_rows = complete_rows.sort_values(by='Score', ascending=False)
    complete_rows['Score'] = complete_rows['Score'].round(3)

    return complete_rows, incomplete_rows

def save_output(scored_df, missing_df, scored_file, missing_file):
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
    parser.add_argument("--w_identity", type=float, default=1/3, help="Weight for identity score (default 1/3)")
    parser.add_argument("--w_distance", type=float, default=1/3, help="Weight for distance score (default 1/3)")
    parser.add_argument("--w_date", type=float, default=1/3, help="Weight for date score (default 1/3)")

    args = parser.parse_args()

    blast_df, taxonomy_df, gbif_df = load_data(args.blast, args.taxonomy, args.gbif)
    merged_df = preprocess_data(blast_df, taxonomy_df, gbif_df)
    scored_df, missing_df = calculate_score(
        merged_df,
        w_identity=args.w_identity,
        w_distance=args.w_distance,
        w_date=args.w_date
    )
    save_output(scored_df, missing_df, args.output, args.missing)

if __name__ == "__main__":
    main()
