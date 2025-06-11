import pandas as pd
from scripts.calculate_score import calculate_score

def test_calculate_score_basic():
    data = {
        'Percent Identity': [100.0, 90.0],
        'Distance_km': [10.0, 100.0],
        'EventDate': ['2025-01-01', '2020-01-01'],
        'Query_ID': ['q1', 'q2'],
        'Accession': ['acc1', 'acc2'],
        'Species': ['species1', 'species2'],
    }
    df = pd.DataFrame(data)

    scored_df, missing_df = calculate_score(df)

    assert not scored_df.empty
    assert 'Score' in scored_df.columns
    assert scored_df['Score'].iloc[0] > scored_df['Score'].iloc[1]
    assert scored_df['identity_score'].iloc[0] == 1.0
    assert 0 <= scored_df['distance_score'].iloc[0] <= 1
    assert 0 <= scored_df['date_score'].iloc[0] <= 1
    assert missing_df.empty  # porque todos os dados estão completos
