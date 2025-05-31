import pandas as pd
from scripts.calculate_score import calculate_score

def test_calculate_score_basic():
    data = {
        'Percent Identity': [100.0, 90.0],
        'Distance_km': [10.0, 100.0],
        'EventDate': ['2025-01-01', '2020-01-01']
    }
    df = pd.DataFrame(data)
    df = calculate_score(df)

    assert 'Score' in df.columns
    assert df['Score'].iloc[0] > df['Score'].iloc[1]
    assert df['identity_score'].iloc[0] == 1.0
    assert 0 <= df['distance_score'].iloc[0] <= 1
    assert 0 <= df['date_score'].iloc[0] <= 1
