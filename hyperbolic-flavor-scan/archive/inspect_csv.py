import pandas as pd
df = pd.read_csv('scan_results.csv')
print("Columns:", df.columns.tolist())
print("First 5 rows:")
print(df.head())
print("Data types:\n", df.dtypes)