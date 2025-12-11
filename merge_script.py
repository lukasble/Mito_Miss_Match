import pandas as pd

breeds = pd.read_csv("breeds_weights.csv", sep=";", decimal=",")
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 200)

load = pd.read_csv("genetic_load_per_sample.tsv", sep="\t")

print("Genetic load rows:", load.shape[0])
print(load.head(), "\n")

print("Breed table rows:", breeds.shape[0])
print(breeds.head(), "\n")

weight_col = breeds.columns[2]
breeds = breeds.rename(columns={weight_col: "weight_value"})
print("Columns in breed table:", breeds.columns.tolist(), "\n")

print("Merging on Sample (genetic load) vs sampleName (breeds)...")
merged = load.merge(breeds, left_on="Sample", right_on="sampleName", how="inner")

print("Merged table shape (rows, columns):", merged.shape)
print("First 10 rows of merged table:")
print(merged.head(10))

out_file = "genetic_load_with_breed.csv"
merged.to_csv(out_file, index=False)
print(f"\nMerged table written to: {out_file}")