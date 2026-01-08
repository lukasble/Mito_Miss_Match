import pandas as pd
import matplotlib.pyplot as plt

# ------------------------
# Paths and directories
# ------------------------
deleterious_tsv = (
    "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/"
    "Sarina/Results/9_deleterious_variants/mito_deleterious_counts.tsv"
)

weights_csv = (
    "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/"
    "breeds_weights.csv"
)

out_png = (
    "/proj/snic2022-6-164/MattC/dog_mitonuclear_conflict_project/"
    "Sarina/Results/9_deleterious_variants/"
    "deleterious_vs_weight_class.png"
)

# ------------------------
# Load deleterious counts
# ------------------------
del_df = pd.read_csv(deleterious_tsv, sep="\t")

# ------------------------
# Defining groups/categories
# ------------------------

# Coyotes to exclude
coyotes = {
    "CLATUS000001", "CLATUS000002", "CLATUS000003", "CLATUS000004"
}

# Wolves
wolves = {
"CLUPAZ000001","CLUPCN000001","CLUPCN000002","CLUPCN000003","CLUPCN000004",
"CLUPCN000005","CLUPCN000006","CLUPCN000007","CLUPCN000009","CLUPCN000010",
"CLUPEA000001","CLUPEU000002","CLUPGR000001","CLUPGR000002","CLUPGR000003",
"CLUPGR000004","CLUPGR000005","CLUPGR000006","CLUPGR000007","CLUPGR000008",
"CLUPGR000009","CLUPGR000010","CLUPGR000011","CLUPGR000012","CLUPIR000001",
"CLUPIR000002","CLUPIR000003","CLUPIR000004","CLUPIR000005","CLUPIR000006",
"CLUPKG000001","CLUPKZ000002","CLUPPT000001","CLUPPT000002","CLUPRU000001",
"CLUPRU000002","CLUPRU000003","CLUPRU000004","CLUPRU000005","CLUPRU000006",
"CLUPRU000007","CLUPRU000008","CLUPRU000009","CLUPRU000010","CLUPRU000011",
"CLUPRU000012","CLUPRU000013","CLUPRU000014","CLUPRU000019","CLUPRU000020",
"CLUPSE000001","CLUPSE000002","CLUPSE000003","CLUPSE000004","CLUPSE000005",
"CLUPSE000006","CLUPTJ000001","CLUPTJ000002","CLUPTJ000003","CLUPTJ000004",
"CLUPTJ000005","CLUPTJ000006","CLUPTJ000007"
}

# ------------------------
# Keep only breed dogs
# ------------------------
def is_breed_dog(sample_id):
    if sample_id in coyotes:
        return False
    if sample_id in wolves:
        return False
    if sample_id.startswith("VILL"):
        return False
    return True

del_df = del_df[del_df["dog_id"].apply(is_breed_dog)]

# ------------------------
# Loading weights
# ------------------------
w_df = pd.read_csv(weights_csv, sep=";")

# Keeping only dog ID and weight 
w_df = w_df.iloc[:, [1, 2]]
w_df.columns = ["dog_id", "weight"]

# Convert weight to float
w_df["weight"] = (
    w_df["weight"]
    .astype(str)
    .str.replace(",", ".", regex=False)
    .astype(float)
)

# ------------------------
# Merging
# ------------------------
df = pd.merge(del_df, w_df, on="dog_id", how="inner")

# ------------------------
# Weight classes
# ------------------------
bins = [0, 5, 10, 20, 40, 80]
labels = ["toy", "small", "medium", "large", "giant"]

df["weight_class"] = pd.cut(
    df["weight"],
    bins=bins,
    labels=labels,
    right=False
)

df = df.dropna(subset=["weight_class"])

# ------------------------
# Plotting the data
# ------------------------
data = [
    df.loc[df["weight_class"] == cls, "deleterious_count"]
    for cls in labels
]

plt.figure()
plt.boxplot(data, tick_labels=labels)
plt.xlabel("Weight class (kg)")
plt.ylabel("Deleterious variant count")
plt.title("Mitochondrial deleterious variants vs dog weight")

plt.tight_layout()
plt.savefig(out_png, dpi=300)
plt.close()

print(f"Plot saved to: {out_png}")
print("Counts per class:")
print(df["weight_class"].value_counts().sort_index())
