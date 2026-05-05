import pandas as pd

original_file = "runlog_results_ML.xlsx"
ml_file = "ml_results_output_extended.xlsx"
output_file = "original_rows_with_ml_risk_probability.xlsx"

original_df = pd.read_excel(original_file)
decision_df = pd.read_excel(ml_file, sheet_name="Decision_Table")
recommended_df = pd.read_excel(ml_file, sheet_name="Recommended_Modulation")

# remove summary rows
original_df = original_df[pd.to_numeric(original_df["SNRdB"], errors="coerce").notna()].copy()

# make merge columns consistent
for df in [original_df, decision_df, recommended_df]:
    df["SNRdB"] = pd.to_numeric(df["SNRdB"], errors="coerce")
    df["NHARQProcesses"] = pd.to_numeric(df["NHARQProcesses"], errors="coerce")

for df in [original_df, decision_df]:
    df["Modulation"] = df["Modulation"].astype(str).str.strip()

recommended_df["Recommended Modulation"] = recommended_df["Recommended Modulation"].astype(str).str.strip()

# bring Decision_Table predicted risk onto every original row
risk_lookup = decision_df[
    [
        "SNRdB",
        "NHARQProcesses",
        "Modulation",
        "Predicted Risk Probability",
        "Recommended Modulation",
        "Is Recommended",
    ]
].copy()

merged_df = original_df.merge(
    risk_lookup,
    on=["SNRdB", "NHARQProcesses", "Modulation"],
    how="left"
)

# add recommended risk probability too
merged_df = merged_df.merge(
    recommended_df[
        [
            "SNRdB",
            "NHARQProcesses",
            "Recommended Risk Probability",
            "Recommendation Reason",
        ]
    ],
    on=["SNRdB", "NHARQProcesses"],
    how="left"
)

# optional baseline selection
merged_df["Baseline Modulation"] = merged_df["SNRdB"].apply(
    lambda x: "16QAM" if x < 10 else "64QAM"
)

merged_df["Is Baseline Selected"] = (
    merged_df["Modulation"] == merged_df["Baseline Modulation"]
).astype(int)

with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    merged_df.to_excel(writer, sheet_name="All_Original_Rows", index=False)

    baseline_selected_df = merged_df[merged_df["Is Baseline Selected"] == 1].copy()
    baseline_selected_df.to_excel(writer, sheet_name="Baseline_Selected_Rows", index=False)

    ml_selected_df = merged_df[merged_df["Is Recommended"] == 1].copy()
    ml_selected_df.to_excel(writer, sheet_name="ML_Selected_Rows", index=False)

print(f"Saved results to {output_file}")