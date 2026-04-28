import pandas as pd

# =========================
# File names
# =========================
baseline_file = "runlog_results_5Kslots.xlsx"
ml_file = "ml_results_output_extended.xlsx"
output_file = "combined_results.xlsx"

# =========================
# Load baseline data
# =========================
baseline = pd.read_excel(baseline_file)

# Remove any summary/median rows
baseline = baseline[pd.to_numeric(baseline["SNRdB"], errors="coerce").notna()].copy()

# Label baseline rows by modulation
baseline["Method"] = baseline["Modulation"]

# =========================
# Load ML-selected data
# =========================
ml = pd.read_excel(ml_file, sheet_name="ML_Selected_Rows")

# Remove any summary/median rows
ml = ml[pd.to_numeric(ml["SNRdB"], errors="coerce").notna()].copy()

# Label ML rows
ml["Method"] = "ML"

# =========================
# Combine datasets
# =========================
combined = pd.concat([baseline, ml], ignore_index=True)

# Optional: sort neatly
combined = combined.sort_values(
    by=["NHARQProcesses", "SNRdB", "Method"]
).reset_index(drop=True)

# =========================
# Save combined file
# =========================
combined.to_excel(output_file, index=False)

print(f"Combined results saved to: {output_file}")
print(f"Baseline rows: {len(baseline)}")
print(f"ML rows: {len(ml)}")
print(f"Total rows: {len(combined)}")