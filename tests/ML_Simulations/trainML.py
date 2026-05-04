import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report

# =========================
# 1. Load dataset
# =========================
excel_file = "runlog_results_ML.xlsx"
df = pd.read_excel(excel_file)

# Remove summary rows like "Median"
df = df[pd.to_numeric(df["SNRdB"], errors="coerce").notna()].copy()
df["SNRdB"] = pd.to_numeric(df["SNRdB"], errors="coerce")
df["NHARQProcesses"] = pd.to_numeric(df["NHARQProcesses"], errors="coerce")

# =========================
# 2. Create overall risk target
# =========================
# Assumes Retransmission Risk and Latency Risk are already 0-to-1 columns
df["OverallRisk"] = (
    (df["Retransmission Risk"] >= 0.75) |
    (df["Latency Risk"] >= 0.75)
).astype(int)

# =========================
# 3. Select features and target
# =========================
feature_cols = [
    "SNRdB",
    "NHARQProcesses",
    "Average Transmissions per TB",
    "Retransmission Risk",
    "Latency Risk",
]

required_cols = feature_cols + [
    "OverallRisk",
    "Modulation",
    "Throughput Efficiency (%)",
    "HARQ Efficiency (%)",
    "Estimated Latency (ms)",
    "Average Throughput (Mbps)",
]

missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise KeyError(f"Missing required columns in Excel file: {missing}")

# Keep only rows with required values
df = df.dropna(subset=required_cols).copy()

X = df[feature_cols]
y = df["OverallRisk"]

# Safety check
if y.nunique() < 2:
    raise ValueError(
        "OverallRisk only has one class. Check your Retransmission Risk / "
        "Latency Risk thresholds or make sure your dataset includes both low-risk "
        "and high-risk modulation rows."
    )

# =========================
# 4. Train/test split
# =========================
X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.30,
    random_state=42,
    stratify=y
)

# =========================
# 5. Model
# =========================
model = LogisticRegression(max_iter=1000)

# =========================
# 6. Train model
# =========================
model.fit(X_train, y_train)

# =========================
# 7. Evaluate model
# =========================
y_pred = model.predict(X_test)

print("=== MODEL PERFORMANCE ===")
print("Accuracy:", round(accuracy_score(y_test, y_pred), 4))
print("\nConfusion Matrix:")
print(confusion_matrix(y_test, y_pred))
print("\nClassification Report:")
print(classification_report(y_test, y_pred, digits=4))

# =========================
# 8. Create candidate decision table
# =========================
decision_rows = []

unique_conditions = (
    df[["SNRdB", "NHARQProcesses", "Modulation"]]
    .drop_duplicates()
    .sort_values(["SNRdB", "NHARQProcesses", "Modulation"])
)

risk_threshold = 0.50

for _, row in unique_conditions.iterrows():
    snr = row["SNRdB"]
    harq = row["NHARQProcesses"]
    mod = row["Modulation"]

    subset = df[
        (df["SNRdB"] == snr) &
        (df["NHARQProcesses"] == harq) &
        (df["Modulation"] == mod)
    ]

    avg_tx_per_tb = subset["Average Transmissions per TB"].mean()
    retx_risk = subset["Retransmission Risk"].mean()
    latency_risk = subset["Latency Risk"].mean()

    test_features = pd.DataFrame({
        "SNRdB": [snr],
        "NHARQProcesses": [harq],
        "Average Transmissions per TB": [avg_tx_per_tb],
        "Retransmission Risk": [retx_risk],
        "Latency Risk": [latency_risk],
    })

    risk_prob = model.predict_proba(test_features[feature_cols])[0, 1]

    decision_rows.append({
        "SNRdB": snr,
        "NHARQProcesses": harq,
        "Modulation Tested": mod,
        "Average Transmissions per TB": round(avg_tx_per_tb, 4),
        "Retransmission Risk": round(retx_risk, 4),
        "Latency Risk": round(latency_risk, 4),
        "Predicted Risk Probability": round(risk_prob, 4),
        "Risk Class": int(risk_prob >= risk_threshold),
        "Mean Throughput Efficiency (%)": subset["Throughput Efficiency (%)"].mean(),
        "Mean HARQ Efficiency (%)": subset["HARQ Efficiency (%)"].mean(),
        "Mean Estimated Latency (ms)": subset["Estimated Latency (ms)"].mean(),
        "Mean Average Throughput (Mbps)": subset["Average Throughput (Mbps)"].mean(),
    })

decision_df = pd.DataFrame(decision_rows)

# =========================
# 9. Select one recommended modulation per SNR/HARQ condition
# =========================
recommendation_rows = []

for (snr, harq), group in decision_df.groupby(["SNRdB", "NHARQProcesses"]):
    safe_group = group[group["Predicted Risk Probability"] <= risk_threshold].copy()

    if not safe_group.empty:
        # Choose the highest-throughput modulation among safe choices
        best_row = safe_group.sort_values(
            ["Mean Average Throughput (Mbps)", "Mean Throughput Efficiency (%)"],
            ascending=False
        ).iloc[0]
        reason = "Highest throughput below risk threshold"
    else:
        # If all modulations are risky, choose the least risky one
        best_row = group.sort_values(
            ["Predicted Risk Probability", "Mean Estimated Latency (ms)"],
            ascending=True
        ).iloc[0]
        reason = "Lowest risk because all candidates exceeded threshold"

    recommendation_rows.append({
        "SNRdB": snr,
        "NHARQProcesses": harq,
        "Recommended Modulation": best_row["Modulation Tested"],
        "Recommended Risk Probability": best_row["Predicted Risk Probability"],
        "Recommendation Reason": reason,
    })

recommendation_df = pd.DataFrame(recommendation_rows)

# Mark the recommended row inside the full decision table
decision_df = decision_df.merge(
    recommendation_df[["SNRdB", "NHARQProcesses", "Recommended Modulation"]],
    on=["SNRdB", "NHARQProcesses"],
    how="left"
)

decision_df["Is Recommended"] = (
    decision_df["Modulation Tested"] == decision_df["Recommended Modulation"]
).astype(int)

print("\n=== CANDIDATE DECISION TABLE ===")
print(decision_df)

print("\n=== RECOMMENDED MODULATION TABLE ===")
print(recommendation_df)

# =========================
# 10. Compare ML choice to actual dataset performance
# =========================
comparison_df = df.merge(
    recommendation_df[["SNRdB", "NHARQProcesses", "Recommended Modulation"]],
    on=["SNRdB", "NHARQProcesses"],
    how="left"
)

ml_selected_df = comparison_df[
    comparison_df["Modulation"] == comparison_df["Recommended Modulation"]
].copy()

# Your parser uses "Final BLER % (post-HARQ)", not "Final BLER %"
# This also falls back to Block Error Rate if needed.
if "Final BLER % (post-HARQ)" in df.columns:
    bler_col = "Final BLER % (post-HARQ)"
elif "Block Error Rate (BLER %)" in df.columns:
    bler_col = "Block Error Rate (BLER %)"
else:
    bler_col = None

metric_cols = [
    "Throughput Efficiency (%)",
    "HARQ Efficiency (%)",
    "Estimated Latency (ms)",
    "Average Throughput (Mbps)",
]

if bler_col is not None:
    metric_cols = [bler_col] + metric_cols

print("\n=== ML-SELECTED PERFORMANCE SUMMARY ===")
print(ml_selected_df[metric_cols].mean(numeric_only=True))

# =========================
# 11. Baseline comparison
# =========================
baseline_df = df.copy()

# Simple fixed-rule baseline:
# below 10 dB -> 16QAM, 10 dB and above -> 64QAM
baseline_df["Baseline Modulation"] = baseline_df["SNRdB"].apply(
    lambda x: "16QAM" if x < 10 else "64QAM"
)

baseline_selected_df = baseline_df[
    baseline_df["Modulation"] == baseline_df["Baseline Modulation"]
].copy()

print("\n=== BASELINE PERFORMANCE SUMMARY ===")
print(baseline_selected_df[metric_cols].mean(numeric_only=True))

# =========================
# 12. Save outputs to Excel
# =========================
with pd.ExcelWriter("ml_results_output_extended.xlsx", engine="openpyxl") as writer:
    decision_df.to_excel(writer, sheet_name="Candidate_Decision_Table", index=False)
    recommendation_df.to_excel(writer, sheet_name="Recommended_Modulation", index=False)
    ml_selected_df.to_excel(writer, sheet_name="ML_Selected_Rows", index=False)
    baseline_selected_df.to_excel(writer, sheet_name="Baseline_Selected_Rows", index=False)

print("\nSaved results to ml_results_output_extended.xlsx")