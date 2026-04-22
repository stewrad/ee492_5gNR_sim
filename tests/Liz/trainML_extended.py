import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report

# =========================
# 1. Load dataset
# =========================
excel_file = "runlog_results_5Kslots.xlsx"   # change if needed
df = pd.read_excel(excel_file)

# Remove summary rows like "Median"
df = df[pd.to_numeric(df["SNRdB"], errors="coerce").notna()].copy()

# =========================
# 2. Create overall risk target
# =========================
df["OverallRisk"] = (
    (df["Retransmission Risk"] == 1) |
    (df["Latency Risk"] == 1)
).astype(int)

# =========================
# 3. Select features and target
# =========================
feature_cols = [
    "SNRdB",
    "NHARQProcesses",
    "Average Transmissions per TB"
]

# Keep only rows with required values
df = df.dropna(subset=feature_cols + ["OverallRisk"]).copy()

X = df[feature_cols]
y = df["OverallRisk"]

# =========================
# 4. Train/test split
# =========================
X_train, X_test, y_train, y_test = train_test_split(
    X, y,
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
# 8. Create ML-based decision table
# =========================
decision_rows = []

unique_conditions = (
    df[["SNRdB", "NHARQProcesses"]]
    .drop_duplicates()
    .sort_values(["SNRdB", "NHARQProcesses"])
)

risk_threshold = 0.50

for _, row in unique_conditions.iterrows():
    snr = row["SNRdB"]
    harq = row["NHARQProcesses"]

    # Pull average HARQ-summary features for this (SNR, NHARQ) condition
    subset = df[
        (df["SNRdB"] == snr) &
        (df["NHARQProcesses"] == harq)
    ]

    avg_tx_per_tb = subset["Average Transmissions per TB"].mean()

    test_features = pd.DataFrame({
        "SNRdB": [snr],
        "NHARQProcesses": [harq],
        "Average Transmissions per TB": [avg_tx_per_tb]
    })

    risk_prob = model.predict_proba(test_features)[0, 1]

    if risk_prob <= risk_threshold:
        recommended_modulation = "64QAM"
    else:
        recommended_modulation = "16QAM"

    decision_rows.append({
        "SNRdB": snr,
        "NHARQProcesses": harq,
        "Predicted Risk Probability": round(risk_prob, 4),
        "Recommended Modulation": recommended_modulation
    })

decision_df = pd.DataFrame(decision_rows)

print("\n=== ML DECISION TABLE ===")
print(decision_df)

# =========================
# 9. Compare ML choice to actual dataset performance
# =========================
comparison_df = df.merge(
    decision_df[["SNRdB", "NHARQProcesses", "Recommended Modulation"]],
    on=["SNRdB", "NHARQProcesses"],
    how="left"
)

ml_selected_df = comparison_df[
    comparison_df["Modulation"] == comparison_df["Recommended Modulation"]
].copy()

print("\n=== ML-SELECTED PERFORMANCE SUMMARY ===")
print(
    ml_selected_df[
        [
            "Final BLER %",
            "Throughput Efficiency (%)",
            "HARQ Efficiency (%)",
            "Estimated Latency (ms)",
            "Average Throughput (Mbps)"
        ]
    ].mean(numeric_only=True)
)

# =========================
# 10. Baseline comparison
# =========================
baseline_df = df.copy()
baseline_df["Baseline Modulation"] = baseline_df["SNRdB"].apply(
    lambda x: "16QAM" if x < 10 else "64QAM"
)

baseline_selected_df = baseline_df[
    baseline_df["Modulation"] == baseline_df["Baseline Modulation"]
].copy()

print("\n=== BASELINE PERFORMANCE SUMMARY ===")
print(
    baseline_selected_df[
        [
            "Final BLER %",
            "Throughput Efficiency (%)",
            "HARQ Efficiency (%)",
            "Estimated Latency (ms)",
            "Average Throughput (Mbps)"
        ]
    ].mean(numeric_only=True)
)

# =========================
# 11. Save outputs to Excel
# =========================
with pd.ExcelWriter("ml_results_output_extended.xlsx", engine="openpyxl") as writer:
    decision_df.to_excel(writer, sheet_name="Decision_Table", index=False)
    ml_selected_df.to_excel(writer, sheet_name="ML_Selected_Rows", index=False)
    baseline_selected_df.to_excel(writer, sheet_name="Baseline_Selected_Rows", index=False)

print("\nSaved results to ml_results_output_extended.xlsx")