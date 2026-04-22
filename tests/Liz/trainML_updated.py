import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report

# =========================
# 1. Load dataset
# =========================
excel_file = "runlog_results_5Kslots.xlsx"   # change if your filename is slightly different
df = pd.read_excel(excel_file)

# =========================
# 2. Create overall risk target
# =========================
# Overall Risk = 1 if either retransmission risk OR latency risk is high
df["OverallRisk"] = (
    (df["Retransmission Risk"] == 1) |
    (df["Latency Risk"] == 1)
).astype(int)

# =========================
# 3. Select features and target
# =========================
# Use only inputs available before/at decision time
X = df[["SNRdB", "NHARQProcesses", "Retransmission Rate (%)", "Modulation"]]
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
# 5. Preprocess + model
# =========================
# One-hot encode Modulation because it's categorical
preprocessor = ColumnTransformer(
    transformers=[
        ("cat", OneHotEncoder(drop="first", handle_unknown="ignore"), ["Modulation"]),
        ("num", "passthrough", ["SNRdB", "NHARQProcesses", "Retransmission Rate (%)"])
    ]
)

model = Pipeline([
    ("preprocessor", preprocessor),
    ("classifier", LogisticRegression(max_iter=1000))
])

# =========================
# 6. Train model
# =========================
model.fit(X_train, y_train)

# =========================
# 7. Evaluate model
# =========================
y_pred = model.predict(X_test)
y_prob = model.predict_proba(X_test)[:, 1]

print("=== MODEL PERFORMANCE ===")
print("Accuracy:", round(accuracy_score(y_test, y_pred), 4))
print("\nConfusion Matrix:")
print(confusion_matrix(y_test, y_pred))
print("\nClassification Report:")
print(classification_report(y_test, y_pred, digits=4))

# =========================
# 8. Create ML-based modulation decision table
# =========================
# For each (SNRdB, NHARQProcesses), score both 16QAM and 64QAM
# Then choose 64QAM if its predicted risk is acceptable; else choose 16QAM.
decision_rows = []

unique_conditions = (
    df[["SNRdB", "NHARQProcesses"]]
    .drop_duplicates()
    .sort_values(["SNRdB", "NHARQProcesses"])
)

# Risk threshold for allowing 64QAM
# You can tighten or loosen this later
risk_threshold = 0.50

for _, row in unique_conditions.iterrows():
    snr = row["SNRdB"]
    harq = row["NHARQProcesses"]

    avg_retx = df[
        (df["SNRdB"] == snr) &
        (df["NHARQProcesses"] == harq)
    ]["Retransmission Rate (%)"].mean()

    test_16 = pd.DataFrame({
        "SNRdB": [snr],
        "NHARQProcesses": [harq],
        "Retransmission Rate (%)": [avg_retx],
        "Modulation": ["16QAM"]
    })

    test_64 = pd.DataFrame({
        "SNRdB": [snr],
        "NHARQProcesses": [harq],
        "Retransmission Rate (%)": [avg_retx],
        "Modulation": ["64QAM"]
    })

    risk_prob_16 = model.predict_proba(test_16)[0, 1]
    risk_prob_64 = model.predict_proba(test_64)[0, 1]

    # Decision rule:
    # If 64QAM predicted risk is acceptable, use it for higher throughput.
    # Otherwise fall back to 16QAM.
    if risk_prob_64 <= risk_threshold:
        recommended_modulation = "64QAM"
    else:
        recommended_modulation = "16QAM"

    decision_rows.append({
        "SNRdB": snr,
        "NHARQProcesses": harq,
        "Predicted Risk Prob (16QAM)": round(risk_prob_16, 4),
        "Predicted Risk Prob (64QAM)": round(risk_prob_64, 4),
        "Recommended Modulation": recommended_modulation
    })

decision_df = pd.DataFrame(decision_rows)

print("\n=== ML DECISION TABLE ===")
print(decision_df)

# =========================
# 9. Compare ML choice to actual dataset performance
# =========================
# Merge ML recommendations back to your actual results so you can compare
# throughput, latency, BLER, etc. for the selected modulation.

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
# Simple baseline rule:
# Use 16QAM for SNR < 10, otherwise 64QAM
# You can change this threshold later if you want.

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
with pd.ExcelWriter("ml_results_output.xlsx", engine="openpyxl") as writer:
    decision_df.to_excel(writer, sheet_name="Decision_Table", index=False)
    ml_selected_df.to_excel(writer, sheet_name="ML_Selected_Rows", index=False)
    baseline_selected_df.to_excel(writer, sheet_name="Baseline_Selected_Rows", index=False)

print("\nSaved results to ml_results_output.xlsx")