import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Streamlit page config
st.set_page_config(page_title="Pharmacokinetics Model", layout="centered")

st.title("💊 Drug Exposure & Elimination")
st.markdown("Simulate drug concentration over time using basic pharmacokinetic models.")

# -------------------------------
# Sidebar Inputs
# -------------------------------
st.sidebar.header("📥 Input Parameters")

# Route of administration
route = st.sidebar.radio("Select Dosing Route", ["IV Bolus", "Oral Dosing"])

# Common inputs
dose = st.sidebar.number_input("Dose (mg)", min_value=0.0, value=500.0, step=50.0)
Vd = st.sidebar.number_input("Volume of Distribution Vd (L)", min_value=0.1, value=50.0)
k = st.sidebar.number_input("Elimination Rate Constant k (1/hr)", min_value=0.01, value=0.2)

# Oral-specific input
ka = None
if route == "Oral Dosing":
    ka = st.sidebar.number_input("Absorption Rate Constant ka (1/hr)", min_value=0.01, value=1.0)

# Multiple dose simulation
multiple_doses = st.sidebar.checkbox("Simulate Multiple Doses", value=False)
tau = 8
n_doses = 1

if multiple_doses:
    tau = st.sidebar.number_input("Dosing Interval τ (hr)", min_value=1, value=8)
    n_doses = st.sidebar.number_input("Number of Doses", min_value=1, value=3)

# Total time to simulate
total_hours = int(n_doses * tau + 24 if multiple_doses else 24)
time = np.linspace(0, total_hours, 1000)

# -------------------------------
# PK Models
# -------------------------------
def iv_bolus(C0, k, t):
    return C0 * np.exp(-k * t)

def oral_dose(D, Vd, ka, k, t):
    # Avoid divide-by-zero
    if ka == k:
        return np.zeros_like(t)
    return (D * ka) / (Vd * (ka - k)) * (np.exp(-k * t) - np.exp(-ka * t))

# -------------------------------
# Simulate Concentration-Time Data
# -------------------------------
concentration = np.zeros_like(time)

for i in range(n_doses):
    t_shift = time - i * tau
    t_shift[t_shift < 0] = 0  # Mask times before the dose
    if route == "IV Bolus":
        C0 = dose / Vd
        concentration += iv_bolus(C0, k, t_shift)
    elif route == "Oral Dosing":
        concentration += oral_dose(dose, Vd, ka, k, t_shift)

# Calculate AUC using trapezoidal rule
AUC = np.trapz(concentration, time)

# -------------------------------
# Output Results
# -------------------------------
st.subheader("📊 Results Summary")

if route == "IV Bolus":
    C0 = dose / Vd
    st.write(f"**Initial Concentration (C₀):** {C0:.2f} mg/L")
elif route == "Oral Dosing":
    Cmax = np.max(concentration)
    Tmax = time[np.argmax(concentration)]
    st.write(f"**Peak Concentration (Cmax):** {Cmax:.2f} mg/L at {Tmax:.1f} hours")

st.write(f"**Area Under Curve (AUC₀–{total_hours}):** {AUC:.2f} mg·hr/L")

# -------------------------------
# Plotting
# -------------------------------
st.subheader("📈 Drug Concentration vs Time")

fig, ax = plt.subplots()
ax.plot(time, concentration, label="Drug Concentration", color='purple')
ax.set_xlabel("Time (hr)")
ax.set_ylabel("Concentration (mg/L)")
ax.set_title("Pharmacokinetic Curve")
ax.grid(True)
ax.legend()
st.pyplot(fig)

# -------------------------------
# Table & Download
# -------------------------------
st.subheader("📋 Concentration Table")
df = pd.DataFrame({
    "Time (hr)": time,
    "Concentration (mg/L)": concentration
})
st.dataframe(df.iloc[::50])  # Show every 50th row

csv = df.to_csv(index=False).encode()
st.download_button("📥 Download CSV", csv, "concentration_data.csv", "text/csv")

# -------------------------------
# Equations & Definitions
# -------------------------------
st.subheader("📐 Mathematical Equations")

st.latex(r"C_{IV}(t) = C_0 \cdot e^{-k t}")
st.latex(r"C_{oral}(t) = \frac{D \cdot k_a}{V_d (k_a - k)} \cdot \left(e^{-k t} - e^{-k_a t}\right)")
st.latex(r"AUC \approx \int_0^T C(t)\, dt")

st.markdown("""
**Where:**
- `C₀` = Initial concentration (mg/L)
- `D` = Dose (mg)
- `Vd` = Volume of distribution (L)
- `k` = Elimination rate constant (1/hr)
- `ka` = Absorption rate constant (1/hr)
- `τ` = Dosing interval (hr)
""")
