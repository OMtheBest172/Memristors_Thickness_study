import pandas as pd
import numpy as np

# =====================================================
# Read COMSOL resistance data
# =====================================================

df = pd.read_csv("fine_90nm_sweep_resistance.csv", skiprows=5)

df.columns = [
    "Thickness",
    "Voltage",
    "Time",
    "Resistance"
]

df = df.apply(pd.to_numeric)

# =====================================================
# Read manually corrected switching times
# =====================================================

switch = pd.read_csv("SwitchingResults.csv")

results = []

# =====================================================
# Loop over all Thickness-Voltage combinations
# =====================================================

for thickness in sorted(df["Thickness"].unique()):

    temp = df[df["Thickness"] == thickness]

    for voltage in sorted(temp["Voltage"].unique()):

        data = temp[temp["Voltage"] == voltage].copy()

        data = data.sort_values("Time")

        t = data["Time"].values
        R = data["Resistance"].values

        # -------------------------------------------------
        # Read manually selected switching time
        # -------------------------------------------------

        row = switch[
            (switch["Thickness (nm)"] == thickness) &
            (switch["Voltage (V)"] == voltage)
        ]

        if len(row) == 0:
            print(f"Missing switching time for {thickness} nm, {voltage} V")
            continue

        switch_time = row["Switching Time (s)"].values[0]

        # -------------------------------------------------
        # Find all points before switching
        # -------------------------------------------------

        mask = t <= switch_time

        t_switch = t[mask]
        R_switch = R[mask]

        if len(t_switch) < 2:
            print(f"Not enough points before switching for {thickness} nm, {voltage} V")
            continue

        # -------------------------------------------------
        # Calculate Power
        # -------------------------------------------------

        power = voltage**2 / R_switch

        # -------------------------------------------------
        # Switching Energy
        # -------------------------------------------------

        switching_energy = np.trapz(power, t_switch)

        # -------------------------------------------------
        # Average Power
        # -------------------------------------------------

        avg_power = switching_energy / switch_time

        # -------------------------------------------------
        # LRS / HRS
        # -------------------------------------------------

        LRS = R[0]
        HRS = R[-1]

        ratio = HRS / LRS

        results.append([
            thickness,
            voltage,
            switch_time,
            switching_energy,
            avg_power,
            LRS,
            HRS,
            ratio
        ])

        print("="*70)
        print(f"Thickness          : {thickness:.0f} nm")
        print(f"Voltage            : {voltage:.1f} V")
        print(f"Switch Time        : {switch_time:.4f} s")
        print(f"Switch Energy      : {switching_energy:.6e} J")
        print(f"Average Power      : {avg_power:.6e} W")
        print(f"LRS                : {LRS:.3f} Ω")
        print(f"HRS                : {HRS:.3f} Ω")
        print(f"HRS/LRS            : {ratio:.3f}")
        print("="*70)

# =====================================================
# Save Results
# =====================================================

results = pd.DataFrame(
    results,
    columns=[
        "Thickness (nm)",
        "Voltage (V)",
        "Switching Time (s)",
        "Switching Energy (J)",
        "Average Power (W)",
        "LRS (Ω)",
        "HRS (Ω)",
        "HRS/LRS"
    ]
)

results.to_csv("SwitchingEnergyResults.csv", index=False)

print("\nSaved as 'SwitchingEnergyResults.csv'")
print(results)