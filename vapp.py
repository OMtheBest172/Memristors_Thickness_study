import pandas as pd
import matplotlib.pyplot as plt

# Read COMSOL CSV (skip metadata)
df = pd.read_csv(
    "ultra_fine_86nm_sweep_resis.csv",
    skiprows=5,
    names=["Thickness", "Voltage", "Time", "Resistance"]
)

# Make sure numeric
df = df.apply(pd.to_numeric)

# Unique thicknesses
thicknesses = sorted(df["Thickness"].unique())

for t in thicknesses:

    plt.figure(figsize=(9,6))

    temp = df[df["Thickness"] == t]

    for v in sorted(temp["Voltage"].unique()):
        d = temp[temp["Voltage"] == v]
        plt.plot(
            d["Time"],
            d["Resistance"],
            label=f"{v:.1f} V"
        )

    plt.title(f"Resistance vs Time ({int(t)} nm WO3)")
    plt.xlabel("Time (s)")
    plt.ylabel("Resistance (Ω)")
    plt.grid(True)
    plt.legend(ncol=2)
    plt.tight_layout()
    plt.show()