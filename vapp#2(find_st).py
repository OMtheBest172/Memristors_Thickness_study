import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------
# Read COMSOL CSV
# ----------------------------------------------------

df = pd.read_csv("vapp.csv", skiprows=5)

df.columns = [
    "Thickness",
    "Voltage",
    "Time",
    "Resistance"
]

df = df.apply(pd.to_numeric)

results = []

# ----------------------------------------------------
# Loop through Thickness
# ----------------------------------------------------

for thickness in sorted(df["Thickness"].unique()):

    temp = df[df["Thickness"] == thickness]

    for voltage in sorted(temp["Voltage"].unique()):

        data = temp[temp["Voltage"] == voltage].copy()

        data = data.sort_values("Time")

        t = data["Time"].values
        R = data["Resistance"].values

        # --------------------------------------------
        # First derivative
        # --------------------------------------------

        dRdt = np.gradient(R, t)

        idx = np.argmax(dRdt)

        switch_time = t[idx]
        switch_R = R[idx]

        # --------------------------------------------
        # HRS/LRS
        # --------------------------------------------

        LRS = R[0]
        HRS = R[-1]
        ratio = HRS / LRS

        # --------------------------------------------
        # Switching Energy
        # E = ∫ V²/R dt
        # --------------------------------------------

        t_switch = t[:idx+1]
        R_switch = R[:idx+1]

        power = voltage**2 / R_switch       # Watts

        switching_energy = np.trapz(
            power,
            t_switch
        )                                   # Joules

        avg_power = switching_energy / switch_time

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

        # ------------------------------------------------
        # Plot
        # ------------------------------------------------

        fig, ax = plt.subplots(1, 2, figsize=(14,5))

        # ---------------- Resistance ----------------

        ax[0].plot(
            t,
            R,
            lw=2,
            color="royalblue"
        )

        ax[0].scatter(
            switch_time,
            switch_R,
            color="red",
            s=80,
            zorder=5,
            label="Switch Point"
        )

        ax[0].axvline(
            switch_time,
            color="red",
            linestyle="--"
        )

        ax[0].set_title(
            f"Resistance\n{int(thickness)} nm, {voltage:.1f} V"
        )

        ax[0].set_xlabel("Time (s)")
        ax[0].set_ylabel("Resistance (Ω)")
        ax[0].grid(True)
        ax[0].legend()

        # ---------------- Derivative ----------------

        ax[1].plot(
            t,
            dRdt,
            lw=2,
            color="darkgreen"
        )

        ax[1].scatter(
            switch_time,
            dRdt[idx],
            color="red",
            s=80,
            zorder=5,
            label="Maximum dR/dt"
        )

        ax[1].axvline(
            switch_time,
            color="red",
            linestyle="--"
        )

        ax[1].set_title("First Derivative")
        ax[1].set_xlabel("Time (s)")
        ax[1].set_ylabel("dR/dt (Ω/s)")
        ax[1].grid(True)
        ax[1].legend()

        plt.suptitle(
            f"{int(thickness)} nm   |   {voltage:.1f} V\n"
            f"Switch Time = {switch_time:.4f} s    "
            f"HRS/LRS = {ratio:.2f}\n"
            f"Switching Energy = {switching_energy:.3e} J    "
            f"Average Power = {avg_power:.3e} W",
            fontsize=12
        )

        plt.tight_layout()

        plt.show()

        # ------------------------------------------------
        # Print
        # ------------------------------------------------

        print("="*70)
        print(f"Thickness           : {thickness:.0f} nm")
        print(f"Voltage             : {voltage:.1f} V")
        print(f"Switch Time         : {switch_time:.5f} s")
        print(f"LRS                 : {LRS:.3f} Ω")
        print(f"HRS                 : {HRS:.3f} Ω")
        print(f"HRS/LRS Ratio       : {ratio:.3f}")
        print(f"Switching Energy    : {switching_energy:.6e} J")
        print(f"Average Power       : {avg_power:.6e} W")
        print("="*70)

# ----------------------------------------------------
# Save Results
# ----------------------------------------------------

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

results.to_csv("SwitchingResults.csv", index=False)

print("\n\nFinal Results:\n")
print(results)