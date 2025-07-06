import numpy as np
import math
import tkinter as tk
from tkinter import messagebox, simpledialog
import matplotlib.pyplot as plt

# File path for the molecular data
MOLECULE_FILE_PATH = r'C:\computational chemistry\water.xyz'

# Load molecular data from the file
def load_molecule_data():
    try:
        data = np.loadtxt(MOLECULE_FILE_PATH, skiprows=2, usecols=[1, 2, 3])
        return data
    except Exception as e:
        messagebox.showerror("Error", f"Unable to load data: {e}")
        return None

# Function to calculate the distance between two atoms
def calculate_distance(coords1, coords2):
    return math.sqrt(sum((c1 - c2) ** 2 for c1, c2 in zip(coords1, coords2)))

# Function to calculate bond potential for a specific equilibrium bond length
def calculate_bond_potential(data):
    try:
        # Select atoms for bond potential calculation
        atom_pair = simpledialog.askstring("Bond Potential", "Enter atom pair (e.g., 1,2):")
        if atom_pair:
            indices = list(map(int, atom_pair.split(',')))
            if len(indices) == 2 and all(1 <= idx <= data.shape[0] for idx in indices):
                coords1 = data[indices[0] - 1]
                coords2 = data[indices[1] - 1]
                actual_distance = calculate_distance(coords1, coords2)

                # Ask for bond stiffness and equilibrium bond length
                bond_stiffness = simpledialog.askfloat("Bond Potential", "Enter bond stiffness constant (k) in N/m:")

                if bond_stiffness is None:
                    messagebox.showerror("Error", "Bond stiffness not provided.")
                    return

                # Validate stiffness constant
                if bond_stiffness <= 0:
                    messagebox.showerror("Error", "Bond stiffness must be a positive value.")
                    return

                # Calculate bond potential using harmonic potential formula
                equilibrium_distance = simpledialog.askfloat(
                    "Bond Potential",
                    "Enter equilibrium bond length (r0) in meters (0.1 to 0.9):"
                )

                if equilibrium_distance is None:
                    messagebox.showerror("Error", "Equilibrium bond length not provided.")
                    return

                # Validate equilibrium bond length range
                if not (0.1 <= equilibrium_distance <= 0.9):
                    messagebox.showerror("Error", "Equilibrium bond length (r0) must be between 0.1 and 0.9 meters.")
                    return

                # Display bond potential calculation
                messagebox.showinfo(
                    "Bond Potential", 
                    f"Equilibrium Bond Length (r0): {equilibrium_distance:.3f} meters\n"
                    f"Actual Distance: {actual_distance:.3f} meters\n"
                    f"Calculated Bond Potential: {0.5 * bond_stiffness * (actual_distance - equilibrium_distance) ** 2:.3e} Joules"
                )
            else:
                messagebox.showerror("Error", "Invalid atom pair. Please select valid indices.")
    except ValueError:
        messagebox.showerror("Error", "Invalid input for bond stiffness or equilibrium bond length.")
    except Exception as e:
        messagebox.showerror("Error", f"An unexpected error occurred: {e}")

# Function to plot bond potential vs. equilibrium bond length
def plot_bond_potential_vs_r0(data):
    try:
        # Select atoms for bond potential calculation
        atom_pair = "1,2"  # Use first two atoms by default for demonstration
        indices = list(map(int, atom_pair.split(',')))
        if len(indices) == 2 and all(1 <= idx <= data.shape[0] for idx in indices):
            coords1 = data[indices[0] - 1]
            coords2 = data[indices[1] - 1]
            actual_distance = calculate_distance(coords1, coords2)

            # Use default bond stiffness constant
            bond_stiffness = 100  # Example stiffness constant in N/m

            # Generate bond potentials for varying equilibrium bond lengths
            r0_values = np.linspace(0.1, 0.9, 100)  # Generate 100 points between 0.1 and 0.9 meters
            bond_potentials = [0.5 * bond_stiffness * (actual_distance - r0) ** 2 for r0 in r0_values]

            # Plot the results
            plt.figure(figsize=(8, 6))
            plt.plot(r0_values, bond_potentials, label=f"Stiffness Constant (k): {bond_stiffness} N/m", color="blue")
            plt.title("Bond Potential vs. Equilibrium Bond Length (r0)", fontsize=14)
            plt.xlabel("Equilibrium Bond Length (r0) [meters]", fontsize=12)
            plt.ylabel("Bond Potential [Joules]", fontsize=12)
            plt.grid(True, linestyle='--', alpha=0.7)
            plt.legend()
            plt.tight_layout()
            plt.show()
        else:
            messagebox.showerror("Error", "Invalid atom pair. Please select valid indices.")
    except Exception as e:
        messagebox.showerror("Error", f"An unexpected error occurred: {e}")

# Display molecule geometry
def show_molecule_geometry(data):
    try:
        geometry = "\n".join([f"Atom {i+1}: {coords.tolist()}" for i, coords in enumerate(data)])
        messagebox.showinfo("Molecule Geometry", f"Molecular Geometry:\n\n{geometry}")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")

# Calculate bond lengths for all atom pairs and display them
def show_bond_lengths(data):
    try:
        bond_lengths = []
        num_atoms = data.shape[0]
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                distance = calculate_distance(data[i], data[j])
                bond_lengths.append(f"Distance between Atom {i+1} and Atom {j+1}: {distance:.3f} meters")
        
        result = "\n".join(bond_lengths)
        messagebox.showinfo("Bond Lengths", f"Calculated Bond Lengths:\n\n{result}")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")

# Display help message
def show_help():
    messagebox.showinfo(
        "Help",
        "This application provides tools to analyze molecular data:\n\n"
        "1. Plot Bond Potential vs. r0: Automatically generates a graph.\n"
        "2. Calculate Bond Potential: Compute bond potential based on stiffness and distance.\n"
        "3. Show Molecule Geometry: Display all atom coordinates.\n"
        "4. Show Bond Lengths: Automatically calculate bond lengths for all atom pairs.\n"
        "5. Quit: Exit the application."
    )

# GUI and main logic
def main():
    data = load_molecule_data()
    if data is None:
        return

    root = tk.Tk()
    root.title("Molecule Analysis")
    root.configure(bg="#e6f7ff")  # Soft blue background

    # Title label
    title_label = tk.Label(root, text="Molecule Analysis", bg="#00509e", fg="white", font=("Arial", 18, "bold"), pady=10)
    title_label.pack(fill="x")

    # Buttons for functionalities
    button_style = {"bg": "#70c1b3", "fg": "black", "font": ("Arial", 12, "bold"), "width": 30}

    tk.Button(root, text="Calculate Bond Potential", command=lambda: calculate_bond_potential(data), **button_style).pack(pady=10)
    tk.Button(root, text="Plot Bond Potential vs. r0", command=lambda: plot_bond_potential_vs_r0(data), **button_style).pack(pady=10)
    tk.Button(root, text="Show Molecule Geometry", command=lambda: show_molecule_geometry(data), **button_style).pack(pady=10)
    tk.Button(root, text="Show Bond Lengths", command=lambda: show_bond_lengths(data), **button_style).pack(pady=10)
    tk.Button(root, text="Help", command=show_help, bg="#ffa41b", fg="black", font=("Arial", 12, "bold"), width=30).pack(pady=10)

    # Quit button
    tk.Button(root, text="Quit", command=root.destroy, bg="#ff4d4d", fg="white", font=("Arial", 12, "bold"), width=30).pack(pady=10)

    root.mainloop()

if __name__ == "__main__":
    main()
