
# Fiber-Reinforced Composite (FRC) Homogenization Tool

This project provides a complete Python-based pipeline for the **mesh generation**, **finite element homogenization**, and **computation of effective mechanical properties** of fiber-reinforced composite materials using periodic representative volume elements (RVEs).

The tool leverages:
- **Gmsh** for mesh generation
- **FEniCS** for finite element analysis
- **NumPy** for data processing

---

## Project Structure

- `FRC.py` – Main script that manages the entire homogenization workflow across multiple RVEs and computes statistical results.
- `gen_mesh.py` – Generates random fiber distributions in a square domain and creates the corresponding mesh.
- `homogenization.py` – Performs the finite element homogenization to calculate the stiffness matrix and effective material properties.

---

## Workflow

1. **Mesh Generation:**  
   Random fibers are distributed within a domain ensuring no overlap (with periodic boundary considerations). The mesh is generated using Gmsh and converted to XDMF format for FEniCS.

2. **Homogenization:**  
   The FEniCS-based solver applies either **uniaxial strain** or **uniaxial stress** to the RVE with periodic boundary conditions and computes the homogenized stiffness matrix.

3. **Results Processing:**  
   After analyzing multiple RVEs, the tool computes the average stiffness matrix, Young's modulus, shear modulus, and Poisson's ratio, along with standard deviations.

---

## Installation

### Requirements
- Python 3.8+
- Gmsh (installed and added to system path)
- FEniCS (recommended via Docker or virtual environment)
- Required Python packages:
  ```bash
  pip install numpy meshio gmsh
  ```

---

## Usage

1. Run the main workflow:
   ```bash
   python FRC.py
   ```
   You will be prompted to input:
   - Number of RVEs
   - Number of fibers
   - Domain size
   - Volume fraction
   - Load type (`strain` or `stress`)

2. The program will:
   - Generate random fiber distributions and mesh them.
   - Perform homogenization for each RVE.
   - Compute average and standard deviation of stiffness matrices.
   - Calculate homogenized Young’s modulus, shear modulus, and Poisson’s ratio.

3. Results will be saved in the `results` directory, and visualization files will be generated in `results/vtk`.

---

## Files and Directories

- `results/output_strain.txt` or `results/output_stress.txt` – Final computed results.
- `results/vtk` – VTK files for displacements, stresses, strains, equivalent strain, and von Mises stress for each case.

---

## Key Features

- Periodic boundary conditions
- Supports both stress-controlled and strain-controlled loading
- Handles multiple RVEs to obtain statistically meaningful homogenized properties
- Generates visualization files compatible with Paraview

---

## Example Input

```text
Enter number of RVEs to analyze: 3
Enter number of fibers: 25
Enter domain size: 10
Enter volume fraction (0-1): 0.5
Enter 'stress' for uniaxial stress and 'strain' for uniaxial strain: strain
```

---

## Notes
- Ensure that Gmsh is properly installed and callable from the terminal.
- FEniCS installation can be tricky on Windows; using Docker or WSL is recommended.
