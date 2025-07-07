import subprocess, numpy as np, json, re, os

def run_mesh_generation(fibers=25, size=10.0, vf=0.5):
    subprocess.run(["python3", "gen_mesh.py", "--fibers", str(fibers), "--size", str(size), "--vf", str(vf)], check=True)

def run_homogenization(load_type="strain", case_id=1):
    result = subprocess.run(["python3", "homogenization.py", "--load_type", load_type, "--case_id", str(case_id)], capture_output=True, text=True)
    return result.stdout

def extract_matrix(output):
    matrix_pattern = r'\[\[.*?\]\]'
    match = re.search(matrix_pattern, output, re.DOTALL)
    if match:
        matrix_str = match.group(0)
        matrix = np.array(json.loads(matrix_str))
        return matrix
    else:
        raise ValueError("Stiffness matrix not found in output.")

try:
    num_rve = int(input("Enter number of RVEs to analyze: "))
    num_fibers = int(input("Enter number of fibers: "))
    domain_size = float(input("Enter domain size: "))
    volume_fraction = float(input("Enter volume fraction (0-1): "))
    load_type = input("Enter 'stress' for uniaxial stress and 'strain' for uniaxial strain: ").strip().lower()

    if load_type not in ["strain", "stress"]:
        raise ValueError("Load type must be either 'strain' or 'stress'.")

except Exception as e:
    print(f"Input Error: {e}")
    exit()

os.makedirs('results', exist_ok=True)
output_file = f'results/output_{load_type}.txt'

with open(output_file, 'w') as f:
    modulus_matrices = []

    for i in range(num_rve):
        f.write(f"\n================ RVE {i+1} of {num_rve} ================\n")
        print(f"\n================ RVE {i+1} of {num_rve} ================")

        run_mesh_generation(num_fibers, domain_size, volume_fraction)

        output = run_homogenization(load_type=load_type,case_id=i+1)
        print(f"Homogenization Output:\n{output}")

        try:
            Chom = extract_matrix(output)
            f.write(f"\nHomogenized Stiffness Matrix (RVE {i+1}):\n{Chom}\n")
            print(f"Homogenized Stiffness Matrix (RVE {i+1}):\n{Chom}\n")
            modulus_matrices.append(Chom)
        except Exception as e:
            print(f"Error extracting matrix: {e}")

    if modulus_matrices:
        modulus_matrices = np.array(modulus_matrices)

        average_matrix = np.mean(modulus_matrices, axis=0)
        std_dev_matrix = np.std(modulus_matrices, axis=0)

        lmbda_avg = average_matrix[0, 1]
        mu_avg = average_matrix[2, 2]

        E_avg = mu_avg * (3 * lmbda_avg + 2 * mu_avg) / (lmbda_avg + mu_avg)
        nu_avg = lmbda_avg / (2 * (lmbda_avg + mu_avg))

        f.write("\n================ FINAL RESULTS ================\n")
        f.write(f"Average Homogenized Stiffness Matrix over {num_rve} RVEs:\n{average_matrix}\n")
        f.write(f"\nStandard Deviation of Homogenized Stiffness Matrix over {num_rve} RVEs:\n{std_dev_matrix}\n")
        f.write(f"\nAverage Homogenized Young's modulus: {E_avg:.3e} Pa\n")
        f.write(f"Average Homogenized Shear modulus: {mu_avg:.3e} Pa\n")
        f.write(f"Average Homogenized Poisson's ratio: {nu_avg:.3f}\n")
        f.write("===============================================\n")

        print("\n================ FINAL RESULTS ================\n")
        print(f"Average Homogenized Stiffness Matrix over {num_rve} RVEs:\n{average_matrix}\n")
        print(f"Standard Deviation of Homogenized Stiffness Matrix over {num_rve} RVEs:\n{std_dev_matrix}\n")
        print(f"Average Homogenized Young's modulus: {E_avg:.3e} Pa")
        print(f"Average Homogenized Shear modulus: {mu_avg:.3e} Pa")
        print(f"Average Homogenized Poisson's ratio: {nu_avg:.3f}")
        print("===============================================\n")
    else:
        print("No valid stiffness matrices were collected.")
