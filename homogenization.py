from fenics import *
import numpy as np, json, os

def homogenization(load_type='strain',case_id=1):
    mesh = Mesh()
    with XDMFFile("mesh/rve.xdmf") as infile:
        infile.read(mesh)

    mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim())
    with XDMFFile("mesh/rve.xdmf") as infile:
        infile.read(mvc, "name_to_read")

    subdomains = cpp.mesh.MeshFunctionSizet(mesh, mvc)

    class PeriodicBoundary(SubDomain):
        def __init__(self, tol=DOLFIN_EPS):
            SubDomain.__init__(self)
            self.tol = tol

        def inside(self, x, on_boundary):
            return (near(x[0], 0.0, self.tol) or near(x[1], 0.0, self.tol)) and on_boundary

        def map(self, x, y):
            if near(x[0], 10.0, self.tol) and near(x[1], 10.0, self.tol):
                y[0] = x[0] - 10.0
                y[1] = x[1] - 10.0
            elif near(x[0], 10.0, self.tol):
                y[0] = x[0] - 10.0
                y[1] = x[1]
            else:
                y[0] = x[0]
                y[1] = x[1] - 10.0

    if load_type == 'strain':
        Ve = VectorElement("CG", mesh.ufl_cell(), 2, dim=2)
        Re = VectorElement("R", mesh.ufl_cell(), 0)
        W = FunctionSpace(mesh, MixedElement([Ve, Re]), constrained_domain=PeriodicBoundary(tol=1e-10))
    else:
        V = VectorFunctionSpace(mesh, "CG", 2, constrained_domain=PeriodicBoundary())
        R = VectorElement("R", mesh.ufl_cell(), 0)
        W = FunctionSpace(mesh, MixedElement([V.ufl_element(), R]))

    v_, lamb_ = TestFunctions(W)
    dv, dlamb = TrialFunctions(W)
    w = Function(W)

    E_matrix, nu_matrix = 4.0E9, 0.34
    E_fiber, nu_fiber = 23.0E9, 0.22
    material_parameters = [(E_matrix, nu_matrix), (E_fiber, nu_fiber)]
    nphases = len(material_parameters)

    def eps(v): return sym(grad(v))

    def sigma(v, i, Eps=None):
        E, nu = material_parameters[i]
        lmbda = E * nu / (1 + nu) / (1 - 2 * nu)
        mu = E / 2 / (1 + nu)
        if load_type == 'strain':
            return lmbda * tr(eps(v) + Eps) * Identity(2) + 2 * mu * (eps(v) + Eps)
        else:
            return lmbda * tr(eps(v)) * Identity(2) + 2 * mu * eps(v)

    dx = Measure('dx', domain=mesh, subdomain_data=subdomains)
    vol = assemble(Constant(1.0) * dx)

    if load_type == 'strain':
        Eps = Constant(((0.0, 0.0), (0.0, 0.0)))
        F = sum([inner(sigma(dv, i, Eps), eps(v_)) * dx(i + 1) for i in range(nphases)])
        a, L = lhs(F), rhs(F)
        a += dot(lamb_, dv) * dx + dot(dlamb, v_) * dx

        def macro_input(i):
            Eps_Voigt = np.zeros((3,))
            Eps_Voigt[i] = 1.0
            return np.array([[Eps_Voigt[0], Eps_Voigt[2] / 2.],
                             [Eps_Voigt[2] / 2., Eps_Voigt[1]]])

        def stress2Voigt(s):
            return as_vector([s[0, 0], s[1, 1], s[0, 1]])

    else:
        Sig = Constant(((0.0, 0.0), (0.0, 0.0)))
        F = sum([inner(sigma(dv, i), eps(v_)) * dx(i + 1) for i in range(nphases)]) - inner(Sig, eps(v_)) * dx
        a, L = lhs(F), rhs(F)
        a += dot(lamb_, dv) * dx + dot(dlamb, v_) * dx

        def macro_input(i):
            Sig_voigt = np.zeros((3,))
            Sig_voigt[i] = 1.0
            return np.array([[Sig_voigt[0], Sig_voigt[2]],
                             [Sig_voigt[2], Sig_voigt[1]]])

        def strain2Voigt(e):
            return as_vector([e[0, 0], e[1, 1], 2 * e[0, 1]])

    C = np.zeros((3, 3))
    S = np.zeros((3, 3))

    cases = ["Exx", "Eyy", "Exy"] if load_type == 'strain' else ["Sxx", "Syy", "Sxy"]

    output_base = f"results/vtk/uniaxial_{load_type}/RVE{case_id}"
    os.makedirs(f"{output_base}", exist_ok=True)
    os.makedirs(f"{output_base}/disp", exist_ok=True)
    os.makedirs(f"{output_base}/strain", exist_ok=True)
    os.makedirs(f"{output_base}/stress", exist_ok=True)
    os.makedirs(f"{output_base}/eq_strain", exist_ok=True)
    os.makedirs(f"{output_base}/von_mises", exist_ok=True)

    for j, case in enumerate(cases):
        print(f"Solving {case} case...")
        if load_type == 'strain':
            Eps.assign(Constant(macro_input(j)))
        else:
            Sig.assign(Constant(macro_input(j)))

        solve(a == L, w)
        v, lamb = w.split()

        v.rename("displacement", "")
        File(f"{output_base}/disp/displacement_{case}.pvd") << v

        strain = project(eps(v), TensorFunctionSpace(mesh, "CG", 1))
        strain.rename("strain", "")
        File(f"{output_base}/strain/strain_{case}.pvd") << strain

        TF = TensorFunctionSpace(mesh, "CG", 1)
        stress = Function(TF)
        for i in range(nphases):
            if load_type == 'strain':
                stress_phase = project(sigma(v, i, Eps), TF, solver_type='cg')
            else:
                stress_phase = project(sigma(v, i), TF, solver_type='cg')

            stress.vector().axpy(1.0, stress_phase.vector())

        stress.rename("stress", "")
        File(f"{output_base}/stress/stress_{case}.pvd") << stress

        strain_eq = project(sqrt(2 / 3 * inner(dev(strain), dev(strain))), FunctionSpace(mesh, "CG", 1))
        strain_eq.rename("equivalent_strain", "")
        File(f"{output_base}/eq_strain/equivalent_strain_{case}.pvd") << strain_eq

        stress_eq = project(sqrt(3 / 2 * inner(dev(stress), dev(stress))), FunctionSpace(mesh, "CG", 1))
        stress_eq.rename("von_mises", "")
        File(f"{output_base}/von_mises/von_mises_{case}.pvd") << stress_eq

        if load_type == 'strain':
            Sigma = np.zeros((3,))
            for k in range(3):
                Sigma[k] = assemble(sum([stress2Voigt(sigma(v, i, Eps))[k] * dx(i + 1) for i in range(nphases)])) / vol
            C[j, :] = Sigma
        else:
            Eps_macro = np.zeros((3,))
            for k in range(3):
                eps_k = sum([strain2Voigt(eps(v))[k] * dx(i + 1) for i in range(nphases)])
                Eps_macro[k] = assemble(eps_k) / vol
            S[j, :] = Eps_macro

    if load_type == 'strain':
        threshold = 1e-1 * np.max(np.abs(C))
        C[np.abs(C) < threshold] = 0
    else:
        C = np.linalg.inv(S)
        threshold = 1e-1 * np.max(np.abs(C))
        C[np.abs(C) < threshold] = 0

    print("Homogenized stiffness matrix:")
    print(json.dumps(C.tolist()))

    lmbda_hom = C[0, 1]
    mu_hom = C[2, 2]
    E_hom = mu_hom * (3 * lmbda_hom + 2 * mu_hom) / (lmbda_hom + mu_hom)
    nu_hom = lmbda_hom / (2 * (lmbda_hom + mu_hom))

    print(f"Homogenized Young's modulus: {E_hom:.3e} Pa")
    print(f"Homogenized Shear modulus: {mu_hom:.3e} Pa")
    print(f"Homogenized Poisson's ratio: {nu_hom:.3f}")
    return C

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--load_type", type=str, default="strain", choices=["strain", "stress"])
    parser.add_argument("--case_id", type=int, default=1)
    args = parser.parse_args()

    homogenization(load_type=args.load_type, case_id=args.case_id) 
