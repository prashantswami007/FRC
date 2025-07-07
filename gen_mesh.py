import gmsh, random, math, meshio

def gen_text_file(NUM_FIBERS_TARGET=25, DOMAIN_SIZE=10.0, VOLUME_FRACTION=0.5):
    OUTPUT_FILE = "mesh/rveinfo.txt"

    area_domain = DOMAIN_SIZE ** 2
    total_fiber_area = area_domain * VOLUME_FRACTION
    approx_radius = math.sqrt(total_fiber_area / (NUM_FIBERS_TARGET * math.pi))
    min_distance_btw_fibers = 0.1 * approx_radius
    MAX_ATTEMPTS_PER_FIBER = 1000

    def intersects(x, y, r, circles):
        for cx, cy, cr in circles:
            if math.hypot(x - cx, y - cy) < (r + cr + min_distance_btw_fibers):
                return True
        return False

    def create_fibers(number):
        fib = []
        for _ in range(number):
            success = False
            for attempt in range(MAX_ATTEMPTS_PER_FIBER):
                r = approx_radius
                x = random.uniform(0, DOMAIN_SIZE)
                y = random.uniform(0, DOMAIN_SIZE)

                images = [(x, y)]
                if x < r:
                    images.append((x + DOMAIN_SIZE, y))
                if x > DOMAIN_SIZE - r:
                    images.append((x - DOMAIN_SIZE, y))
                if y < r:
                    images.append((x, y + DOMAIN_SIZE))
                if y > DOMAIN_SIZE - r:
                    images.append((x, y - DOMAIN_SIZE))

                if x < r:
                    if y < r:
                        images.append((x + DOMAIN_SIZE, y + DOMAIN_SIZE))
                    if y > DOMAIN_SIZE - r:
                        images.append((x + DOMAIN_SIZE, y - DOMAIN_SIZE))
                if x > DOMAIN_SIZE - r:
                    if y < r:
                        images.append((x - DOMAIN_SIZE, y + DOMAIN_SIZE))
                    if y > DOMAIN_SIZE - r:
                        images.append((x - DOMAIN_SIZE, y - DOMAIN_SIZE))

                if all(not intersects(ix, iy, r, fib) for ix, iy in images):
                    for ix, iy in images:
                        fib.append((ix, iy, r))
                    success = True
                    break
            if not success:
                return None
        return fib

    n = 0
    while True:
        n += 1
        fibers = []
        fib = create_fibers(NUM_FIBERS_TARGET)
        if fib is not None:
            fibers.extend(fib)
            print(f"Successfully placed {NUM_FIBERS_TARGET} physical fibers.")
            break
        else:
            print("Failed to place all fibers. Restarting from scratch...", n)

    with open(OUTPUT_FILE, "w") as f:
        for x, y, r in fibers:
            f.write(f"Centroid: ({x:.6f}, {y:.6f})\n")
            f.write(f"Circumradius: {r:.6f}\n\n")

    print(f"Wrote {len(fibers)} total (including images) to {OUTPUT_FILE}")

def read_fibers(txt_path):
    fibers = []
    with open(txt_path, 'r') as f:
        lines = f.readlines()
    for i in range(0, len(lines), 3):
        line_c = lines[i].strip()
        line_r = lines[i + 1].strip()
        x, y = map(float, line_c.replace('Centroid:', '').strip('() ').split(','))
        r = float(line_r.replace('Circumradius:', '').strip())
        fibers.append((x, y, r))
    return fibers


def create_circle(x, y, r, num_segments=20):
    points = []
    for i in range(num_segments):
        theta = 2 * math.pi * i / num_segments
        px = x + r * math.cos(theta)
        py = y + r * math.sin(theta)
        pt = gmsh.model.occ.addPoint(px, py, 0)
        points.append(pt)

    lines = []
    for i in range(num_segments):
        start = points[i]
        end = points[(i + 1) % num_segments]
        l = gmsh.model.occ.addLine(start, end)
        lines.append(l)

    loop = gmsh.model.occ.addCurveLoop(lines)
    surf = gmsh.model.occ.addPlaneSurface([loop])
    return (2, surf)



def create_mesh(fibers, domain_size=10.0):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("RVE")

    square_tag = gmsh.model.occ.addRectangle(0, 0, 0, domain_size, domain_size)
    square_entity = (2, square_tag)

    trimmed_fibers = []
    for x, y, r in fibers:
        disk = create_circle(x, y, r, num_segments=20)
        intersected, _ = gmsh.model.occ.intersect([disk], [square_entity], removeObject=True, removeTool=False)
        trimmed_fibers.extend(intersected)

    gmsh.model.occ.synchronize()
    matrix_cut, _ = gmsh.model.occ.cut([square_entity], trimmed_fibers, removeObject=True, removeTool=False)
    gmsh.model.occ.synchronize()

    all_surfaces = gmsh.model.occ.getEntities(dim=2)
    trimmed_fiber_tags = [tag for dim, tag in trimmed_fibers if dim == 2]
    matrix_tags = [tag for dim, tag in all_surfaces if dim == 2 and tag not in trimmed_fiber_tags]

    if matrix_tags:
        gmsh.model.addPhysicalGroup(2, matrix_tags, tag=1)
        gmsh.model.setPhysicalName(2, 1, "Matrix")

    if trimmed_fiber_tags:
        gmsh.model.addPhysicalGroup(2, trimmed_fiber_tags, tag=2)
        gmsh.model.setPhysicalName(2, 2, "Fibers")

    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.2)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.4)

    gmsh.model.mesh.generate(2)
    gmsh.write("mesh/rve.msh")
    gmsh.finalize()


def convert_to_XDMF():
    msh = meshio.read("mesh/rve.msh")

    if "triangle" in msh.cells_dict:
        cell_type = "triangle"
    else:
        raise ValueError("Mesh must contain triangle elements.")

    cells = [(cell_type, msh.cells_dict[cell_type])]
    physical_tags = msh.cell_data_dict["gmsh:physical"][cell_type]

    
    meshio.write("mesh/rve.xdmf", meshio.Mesh(
        points=msh.points[:, :2],
        cells=cells,
        cell_data={"name_to_read": [physical_tags]}
    ))

    print("Mesh saved in XDMF.")




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--fibers", type=int, default=25)
    parser.add_argument("--size", type=float, default=10.0)
    parser.add_argument("--vf", type=float, default=0.5)
    args = parser.parse_args()

    gen_text_file(NUM_FIBERS_TARGET=args.fibers, DOMAIN_SIZE=args.size, VOLUME_FRACTION=args.vf)
    fibers = read_fibers("mesh/rveinfo.txt")
    create_mesh(fibers, domain_size=args.size)
    convert_to_XDMF()

