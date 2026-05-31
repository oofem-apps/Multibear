import gmsh
import warnings


def generate_concrete_mesh(globvar):

    By = globvar.By
    Bx = By / 2
    Bz = globvar.Bz
    Plate_By = globvar.Bp_y
    Plate_Bx = globvar.Bp /2
    cover = globvar.cover



    # ignore cover
    #if globvar.flag_ignore_cover == 2:
    #    cover = 0

    cover = 0 #always ignor cover

    Bx = Bx - cover
    By = By - 2*cover


    Elem_size_YZ = globvar.elem_size_YZ
    Elem_size_X = globvar.elem_size_X

    Elem_size_plate = Elem_size_YZ
    # 5913, nelem: 6334

    # check if elem_size is doable
    if Elem_size_plate >= (Bz - 2 * Plate_By) / 4:
        Elem_size_plate = round( ((Bz - 2 * Plate_By) / 4) * 0.95 , 5 )
        warnings.warn("⚠️ Element size YZ too big ⚠️ ")
        print(f"New Element size near plate: {Elem_size_plate}")


    N_X = round(Bx / Elem_size_X)
    #plate_Bx to Bx ratio
    Bx_ratio = Bx / Plate_Bx
    N_X1 = round(N_X / Bx_ratio)
    N_X2 = N_X - N_X1


    gmsh.initialize()

    # Mesh settings
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", 1.0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", min(Elem_size_YZ / 1.5 , Elem_size_plate))
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", Elem_size_YZ)

    # 1) Points
    coords = {
        'A': (Plate_Bx, 0 + cover, Bz),
        'B': (Plate_Bx, By + cover, Bz),
        'C': (Plate_Bx, By + cover, 0),
        'D': (Plate_Bx, 0 + cover, 0),
        'E': (Plate_Bx, Plate_By, Bz),
        'F': (Plate_Bx, Plate_By + Elem_size_plate, Bz),
        'G': (Plate_Bx, 0 + cover, Bz - 2 * (Plate_By + Elem_size_plate)),
        'H': (Plate_Bx, 0 + cover, Bz - 2 * Plate_By),
        'I': (Plate_Bx, Plate_By + Elem_size_plate , 0)
    }
    p = {k: gmsh.model.occ.addPoint(*v, 0.1) for k, v in coords.items()}

    # 2) All boundary & partition lines
    lines = {}
    # outer boundary
    lines['AE'] = gmsh.model.occ.addLine(p['A'], p['E'])
    lines['EF'] = gmsh.model.occ.addLine(p['E'], p['F'])
    lines['FB'] = gmsh.model.occ.addLine(p['F'], p['B'])
    lines['BC'] = gmsh.model.occ.addLine(p['B'], p['C'])
    lines['CI'] = gmsh.model.occ.addLine(p['C'], p['I'])
    lines['ID'] = gmsh.model.occ.addLine(p['I'], p['D'])
    lines['DG'] = gmsh.model.occ.addLine(p['D'], p['G'])
    lines['GH'] = gmsh.model.occ.addLine(p['G'], p['H'])
    lines['HA'] = gmsh.model.occ.addLine(p['H'], p['A'])
    # interior splitting lines
    lines['EH'] = gmsh.model.occ.addLine(p['E'], p['H'])
    lines['FG'] = gmsh.model.occ.addLine(p['F'], p['G'])
    lines['FI'] = gmsh.model.occ.addLine(p['F'], p['I'])

    # 3) Big outer surface
    outer_loop = gmsh.model.occ.addCurveLoop([
        lines['AE'], lines['EF'], lines['FB'], lines['BC'],
        lines['CI'], lines['ID'], lines['DG'], lines['GH'], lines['HA']
    ])
    outer_surf = gmsh.model.occ.addPlaneSurface([outer_loop])

    # 4) Fragment it by all interior lines
    entities_to_split = [(2, outer_surf)]
    cutting_curves = [(1, lines[l]) for l in ['EH', 'FG', 'FI']]
    subs, cuts = gmsh.model.occ.fragment(entities_to_split, cutting_curves)

    # Sync geometry
    gmsh.model.occ.synchronize()

    n_y = round((By - (Plate_By + Elem_size_YZ )) / Elem_size_YZ)
    n_z = round(Bz / Elem_size_YZ)

    if n_y % 2 != 0:
        n_y += 1
    if n_z % 2 != 0:
        n_z += 1

    gmsh.model.mesh.setTransfiniteCurve(19, n_y + 1) #FB
    gmsh.model.mesh.setTransfiniteCurve(17, n_y + 1) #CI
    gmsh.model.mesh.setTransfiniteCurve(12, n_z + 1) #FI
    gmsh.model.mesh.setTransfiniteCurve(18, n_z + 1) #BC

    gmsh.model.mesh.setTransfiniteSurface(3)

    n_p = round(Elem_size_YZ / (Elem_size_YZ / 2))
    n_t = round(((5**0.5)* Plate_By) / (Elem_size_YZ / 2))

    if n_p % 2 != 0:
        n_p += 1
    if n_t % 2 != 0:
        n_t += 1

    gmsh.model.mesh.setTransfiniteCurve(16, n_p + 1) #EF
    gmsh.model.mesh.setTransfiniteCurve(15, n_p + 1) #GH
    gmsh.model.mesh.setTransfiniteCurve(10, n_t + 1) #EH
    gmsh.model.mesh.setTransfiniteCurve(11, n_t + 1) #FG

    gmsh.model.mesh.setTransfiniteSurface(2)

    # Setting quad mesh for all resulting surfaces
    selected_surfaces = [2, 3, 4]
    for surface_id in selected_surfaces:
        gmsh.model.mesh.setRecombine(2, surface_id)

    #  Define size field to refine mesh near line EH
    eh_id = lines['EH']

    # Field 1: Distance from line EH
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "EdgesList", [eh_id])

    # Field 2: Threshold to control element size based on distance
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", Elem_size_YZ / 2)  # Smallest elements near the line
    gmsh.model.mesh.field.setNumber(2, "SizeMax", Elem_size_YZ)      # Larger elements further from the line
    gmsh.model.mesh.field.setNumber(2, "DistMin", Elem_size_YZ * 1)
    gmsh.model.mesh.field.setNumber(2, "DistMax", Elem_size_YZ * 2)

    # Set Field 2 as the background mesh size controller
    gmsh.model.mesh.field.setAsBackgroundMesh(2)

    #  Create 3D model by extrusion
    volumes_1 = gmsh.model.occ.extrude([(2, surf) for dim, surf in subs if dim == 2],
                                    -Plate_Bx, 0, 0,
                                    numElements=[N_X1],
                                    recombine=True
                                    )

    volumes_2 = gmsh.model.occ.extrude([(2, surf) for dim, surf in subs if dim == 2],
                                    Bx - Plate_Bx, 0, 0,
                                    numElements=[N_X2 * 2],
                                    recombine=True
                                    )

    # Sync after extrusion
    gmsh.model.occ.synchronize()

    # Obtaining ID volumes from extrusion (3 dimensions only)
    volumes_1_ids = [vol[1] for vol in volumes_1 if vol[0] == 3]
    volumes_2_ids = [vol[1] for vol in volumes_2 if vol[0] == 3]

    # Fuze 2 volumes
    volumes_fused, _ = gmsh.model.occ.fuse(
        [(3, vol_id) for vol_id in volumes_1_ids],
        [(3, vol_id) for vol_id in volumes_2_ids]
    )

    # physical groups
    for i, (dim, vol_id, *_) in enumerate(volumes_fused, start=1):
        if dim == 3:
            gmsh.model.addPhysicalGroup(3, [vol_id], tag=i)
            gmsh.model.setPhysicalName(3, i, f"fused_volume{i}")


    #gmsh.model.occ.removeAllDuplicates()

    #  Generate 3D mesh
    gmsh.model.mesh.generate(3)

    # Show mesh and end function
    if globvar.flag_show_mesh:
        # show mesh
        gmsh.fltk.run()
        # close gmsh
        gmsh.finalize()
        return

    ############
    # OOFEM INPUT
    ############

    node_tags, coords, _        = gmsh.model.mesh.getNodes()
    coords = list(zip(coords[0::3], coords[1::3], coords[2::3]))
    elem_types, elem_tags, elem_nodes = gmsh.model.mesh.getElements()

    ### Exporting ###

    # nodes
    nodes_lines_concrete = []
    for tag, (x, y, z) in zip(node_tags, coords):
        if z == Bz and x <= Plate_Bx and y <= Plate_By:
            line = f"slaveNode {tag}   coords 3  {x:.6f} {y:.6f} {z:.6f}    dofType 3 2 2 2   masterDofMan 2 100000 100001    weights 2 1. {(-1+2*(y/Plate_By)):.6f}"
        else:
            line = f"node {tag}   coords 3  {x:.6f} {y:.6f} {z:.6f}"

            if z == 0:
                line += "   bc 3 1 1 1"
            elif x == 0:
                line += "   bc 3 1 0 0"

        line += "\n"
        nodes_lines_concrete.append(line)

    # elemnts
    elem_lines_concrete = []
    elem_counter = 1
    for etype, tags, nodes in zip(elem_types, elem_tags, elem_nodes):

        if etype == 5:
            # square bricks
            nnode = 8
            for i, t in enumerate(tags):
                n = nodes[i * nnode: (i + 1) * nnode]
                # rearenge nodes
                nds = [n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7]]

                nds_str = " ".join(str(n) for n in nds)
                elem_lines_concrete.append(f"LSpaceBB {elem_counter}   nodes 8 {nds_str}    crossSect 1   mat 1   nlgeo 1\n")
                elem_counter += 1

        elif etype == 6:
            # triangle bricks
            nnode = 6
            for i, t in enumerate(tags):
                n = nodes[i * nnode: (i + 1) * nnode]  # [n1, n2, n3, n4, n5, n6]
                # duplicate nodes
                nds8 = [n[0], n[1], n[2], n[2], n[3], n[4], n[5], n[5]]  # [1, 2, 3, 3, 4, 5, 6, 6]

                nds_str = " ".join(str(x) for x in nds8)
                elem_lines_concrete.append(f"LSpace {elem_counter} nodes 8 {nds_str}   crossSect 1   mat 1   nlgeo 1\n")
                elem_counter += 1

        else:
            continue  # skip non use elements

    # count component size
    ndofman_concrete = len(nodes_lines_concrete)
    nelem_concrete = len(elem_lines_concrete)

    # close gmsh
    gmsh.finalize()

    return nodes_lines_concrete, elem_lines_concrete, ndofman_concrete, nelem_concrete

