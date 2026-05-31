import gmsh

def generate_rebars_mesh(globvar):

    # variables from globvars
    topology_layout = globvar.topology_layout

    By = globvar.By
    Bx = By / 2
    Bz = globvar.Bz

    DL = (globvar.rebars_CS[globvar.DL]).diam
    DT = (globvar.rebars_CS[globvar.DT]).diam
    DS = (globvar.rebars_CS[globvar.DS]).diam
    dS = globvar.dS
    corner_radius = globvar.corner_radius
    cover = globvar.cover
    Hz = globvar.s

    elem_size = min(globvar.elem_size_YZ, globvar.elem_size_X)

    rebar = []
    small_rebar = []

    gmsh.initialize()
    gmsh.model.add("rebar")

    ### SPIRALS ###
    if topology_layout == 0:

        # points
        p1 = gmsh.model.occ.addPoint(0, cover + DL / 2, Bz - Hz/2)  # start
        p2 = gmsh.model.occ.addPoint(0, By - (cover + DL / 2),  Bz - Hz/2)  # end

        c1 = gmsh.model.occ.addPoint(0, By / 2,  Bz - Hz/2)  # center of circle

        # circle
        a1 = gmsh.model.occ.addCircleArc(p2, c1, p1)

        # sync
        gmsh.model.occ.synchronize()

        # remove center point
        gmsh.model.occ.remove([(0, c1)], recursive=False)

        rebar = [(1, a1)]

        ### multipal spirals ###
        if DS != 0:
            ## points
            s_p11 = gmsh.model.occ.addPoint(Bx - (cover + dS / 2 + DS / 2), cover + DS / 2, Bz - Hz)
            s_p12 = gmsh.model.occ.addPoint(Bx - (cover + dS / 2 + DS / 2), cover + dS + DS / 2, Bz - Hz)

            s_p21 = gmsh.model.occ.addPoint(Bx - (cover + dS / 2 + DS / 2), By - (cover + DS / 2), Bz - Hz)
            s_p22 = gmsh.model.occ.addPoint(Bx - (cover + dS / 2 + DS / 2), By - (cover + dS + DS / 2), Bz - Hz)

            # center of spirals
            s_c1 = gmsh.model.occ.addPoint(Bx - (cover + dS / 2 + DS / 2), cover + dS / 2 + DS / 2, Bz - Hz)
            s_c2 = gmsh.model.occ.addPoint(Bx - (cover + dS / 2 + DS / 2), By - (cover + dS / 2 + DS / 2), Bz - Hz)

            ## make circle
            s_a11 = gmsh.model.occ.addCircleArc(s_p11, s_c1, s_p12)
            s_a12 = gmsh.model.occ.addCircleArc(s_p12, s_c1, s_p11)

            s_a21 = gmsh.model.occ.addCircleArc(s_p21, s_c2, s_p22)
            s_a22 = gmsh.model.occ.addCircleArc(s_p22, s_c2, s_p21)

            # sync
            gmsh.model.occ.synchronize()

            # remove center points
            gmsh.model.occ.remove([(0, s_c1), (0, s_c2)], recursive=False)

            small_rebar = [(1, s_a11), (1, s_a12), (1, s_a21), (1, s_a22)]

    ### TIES ###
    elif topology_layout == 1:

        cover = cover + DT/2

        # points
        p1 = gmsh.model.occ.addPoint(0, cover,  Bz - Hz/2)
        p2 = gmsh.model.occ.addPoint(Bx -(cover + corner_radius), cover,  Bz - Hz/2)
        p3 = gmsh.model.occ.addPoint(Bx - cover, cover + corner_radius,  Bz - Hz/2)
        p4 = gmsh.model.occ.addPoint(Bx - cover, By -(cover + corner_radius),  Bz - Hz/2)
        p5 = gmsh.model.occ.addPoint(Bx -(cover + corner_radius), By - cover,  Bz - Hz/2)
        p6 = gmsh.model.occ.addPoint(0, By - cover,  Bz - Hz/2)

        # center of corner circles
        c1 = gmsh.model.occ.addPoint(Bx -(cover + corner_radius), cover + corner_radius,  Bz - Hz/2)
        c2 = gmsh.model.occ.addPoint(Bx -(cover + corner_radius), By -(cover + corner_radius),  Bz - Hz/2)

        # curves
        l1 = gmsh.model.occ.addLine(p1, p2)
        a1 = gmsh.model.occ.addCircleArc(p2, c1, p3)
        l2 = gmsh.model.occ.addLine(p3, p4)
        a2 = gmsh.model.occ.addCircleArc(p4, c2, p5)
        l3 = gmsh.model.occ.addLine(p5, p6)

        gmsh.model.occ.remove([(0, c1), (0, c2)], recursive=False)

        # sync
        gmsh.model.occ.synchronize()

        # Base tie for copy
        rebar = [(1, l1), (1, a1), (1, l2), (1, a2), (1, l3)]

        DS = 0

    else:
        raise ValueError("Unknown topology")

    # Copy rebars
    for i in range(1, int(Bz / Hz)):
        new_entities = gmsh.model.occ.copy(rebar)
        gmsh.model.occ.translate(new_entities, 0, 0, -(i * Hz))

    if DS != 0:
        for i in range(1, int(Bz / Hz) - 1):
            new_entities_small = gmsh.model.occ.copy(small_rebar)
            gmsh.model.occ.translate(new_entities_small, 0, 0, -(i * Hz))

    # generate outlines for showing mesh
    if globvar.flag_show_mesh:

        # reference points
        rp1 = gmsh.model.occ.addPoint(0, 0, 0)
        rp2 = gmsh.model.occ.addPoint(Bx, 0, 0)
        rp3 = gmsh.model.occ.addPoint(Bx, By, 0)
        rp4 = gmsh.model.occ.addPoint(0, By, 0)
        rp5 = gmsh.model.occ.addPoint(0, 0, Bz)
        rp6 = gmsh.model.occ.addPoint(Bx, 0, Bz)
        rp7 = gmsh.model.occ.addPoint(Bx, By, Bz)
        rp8 = gmsh.model.occ.addPoint(0, By, Bz)

        # reference lines
        rl1 = gmsh.model.occ.addLine(rp1, rp2)
        rl2 = gmsh.model.occ.addLine(rp2, rp3)
        rl3 = gmsh.model.occ.addLine(rp3, rp4)
        rl4 = gmsh.model.occ.addLine(rp4, rp1)
        rl5 = gmsh.model.occ.addLine(rp5, rp6)
        rl6 = gmsh.model.occ.addLine(rp6, rp7)
        rl7 = gmsh.model.occ.addLine(rp7, rp8)
        rl8 = gmsh.model.occ.addLine(rp8, rp5)
        rl9 = gmsh.model.occ.addLine(rp1, rp5)
        rl10 = gmsh.model.occ.addLine(rp2, rp6)
        rl11 = gmsh.model.occ.addLine(rp3, rp7)
        rl12 = gmsh.model.occ.addLine(rp4, rp8)


    # sync
    gmsh.model.occ.synchronize()

    # Mesh density
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), elem_size / 2)

    # generate 1d mesh
    gmsh.model.mesh.generate(1)

    # show mesh
    if globvar.flag_show_mesh:
        gmsh.fltk.run()
        # end gmsh
        gmsh.finalize()
        return

    #############
    #OOOFEM INPUT FILE
    #############

    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    coords = list(zip(coords[0::3], coords[1::3], coords[2::3]))
    elem_types, elem_tags, elem_nodes = gmsh.model.mesh.getElements()

    nodes_lines_rebars = []
    ndofman_concrete = int(globvar.ndofman_concrete)
    border_nodes = []
    z_1 = Bz
    for tag, (x, y, z) in zip(node_tags, coords):
        line = f"hangingNode {int(tag + ndofman_concrete)}   coords 3  {x:.6f} {y:.6f} {z:.6f}     doftype 3 2 2 2\n"
        nodes_lines_rebars.append(line)

        if z > z_1 :
            border_nodes.append(tag)
        z_1 = z


    elem_lines_rebars = []
    elem_counter = int(globvar.nelem_concrete + 1)

    border_nodes.pop(-1)
    while len(border_nodes) > 1:
        border_nodes.pop()

    crossSect = 2
    for etype, tags, nodes in zip(elem_types, elem_tags, elem_nodes):

        if etype == 1:
            nnode = 2
            for i, t in enumerate(tags):
                n = nodes[i * nnode: (i + 1) * nnode]

                if min(n) in border_nodes:
                    crossSect = 3

                nds_str = " ".join(str(int(nj + ndofman_concrete)) for nj in n)
                elem_lines_rebars.append(f"Truss3d {elem_counter} nodes 2 {nds_str}     crossSect {crossSect}   mat 2\n")
                elem_counter += 1


    # count component size
    ndofman_rebars = len(nodes_lines_rebars)
    nelem_rebars = len(elem_lines_rebars)

    # end gmsh
    gmsh.finalize()
    
    return nodes_lines_rebars, elem_lines_rebars, ndofman_rebars, nelem_rebars


