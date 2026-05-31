
def write_oofem_input(globvar):

    Plate_By = globvar.Bp_y

    ndofman = int( globvar.ndofman_concrete + globvar.ndofman_rebar ) + 2
    nelem = int( globvar.nelem_concrete + globvar.nelem_rebar )

    DS = (globvar.rebars_CS[globvar.DS]).diam
    ncrosssect = 2
    crosssect_3 =""
    min_step_length = globvar.min_step_length

    project_name = f"{globvar.project_name}"

    if globvar.topology_layout == 0 and DS != 0:
        ncrosssect += 1
        crosssect_3 = f"SimpleCS 3 area {globvar.rebars_CS[globvar.DS].area :.8f}\n"

    if globvar.load_plate_length <= 0.3:
        min_step_length = min_step_length / 2


    #### Writing file ###
    file_path_name = globvar.file_path / f"{project_name}.in"
    with open(file_path_name, "w") as f:

        ### output file name
        f.write(f"{project_name}.out\n")

        ### extractor record
        #f.write("#%BEGIN_CHECK%\n#TIME\n#LOADLEVEL\n#DOFMAN number 100000 dof 3 unknown d\n#%END_CHECK%\n")

        ### job description
        f.write(f"oofem input file for multibear    plate_By / By ratio: {str(globvar.file_path)[-5:]}\n")

        ### output file info
        f.write(f"# Bx (half use for FEM) = {globvar.Bx / 2} \n# By = {globvar.By} \n# Bz = {globvar.Bz} \n# Plate_Bx = {globvar.Bp / 2} \n# Plate_By = {Plate_By} \n# Elemsize_X = {globvar.elem_size_X} \n# Elemsize_YZ = {globvar.elem_size_YZ}\n")

        f.write("\n")

        ### analysis record
        f.write(f"NonLinearStatic nsteps {globvar.n_steps} stiffmode 2 rtolv 1e-4 stepLength 1.e-4 minStepLength {min_step_length} hpc 2 100000 3 hpcw 1 -1. Psi 0.0 MinIter 3 MaxIter 200 ReqIterations 20 lstype 4 smtype 8 initialguess 1 maxrestarts 0 renumber 1 nmodules 2\n")

        ### export module record
        f.write("VTKXML primvars 1 1 vars 3 1 2 4 cellvars 2 13 27 tstep_step 10 stype 0 regionsets 1 2\n")
        f.write("VTKXML primvars 1 1 vars 2 1 4 cellvars 1 27 tstep_step 10 stype 0 regionsets 1 3\n")

        f.write("\n")

        ### domain record
        f.write("domain 3d\n")

        ### output manager record
        f.write("outputmanager tstep_all dofman_output {100000}\n")

        ### components size record
        f.write(f"ndofman {ndofman} nelem {nelem} ncrosssect {ncrosssect} nmat 2 nbc 2 nic 0 nltf 2 nset 3  \n")

        #### nodes
        f.write("\n# NODES\n")

        f.write("\n# Concrete\n")
        f.writelines(globvar.nodes_concrete)

        f.write("\n# Rebars\n")
        f.writelines(globvar.nodes_rebar)

        f.write("\n# Master nodes\n")
        f.write("node 100000 coords 3 0. 0. 0. bc 3 1 0 0 load 1 2 \n")
        f.write("node 100001 coords 3 0. 0. 0. bc 3 1 1 0 \n")


        #### elements
        f.write("\n# ELEMENTS\n")

        f.write("\n# Concrete\n")
        f.writelines(globvar.elements_concrete)

        f.write("\n# Rebars\n")
        f.writelines(globvar.elements_rebar)

        f.write("\n")

        ### cross section records(s)
        f.write("SimpleCS 1 thick 1.\n")
        f.write(f"SimpleCS 2 area {globvar.rebars_CS[globvar.DL].area :.8f}\n")
        f.write(crosssect_3)

        f.write("\n")

        ### material type record(s)
        for mat in globvar.materials:
            mat_type = mat.give_material_type()
            if (mat_type == "concrete"):
                Ec = mat.Eci
                fcm = mat.fcm
                fctm = mat.fctm
                nu = mat.nu
                wf = mat.wf
                ecc = mat.ecc
                kinit = mat.kinit
                Hp = mat.Hp
                dilation = mat.dilation
                Ahard = mat.Ahard
                Bhard = mat.Bhard
                Chard = mat.Chard
                Dhard = mat.Dhard
                Asoft = mat.Asoft
                efc = mat.efc
                wf1 = mat.wf1
                ft1 = mat.ft1

        f.write(f"Con2DPM 1 d 0. E {(Ec):.3f} n {(nu):.3f} talpha 0. fc {(fcm):.3f} ft {(fctm):.3f} hp {(Hp):.3f} Asoft {(Asoft):.3f} dilation {(dilation):.3f} wf {(wf):.4e} Ahard {(Ahard):.4e} Bhard {(Bhard):4e} Chard {(Chard):.4e} Dhard {(Dhard):.4e} wf1 {(wf1):.3f} ft1 {(ft1):.3f} ecc {(ecc):.3f} kinit {(kinit):.3f} efc {(efc):.4e}\n")
        f.write(f"MisesMat 2 d 0. talpha 0. E 200.e3 n 0.3 sig0 {globvar.fy_lat}  H 0. omega_crit 1. a 0.\n")

        f.write("\n")

        ### load, boundary conditions record(s)
        f.write("BoundaryCondition  1 loadTimeFunction 1 prescribedvalue 0.\n")
        f.write("NodalLoad 2 loadTimeFunction 1 Components 3 0. 0. -0.5 ")
        f.write("\n")

        ### time function(s)
        f.write("ConstantFunction 1 f(t) 1.0\n")
        f.write("PiecewiseLinFunction 2 nPoints 2 t 2 0. 1. f(t) 2 0. 1.\n")

        f.write("\n")

        ### set record(s)
        f.write("Set 1 allElements\n")
        f.write(f"Set 2 elementranges {{(1 {int( globvar.nelem_concrete)} )}}\n")
        f.write(f"Set 3 elementranges {{( {int( globvar.nelem_concrete) + 1 } {int( globvar.nelem_concrete) + int( globvar.nelem_rebar)} )}}\n")

    #################

    print(f"\nndofman: {ndofman}, nelem: {nelem}\n")
