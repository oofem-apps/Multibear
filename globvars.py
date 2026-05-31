from cross_sections import rebar_CS
from PySide2.QtCore import Qt
from datetime import date
from tasks import task, Task_status
import math

class globvars:
    
    def __init__(self):
    
        ## definition of global variables and containers etc...
        self.base_name = "multibear_oofem"

        self.program_name = "Multibear"

        self.project_name = "FEM_multibear"
     
        self.filename_in = self.base_name + ".in"
        self.filename_out = self.base_name + ".out"

        self.file_path =""

        self.oofem_folder = []
        self.cpu_nr =  []

        self.loading = None

        ### pre-defined cross-sections
        self.cross_section_types = ["5S4", "6S4", "user-square", "user-rectangle"]
        self.cross_section_type = self.cross_section_types[0]

        #Layout
        self.topology_layout = 0
        self.id_plot = 0

        # TOPOLOGY:

        self.Bx = 0.4
        self.By = self.Bx
        self.Bz = round(self.Bx * 1.5, 3)
        self.s = self.Bx / 10.
        self.cover = 0.02
        self.max_aspect_ratio = 2.
        self.global_div_length = 10.e-3

        self.elem_size_X =  0.035
        self.elem_size_YZ = 0.035


        self.load_plate_width = 0.5
        self.Bp = self.load_plate_width * self.Bx
        self.load_plate_length = 0.5
        self.Bp_y = self.load_plate_length * self.By

        # DUMMY NODES

        #self.master_node_force = 1000000   #remove
        #self.master_node_shortening = self.master_node_force+1 #remove
        #self.master_node_bending_x = self.master_node_force+2  #remove
        #self.master_node_bending_y = self.master_node_force+3  #remove

 
        # global numbers
        self.ndofman = 0
        self.nelem = 0
        self.n_steps = 250 #max number of steps, usually 50-100 is enough
        self.estimate = 1.
        self.min_step_length = 4.e-5

        # containers - mesh
        self.master_nodes = []
        self.hanging_nodes = []
        self.brick_elements = []
        self.truss_elements = []

        self.materials = []
        self.rebars = []
        self.cross_sections = []

        self.bcs = []
        self.ltfs = []
        self.sets = []
                
        # MATERIALS:
        
        self.material_type = ("concrete", "rebar", "tendon") 

        self.rebarmat = ['fy 500 MPa, eps = infty', # 0 
                         'SD420 fy 473 MPa, fu = 668 MPa, eps = 10%', # 1
                         'SD420 fy 497 MPa, fu = 723 MPa, eps = 10%', # 2
                         'B500A fy 500 MPa, fu = 525 MPa, eps = 2.5%', # 3
                         'B500B fy 500 MPa, fu = 540 MPa, eps = 5.0%', # 4
                         'B500C fy 500 MPa, fu = 575 MPa, eps = 7.5%'] # 5

        self.tendonmat = ['1800 MPa, todo, todo'] # 0

        self.fcm = 'C30'
        self.fy_vert = 500.
        self.fy_lat = 500.
        
        self.Es = 210.e3
               
        delta_f = 8.
        self.concretes = {}
        self.concretes['C20'] = 20.+delta_f
        self.concretes['C25'] = 25.+delta_f
        self.concretes['C30'] = 30.+delta_f
        self.concretes['C35'] = 35.+delta_f
        self.concretes['C40'] = 40.+delta_f
        self.concretes['C45'] = 45.+delta_f
        self.concretes['C50'] = 50.+delta_f
        
        # REINFORCEMENT
        self.minimum_radius = 0.10/2.
        self.rho_lat = 0.
        self.rho_lat_tie = 0.
        self.rho_vert = 0.
        self.n_v_bars = 20 #multiple of 4
        
        self.rebars_CS = {}
        self.rebars_CS['none'] = rebar_CS(diam = 0. ) #careful with this one!

        self.rebars_CS['D6'] = rebar_CS(diam = 6.e-3)
        self.rebars_CS['D8'] = rebar_CS(diam = 8.e-3)
        self.rebars_CS['D10'] = rebar_CS(diam = 10.e-3)
        self.rebars_CS['D12'] = rebar_CS(diam = 12.e-3)
        self.rebars_CS['D14'] = rebar_CS(diam = 14.e-3)
        self.rebars_CS['D16'] = rebar_CS(diam = 16.e-3)
        self.rebars_CS['D18'] = rebar_CS(diam = 18.e-3)
        self.rebars_CS['D20'] = rebar_CS(diam = 20.e-3)
        self.rebars_CS['D22'] = rebar_CS(diam = 22.e-3)
        self.rebars_CS['D25'] = rebar_CS(diam = 25.e-3)
        self.rebars_CS['D28'] = rebar_CS(diam = 28.e-3)
        self.rebars_CS['D32'] = rebar_CS(diam = 32.e-3)
        
        self.rebars_CS['#3'] = rebar_CS(diam = 9.525e-3, area = 71.e-6)
        self.rebars_CS['#4'] = rebar_CS(diam = 12.7e-3, area = 129.e-6)
        self.rebars_CS['#5'] = rebar_CS(diam = 15.875e-3, area = 200.e-6)
        self.rebars_CS['#6'] = rebar_CS(diam = 19.05e-3, area = 284.e-6)
        self.rebars_CS['#7'] = rebar_CS(diam = 22.225e-3, area = 387.e-6)
        self.rebars_CS['#8'] = rebar_CS(diam = 25.4e-3, area = 509.e-6)
        self.rebars_CS['#9'] = rebar_CS(diam = 28.65e-3, area = 645.e-6)
        self.rebars_CS['#10'] = rebar_CS(diam = 32.26e-3, area = 819.e-6)
        self.rebars_CS['#11'] = rebar_CS(diam = 35.81e-3, area = 1006.e-6)
        self.rebars_CS['#14'] = rebar_CS(diam = 43.e-3, area = 1452.e-6)
        self.rebars_CS['#18'] = rebar_CS(diam = 57.33e-3, area = 2581.e-6)


        # REBAR DIAMETERS
        self.DS = '#3'
        self.DL = '#3'
        self.DT = '#4'
        self.DV = 'none'

        self.DT_equivalent = 0.0133
        self.DL_equivalent = 0.009056

        # SPIRAL DIAMETERS
        self.dL = self.Bx - 2.*self.cover - (self.rebars_CS[self.DL]).diam
        self.dS = (self.Bx - 2.*self.cover) / 3

        # TIE DIAMETERS
        self.dT = self.Bx - 2.*self.cover  - (self.rebars_CS[self.DT]).diam
        self.corner_radius = 4.5 * (self.rebars_CS[self.DT]).diam
        # FLAGS - MISC

        self.flag_show_Bp = True

        self.debug_flag = False
        
        self.flag_problem_changed = True
        self.flag_output_generated = False
        self.flag_analyses_run = False
        
        self.flag_oofem_selected = False

        self.flag_loading_selected = False
        self.flag_loading_display = Qt.Unchecked
        
        # FLAGS - GUI STATES
        self.flag_id_steel = Qt.Unchecked
        self.flag_id_ACI = Qt.Unchecked
        self.flag_id_ASCE = Qt.Unchecked
        self.flag_id_concrete = Qt.Checked
        self.flag_id_CTU = Qt.Checked
        self.flag_id_CTU_concrete = Qt.Unchecked
        self.flag_id_FEM = Qt.Checked
        self.flag_id_FEM_paths = Qt.Checked

        # neutral axis - for debugging and interactive drawing
        self.flag_active_neutral_axis = Qt.Unchecked
        self.c_neutral_axis = self.By/2.

        self.flag_ignore_cover = Qt.Checked
        self.flag_show_mesh = False

        ### TASKS
        self.ecc_nr = 9
        self.tasks = []
        self.ecc_min = 1 / 6
        self.ecc_max = 2 / 3
        for i in range(self.ecc_nr):

            if ( i % 2 ):
                status = Task_status.UNSELECTED
            else:
                status = Task_status.SELECTED

            eccentricity_normalized = (1, self.ecc_min + (self.ecc_max - self.ecc_min) * i / (self.ecc_nr - 1))

            t = task(status, eccentricity_normalized)
                
            self.tasks.append(t)

        t = task(Task_status.PREDEFINED, [math.nan,math.nan])
        self.tasks.append(t)

        ### OOFEM nodes and elements lists
        self.nodes_concrete = []
        self.nodes_rebar = []
        self.elements_concrete = []
        self.elements_rebar = []

        self.ndofman_concrete = 0
        self.ndofman_rebar = 0
        self.nelem_concrete = 0
        self.nelem_rebar = 0



