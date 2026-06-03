import sys
import signal
import random
import matplotlib
import numpy as np
import math
from pathlib import Path
from time import time, sleep
import os
# to find the number of cores
import multiprocessing
import shutil
import pandas as pd


from nodes import node
from tasks import Task_status, task

matplotlib.use('Qt5Agg')

from PySide2.QtWidgets import QMainWindow, QApplication, QVBoxLayout, QGridLayout, QWidget, QLabel, \
    QPushButton, QLineEdit, QMessageBox, QDoubleSpinBox, QComboBox, QCheckBox, \
    QPlainTextEdit, QProgressBar, QFileDialog, QStyle, QSplashScreen, QStackedLayout,  QFrame

from PySide2.QtCore import Qt

from PySide2.QtCore import QObject, QThread, Signal
from PySide2.QtGui import QGuiApplication, QPixmap, QIcon, QPainter

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.transforms import Affine2D

import warnings
warnings.simplefilter('always', UserWarning)

import logging

from globvars import globvars

from rebars import rebar, hoop, spiral, tie
from cross_sections import cross_section
from elements import element
from materials import concrete_mat, rebar_mat, tendon_mat, window_CDPM2, window_Mises


from concrete_mesh import generate_concrete_mesh
from rebar_mesh import generate_rebars_mesh
from oofem_input import write_oofem_input
from diagram_functions import ACI, ASCE, id_fcm, AS_plain_concrete, AS_TIE, AS_MSR


class Worker(QObject):

    finished = Signal()
    progress = Signal(object, object, object)
    result = Signal(object)
    fail = Signal(object)

    def __init__(self, task, oofem_folder, n_cores, project_name):
        super().__init__()

    
        self.task = task
        self.load = []
        
        self.oofem_folder = oofem_folder
        self.n_cores = n_cores
        self.project_name = project_name
 
    #@Slot()  # QtCore.Slot
    def run(self):

        current_path = os.getcwd()

        problem_path = self.task.file_path
        os.chdir(problem_path)
        sys.path.append(self.oofem_folder)
        # the import needs to be in the "run" function, the imports from the above do not work
        try:
            import oofempy
            print("oofempy succsesfully imported")
        except ImportError:
            message = "oofempy import failed"
            warnings.warn(message)
            self.fail.emit(message)
            os.chdir(current_path)
            return

        #oofempy.init(logLevel=3, numberOfThreads=self.n_cores)
        oofempy.init(logLevel=2, numberOfThreads=self.n_cores) 
        self.task.status = Task_status.PROGRESS
        # TODO: task should have its name
        dr = oofempy.OOFEMTXTDataReader(f"{self.project_name}.in")
        problem = oofempy.InstanciateProblem(dr, oofempy.problemMode.processor, False, None, False)
        domain = problem.giveDomain(1)

        problem.init()
        problem.checkProblemConsistency()
        activeMStep = problem.giveMetaStep(1)
        problem.initMetaStepAttributes(activeMStep);

        load_level = 0.
        max_load_level = 0.

        n_steps = problem.giveNumberOfSteps()

        for timeStep in range(n_steps):
            if (load_level == max_load_level):

                problem.preInitializeNextStep()
                problem.giveNextStep()
                currentStep = problem.giveCurrentStep()
                problem.initializeYourself( currentStep )
                problem.solveYourselfAt( currentStep )
                problem.updateYourself( currentStep )
                load_level = problem.giveLoadLevel()
                
                if (load_level > max_load_level):
                    max_load_level = load_level
                
                problem.terminate( currentStep )
                self.progress.emit(timeStep, load_level, task)
                self.load.append(load_level)

            else:
                break


        problem.terminateAnalysis()

        self.task.load_level = self.load
        self.task.max_load = max_load_level

        self.result.emit(self.task)
        print ("FINISHED TASK")
        self.task.status = Task_status.COMPLETED
        self.finished.emit()
        
        os.chdir(current_path)


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

        
class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.w_cdpm2 = None

        self.globvar = globvars()

        self.setWindowTitle(self.globvar.program_name)


        # all ids to be updated
        self.need_update_ids_flag = True

        # init info labels
        self.label_info_topology = QLabel()
        self.label_info_topology_tie = QLabel()


        # MATERIALS
        # material: CONCRETE
        mat =  concrete_mat(globvar = self.globvar, nr = 1)
        self.globvar.materials.append(mat)

        #  material: LATERAL REINFORCEMENT
        mat =  rebar_mat(globvar = self.globvar, nr = 2, mat = self.globvar.rebarmat[0])
        self.globvar.materials.append(mat)

        #  material: VERTICAL REINFORCEMENT
        mat =  rebar_mat(globvar = self.globvar, nr = 3, mat = self.globvar.rebarmat[0])
        self.globvar.materials.append(mat)

        self.cdpm2 = self.globvar.materials[0]
        self.mises_lat = self.globvar.materials[1]
        self.mises_vert = self.globvar.materials[2]


        ### LATERAL REINFORCEMENT - topology
        self.vec = [0., 0., 1.]

        rebar_DS = (self.globvar.rebars_CS[self.globvar.DS]).diam
        rebar_DL = (self.globvar.rebars_CS[self.globvar.DL]).diam

        # positions of small spirals to fit at the corners
        auxX = self.globvar.Bx/2. - self.globvar.cover - self.globvar.dS/2 - rebar_DS/2.
        auxY = self.globvar.By/2. - self.globvar.cover - self.globvar.dS/2 - rebar_DS/2.

        ### large spiral
        center  = [0., 0., 0.]
        start = [ center[0], center[1] - self.globvar.dL/2., center[2] ]

        spiral_L = spiral(globvar = self.globvar, XYZ = center, cs_nr = 2, cs_tag = self.globvar.DL, cs_mat = 2, vec = self.vec, pitch = self.globvar.s, axial_length = self.globvar.s, start = start, tag ='L1')
        #spiral_L.check_topology_consistency(self.globvar)
        self.globvar.rebars.append(spiral_L)

        ### small spirals:
        # S1
        center = [auxX, -auxY, 0.]
        start = [ center[0], center[1] - self.globvar.dS/2., center[2] ]

        spiral_S = spiral(globvar = self.globvar, XYZ = center, cs_nr = 3, cs_tag = self.globvar.DS, cs_mat = 2, vec = self.vec, pitch = self.globvar.s, axial_length = self.globvar.s, start = start, tag ='S1')
        #spiral_S.check_topology_consistency(self.globvar)
        self.globvar.rebars.append(spiral_S)

        # S2
        center = [auxX, auxY, 0.]
        start = [ center[0], center[1] - self.globvar.dS/2., center[2] ]

        spiral_S = spiral(globvar = self.globvar, XYZ = center, cs_nr = 3, cs_tag = self.globvar.DS, cs_mat = 2, vec = self.vec, pitch = self.globvar.s, axial_length = self.globvar.s, start = start, tag ='S2')
        #spiral_S.check_topology_consistency(self.globvar)
        self.globvar.rebars.append(spiral_S)

        # S3
        center = [-auxX, auxY, 0.]
        start = [ center[0], center[1] - self.globvar.dS/2., center[2] ]

        spiral_S = spiral(globvar = self.globvar, XYZ = center, cs_nr = 3, cs_tag = self.globvar.DS, cs_mat = 2, vec = self.vec, pitch = self.globvar.s, axial_length = self.globvar.s, start = start, tag ='S3')
        #spiral_S.check_topology_consistency(self.globvar)
        self.globvar.rebars.append(spiral_S)

        # S4
        center = [-auxX, -auxY, 0.]
        start = [ center[0], center[1] - self.globvar.dS/2., center[2] ]

        spiral_S = spiral(globvar = self.globvar, XYZ = center, cs_nr = 3, cs_tag = self.globvar.DS, cs_mat = 2, vec = self.vec, pitch = self.globvar.s, axial_length = self.globvar.s, start = start, tag ='S4')
        #spiral_S.check_topology_consistency(self.globvar)
        self.globvar.rebars.append(spiral_S)

        ### Large tie
        center = [0., 0., 0.]
        vec = [0., 0., 1.]

        tie_L = tie( globvar=self.globvar, XYZ=center, cs_nr=2, cs_tag=self.globvar.DL, cs_mat=2, width=self.globvar.dT, height=self.globvar.dT, vec=vec, tag='TL')

        self.globvar.rebars.append(tie_L)

        ## VERTICAL REINFORCEMENT

        # generation of vertical reinforcement for the first time
        self.update_vertical_rebars()


        #########################
        # TOPOLOGY etc.
        #########################


        ### MSR

        # column breadth
        self.widget_Bx = QDoubleSpinBox()
        self.widget_Bx.setMinimum(0.2)
        self.widget_Bx.setMaximum(1.0)
        self.widget_Bx.setDecimals(3)
        self.widget_Bx.setValue( self.globvar.Bx )
        self.widget_Bx.setSingleStep(50./1000.)
        self.widget_Bx.valueChanged.connect(self.set_Bx)
        self.widget_Bx.valueChanged.connect(self.update_topology_plot)

        # spiral pitch
        self.widget_s = QDoubleSpinBox()
        self.widget_s.setMinimum(0.03)
        self.widget_s.setMaximum(self.globvar.Bz / 2)
        self.widget_s.setDecimals(3)
        self.widget_s.setValue(self.globvar.s)

        self.widget_s.setSingleStep(10. / 1000.)
        self.widget_s.valueChanged.connect(self.set_s)
        self.widget_s.valueChanged.connect(self.update_topology_plot)


        # cover of reinforcement
        self.widget_cover = QDoubleSpinBox()
        self.widget_cover.setMinimum(0.)
        self.widget_cover.setMaximum(0.1)
        self.widget_cover.setDecimals(3)
        self.widget_cover.setValue( self.globvar.cover )

        self.widget_cover.setSingleStep(5./1000.)
        self.widget_cover.valueChanged.connect(self.set_cover)
        self.widget_cover.valueChanged.connect(self.update_topology_plot)


        # axial diameter of small spiral
        self.widget_dS = QDoubleSpinBox()

        self.widget_dS.setMinimum(self.globvar.minimum_radius*2.)
        self.widget_dS.setMaximum(self.globvar.Bx /2. - self.globvar.cover - (self.globvar.rebars_CS[self.globvar.DS]).diam)
        self.widget_dS.setValue( self.globvar.dS )

        self.widget_dS.setSingleStep(10./1000.)
        self.widget_dS.valueChanged.connect(self.set_dS)
        self.widget_dS.valueChanged.connect(self.update_topology_plot)

        # rebar diameter of small spirals
        self.widget_DS = QComboBox()
        self.widget_DS.addItems(self.globvar.rebars_CS.keys())
        self.widget_DS.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DS))
        self.widget_DS.currentTextChanged.connect(self.set_DS)
        self.widget_DS.currentTextChanged.connect(self.update_topology_plot)

        # rebar diameter of large spiral(s)
        self.widget_DL = QComboBox()
        self.widget_DL.addItems(self.globvar.rebars_CS.keys())
        self.widget_DL.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DL))
        self.widget_DL.currentTextChanged.connect(self.set_DL)
        self.widget_DL.currentTextChanged.connect(self.update_topology_plot)

        # rebar diameter of vertical rebars
        self.widget_DV = QComboBox()
        self.widget_DV.addItems(self.globvar.rebars_CS.keys())
        self.widget_DV.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DV))
        self.widget_DV.currentTextChanged.connect(self.set_DV)
        self.widget_DV.currentTextChanged.connect(self.update_topology_plot)

        # load plate width
        self.widget_load_plate = QDoubleSpinBox()
        self.widget_load_plate.setMinimum(0.25)
        self.widget_load_plate.setMaximum(1.0)
        self.widget_load_plate.setDecimals(2)
        self.widget_load_plate.setValue(self.globvar.load_plate_width)

        self.widget_load_plate.setSingleStep(5. / 100.)
        self.widget_load_plate.valueChanged.connect(self.set_load_plate_width)
        self.widget_load_plate.valueChanged.connect(self.update_topology_plot)

        self.widget_show_load_plate = QCheckBox()
        self.widget_show_load_plate.setChecked(self.globvar.flag_show_Bp)
        self.widget_show_load_plate.toggled.connect(self.set_show_Bp)


        ### TIE

        # column breadth
        self.widget_Bx_TIE = QDoubleSpinBox()
        self.widget_Bx_TIE.setMinimum(0.2)
        self.widget_Bx_TIE.setMaximum(1.0)
        self.widget_Bx_TIE.setDecimals(3)
        self.widget_Bx_TIE.setValue(self.globvar.Bx)

        self.widget_Bx_TIE.setSingleStep(50. / 1000.)
        self.widget_Bx_TIE.valueChanged.connect(self.set_Bx)
        self.widget_Bx_TIE.valueChanged.connect(self.update_topology_plot)


        # ties spacing
        self.widget_s_TIE = QDoubleSpinBox()
        self.widget_s_TIE.setMinimum(0.03)
        self.widget_s_TIE.setMaximum(self.globvar.Bz / 2)
        self.widget_s_TIE.setDecimals(3)
        self.widget_s_TIE.setValue(self.globvar.s)

        self.widget_s_TIE.setSingleStep(10. / 1000.)
        self.widget_s_TIE.valueChanged.connect(self.set_s)
        self.widget_s_TIE.valueChanged.connect(self.update_topology_plot)

        # ties reinforcement cover
        self.widget_cover_TIE = QDoubleSpinBox()
        self.widget_cover_TIE.setMinimum(0.)
        self.widget_cover_TIE.setMaximum(0.1)
        self.widget_cover_TIE.setDecimals(3)
        self.widget_cover_TIE.setValue(self.globvar.cover)

        self.widget_cover_TIE.setSingleStep(5. / 1000.)
        self.widget_cover_TIE.valueChanged.connect(self.set_cover)
        self.widget_cover_TIE.valueChanged.connect(self.update_topology_plot)

        # number of vertical rebars
        self.widget_n_v_bars = QDoubleSpinBox()
        self.widget_n_v_bars.setMinimum(4)
        self.widget_n_v_bars.setDecimals(0)
        self.widget_n_v_bars.setValue(self.globvar.n_v_bars)

        self.widget_n_v_bars.setSingleStep(4)
        self.widget_n_v_bars.valueChanged.connect(self.set_n_v_bars)
        self.widget_n_v_bars.valueChanged.connect(self.update_topology_plot)

        # rebar diameter of vertical rebars
        self.widget_DV_TIE = QComboBox()
        self.widget_DV_TIE.addItems(self.globvar.rebars_CS.keys())
        self.widget_DV_TIE.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DV))
        self.widget_DV_TIE.currentTextChanged.connect(self.set_DV)
        self.widget_DV_TIE.currentTextChanged.connect(self.update_topology_plot)

        # rebar diameter of large tie(s)
        self.widget_DT = QComboBox()
        self.widget_DT.addItems(self.globvar.rebars_CS.keys())
        self.widget_DT.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DT))
        self.widget_DT.currentTextChanged.connect(self.set_DT)
        self.widget_DT.currentTextChanged.connect(self.update_topology_plot)

        # load plate width
        self.widget_load_plate_tie = QDoubleSpinBox()
        self.widget_load_plate_tie.setMinimum(0.25)
        self.widget_load_plate_tie.setMaximum(1.0)
        self.widget_load_plate_tie.setDecimals(2)
        self.widget_load_plate_tie.setValue(self.globvar.load_plate_width)

        self.widget_load_plate_tie.setSingleStep(5. / 100.)
        self.widget_load_plate_tie.valueChanged.connect(self.set_load_plate_width)
        self.widget_load_plate_tie.valueChanged.connect(self.update_topology_plot)

        self.widget_show_load_plate_tie = QCheckBox()
        self.widget_show_load_plate_tie.setChecked(self.globvar.flag_show_Bp)
        self.widget_show_load_plate_tie.toggled.connect(self.set_show_Bp)

        #########################
        # TOPOLOGY & DEFINITION CANVAS
        #########################

        self.topology_canvas = MplCanvas(self, width=5, height=4, dpi=100)


       # MSR layout

        topology_specs_layout_msr = QGridLayout()

        row = 0
        topology_specs_layout_msr.addWidget(QLabel("Geometry"), row, 0, 1, 1)

        row += 1
        label_Bx = QLabel("Bx = ")
        label_Bx_unit = QLabel("[m]")
        topology_specs_layout_msr.addWidget(label_Bx, row, 0, 1, 1, Qt.AlignRight)
        topology_specs_layout_msr.addWidget(self.widget_Bx, row, 1)
        topology_specs_layout_msr.addWidget(label_Bx_unit, row, 2)

        row += 1
        label_cover = QLabel("c = ")
        label_cover_unit = QLabel("[m]")
        topology_specs_layout_msr.addWidget(label_cover, row, 0, 1, 1, Qt.AlignRight)
        topology_specs_layout_msr.addWidget(self.widget_cover, row, 1)
        topology_specs_layout_msr.addWidget(label_cover_unit, row, 2)

        row += 1
        label_H = QLabel("s = ")
        label_H_unit = QLabel("[m]")
        topology_specs_layout_msr.addWidget(label_H, row, 0, 1, 1, Qt.AlignRight)
        topology_specs_layout_msr.addWidget(self.widget_s, row, 1)
        topology_specs_layout_msr.addWidget(label_H_unit, row, 2)

        row += 1
        label_dS = QLabel("dS = ")
        label_dS_unit = QLabel("[m]")
        topology_specs_layout_msr.addWidget(label_dS, row, 0, 1, 1, Qt.AlignRight)
        topology_specs_layout_msr.addWidget(self.widget_dS, row, 1)
        topology_specs_layout_msr.addWidget(label_dS_unit, row, 2)


        row = 0
        topology_specs_layout_msr.addWidget(QLabel("Reinforcement"), row, 3, 1, 3)

        row += 1
        label_DS = QLabel("DS = ")
        label_DS_unit = QLabel("[m]")
        topology_specs_layout_msr.addWidget(label_DS, row, 3, 1, 1, Qt.AlignRight)
        topology_specs_layout_msr.addWidget(self.widget_DS, row, 4)
        topology_specs_layout_msr.addWidget(label_DS_unit, row, 5)

        row += 1
        label_DL = QLabel("DL = ")
        label_DL_unit = QLabel("[m]")
        topology_specs_layout_msr.addWidget(label_DL, row, 3, 1, 1, Qt.AlignRight)
        topology_specs_layout_msr.addWidget(self.widget_DL, row, 4)
        topology_specs_layout_msr.addWidget(label_DL_unit, row, 5)

        row += 1
        label_DV = QLabel("DV = ")
        label_DV_unit = QLabel("[m]")
        topology_specs_layout_msr.addWidget(label_DV, row, 3, 1, 1, Qt.AlignRight)
        topology_specs_layout_msr.addWidget(self.widget_DV, row, 4)
        topology_specs_layout_msr.addWidget(label_DV_unit, row, 5)


        row = 5

        label_load_plate = QLabel("Load plate width / Bx")
        self.label_load_plate_real = QLabel()
        label_show_load_plate = QLabel("Show")
        topology_specs_layout_msr.addWidget(label_load_plate, row, 0, 1, 3,  Qt.AlignLeft)
        topology_specs_layout_msr.addWidget(self.widget_load_plate, row, 2, Qt.AlignRight)
        topology_specs_layout_msr.addWidget(self.label_load_plate_real, row, 3, 1, 2, Qt.AlignLeft)
        topology_specs_layout_msr.addWidget(label_show_load_plate, row, 4,Qt.AlignRight)
        topology_specs_layout_msr.addWidget(self.widget_show_load_plate, row, 5, 1, 1, Qt.AlignLeft)


        row += 1
        topology_specs_layout_msr.addWidget(self.label_info_topology, row, 0, 1, 6, Qt.AlignLeft)

        # TIE

        topology_specs_layout_tie = QGridLayout()

        row = 0
        self.s_min_label = QLabel("")
        topology_specs_layout_tie.addWidget(self.s_min_label, row, 0, 1, 6)
        self.s_min_label.setAlignment(Qt.AlignCenter)

        row = 1
        topology_specs_layout_tie.addWidget(QLabel("Geometry"), row, 0, 1, 1)

        row += 1
        label_Bx = QLabel("Bx = ")
        label_Bx_unit = QLabel("[m]")
        topology_specs_layout_tie.addWidget(label_Bx, row, 0, 1, 1, Qt.AlignRight)
        topology_specs_layout_tie.addWidget(self.widget_Bx_TIE, row, 1)
        topology_specs_layout_tie.addWidget(label_Bx_unit, row, 2)

        row += 1
        label_cover = QLabel("c = ")
        label_cover_unit = QLabel("[m]")
        topology_specs_layout_tie.addWidget(label_cover, row, 0, 1, 1, Qt.AlignRight)
        topology_specs_layout_tie.addWidget(self.widget_cover_TIE, row, 1)
        topology_specs_layout_tie.addWidget(label_cover_unit, row, 2)

        row += 1
        label_H = QLabel("s = ")
        label_H_unit = QLabel("[m]")
        topology_specs_layout_tie.addWidget(label_H, row, 0, 1, 1, Qt.AlignRight)
        topology_specs_layout_tie.addWidget(self.widget_s_TIE, row, 1)
        topology_specs_layout_tie.addWidget(label_H_unit, row, 2)

        row = 1
        topology_specs_layout_tie.addWidget(QLabel("Reinforcement"), row, 3, 1, 3)

        row += 1
        label_DT = QLabel("DT = ")
        label_DT_unit = QLabel("[m]")
        topology_specs_layout_tie.addWidget(label_DT, row, 3, 1, 1, Qt.AlignRight)
        topology_specs_layout_tie.addWidget(self.widget_DT, row, 4)
        topology_specs_layout_tie.addWidget(label_DT_unit, row, 5)

        row += 1
        label_DV = QLabel("DV = ")
        label_DV_unit = QLabel("[m]")
        topology_specs_layout_tie.addWidget(label_DV, row, 3, 1, 1, Qt.AlignRight)
        topology_specs_layout_tie.addWidget(self.widget_DV_TIE, row, 4)
        topology_specs_layout_tie.addWidget(label_DV_unit, row, 5)

        row += 1
        label_n_v_bars = QLabel("# V-Bars= ")
        label_n_v_bars_unit = QLabel("[-]")
        topology_specs_layout_tie.addWidget(label_n_v_bars, row, 3, 1, 1, Qt.AlignRight)
        topology_specs_layout_tie.addWidget(self.widget_n_v_bars, row, 4)
        topology_specs_layout_tie.addWidget(label_n_v_bars_unit, row, 5, 1, 1, Qt.AlignLeft)

        row += 1

        label_load_plate_tie = QLabel("Load plate width / Bx")
        self.label_load_plate_real_tie = QLabel()
        label_show_load_plate_tie = QLabel("Show")
        topology_specs_layout_tie.addWidget(label_load_plate_tie, row, 0, 1, 3, Qt.AlignLeft)
        topology_specs_layout_tie.addWidget(self.widget_load_plate_tie, row, 2, Qt.AlignRight)
        topology_specs_layout_tie.addWidget(self.label_load_plate_real_tie, row, 3, 1, 2, Qt.AlignLeft)
        topology_specs_layout_tie.addWidget(label_show_load_plate_tie, row, 4, Qt.AlignRight)
        topology_specs_layout_tie.addWidget(self.widget_show_load_plate_tie, row, 5, 1, 1, Qt.AlignLeft)



        row = 5
        row += 1
        topology_specs_layout_tie.addWidget(self.label_info_topology_tie, row, 0, 1, 6, Qt.AlignLeft)


        ### MAIN TOPOLOGY LAYOUT

        # layout switcher
        self.topology_switch = QComboBox()
        self.topology_switch.addItems(["MSR", "TIEs"])
        self.topology_switch.setCurrentIndex(self.globvar.topology_layout)
        self.topology_switch.currentIndexChanged.connect(self.topology_set)
        self.topology_switch.currentIndexChanged.connect(self.update_topology_plot)


        # turning layouts to widgets
        t_layout_1 = QWidget()
        t_layout_1.setLayout(topology_specs_layout_msr)

        t_layout_2 = QWidget()
        t_layout_2.setLayout(topology_specs_layout_tie)

        #stacket layout
        self.topology_specs_layout = QStackedLayout()
        self.topology_specs_layout.addWidget(t_layout_1)
        self.topology_specs_layout.addWidget(t_layout_2)

        self.topology_specs_layout_widget = QWidget()
        self.topology_specs_layout_widget.setLayout(self.topology_specs_layout)

        #main layout
        topology_layout = QVBoxLayout()
        topology_layout.addWidget(self.topology_switch)
        #topology_layout.addWidget(topology_switch_layout_widget)
        topology_layout.addWidget(self.topology_specs_layout_widget)

        #########################
        # INTERACTION DIAGRAMS etc.
        #########################

        ### basic material properties
        # concrete
        self.widget_fcm = QComboBox()
        self.widget_fcm.addItems(self.globvar.concretes.keys())
        self.widget_fcm.setCurrentIndex(list(self.globvar.concretes.keys()).index(self.globvar.fcm))
        self.widget_fcm.currentTextChanged.connect(self.set_fcm)
        self.widget_fcm.currentTextChanged.connect(self.update_diagram_plot)

        button_CDPM2 = QPushButton(" Concrete properties ")
        button_CDPM2.clicked.connect(self.set_CDPM2)

        # steel - vertical reinforcement
        self.widget_fy_vert = QDoubleSpinBox()
        self.widget_fy_vert.setMinimum(200)
        self.widget_fy_vert.setMaximum(700)
        self.widget_fy_vert.setDecimals(0)
        self.widget_fy_vert.setValue( self.globvar.fy_vert )

        self.widget_fy_vert.setSingleStep(20.)
        self.widget_fy_vert.valueChanged.connect(self.set_fy_vert)
        self.widget_fy_vert.valueChanged.connect(self.update_diagram_plot)

        # steel - lateral reinforcement
        self.widget_fy_lat = QDoubleSpinBox()
        self.widget_fy_lat.setMinimum(200)
        self.widget_fy_lat.setMaximum(700)
        self.widget_fy_lat.setDecimals(0)
        self.widget_fy_lat.setValue( self.globvar.fy_lat )

        self.widget_fy_lat.setSingleStep(20.)
        self.widget_fy_lat.valueChanged.connect(self.set_fy_lat)
        self.widget_fy_lat.valueChanged.connect(self.update_diagram_plot)

        ### checkboxes for interaction diagrams
        checkbox_id_steel = QCheckBox()
        checkbox_id_steel.setCheckState(self.globvar.flag_id_steel)
        checkbox_id_steel.stateChanged.connect(self.set_id_steel)
        checkbox_id_steel.stateChanged.connect(self.update_diagram_plot)

        checkbox_id_ACI = QCheckBox()
        checkbox_id_ACI.setCheckState(self.globvar.flag_id_ACI)
        checkbox_id_ACI.stateChanged.connect(self.set_id_ACI)
        checkbox_id_ACI.stateChanged.connect(self.update_diagram_plot)

        checkbox_id_ASCE = QCheckBox()
        checkbox_id_ASCE.setCheckState(self.globvar.flag_id_ASCE)
        checkbox_id_ASCE.stateChanged.connect(self.set_id_ASCE)
        checkbox_id_ASCE.stateChanged.connect(self.update_diagram_plot)

        checkbox_id_CTU = QCheckBox()
        checkbox_id_CTU.setCheckState(self.globvar.flag_id_CTU)
        checkbox_id_CTU.stateChanged.connect(self.set_id_CTU)
        checkbox_id_CTU.stateChanged.connect(self.update_diagram_plot)

        checkbox_id_CTU_concrete = QCheckBox()
        checkbox_id_CTU_concrete.setCheckState(self.globvar.flag_id_CTU_concrete)
        checkbox_id_CTU_concrete.stateChanged.connect(self.set_id_CTU_concrete)
        checkbox_id_CTU_concrete.stateChanged.connect(self.update_diagram_plot)

        checkbox_id_FEM = QCheckBox()
        checkbox_id_FEM.setCheckState(self.globvar.flag_id_FEM)
        checkbox_id_FEM.stateChanged.connect(self.set_id_FEM)
        checkbox_id_FEM.stateChanged.connect(self.update_diagram_plot)

        self.checkbox_id_FEM_paths = QCheckBox()
        self.checkbox_id_FEM_paths.setCheckState(self.globvar.flag_id_FEM_paths)
        self.checkbox_id_FEM_paths.stateChanged.connect(self.set_id_FEM_paths)
        self.checkbox_id_FEM_paths.stateChanged.connect(self.update_diagram_plot)

        # solution for neutral axis visible
        checkbox_active_neutral_axis = QCheckBox()
        checkbox_active_neutral_axis.setCheckState(self.globvar.flag_active_neutral_axis)
        checkbox_active_neutral_axis.stateChanged.connect(self.set_active_neutral_axis)
        checkbox_active_neutral_axis.stateChanged.connect(self.update_topology_plot)
        checkbox_active_neutral_axis.stateChanged.connect(self.update_diagram_plot)

        # display loading points
        self.checkbox_show_loading = QCheckBox()
        self.checkbox_show_loading.setCheckState(self.globvar.flag_loading_display)
        self.checkbox_show_loading.stateChanged.connect(self.set_loading_display)
        self.checkbox_show_loading.stateChanged.connect(self.update_diagram_plot)

        #########################
        # INTERACTION DIAGRAM CANVAS
        #########################

        self.diagram_canvas = MplCanvas(self, width=6, height=4, dpi=100)

        id_specs_layout = QGridLayout()

        row = 0

        row += 1
        self.id_switcher = QComboBox()
        self.id_switcher.addItems([" N / M diagram ", " fb / d diagram "])
        self.id_switcher.setCurrentIndex(self.globvar.id_plot)
        self.id_switcher.currentIndexChanged.connect(self.set_id_plot)
        #id_specs_layout.addWidget(self.id_switcher, row, 2, 1, 2, Qt.AlignRight)
        id_specs_layout.addWidget(self.id_switcher, row, 1, 1, 3)

        row += 1
        label_id_materials = QLabel("Materials")
        font_id = label_id_materials.font()
        font_id.setPointSize(13)
        label_id_materials.setFont(font_id)
        id_specs_layout.addWidget(label_id_materials, row, 0, 1, 4, Qt.AlignLeft)

        row += 1
        label_fy_lat = QLabel("Yield stress of reinforcement")
        label_fy_lat_unit = QLabel("[MPa]")
        id_specs_layout.addWidget(label_fy_lat, row, 0, 1, 3)
        id_specs_layout.addWidget(self.widget_fy_lat, row, 2, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(label_fy_lat_unit, row, 3, Qt.AlignLeft)

        row += 1
        label_fcm = QLabel("Concrete class (characteristic strength)")
        label_fcm_unit = QLabel("[MPa]")
        id_specs_layout.addWidget(label_fcm, row, 0, 1, 3, Qt.AlignLeft)
        id_specs_layout.addWidget(self.widget_fcm, row, 2, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(label_fcm_unit, row, 3, Qt.AlignLeft)

        row+= 1
        id_specs_layout.addWidget(button_CDPM2, row, 2, 1, 2, Qt.AlignRight)

        row += 1
        id_line =  QFrame()
        id_line.setFrameShape(QFrame.HLine)
        id_line.setStyleSheet("color: grey")
        id_specs_layout.addWidget(id_line, row, 0, 1, 4)

        row += 1
        label_id_results = QLabel("Results")
        label_id_results.setFont(font_id)
        id_specs_layout.addWidget(label_id_results, row, 0, 1, 4, Qt.AlignLeft)

        row += 1
        label_id_FEM = QLabel("FEM: ")
        font_id_results = label_id_FEM.font()
        font_id_results.setPointSize(11)
        label_id_FEM.setFont(font_id_results)
        id_specs_layout.addWidget(label_id_FEM, row, 0, 1, 2, Qt.AlignLeft)

        row += 1
        id_specs_layout.addWidget(QLabel("Diagram"), row, 0, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(checkbox_id_FEM, row, 1)
        id_specs_layout.addWidget(QLabel("Loading path"), row, 2, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(self.checkbox_id_FEM_paths, row, 3)

        row += 1
        label_id_standards = QLabel("Standards: ")
        label_id_standards.setFont(font_id_results)
        id_specs_layout.addWidget(label_id_standards, row, 0, 1, 2, Qt.AlignLeft)

        row += 1
        id_specs_layout.addWidget(QLabel("ACI"), row, 0, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(checkbox_id_ACI, row, 1)

        id_specs_layout.addWidget(QLabel("ASCE"), row, 2, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(checkbox_id_ASCE, row, 3)

        #row += 1
        #id_specs_layout.addWidget(QLabel("Fcm"), row, 0, 1, 1, Qt.AlignRight)
        #id_specs_layout.addWidget(checkbox_id_concrete, row, 1)

        row += 1
        label_id_CTU = QLabel("Analytical Solution: ")
        label_id_CTU.setFont(font_id_results)
        id_specs_layout.addWidget(label_id_CTU, row, 0, 1, 1, Qt.AlignLeft)
        row += 1
        id_specs_layout.addWidget(QLabel("Plain\nconcrete"), row, 0, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(checkbox_id_CTU_concrete, row, 1)

        id_specs_layout.addWidget(QLabel("Reinforced\nconcrete"), row, 2, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(checkbox_id_CTU, row, 3)

        #########################
        # FEM etc.
        #########################

        # finite element mesh
        self.widget_dx = QDoubleSpinBox()
        self.widget_dx.setMinimum( self.globvar.Bx/200. )
        self.widget_dx.setMaximum( self.globvar.Bx/8. )
        self.widget_dx.setDecimals(4)
        self.widget_dx.setSingleStep(5./1000.)
        self.widget_dx.setValue( self.globvar.elem_size_X )
        self.widget_dx.valueChanged.connect( self.set_dx )


        self.widget_dy = QDoubleSpinBox()
        self.widget_dy.setMinimum( self.globvar.By/200. )
        self.widget_dy.setMaximum( self.globvar.By/8. )
        self.widget_dy.setDecimals(4)
        self.widget_dy.setSingleStep(5./1000.)
        self.widget_dy.setValue( self.globvar.elem_size_YZ )
        self.widget_dy.valueChanged.connect( self.set_dy )

        self.checkbox_ignore_cover = QCheckBox()
        self.checkbox_ignore_cover.setCheckState(self.globvar.flag_ignore_cover)
        self.checkbox_ignore_cover.stateChanged.connect(self.set_ignore_cover)


        self.show_concrete_button = QPushButton("Concrete")
        self.show_concrete_button.clicked.connect(self.show_concrete_mesh)

        self.show_rebars_button = QPushButton("Rebars")
        self.show_rebars_button.clicked.connect(self.show_rebars_mesh)

        self.button_generate_oofem_input = QPushButton("Generate oofem input(s)")
        self.button_generate_oofem_input.clicked.connect(self.generate_inputs)

        button_run_oofem = QPushButton("Run FEM analysis")
        button_run_oofem.clicked.connect(self.run_oofem_problems )

        #########################
        # FEM MODEL LAYOUT
        #########################

        fem_model_layout = QGridLayout()

        row = 0

        row += 1
        label_elem_size = QLabel("Elements size:")
        fem_model_layout.addWidget(label_elem_size, row, 0, 1, 3, Qt.AlignLeft)

        #row += 1
        label_dx = QLabel("dx ≈ ")
        label_dx_unit = QLabel("[m]")
        fem_model_layout.addWidget(label_dx, row, 1, 1, 1, Qt.AlignRight)
        fem_model_layout.addWidget(self.widget_dx, row, 2)
        fem_model_layout.addWidget(label_dx_unit, row, 3, Qt.AlignLeft)

        row += 1
        label_dy = QLabel("dy, dz ≈ ")
        label_dy_unit = QLabel("[m]")
        fem_model_layout.addWidget(label_dy, row, 1, 1, 1, Qt.AlignRight)
        fem_model_layout.addWidget(self.widget_dy, row, 2)
        fem_model_layout.addWidget(label_dy_unit, row, 3, Qt.AlignLeft)

        '''
        row += 1
        label_ignore_cover = QLabel("Ignore concrete cover")
        fem_model_layout.addWidget(label_ignore_cover, row, 0, 1, 2, Qt.AlignRight)
        fem_model_layout.addWidget(self.checkbox_ignore_cover, row, 2)
        '''

        row += 1
        fem_model_layout.addWidget(QLabel("Generate and show mesh:"), row, 0, 1, 2, Qt.AlignLeft)
        fem_model_layout.addWidget(self.show_concrete_button, row, 2, 1, 1,  Qt.AlignRight)
        fem_model_layout.addWidget(self.show_rebars_button, row, 3, 1, 1,  Qt.AlignLeft)

        row += 1
        label_project_name = QLabel("Project name:")
        self.lineEditProject = QLineEdit()
        self.lineEditProject.setMaxLength(50)
        self.lineEditProject.setText(self.globvar.project_name)
        self.lineEditProject.textChanged.connect(self.set_project_name)

        fem_model_layout.addWidget(label_project_name, row, 0)
        fem_model_layout.addWidget(self.lineEditProject, row, 1, 1, 3)

        row += 1
        lineEditOOFEM = QLineEdit()
        lineEditOOFEM.setMaxLength(200)
        lineEditOOFEM.setPlaceholderText("specify OOFEM folder")
        lineEditOOFEM.textChanged.connect(self.set_oofem_folder)

        button_oofem_folder = QPushButton("OOFEM folder")
        button_oofem_folder.clicked.connect(self.select_oofem_folder)
        button_oofem_folder.clicked.connect( lambda: lineEditOOFEM.setText(self.globvar.oofem_folder) )

        fem_model_layout.addWidget(button_oofem_folder, row, 0)
        fem_model_layout.addWidget(lineEditOOFEM, row, 1, 1, 3)

        #########################
        # FEM CHECKBOXES LAYOUT
        #########################

        self.fem_checkboxes_layout = QGridLayout()

        self.fem_checkboxes = []

        row = 0
        label_checkboxes = QLabel("Normalized loading plate depth selection:")
        self.fem_checkboxes_layout.addWidget(label_checkboxes, row, 0, 1, self.globvar.ecc_nr)


        row += 1
        for i in range(self.globvar.ecc_nr):
            label_i = QLabel(f"{( self.globvar.tasks[i].eccentricity_normalized[1] ):.3f}" )
            self.fem_checkboxes_layout.addWidget(label_i, row, i)

        row += 1
        for i in range(self.globvar.ecc_nr):

            if (self.globvar.tasks[i].status == Task_status.PREDEFINED ):
                continue

            checkbox = QCheckBox()
            if ( self.globvar.tasks[i].status == Task_status.UNSELECTED):
                checkbox.setCheckState( Qt.Unchecked)
            else:
                checkbox.setCheckState( Qt.Checked)

            checkbox.stateChanged.connect( self.update_ecc_selection )
            self.fem_checkboxes.append(checkbox)


        for checkbox in self.fem_checkboxes:
            index = self.fem_checkboxes.index(checkbox)
            #print(index)
            self.fem_checkboxes_layout.addWidget(checkbox, row, index, Qt.AlignCenter)

        #########################
        # FEM ANALYSIS LAYOUT
        #########################

        self.fem_analysis_layout = QGridLayout()

        row = 0
        self.fem_analysis_layout.addWidget(self.button_generate_oofem_input, row, 0, 1, 4)

        row += 1
        self.fem_analysis_layout.addWidget(button_run_oofem, row, 0, 1, 1)

        cpu_label= QLabel("Threads # (less than " + str(multiprocessing.cpu_count()) + "):" )
        cpu_label.setAlignment(Qt.AlignVCenter)

        cpu_nr = QDoubleSpinBox()
        cpu_nr.setDecimals(0)
        recommended_cpu = math.ceil(multiprocessing.cpu_count()*3./4.)
        self.set_cpu_nr( recommended_cpu )
        cpu_nr.setValue( recommended_cpu )
        cpu_nr.setMinimum( 1 )
        cpu_nr.setMaximum( multiprocessing.cpu_count() )

        cpu_nr.valueChanged.connect(self.set_cpu_nr)

        self.fem_analysis_layout.addWidget(cpu_label, row, 1, 1, 2)
        self.fem_analysis_layout.addWidget(cpu_nr, row, 3)

        row += 1

        label_console = QLabel("Analysis results:")
        self.fem_analysis_layout.addWidget(label_console, row, 0, 1, 4)

        row += 1
        self.console = QPlainTextEdit()
        self.console.setReadOnly(True)
        self.fem_analysis_layout.addWidget(self.console, row, 0, 1, 4)

        row += 1
        progress_label = QLabel("Load vs. estimate:")
        self.fem_analysis_layout.addWidget(progress_label, row, 0, 1, 1)
        self.progressBar = QProgressBar()
        self.progressBar.setMinimum(0)
        self.progressBar.setMaximum(200)
        self.progressBar.setFormat("%v%")
        self.fem_analysis_layout.addWidget(self.progressBar, row, 1, 1, 3)

        fem_group_layout = QVBoxLayout()
        fem_group_layout.addLayout(fem_model_layout)
        fem_group_layout.addLayout(self.fem_checkboxes_layout)
        fem_group_layout.addLayout(self.fem_analysis_layout)

        #########################
        # MAIN LAYOUT
        #########################

        main_layout = QGridLayout()

        label_L = QLabel("Topology Specification")
        label_L.setAlignment(Qt.AlignCenter)
        font_heading = label_L.font()
        font_heading.setPointSize(16)
        font_heading.setBold(True)
        label_L.setFont(font_heading)
        main_layout.addWidget(label_L, 0, 0)
        main_layout.addWidget(self.topology_canvas, 1, 0)
        main_layout.addLayout(topology_layout, 2, 0)

        label_C = QLabel("Materials and Results")
        label_C.setAlignment(Qt.AlignCenter)
        label_C.setFont(font_heading)

        main_layout.addWidget(label_C, 0, 1)
        main_layout.addWidget(self.diagram_canvas, 1, 1)
        main_layout.addLayout(id_specs_layout, 2, 1)

        label_R = QLabel("FEM Definition")
        label_R.setAlignment(Qt.AlignCenter)
        label_R.setFont(font_heading)
        main_layout.addWidget(label_R, 0, 2)
        main_layout.addLayout(fem_group_layout, 1, 2, 2, 1)

        self.update_topology_plot()
        self.update_diagram_plot()

        # TODO: improve
        self.globvar.sets = max( len(self.globvar.materials), len(self.globvar.cross_sections) )

        widget = QWidget()
        widget.setLayout(main_layout)

        self.setCentralWidget(widget)


        # prevent maximization and rescaling
        #self.setFixedSize(1280, 720)

        #https://stackoverflow.com/questions/60815433/deprecationwarning-qdesktopwidget-availablegeometryint-screen-const-is-deprec
        '''
        self.setGeometry(
            QStyle.alignedRect(
            Qt.LeftToRight,
            Qt.AlignCenter,
            self.size(),
            QGuiApplication.primaryScreen().availableGeometry(),
            ),
        )
        '''

        #self.show()


    # closes all windows when the main is closed
    def closeEvent(self, event):
        for window in QApplication.topLevelWidgets():
            window.close()
        #if self.w:
        #    self.w.close()

    def define_CDPM2(self, cdpm2):
        self.w_cdpm2 = window_CDPM2(cdpm2)
        self.w_cdpm2.show()
        self.w_cdpm2.setGeometry(
            QStyle.alignedRect(
            Qt.LeftToRight,
            Qt.AlignCenter,
            self.w_cdpm2.size(),
            QGuiApplication.primaryScreen().availableGeometry(),
            ),
        )
    def draw_segment(self, p0, p1, diameter, axes, color=(1.0, 0.5, 0.0, 1.0 )): #Draws a rectangular rebar segment between p0 and p1.
        """parameters:
            p0, p1 : list or array
                2D coordinates of segment endpoints
            diameter : float
                Thickness of the rebar (height of rectangle)
            axes : matplotlib axes
                Target axes for drawing
            color : RGBA tuple
                Fill color of the rectangle
        """
        p0 = np.array(p0[:2])   #take only x and y for 2D
        p1 = np.array(p1[:2])

        # Segment center
        center = (p0 + p1) / 2

        # Length and angle of the segment
        dx, dy = p1 - p0
        length = np.sqrt(dx ** 2 + dy ** 2)
        angle_rad = np.arctan2(dy, dx)

        # Create a horizontal rectangle centered at origin
        rect = Rectangle(
            (-length / 2, -diameter / 2),
            width=length,
            height=diameter,
            linewidth=1.0,
            edgecolor=(0.0 , 0.0 , 0.0 , 1.0),
            facecolor=color
        )

        # Apply rotation and translation
        transform = (
                Affine2D()
                .rotate(angle_rad)
                .translate(center[0], center[1])
                + axes.transData
        )

        rect.set_transform(transform)
        axes.add_patch(rect)


    def define_Mises(self, mises):
        self.w_mises = window_Mises(mises)
        self.w_mises.show()
        self.w_mises.setGeometry(
            QStyle.alignedRect(
            Qt.LeftToRight,
            Qt.AlignCenter,
            self.w_mises.size(),
            QGuiApplication.primaryScreen().availableGeometry(),
            ),
        )



    def warning_dialog_ok(self,message):
         button = QMessageBox.warning(self, "Warning", message, buttons=QMessageBox.Ok, defaultButton=QMessageBox.Ok)


    def info_dialog(self,message):
         button = QMessageBox.information(self, "Info", message, buttons=QMessageBox.Ok, defaultButton=QMessageBox.Ok)



    def update_topology_information(self):

        rho_vert = self.compute_vertical_reinforcement_ratio()
        self.globvar.rho_vert = rho_vert

        rho_lat = self.compute_lateral_reinforcement_ratio()
        self.globvar.rho_lat = rho_lat

        rho_lat_tie = self.compute_lateral_reinforcement_ratio_tie()
        self.globvar.rho_lat_tie = rho_lat_tie

        Aconf = self.compute_total_confined_area()
        Aconf_ties = ( self.globvar.Bx - self.globvar.cover - self.globvar.rebars_CS[self.globvar.DT].diam / 2 ) **2
        Aconf_ties = Aconf_ties / 3
        A = self.globvar.Bx * self.globvar.By

        sigL_L = 0.
        sigL_S = 0.

        ke_L = 0.
        ke_S = 0.

        as_v_1 = math.pi * ((self.globvar.rebars_CS[self.globvar.DV]).diam **2) /4
        if as_v_1 == 0:
            min_v_r = 0.
            max_v_r = 0.
        else:
            min_v_r = math.ceil((0.01 * A)/as_v_1)
            max_v_r = math.floor((0.06 * A)/as_v_1)

        for rebar in self.globvar.rebars:
            if ( rebar.give_rebar_type() == 'spiral' and rebar.tag[0] == 'S' ):
                sigL_S = rebar.compute_confinement( self.globvar, self.globvar.fy_lat )
                ke_S = rebar.compute_confinement_effectiveness()
                break

        for rebar in self.globvar.rebars:
            if ( rebar.give_rebar_type() == 'spiral' and rebar.tag[0] == 'L' ):
                sigL_L = rebar.compute_confinement( self.globvar, self.globvar.fy_lat )
                ke_L = rebar.compute_confinement_effectiveness()
                break

        self.label_info_topology.setText(f"Reinforcement ratios: {(rho_lat*100.):.3f}% (lateral) \nConfinement: sigL_L = {(sigL_L):.3f} MPa, sigL_S = {(sigL_S):.3f} MPa \nEffectiveness: ke_L = {(ke_L):.3f}, ke_S = {(ke_S):.3f}\n"
                                         f"Sigma eff.: sigL_eff_L = {(sigL_L * ke_L):.3f} MPa, sigL_eff_S = {(sigL_S * ke_S):.3f} MPa \nArea = {(A):.3f}m^2, A_eff = {(Aconf):.3f}m^2, A_eff/A = {(Aconf/A):.3f} [-]")
        self.label_info_topology_tie.setText(f"Reinforcement ratios: {(rho_lat_tie * 100.):.3f}% (lateral) \nACI requirements (1% ~ 6%): {(min_v_r)} ~ {(max_v_r)} # of vertical bars \nArea = {(A):.3f} m^2, A_eff = {(Aconf_ties):.3f} m^2, A_eff/A = {(Aconf_ties / A):.3f} [-]")
        font_info = self.label_info_topology.font()
        font_info.setPointSize(11)
        self.label_info_topology.setFont(font_info)
        self.label_info_topology_tie.setFont(font_info)

        Bp = self.globvar.Bx * round(self.globvar.load_plate_width, 2)
        self.label_load_plate_real.setText(f"[-] ≈ {round(Bp, 3)} m")
        self.label_load_plate_real_tie.setText(f"[-] ≈ {round(Bp, 3)} m")

    def update_topology_plot(self):

        self.update_topology_information()

        logger.info("Topology updated")

        # all ids to be updated
        self.need_update_ids_flag = True

        self.topology_canvas.axes.cla()  # Clear the canvas.

        edge_color_opacity = 1 # 0<val<1
        face_color_opacity = 1. # 0<val<1

        face_color_opacity_V = 1.

        self.topology_canvas.axes.add_patch(matplotlib.patches.Rectangle((-self.globvar.Bx/2., -self.globvar.By/2.), self.globvar.Bx, self.globvar.By, edgecolor=(0, 0, 0, edge_color_opacity), facecolor=(192./255.,192./255.,192./255., face_color_opacity), linewidth=2))


        # draw load plate
        if (self.globvar.flag_show_Bp):

            width = self.globvar.Bp
            height = self.globvar.Bx
            x = -width / 2
            y = -height / 2

            hatch_rect = Rectangle((x, y), width, height, facecolor='none', edgecolor='red', hatch='///', linewidth=1.0, zorder=5)
            self.topology_canvas.axes.add_patch(hatch_rect)


        for rebar in self.globvar.rebars:

            reb_type = rebar.give_rebar_type()
            if (reb_type == 'rebar'):
                # get first coordinate, expecting vertical orientation
                XYZ = rebar.XYZ[0]
                xy = [XYZ[0], XYZ[1]]
                radius = rebar.give_diameter(self.globvar)/2.

                circle_path = matplotlib.path.Path.circle(xy, radius, readonly=False)
                circle_patch = matplotlib.patches.PathPatch(circle_path, edgecolor=(0, 0, 0, edge_color_opacity), facecolor=(0.0, 1.0, 0.0, face_color_opacity_V), linewidth=1.)
                self.topology_canvas.axes.add_patch(circle_patch)


            elif (reb_type in ['hoop', 'spiral']) and rebar.tag[0] == 'S' :
                if self.globvar.topology_layout == 0:
                    xy = [rebar.XYZ[0],rebar.XYZ[1]]
                    annulus_patch = matplotlib.patches.Annulus(xy, rebar.radius + rebar.give_diameter(self.globvar)/2., rebar.give_diameter(self.globvar), edgecolor=(0, 0, 0, edge_color_opacity), facecolor=(1., 0.0, 0.0, face_color_opacity), linewidth=1)
                    self.topology_canvas.axes.add_patch(annulus_patch)

            elif (reb_type in ['hoop', 'spiral']) and rebar.tag[0] == 'L':
                if self.globvar.topology_layout == 0:
                    xy = [rebar.XYZ[0],rebar.XYZ[1]]
                    annulus_patch = matplotlib.patches.Annulus(xy, rebar.radius + rebar.give_diameter(self.globvar)/2., rebar.give_diameter(self.globvar), edgecolor=(0, 0, 0, edge_color_opacity), facecolor=(30./255.,144./255.,255./255., face_color_opacity), linewidth=1)
                    self.topology_canvas.axes.add_patch(annulus_patch)

            elif reb_type == 'tie' :
                if self.globvar.topology_layout == 1:
                    # draw the tie from rebar_points
                    points = rebar.rebar_points

                    if len(points) < 2:
                        print(f"⚠️ Not enough points to draw tie: tag={rebar.tag}")
                        continue

                    diam = self.globvar.rebars_CS[self.globvar.DT].diam

                    for i in range(len(points) - 1):
                        self.draw_segment(points[i], points[i + 1], diam, self.topology_canvas.axes)


            else:
                print(f"⚠️ Unknown rebar type: {reb_type}, tag: {rebar.tag}")
                raise ValueError


        #self.topology_canvas.axes.plot(self.xdata, self.ydata, color = 'chartreuse')
        # Trigger the topology_canvas to update and redraw.
        self.topology_canvas.axes.set_aspect('equal')
        self.topology_canvas.axes.autoscale_view()

        # remove border lines
        self.topology_canvas.axes.spines['top'].set_visible(False)
        self.topology_canvas.axes.spines['right'].set_visible(False)
        self.topology_canvas.axes.spines['bottom'].set_visible(True)
        self.topology_canvas.axes.spines['left'].set_visible(True)

        spines_size = 1.6
        self.topology_canvas.axes.spines['bottom'].set_linewidth(spines_size)
        self.topology_canvas.axes.spines['left'].set_linewidth(spines_size)

        self.topology_canvas.axes.xaxis.tick_bottom()
        self.topology_canvas.axes.yaxis.tick_left()

        #arrows
        self.topology_canvas.axes.plot(1, 0, ">k", transform=self.topology_canvas.axes.transAxes, clip_on=False,
                                       markersize=6)
        self.topology_canvas.axes.plot(0, 1, "^k", transform=self.topology_canvas.axes.transAxes, clip_on=False,
                                       markersize=6)
        # xy label
        self.topology_canvas.axes.text(1.04, 0, 'x', transform=self.topology_canvas.axes.transAxes,
                                       fontsize=11, fontweight='bold', va='center', ha='left', clip_on=False)
        self.topology_canvas.axes.text(0, 1.04, 'y', transform=self.topology_canvas.axes.transAxes,
                                       fontsize=11, fontweight='bold', va='bottom', ha='center', clip_on=False)

        self.topology_canvas.draw()
        self.update_diagram_plot()



    def update_diagram_plot(self):

        logger.info("Diagram updated")

        if ( self.globvar.id_plot == 0 ):

            self.update_diagram_plot_N_M()

        elif ( self.globvar.id_plot == 1 ):

            self.update_diagram_plot_fb_d()

        else:
            self.warning_dialog_ok("Unknown canvas id")
            warnings.warn("Unknown canvas id")
            raise ValueError



    def update_diagram_plot_N_M(self):

        self.diagram_canvas.axes.cla()  # Clear the canvas.

        self.diagram_canvas.axes.set_ylabel('N [MN]', fontsize=9)
        self.diagram_canvas.axes.set_xlabel('M [MNm]', fontsize=9)

        self.diagram_canvas.axes.grid(True)

        self.diagram_canvas.axes.autoscale_view()
        self.diagram_canvas.axes.spines['top'].set_visible(False)
        self.diagram_canvas.axes.spines['right'].set_visible(False)

        ### functions to draw the interaction diagrams
        Bp_y_array = np.linspace(0, self.globvar.By, 100)

        #ACI
        if ( self.globvar.flag_id_ACI):
            N_max_ACI, M_max_ACI, f_b, d = ACI(Bp_y_array, self.globvar)
            self.diagram_canvas.axes.plot(M_max_ACI, N_max_ACI, label='ACI', color= 'green', linestyle='solid')

        #ASCE
        if ( self.globvar.flag_id_ASCE):
            N_max_ASCE, M_max_ASCE, f_b, d = ASCE(self.globvar)
            self.diagram_canvas.axes.plot(M_max_ASCE, N_max_ASCE, label='ASCE', color='green', linestyle='none', marker='o',markersize=5)

        #concrete
        if ( self.globvar.flag_id_concrete):
            N_max_concrete, M_max_concrete, f_b, d = id_fcm(Bp_y_array, self.globvar)
            self.diagram_canvas.axes.plot(M_max_concrete, N_max_concrete, label='fcm', color='purple', linestyle='solid')

        #CTU total
        if (self.globvar.flag_id_CTU):

            N_max_AS_MSR, M_max_AS_MSR, f_b_MSR, d_MSR = AS_MSR(Bp_y_array, self.globvar)
            self.diagram_canvas.axes.plot(M_max_AS_MSR, N_max_AS_MSR, label='MSR', color='blue', linestyle='solid')

            N_max_AS_TIE, M_max_AS_TIE, f_b_TIE, d_TIE = AS_TIE(Bp_y_array, self.globvar)
            self.diagram_canvas.axes.plot(M_max_AS_TIE, N_max_AS_TIE, label='Ties', color='blue', linestyle='--')

        if (self.globvar.flag_id_CTU_concrete):
            N_max, M_max, f_b, d = AS_plain_concrete(Bp_y_array, self.globvar)
            self.diagram_canvas.axes.plot(M_max, N_max, label='plain concrete', color='red', linestyle='--')

        ## FEM
        if self.globvar.flag_id_FEM:

            max_M = []
            max_N = []

            for task in self.globvar.tasks:
                #if task.status in (Task_status.COMPLETED, Task_status.PREDEFINED):
                if task.status == Task_status.COMPLETED:
                    max_M.append(task.max_MN[0])
                    max_N.append(task.max_MN[1])

            if self.globvar.flag_id_FEM_paths:
                for task in self.globvar.tasks:
                    if task.status == Task_status.COMPLETED:

                        aux_M = [0.0] + [abs(m) for m in task.M]
                        aux_N = [0.0] + [abs(n) for n in task.N]

                        self.diagram_canvas.axes.plot(
                            aux_M, aux_N,
                            color=(192. / 255., 192. / 255., 192. / 255.),
                            linestyle='solid', linewidth=0.5, marker='x',
                            markeredgecolor='black', markeredgewidth=0.5, markersize=2
                        )
            '''
            # Envelope curve
            if len(max_M) >= 3:
                self.diagram_canvas.axes.plot(max_M, max_N, '-or', markeredgewidth=1.5, markersize=6, label="FEM",
                                              markeredgecolor='black')
            elif len(max_M) > 1:
                self.diagram_canvas.axes.plot(max_M, max_N, 'or', markeredgewidth=1.5, markersize=6, label="FEM",
                                              markeredgecolor='black')
            '''

            if len(max_M) >= 1:
                self.diagram_canvas.axes.plot(max_M, max_N, 'or', markeredgewidth=1.5, markersize=6, label="FEM",
                                            markeredgecolor='black')


        self.diagram_canvas.axes.legend(loc='lower right')
        self.diagram_canvas.figure.tight_layout()
        self.diagram_canvas.draw()

        self.need_update_ids_flag = False


    def update_diagram_plot_fb_d(self):

        self.diagram_canvas.axes.cla()  # Clear the canvas.

        self.diagram_canvas.axes.set_xlabel('d [m]')
        self.diagram_canvas.axes.set_ylabel('fb [MPa]')
        #self.diagram_canvas.axes.invert_yaxis()

        self.diagram_canvas.axes.grid(True)

        self.diagram_canvas.axes.autoscale_view()
        self.diagram_canvas.axes.spines['top'].set_visible(False)
        self.diagram_canvas.axes.spines['right'].set_visible(True)

        ### functions to draw the interaction diagrams
        Bp_y_array = np.linspace(0, self.globvar.By, 100)

        #ACI
        if ( self.globvar.flag_id_ACI):
            N_max_ACI, M_max_ACI, f_b, d = ACI(Bp_y_array, self.globvar)
            self.diagram_canvas.axes.plot(d, f_b, label='ACI', color= 'green', linestyle='solid')

        #ASCE
        if ( self.globvar.flag_id_ASCE):
            N_max_ASCE, M_max_ASCE, f_b, d = ASCE(self.globvar)
            self.diagram_canvas.axes.plot(d, f_b, label='ASCE', color='green', marker='o',markersize=5,linestyle='none')

        #concrete
        if ( self.globvar.flag_id_concrete):
            N_max_concrete, M_max_concrete, f_b, d = id_fcm(Bp_y_array, self.globvar)
            self.diagram_canvas.axes.plot(d, f_b, label='fcm', color='purple', linestyle='solid')

        #CTU total
        if (self.globvar.flag_id_CTU):

            N_max_AS_MSR, M_max_AS_MSR, f_b_MSR, d_MSR = AS_MSR(Bp_y_array, self.globvar)
            self.diagram_canvas.axes.plot(d_MSR, f_b_MSR, label='MSR', color='blue', linestyle='solid')

            N_max_AS_TIE, M_max_AS_TIE, f_b_TIE, d_TIE = AS_TIE(Bp_y_array, self.globvar)
            self.diagram_canvas.axes.plot(d_TIE, f_b_TIE, label='Ties', color='blue', linestyle='--')


        if (self.globvar.flag_id_CTU_concrete):
            N_max, M_max_, f_b, d, = AS_plain_concrete(Bp_y_array, self.globvar)
            self.diagram_canvas.axes.plot(d, f_b, label='plain concrete', color='red', linestyle='--')

        ## FEM
        if self.globvar.flag_id_FEM:

            max_fb = []
            max_d = []

            for task in self.globvar.tasks:
                # if task.status in (Task_status.COMPLETED, Task_status.PREDEFINED):
                if task.status == Task_status.COMPLETED:
                    max_fb.append(task.max_fb_d[0])
                    max_d.append(task.max_fb_d[1])

            if self.globvar.flag_id_FEM_paths:
                for task in self.globvar.tasks:
                    if task.status == Task_status.COMPLETED:

                        aux_d = [abs(m) for m in task.d]
                        aux_fb = [abs(n) for n in task.fb]

                        aux_d.insert(0,task.d[0])
                        aux_fb.insert(0,0.)


                        self.diagram_canvas.axes.plot(
                            aux_d, aux_fb,
                            color=(192. / 255., 192. / 255., 192. / 255.),
                            linestyle='solid', linewidth=0.5, marker='x',
                            markeredgecolor='black', markeredgewidth=0.5, markersize=2
                        )

            '''
            if len(max_d) >= 3:
                self.diagram_canvas.axes.plot(max_d, max_fb, '-or', markeredgewidth=1.5, markersize=6, label="FEM",
                                              markeredgecolor='black')
            elif len(max_d) >= 1:
                self.diagram_canvas.axes.plot(max_d, max_fb, 'or', markeredgewidth=1.5, markersize=6, label="FEM",
                                              markeredgecolor='black')
            '''
            if len(max_d) >= 1:
                self.diagram_canvas.axes.plot(max_d, max_fb, 'or', markeredgewidth=1.5, markersize=6, label="FEM",
                                              markeredgecolor='black')

        self.diagram_canvas.axes.legend(loc='lower right')
        self.diagram_canvas.figure.tight_layout()
        self.diagram_canvas.draw()

        self.need_update_ids_flag = False

    def update_large_ties(self):
        self.globvar.dT = self.globvar.Bx - 2. * self.globvar.cover - (self.globvar.rebars_CS[self.globvar.DT]).diam
        for reb in self.globvar.rebars[:]:
            if reb.give_rebar_type() == 'tie' and reb.tag.startswith('TL'):
                # Reuse original position
                center = reb.XYZ
                vec = [0., 0., 1.]

                diameter = self.globvar.rebars_CS[self.globvar.DT].diam
                cs_tag = self.globvar.DT

                # definition of ben diameters
                if cs_tag.startswith('D'):  # EC
                    if diameter <= 0.016:
                        corner_radius = 4.5 * diameter
                    else:
                        corner_radius = 7.5 * diameter
                elif cs_tag.startswith('#'):  #ACI
                    if diameter <= 25.4e-3:  #less than #8
                        corner_radius = 4.5 * diameter
                    elif diameter <= 35.81e-3:  #less than #11
                        corner_radius = 6.5 * diameter
                    else:       #14 - #18
                        corner_radius = 8.5 * diameter
                else:
                    raise ValueError(f"Uknown format: {cs_tag}")

                # Create update tie
                new_tie = tie(
                    globvar=self.globvar,
                    XYZ=center,
                    cs_nr=reb.cs_nr,
                    cs_tag=self.globvar.DT,
                    cs_mat=reb.give_mat_nr(self.globvar),
                    width=self.globvar.dT,
                    height=self.globvar.dT,
                    vec=vec,
                    tag=reb.tag,
                    corner_radius=corner_radius
                )
                self.globvar.corner_radius = corner_radius
                # Replace old tie in the list
                self.globvar.rebars[self.globvar.rebars.index(reb)] = new_tie

    def update_small_spirals(self):

        for reb in self.globvar.rebars[:]:
            if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'S') ):

                dS = reb.radius * 2.
                rebar_DS = reb.give_diameter(self.globvar)

                auxX = self.globvar.Bx/2. - self.globvar.cover - dS/2  - rebar_DS/2.
                auxY = self.globvar.By/2. - self.globvar.cover - dS/2. - rebar_DS/2.

                # position of center and start based on previous position
                center = [auxX * np.sign(reb.XYZ[0]), auxY * np.sign(reb.XYZ[1]), reb.XYZ[2]]
                start = [ center[0], center[1] - dS/2., center[2] ]

                spiral_S = spiral(globvar = self.globvar, XYZ = center, cs_nr = reb.cs_nr, cs_tag = reb.give_cs_tag(self.globvar), cs_mat = reb.give_mat_nr(self.globvar), vec = self.vec, pitch = self.globvar.s, axial_length = self.globvar.s, start = start, tag = reb.tag)

                #spiral_S.check_topology_consistency(self.globvar)

                # switch the items in the container
                self.globvar.rebars[ self.globvar.rebars.index(reb) ] = spiral_S


    def update_large_spirals(self):

        for reb in self.globvar.rebars[:]:
            if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'L') ):

                dL = min(self.globvar.Bx, self.globvar.By ) - 2.*self.globvar.cover
                dL -= reb.give_diameter(self.globvar)

                # position of center and start based on previous position
                center = reb.XYZ
                start = [ center[0], center[1] - dL/2., center[2] ]

                spiral_L = spiral(globvar = self.globvar, XYZ = center, cs_nr = reb.cs_nr, cs_tag = reb.give_cs_tag(self.globvar), cs_mat = reb.give_mat_nr(self.globvar), vec = self.vec, pitch = self.globvar.s, axial_length = self.globvar.s, start = start, tag = reb.tag)
                #spiral_L.check_topology_consistency(self.globvar)

                # switch the items in the container
                self.globvar.rebars[ self.globvar.rebars.index(reb) ] = spiral_L


    def update_vertical_rebars(self):
        self.remove_vertical_rebars()
        self.create_vertical_rebars()


    ####
    # VARIABLES FOR TOPOLOGY
    ####

    def compute_equivalent_ties(self):

        vol_concrete = self.globvar.Bx * self.globvar.By * self.globvar.s
        rho_msr = self.compute_lateral_reinforcement_ratio()

        target_vol_steel = rho_msr * vol_concrete

        DT_current = 0.01
        tolerance = 1e-3
        max_iter = 20

        for i in range(max_iter):

            # corner radius
            if DT_current <= 0.016:
                r_corner = 4.5 * DT_current
            else:
                r_corner = 7.5 * DT_current

            # tie lenght
            l_tie = (
                    r_corner * 2 * math.pi +
                    4 * (self.globvar.Bx - 2 * r_corner - 2 * self.globvar.cover)
            )

            if l_tie <= 0:
                return 0.0

            # CS area
            CS_area = target_vol_steel / l_tie

            # new diam
            DT_new = math.sqrt((4 * CS_area) / math.pi)

            # check if tolerance is satisfied
            if abs(DT_new - DT_current) < tolerance:
                DT_current = DT_new
                break

            # next iter
            DT_current = DT_new

        # change
        self.globvar.DT_equivalent = DT_current

    def compute_equivalent_MSR(self):

        vol_concrete = self.globvar.Bx * self.globvar.By * self.globvar.s
        rho_tie = self.compute_lateral_reinforcement_ratio_tie()

        target_vol_steel = vol_concrete * rho_tie

        large_spiral_d = self.globvar.Bx - 2 * self.globvar.cover
        small_spiral_d = large_spiral_d / 3

        spirals_lenght = large_spiral_d * math.pi + 4 * (small_spiral_d * math.pi)

        CS_area = target_vol_steel / spirals_lenght
        DL_eqq = math.sqrt((4 * CS_area) / math.pi)

        self.globvar.DL_equivalent = DL_eqq


    # set cross-section width & update
    def set_Bx(self, s):

        try:
            Bx = float(s)
        except ValueError:
            Bx_default = self.globvar.Bx
            message = "unsupported value, using Bx = " + str(Bx_default)
            self.warning_dialog_ok(message)
            warnings.warn(message)
            Bx = Bx_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.Bx == Bx ):
            return

        # change
        if ( self.change_model_definition() ):

            self.globvar.Bx = Bx
            self.globvar.By = Bx
            self.globvar.Bz = round(Bx * 1.5, 3)
            self.globvar.Bp = Bx * self.globvar.load_plate_width
            self.globvar.dS = (Bx - 2 * self.globvar.cover) / 3
            if self.globvar.dS < self.globvar.minimum_radius * 2:
                self.globvar.dS = self.globvar.minimum_radius * 2

            logger.info(f"Setting Bx = {(self.globvar.Bx):.3f}")
            logger.info(f"Setting By = {(self.globvar.By):.3f}")
            logger.info(f"Setting Bz = {(self.globvar.Bz):.3f}")
            logger.info(f"Setting Bp = {(self.globvar.Bp):.3f}")

            self.update_small_spirals()
            self.update_large_spirals()
            self.update_large_ties()
            self.update_vertical_rebars()

            # determine minimum element size - default values for computational efficiency
            def_elem_size_horizontal = min(self.globvar.Bx, self.globvar.By)/15.

            # update bounds for FE mesh
            self.widget_dx.setMinimum( self.globvar.Bx/200. )
            self.widget_dx.setMaximum( self.globvar.Bx/8. )
            self.widget_dx.setValue( def_elem_size_horizontal )

            self.widget_dy.setMinimum( self.globvar.By/200. )
            self.widget_dy.setMaximum( self.globvar.By/8. )
            self.widget_dy.setValue( def_elem_size_horizontal )


        # sync widgets
        self.widget_Bx.blockSignals(True)
        self.widget_Bx.setValue(self.globvar.Bx)
        self.widget_Bx.blockSignals(False)
        self.widget_Bx_TIE.blockSignals(True)
        self.widget_Bx_TIE.setValue(self.globvar.Bx)
        self.widget_Bx_TIE.blockSignals(False)
        self.widget_dS.setMaximum(max(self.globvar.Bx / 2. - self.globvar.cover - (self.globvar.rebars_CS[self.globvar.DS]).diam, self.globvar.minimum_radius * 2))
        self.widget_dS.setValue(self.globvar.dS)



    # set spiral pitch & update
    def set_s(self, s):

        try:
            s = float(s)
        except ValueError:
            s_default = self.globvar.s
            self.warning_dialog_ok("unsupported value, using s = " + str(s_default))
            warnings.warn("unsupported value, using s = " + str(s_default))
            s = s_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.s == s ):
            return

        # change
        if ( self.change_model_definition() ):

            self.globvar.s = s
            logger.info(f"Setting s = {(self.globvar.s):.3f}")

            self.update_small_spirals()
            self.update_large_spirals()
            self.update_large_ties()
            self.update_vertical_rebars()

        # widget sync
        self.widget_s.blockSignals(True)
        self.widget_s.setValue(self.globvar.s)
        self.widget_s.blockSignals(False)
        self.widget_s_TIE.blockSignals(True)
        self.widget_s_TIE.setValue(self.globvar.s)
        self.widget_s_TIE.blockSignals(False)


    # set concrete cover and update spirals topology
    def set_cover(self, s):

        try:
            cover = float(s)
        except ValueError:
            cover_default = self.globvar.cover
            self.warning_dialog_ok("unsupported value, using cover = " + str(cover_default) )
            warnings.warn("unsupported value, using cover = " + str(cover_default))
            cover = cover_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.cover == cover ):
            return

        # change
        if ( self.change_model_definition() ):

            self.globvar.cover = cover
            self.globvar.dS = (self.globvar.Bx - 2 * cover) / 3
            logger.info(f"Setting cover = {(self.globvar.cover):.3f}")


            self.update_small_spirals()
            self.update_large_spirals()
            self.update_large_ties()
            self.update_vertical_rebars()

        # widget sync
        self.widget_cover.blockSignals(True)
        self.widget_cover.setValue(self.globvar.cover)
        self.widget_cover.blockSignals(False)
        self.widget_cover_TIE.blockSignals(True)
        self.widget_cover_TIE.setValue(self.globvar.cover)
        self.widget_cover_TIE.blockSignals(False)
        self.widget_dS.setMaximum(max(self.globvar.Bx / 2. - self.globvar.cover - (self.globvar.rebars_CS[self.globvar.DS]).diam, self.globvar.minimum_radius * 2))
        self.widget_dS.setValue(self.globvar.dS)


    def set_n_v_bars(self, s):

        try:
            n_v_bars = int(s)
        except ValueError:
            n_v_bars_default = self.globvar.s
            self.warning_dialog_ok("unsupported value, using # of vertical bars = " + str(n_v_bars_default) )
            warnings.warn("unsupported value, using # of vertical bars = " + str(n_v_bars_default))
            n_v_bars = n_v_bars_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.n_v_bars == n_v_bars ):
            return

        # change
        if ( self.change_model_definition() ):

            self.globvar.n_v_bars = n_v_bars
            logger.info(f"Setting n_v_bars = {(self.globvar.n_v_bars):.3f}")


            self.update_vertical_rebars()

        else:
            self.widget_n_v_bars.setValue(self.globvar.n_v_bars)

    # set load plate width
    def set_load_plate_width (self, s):

        try:
            load_plate_width = float(s)
        except ValueError:
            load_plate_width_default = self.globvar.load_plate_width
            self.warning_dialog_ok("unsupported value, using load plate width = " + str(load_plate_width_default))
            warnings.warn("unsupported value, using load plate width = " + str(load_plate_width_default))
            load_plate_width = load_plate_width_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.load_plate_width == load_plate_width ):
            return

            # change
        if (self.change_model_definition()):

            self.globvar.load_plate_width = load_plate_width
            self.globvar.Bp = load_plate_width * self.globvar.Bx
            logger.info(f"Setting load plate width = {(self.globvar.load_plate_width):.3f}")
            logger.info(f"Setting Bp = {(self.globvar.Bp):.3f}")

            """
            self.update_small_spirals()
            self.update_large_spirals()
            self.update_large_ties()
            self.update_vertical_rebars()
            self.diagram_MSR.update_MSR_intersections(self.globvar)
            """
        # widget sync
        self.widget_load_plate.blockSignals(True)
        self.widget_load_plate.setValue(self.globvar.load_plate_width)
        self.widget_load_plate.blockSignals(False)
        self.widget_load_plate_tie.blockSignals(True)
        self.widget_load_plate_tie.setValue(self.globvar.load_plate_width)
        self.widget_load_plate_tie.blockSignals(False)



    # set axial diameter to all small spirals
    def set_dS(self, s):

        try:
            dS = float(s)
        except ValueError:
            dS_default = self.globvar.rebars[0].radius * 0.3 * 2.
            self.warning_dialog_ok("unsupported value, using dS = " + str(dS_default) )
            warnings.warn("unsupported value, using dS = " + str(dS_default))
            dS = dS_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.dS == dS ):
            return

        # change
        if ( self.change_model_definition() ):

            self.globvar.dS = dS
            logger.info(f"Setting dS = {(self.globvar.dS):.3f}")

            for reb in self.globvar.rebars:
                if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'S') ):

                    reb.radius = dS/2.

            self.update_small_spirals()
            self.update_vertical_rebars()

        else:
            self.widget_dS.setValue( self.globvar.dS )



    # set rebar diameter to all small spirals
    def set_DS(self, s):

        try:
            rebar_DS = (self.globvar.rebars_CS[s]).diam

        except ValueError:
            rebar_DS_default = (self.globvar.rebars_CS[self.globvar.DS]).diam
            self.warning_dialog_ok("unsupported value, using DS = " + str(rebar_DS_default) )
            warnings.warn("unsupported value, using DS = " + str(rebar_DS_default))
            rebar_DS = rebar_DS_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.DS == s ):
            return

        # change
        if ( self.change_model_definition() ):

            self.globvar.DS = s
            logger.info("Setting DS = " + self.globvar.DS)
            area = self.globvar.rebars_CS[self.globvar.DS].area

            # update cross-section first - otherwise a new one is created in "update_xxx_spirals"
            for reb in self.globvar.rebars:
                if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'S') ):
                    cs = reb.find_cross_section(self.globvar, reb.cs_nr)
                    cs.set_properties_from_cs_tag(self.globvar, self.globvar.DS)
                    break

            self.update_small_spirals()
            self.update_vertical_rebars()

        else:
            self.widget_DS.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DS))

        self.widget_dS.setMaximum(max(self.globvar.Bx / 2. - self.globvar.cover - (self.globvar.rebars_CS[self.globvar.DS]).diam, self.globvar.minimum_radius * 2))

        self.compute_equivalent_ties()

    # set axial diameter to large spiral(s)
    def set_dL(self, s):

        try:
            dL = float(s)
        except ValueError:
            dL_default = min(self.globvar.Bx, self.globvar.By )/2. - 2.*self.globvar.cover
            dL_default -= (self.globvar.rebars_CS[self.globvar.DL]).diam

            self.warning_dialog_ok("unsupported value, using dL = " + str(dL_default) )
            warnings.warn("unsupported value, using dL = " + str(dL_default))
            dL = dL_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.dL == dL ):
            return

        # change
        if ( self.change_model_definition() ):

            self.globvar.dL = dL
            logger.info(f"Setting dL = {(self.globvar.dL):.3f}")

            for reb in self.globvar.rebars:
                if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'L') ):

                    reb.radius = dL/2.

            self.update_large_spirals()
            self.update_large_ties()
            self.update_vertical_rebars()


        else:
            self.widget_dL.setValue( self.globvar.dL )


    # set rebar diameter to large spiral
    def set_DL(self, s):

        try:
            rebar_DL = (self.globvar.rebars_CS[s]).diam

            if (s == "none"):   #exclude none value
                rebar_DL_default = (self.globvar.rebars_CS[self.globvar.DL]).diam
                self.warning_dialog_ok("Main reinforcement must be define, using DL= " + str(rebar_DL_default) )
                rebar_DL = rebar_DL_default
                self.widget_DL.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DL))
                return

        except ValueError:
            rebar_DL_default = (self.globvar.rebars_CS[self.globvar.DL]).diam
            self.warning_dialog_ok("unsupported value, using DL  = " + str(rebar_DL_default) )
            warnings.warn("unsupported value, using DL  = " + str(rebar_DL_default))
            rebar_DL = rebar_DL_default
            self.widget_DL.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DL))


        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.DL == s ):
            return

        # change
        if ( self.change_model_definition() ):
            self.globvar.DL = s
            logger.info(f"Setting DL = " + self.globvar.DL)

            area = self.globvar.rebars_CS[self.globvar.DL].area

            for reb in self.globvar.rebars:
                if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'L') ):
                    cs = reb.find_cross_section(self.globvar, reb.cs_nr)
                    cs.set_properties_from_cs_tag(self.globvar, self.globvar.DL)
                    break



            self.update_large_spirals()
            #self.update_large_ties()
            #self.update_vertical_rebars()
            #self.diagram_MSR.update_MSR_intersections(self.globvar)

            self.set_DS(s)
            # sync widgets
            self.widget_DL.blockSignals(True)
            self.widget_DL.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DL))
            self.widget_DL.blockSignals(False)
            self.widget_DS.blockSignals(True)
            self.widget_DS.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DS))
            self.widget_DS.blockSignals(False)
        else:
            self.widget_DL.blockSignals(True)
            self.widget_DL.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DL))
            self.widget_DL.blockSignals(False)
            self.widget_DS.blockSignals(True)
            self.widget_DS.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DS))
            self.widget_DS.blockSignals(False)

    # set rebar diameter to large ties
    def set_DT(self, s):

        try:
            rebar_DT = (self.globvar.rebars_CS[s]).diam

            if (s == "none"):  # exclude none value
                rebar_DT_default = (self.globvar.rebars_CS[self.globvar.DT]).diam
                self.warning_dialog_ok(
                    "Main reinforcement must be define, using DT = " + str(rebar_DT_default))
                rebar_DT = rebar_DT_default
                self.widget_DT.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DT))
                return

        except ValueError:
            rebar_DT_default = (self.globvar.rebars_CS[self.globvar.DT]).diam
            self.warning_dialog_ok("unsupported value, using DL (same as DT) = " + str(rebar_DT_default))
            warnings.warn("unsupported value, using DL (same as DT) = " + str(rebar_DT_default))
            rebar_DT = rebar_DT_default
            self.widget_DT.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DT))

        # not to show the dialog again if the value has been set to its former setting
        if (self.globvar.DT == s):
            return

        # change
        if (self.change_model_definition()):
            self.globvar.DT = s
            logger.info(f"Setting DT = " + self.globvar.DT)

            area = self.globvar.rebars_CS[self.globvar.DT].area

            self.update_large_ties()
            self.update_vertical_rebars()

        else:
            self.widget_DT.blockSignals(True)
            self.widget_DT.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DT))
            self.widget_DT.blockSignals(False)

        self.compute_equivalent_MSR()

    def set_DV(self, s):

        try:
            rebar_DV = (self.globvar.rebars_CS[s]).diam
        except ValueError:
            rebar_DV_default = (self.globvar.rebars_CS[self.globvar.DV]).diam
            self.warning_dialog_ok("unsupported value, using DV = " + str(rebar_DV_default) )
            warnings.warn("unsupported value, using DV = " + str(rebar_DV_default))
            rebar_DV = rebar_DV_default


        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.DV == s ):
            return

        # change
        if ( self.change_model_definition() ):

            self.globvar.DV = s
            logger.info(f"Setting DV = " + self.globvar.DV)

            for reb in self.globvar.rebars:
                if ( (reb.give_rebar_type() == 'rebar') and (reb.tag[0] == 'V') ):
                    cs = reb.find_cross_section(self.globvar, reb.cs_nr)
                    cs.set_properties_from_cs_tag(self.globvar, self.globvar.DV)
                    break
                    # cross-section needs to be updated!

            self.update_vertical_rebars()

            # update cross-section
            for reb in self.globvar.rebars:
                if ( (reb.give_rebar_type() == 'rebar') and (reb.tag[0] == 'V') ):
                    #cs = reb.find_cross_section(self.globvar, reb.cs)
                    #cs.area = reb.area
                    # TODO
                    break

        else:
            self.widget_DV.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DV))

        # sync widgets
        self.widget_DV.blockSignals(True)
        self.widget_DV.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DV))
        self.widget_DV.blockSignals(False)
        self.widget_DV_TIE.blockSignals(True)
        self.widget_DV_TIE.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DV))
        self.widget_DV_TIE.blockSignals(False)

    def topology_set(self):
        select_layout = self.topology_switch.currentIndex()
        self.topology_specs_layout.setCurrentIndex(select_layout)
        self.globvar.topology_layout = select_layout
        self.update_vertical_rebars()

    ####
    # VARIABLE FOR INTERACTION DIAGRAM
    ####

    def set_id_plot(self):
        select_id = self.id_switcher.currentIndex()
        self.globvar.id_plot = select_id
        self.update_diagram_plot()

    # set concrete grade
    def set_fcm(self, s):

        # selection from catalogue
        try:
            fcm_val = float(self.globvar.concretes[s] )

        except ValueError:
            fcm_default = (self.globvar.concretes[self.globvar.fcm])
            self.warning_dialog_ok("unsupported value, using fcm = " + str(fcm_default) )
            warnings.warn("unsupported value, using fcm = " + str(fcm_default))
            fcm_val = fcm_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.fcm == s ):
            return

        # change
        if ( self.change_model_definition() ):
            self.globvar.fcm = s
            logger.info(f"Setting fcm = " + self.globvar.fcm)
            self.cdpm2.fcm = fcm_val
            self.cdpm2.predict_concrete_parameters()

        else:
            self.widget_fcm.setCurrentIndex(list(self.globvar.concretes.keys()).index(self.globvar.fcm))


    # set steel grade for vertical (longitudinal) reinforcement
    def set_fy_vert(self, s):

        try:
            fy = float(s)
        except ValueError:
            fy_default = self.globvar.fy_vert
            self.warning_dialog_ok("unsupported value, using fy = " + str(fy_default) )
            warnings.warn("unsupported value, using fy = " + str(fy_default))
            fy = fy_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.fy_vert == fy ):
            return

        # change
        if ( self.change_model_definition() ):

            self.globvar.fy_vert = fy
            logger.info(f"Setting fy_vert = {(self.globvar.fy_vert):.3f}")
            self.mises_vert.sig_0 = fy

        else:
            self.widget_fy_vert.setValue( self.globvar.fy_vert )



    # set steel grade for lateral (transverse) reinforcement
    def set_fy_lat(self, s):

        try:
            fy = float(s)
        except ValueError:
            fy_default = self.globvar.fy_lat
            self.warning_dialog_ok("unsupported value, using fy = " + str(fy_default) )
            warnings.warn("unsupported value, using fy = " + str(fy_default))
            fy = fy_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.fy_lat == fy ):
            return

        # change
        if ( self.change_model_definition() ):

            self.globvar.fy_lat = fy
            logger.info(f"Setting fy_lat = {(self.globvar.fy_lat):.3f}")
            self.mises_lat.sig_0 = fy

        else:
            self.widget_fy_lat.setValue( self.globvar.fy_lat )


    def set_dx(self,s):

        try:
            dx = float(s)
        except ValueError:
            dx_default = self.globvar.elem_size_X
            self.warning_dialog_ok("unsupported value, using = " + str(dx_default) )
            warnings.warn("unsupported value, using = " + str(dx_default))
            dx = dx_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.elem_size_X == dx ):
            return

        # change
        if ( self.change_model_definition() ):

            self.globvar.elem_size_X = dx

        else:
            self.widget_dx.setValue( self.globvar.elem_size_X )
            logger.info(f"Setting elem_size_X = {(self.globvar.elem_size_X):.3f}")

    def set_dy(self,s):

        try:
            dy = float(s)
        except ValueError:
            dy_default = self.globvar.elem_size_YZ
            self.warning_dialog_ok("unsupported value, using = " + str(dy_default) )
            warnings.warn("unsupported value, using = " + str(dy_default))
            dy = dy_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.elem_size_YZ == dy ):
            return

        # change
        if ( self.change_model_definition() ):

            self.globvar.elem_size_YZ = dy

        else:
            self.widget_dy.setValue( self.globvar.elem_size_YZ )
            logger.info(f"Setting elem_size_YZ = {(self.globvar.elem_size_YZ):.3f}")

    def show_concrete_mesh(self):
        # turn off buttons and turn on viewer
        self.show_concrete_button.setEnabled(False)
        self.show_rebars_button.setEnabled(False)
        self.button_generate_oofem_input.setEnabled(False)
        self.globvar.flag_show_mesh = True

        # generate and show mesh
        generate_concrete_mesh(self.globvar)

        # default
        self.globvar.flag_show_mesh = False
        self.show_concrete_button.setEnabled(True)
        self.show_rebars_button.setEnabled(True)
        self.button_generate_oofem_input.setEnabled(True)
        print("heloo")

    def show_rebars_mesh(self):
        # turn off buttons and turn on viewer
        self.show_concrete_button.setEnabled(False)
        self.show_rebars_button.setEnabled(False)
        self.button_generate_oofem_input.setEnabled(False)
        self.globvar.flag_show_mesh = True

        # generate and show mesh
        generate_rebars_mesh(self.globvar)

        # default
        self.globvar.flag_show_mesh = False
        self.show_concrete_button.setEnabled(True)
        self.show_rebars_button.setEnabled(True)
        self.button_generate_oofem_input.setEnabled(True)

    def set_CDPM2(self):
        # change
        if ( self.change_model_definition() ):
            self.define_CDPM2(self.cdpm2)


    def set_mises_vert(self):
        # change
        if ( self.change_model_definition() ):
            self.define_Mises(self.mises_vert)


    def set_mises_lat(self):
        # change
        if ( self.change_model_definition() ):
            self.define_Mises(self.mises_lat)



    def set_id_steel(self,s):
        self.globvar.flag_id_steel = s
    def set_id_ACI(self,s):
        self.globvar.flag_id_ACI = s

    def set_id_ASCE(self,s):
        self.globvar.flag_id_ASCE = s

    def set_id_CTU(self,s):
        self.globvar.flag_id_CTU = s

    def set_id_CTU_concrete(self,s):
        self.globvar.flag_id_CTU_concrete = s

    def set_id_FEM(self,s):
        self.globvar.flag_id_FEM = s

        if ( not self.globvar.flag_id_FEM ):
            self.checkbox_id_FEM_paths.setCheckState(Qt.Unchecked)


    def set_id_FEM_paths(self,s):
        self.globvar.flag_id_FEM_paths = s

    def set_active_neutral_axis(self,s):
        self.globvar.flag_active_neutral_axis = s

    def set_loading_display(self,s):
        self.globvar.flag_loading_display = s


    # set steel grade for lateral (transverse) reinforcement
    def set_c_neutral_axis(self, s):

        try:
            c = float(s)
        except ValueError:
            c_default = self.globvar.c_neutral_axis
            self.warning_dialog_ok("unsupported value, using c = " + str(c_default) )
            warnings.warn("unsupported value, c = " + str(c_default))
            c = c_default

        self.globvar.c_neutral_axis = c

    def set_show_Bp(self, s):
        self.globvar.flag_show_Bp = s
        self.update_topology_plot()
        self.widget_show_load_plate.blockSignals(True)
        self.widget_show_load_plate.setChecked(s)
        self.widget_show_load_plate.blockSignals(False)
        self.widget_show_load_plate_tie.blockSignals(True)
        self.widget_show_load_plate_tie.setChecked(s)
        self.widget_show_load_plate_tie.blockSignals(False)

    def set_ignore_cover(self,s):

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.flag_ignore_cover == s ):
            return

        if ( self.change_model_definition() ):
            self.globvar.flag_ignore_cover = s
        else:
            self.checkbox_ignore_cover.setCheckState(self.globvar.flag_ignore_cover)

        print(f"cover_check : {self.globvar.flag_ignore_cover}")

    def change_model_definition(self):
        # if the definition has been changed previously, do nothing
        if (self.globvar.flag_problem_changed):
            self.need_update_ids_flag = True
            return True

        # otherwise ask whether the content should be deleted
        else:
            button = QMessageBox.question(self, "Change in problem definition detected", " Delete computed results/generated inputs and proceed?")

        if button == QMessageBox.Yes:
            # change flags
            self.globvar.flag_problem_changed = True
            self.globvar.flag_output_generated = False

            # project directory will be cleared in a sequel upon input creation
            for task in self.globvar.tasks:
                task.reset()

            self.update_ecc_selection()

            self.console.clear()

            self.progressBar.setValue(0)

            self.need_update_ids_flag = True
            return True
        else:
            return False



    def remove_vertical_rebars(self):

        # removing all vertical rebars, iterating over a copy of a list
        for reb in self.globvar.rebars[:]:
            if ( (reb.give_rebar_type() == 'rebar') and (reb.tag[0] == 'V') ):
                self.globvar.rebars.remove(reb)

    def create_vertical_rebars(self):
        if (self.globvar.topology_layout == 0):
            self.create_vertical_rebars_spirals()

        elif (self.globvar.topology_layout == 1):
            self.create_vertical_rebars_ties()

        else:
            self.warning_dialog_ok("unknow layout")
            warnings.warn("unknow layout")

    def create_vertical_rebars_ties(self):
        DV = self.globvar.rebars_CS[self.globvar.DV].diam
        H = self.globvar.s
        n_total = self.globvar.n_v_bars

        def sample_along_polyline(points, num_samples):
            points = np.array(points)
            seg_lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
            total_length = np.sum(seg_lengths)
            distances = np.linspace(0, total_length, num_samples, endpoint=False)

            result = []
            cum_length = 0.0
            seg_idx = 0
            for d in distances:
                while seg_idx < len(seg_lengths) - 1 and cum_length + seg_lengths[seg_idx] < d:
                    cum_length += seg_lengths[seg_idx]
                    seg_idx += 1
                seg_start = points[seg_idx]
                seg_end = points[seg_idx + 1]
                seg_length = seg_lengths[seg_idx]
                t = (d - cum_length) / seg_length if seg_length > 0 else 0
                pt = (1 - t) * seg_start + t * seg_end
                result.append(pt)
            return result

        for reb in self.globvar.rebars:
            if reb.give_rebar_type() == 'tie' and reb.tag.startswith('TL'):
                cs_tag = reb.give_cs_tag(self.globvar)
                cs_mat = reb.give_mat_nr(self.globvar)

                shrink = DV + self.globvar.rebars_CS[self.globvar.DL].diam
                scale = (reb.width - shrink) / reb.width
                scaled_r = reb.r * scale

                small_tie = tie(
                    globvar=self.globvar,
                    XYZ=reb.XYZ,
                    cs_nr=reb.cs_nr,
                    cs_tag=cs_tag,
                    cs_mat=cs_mat,
                    width=reb.width - shrink,
                    height=reb.height - shrink,
                    vec=reb.vec,
                    corner_radius=scaled_r,
                    div_length=reb.div_length,
                    tag='TEMP',
                    connect=False
                )

                points = list(small_tie.rebar_points)

                # find the corner rebar
                corner_idx = max(range(len(points)), key=lambda i: points[i][0] + points[i][1])
                points_rotated = points[corner_idx:] + points[:corner_idx]

                sampled_pts = sample_along_polyline(points_rotated, n_total)

                for i, pt in enumerate(sampled_pts):
                    origin = [pt[0], pt[1], 0.0]
                    end = [pt[0], pt[1], H]
                    tag = f'V_TL{i + 1}_{reb.tag}'
                    new_rebar = rebar(
                        globvar=self.globvar,
                        XYZ=[origin, end],
                        cs_nr=4,
                        cs_tag=self.globvar.DV,
                        cs_mat=3,
                        tag=tag
                    )
                    self.globvar.rebars.append(new_rebar)

            for reb in self.globvar.rebars:
                if reb.give_rebar_type() == 'tie' and reb.tag.startswith("TL"):
                    lenght = reb.compute_rebar_length()
                    break

        s_v_bars = lenght / n_total - self.globvar.rebars_CS[self.globvar.DV].diam
        s_min = max(1.5 * (self.globvar.rebars_CS[self.globvar.DV]).diam, 0.0381)


        if s_v_bars < s_min:
            self.s_min_label.setText(f"Vertical reinforcement does not comply with the ACI\n [s_min = {s_min:.3f} m]")
            self.s_min_label.setStyleSheet("QLabel { color: red; background-color: yellow; font-weight: bold}")
        else:
            self.s_min_label.setStyleSheet("color: transparent; background-color: white")
            self.s_min_label.setText("")



    def create_vertical_rebars_spirals(self):

        DV = self.globvar.rebars_CS[self.globvar.DV].diam

        #cs = self.find_cross_section(globvar, self.cs)
        #cs.D =
        #cs.A =

        # large spirals - vertical quarters - inner
        for reb in self.globvar.rebars:
            if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'L') ):

                # start with axial radius, subtract spiral rebar radius and radius of vertical rebar
                diam = reb.give_diameter(self.globvar)
                plan_radius = reb.radius - diam/2. - DV/2.

                for i in range(4):
                    if (i == 0): # X+
                        topology_origin = [reb.XYZ[0] + plan_radius, reb.XYZ[1], 0.]
                    elif (i == 1): # X-
                        topology_origin = [reb.XYZ[0] - plan_radius, reb.XYZ[1], 0.]
                    elif (i == 2): # Y+
                        topology_origin = [reb.XYZ[0], reb.XYZ[1] + plan_radius, 0.]
                    elif (i == 3): # Y-
                        topology_origin = [reb.XYZ[0], reb.XYZ[1] - plan_radius, 0.]

                    topology_end = [topology_origin[0], topology_origin[1], self.globvar.s]
                    topology = [topology_origin, topology_end]

                    rebar_V = rebar(globvar = self.globvar, XYZ = topology, cs_nr = 4, cs_tag = self.globvar.DV, cs_mat = 3, tag = 'V_center_'+reb.tag)
                    #rebar_V.check_topology_consistency(self.globvar)
                    self.globvar.rebars.append(rebar_V)


        # large spirals - diagonals - inner
        for reb in self.globvar.rebars:
            if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'L') ):

                # start with axial radius, subtract spiral rebar radius and radius of vertical rebar
                diam = reb.give_diameter(self.globvar)
                plan_radius = reb.radius - diam/2. - DV/2.
                radius_projection = plan_radius/math.sqrt(2)

                for i in range(4):
                    if (i == 0): # X+, Y+
                        topology_origin = [reb.XYZ[0] + radius_projection, reb.XYZ[1] + radius_projection, 0.]
                    elif (i == 1): # X-, Y+
                        topology_origin = [reb.XYZ[0] - radius_projection, reb.XYZ[1] + radius_projection, 0.]
                    elif (i == 2): # X-, Y-
                        topology_origin = [reb.XYZ[0] - radius_projection, reb.XYZ[1] - radius_projection, 0.]
                    elif (i == 3): # X+, Y-
                        topology_origin = [reb.XYZ[0] + radius_projection, reb.XYZ[1] - radius_projection, 0.]

                    topology_end = [topology_origin[0], topology_origin[1], self.globvar.s]
                    topology = [topology_origin, topology_end]

                    rebar_V = rebar(globvar = self.globvar, XYZ = topology, cs_nr = 4, cs_tag = self.globvar.DV, cs_mat = 3, tag = 'V_diag_'+reb.tag)
                    #rebar_V.check_topology_consistency(globvar = self.globvar)
                    self.globvar.rebars.append(rebar_V)


        # small spirals - diagonals - inner, outer
        # flag_inner_S = True
        flag_inner_S = False
        flag_outer_S = True

        for reb in self.globvar.rebars:
            if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'S') ):
                # start with axial radius
                plan_radius = reb.radius
                # math.sqrt( (reb.start[0]-reb.XYZ[0])**2 + (reb.start[1]-reb.XYZ[1])**2 )
                plan_radius -= reb.give_diameter(self.globvar)/2. # inner radius
                plan_radius -= DV/2. # to axis of vertical rebar

                radius_projection = plan_radius/math.sqrt(2)

                for i in range(2):
                    # spiral center X+ Y+

                    if ( reb.XYZ[0] > 0.):
                        aux_X =  1.
                    else:
                        aux_X = -1.

                    if ( reb.XYZ[1] > 0.):
                        aux_Y =  1.
                    else:
                        aux_Y = -1.

                    if ( (i == 0) and flag_inner_S ):
                        topology_origin = [reb.XYZ[0] - aux_X*radius_projection, reb.XYZ[1] - aux_Y*radius_projection, 0.]

                    elif ( (i == 1) and flag_outer_S ):
                        topology_origin = [reb.XYZ[0] + aux_X*radius_projection, reb.XYZ[1] + aux_Y*radius_projection, 0.]

                    topology_end = [topology_origin[0], topology_origin[1], self.globvar.s]
                    topology = [topology_origin, topology_end]

                    rebar_V = rebar(globvar = self.globvar, XYZ = topology, cs_nr = 4, cs_tag = self.globvar.DV, cs_mat = 3, tag = 'V_'+reb.tag)
                    #rebar_V.check_topology_consistency(self.globvar)
                    self.globvar.rebars.append(rebar_V)



        # large vs. small spirals - inner or outer intersections
        # eg. S1 vs. L1 etc.
        # flag_inner_L_vs_S = True
        flag_inner_L_vs_S = False

        for reb_L in self.globvar.rebars:
            # outer loop over all large spirals and find all potential intersecting small spirals
            if ( (reb_L.give_rebar_type() == 'spiral') and (reb_L.tag[0] == 'L') ):
                plan_radius_L = reb_L.radius
                #math.sqrt( (reb_L.start[0]-reb_L.XYZ[0])**2 + (reb_L.start[1]-reb_L.XYZ[1])**2 )

                if flag_inner_L_vs_S:
                    r_L_to_V = plan_radius_L - reb_L.give_diameter(self.globvar)/2. - DV/2.
                else:
                    r_L_to_V = plan_radius_L + reb_L.give_diameter(self.globvar)/2. + DV/2.

                # inner loop over inner spirals
                for reb_S in self.globvar.rebars:
                    # outer loop over all large spirals and find all potential intersecting small spirals
                    if ( (reb_S.give_rebar_type() == 'spiral') and (reb_S.tag[0] == 'S') ):
                        plan_radius_S = reb_S.radius
                        # plan_radius_S = math.sqrt( (reb_S.start[0]-reb_S.XYZ[0])**2 + (reb_S.start[1]-reb_S.XYZ[1])**2)
                        # calculate distance center_L to center_S
                        center_to_center = math.sqrt( (reb_L.XYZ[0]-reb_S.XYZ[0])**2 + (reb_L.XYZ[1]-reb_S.XYZ[1])**2 )

                        if (center_to_center < plan_radius_L + plan_radius_S - reb_L.give_diameter(self.globvar)/2. - reb_S.give_diameter(self.globvar)/2. - DV):
                            r_S_to_V = plan_radius_S - reb_S.give_diameter(self.globvar)/2. - DV/2.

                            DD = center_to_center/math.sqrt(2.)

                            alpha = math.acos ( (2*DD**2 + r_L_to_V**2 - r_S_to_V**2 )/ (2.*math.sqrt(2.) * DD * r_L_to_V) );
                            beta = math.pi/4.-alpha;

                            D_sin = r_L_to_V*math.sin(beta);
                            D_cos = r_L_to_V*math.cos(beta);

                            if ( reb_S.XYZ[0] > reb_L.XYZ[0] ):
                                aux_X = 1.
                            else:
                                aux_X = -1.

                            if ( reb_S.XYZ[1] > reb_L.XYZ[1] ):
                                aux_Y = 1.
                            else:
                                aux_Y = -1.

                            ## nr 1
                            topology_origin = [reb_L.XYZ[0] + aux_X*D_sin, reb_L.XYZ[1] + aux_Y*D_cos, 0.]
                            topology_end = [topology_origin[0], topology_origin[1], self.globvar.s]
                            topology = [topology_origin, topology_end]

                            rebar_V = rebar(globvar = self.globvar, XYZ = topology, cs_nr = 4, cs_tag = self.globvar.DV, cs_mat = 3, tag = 'V_' + reb_L.tag + reb_S.tag)
                            #rebar_V.check_topology_consistency(self.globvar)
                            self.globvar.rebars.append(rebar_V)

                            ## nr 2
                            topology_origin = [reb_L.XYZ[0] + aux_X*D_cos, reb_L.XYZ[1] + aux_Y*D_sin, 0.]
                            topology_end = [topology_origin[0], topology_origin[1], self.globvar.s]
                            topology = [topology_origin, topology_end]

                            rebar_V = rebar(globvar = self.globvar, XYZ = topology, cs_nr = 4, cs_tag = self.globvar.DV, cs_mat = 3, tag = 'V_' + reb_L.tag + reb_S.tag)
                            #rebar_V.check_topology_consistency(self.globvar)
                            self.globvar.rebars.append(rebar_V)

    def annotate_task_results(self,task):

        text = (f"e\u0302={(task.eccentricity_normalized[1]):.2f}\n[{(task.max_MN[0]):.2f}, {(task.max_MN[1]):.2f}]" )

        dN = -2.
        dM = -0.25

        self.diagram_canvas.axes.annotate(
        text,
        fontsize=8,
        xy=(task.max_MN[0], task.max_MN[1]), xycoords='data',
        xytext=(task.max_MN[0]+dM,task.max_MN[1]+dN),
        arrowprops=dict(arrowstyle="->",
        connectionstyle="arc3,rad=.2"))


    def compute_lateral_reinforcement_ratio(self):

        vol_concrete = self.globvar.Bx * self.globvar.By * self.globvar.s
        vol_steel = 0.

        for reb in self.globvar.rebars:
            if ( reb.give_rebar_type() == 'spiral') :

                vol_steel += reb.compute_rebar_volume(self.globvar)

        rho = vol_steel / vol_concrete

        return rho

    def compute_lateral_reinforcement_ratio_tie(self):

        vol_concrete = self.globvar.Bx * self.globvar.By * self.globvar.s

        vol_steel = ( (self.globvar.corner_radius * 2 * math.pi +
                     4 * (self.globvar.Bx - 2*self.globvar.corner_radius - 2*self.globvar.cover))
            * ( (math.pi * self.globvar.rebars_CS[self.globvar.DT].diam**2) / 4 )
        )

        rho = vol_steel / vol_concrete

        return rho


    def compute_vertical_reinforcement_ratio(self):
        vol_concrete = self.globvar.Bx * self.globvar.By * self.globvar.s
        vol_steel = 0.

        for reb in self.globvar.rebars:
            if (reb.give_rebar_type() == 'rebar' and
                reb.tag is not None and
                reb.tag[0] == 'V'):
                vol_steel += reb.compute_rebar_volume(self.globvar)

        rho = vol_steel / vol_concrete
        return rho

    def generate_mesh(self):

        ### generate concrete mesh
        nodes_c, elem_c, ndofman_c , nelem_c = generate_concrete_mesh(self.globvar)

        # save variables
        self.globvar.nodes_concrete = nodes_c
        self.globvar.elements_concrete = elem_c
        self.globvar.ndofman_concrete = ndofman_c
        self.globvar.nelem_concrete = nelem_c

        ### generate rebar mesh
        nodes_r , elem_r, ndofman_r, nelem_r = generate_rebars_mesh(self.globvar)

        #save variables
        self.globvar.nodes_rebar = nodes_r
        self.globvar.elements_rebar = elem_r
        self.globvar.ndofman_rebar = ndofman_r
        self.globvar.nelem_rebar = nelem_r

        write_oofem_input(self.globvar)

    def generate_inputs(self):

        # create project path
        project = Path(self.globvar.project_name)

        # check if folder already exists
        if project.exists() and self.globvar.flag_problem_changed :
            # ask if user wants to delete old files
            question = QMessageBox.question(self, "Directory Already Exists",
                                            f"Folder '{project}' already exists.\nOverwrite?")
            if question != QMessageBox.Yes:
                logger.info("Input(s) generation cancelled")
                return

            # delete old files
            shutil.rmtree(project)

        # create new project folder
        project.mkdir(parents=True, exist_ok=True)

        self.Bp_y_values = []

        for task in self.globvar.tasks:

            if task.status == Task_status.SELECTED:
                load_plate_length = task.eccentricity_normalized[1]

                # modify actual eccentricity
                task.eccentricity_actual = [ task.eccentricity_normalized[0] * self.globvar.Bx, task.eccentricity_normalized[1] * self.globvar.By ]
                self.Bp_y_values.append(load_plate_length)
                self.globvar.load_plate_length = load_plate_length
                self.globvar.Bp_y = load_plate_length * self.globvar.By

                # create path
                folder_name = f"{load_plate_length:.3f}"
                task_folder = project / folder_name
                self.globvar.file_path = task_folder
                task.file_path = task_folder

                if task_folder.exists():
                    logger.info(f"Input for configuration {folder_name} already exist. Skip")
                    continue

                task_folder.mkdir(parents=True, exist_ok=True)

                logger.info(f"Generating oofem input, normalized loading plate depth: {load_plate_length:.3f}")

                self.generate_mesh()

        # CHANGE FLAG
        self.globvar.flag_output_generated = True
        self.globvar.flag_problem_changed = False

    def compute_total_confined_area(self):

        A_conf = 0.
        # summ all spiral areas
        for reb in self.globvar.rebars:
            if ( reb.give_rebar_type() == 'spiral') :
                A_conf += self.compute_circle_area(reb)

                # subtract all multiple-confined areas
                for reb_intersect in self.globvar.rebars:
                    if ( reb_intersect.give_rebar_type() == 'spiral'):
                        # to subtract the intersection only once
                        if (self.globvar.rebars.index(reb_intersect) > self.globvar.rebars.index(reb) ):

                            A_conf -= self.compute_intersecting_area_of_two_circles(reb,reb_intersect)

        return A_conf


    def compute_circle_area(self,circ):
        # compute area - axial-wise

        A_circle = math.pi * circ.radius**2

        return A_circle


    def compute_intersecting_area_of_two_circles(self,circ_1,circ_2):
        # assuming circles in horizontal plane
        # computed axial-wise
        # https://mathworld.wolfram.com/Circle-CircleIntersection.html

        # distance between centers
        d = math.sqrt( (circ_1.XYZ[0]-circ_2.XYZ[0])**2 + (circ_1.XYZ[1]-circ_2.XYZ[1])**2 )
        # radius of first and second circles
        R1 = circ_1.radius
        R2 = circ_2.radius

        A_intersect = 0.

        if ( d < R1+R2 ):
            d1 = (d**2 - R2**2 + R1**2) / (2.*d)
            d2 = d - d1

            A_intersect += R1**2 * math.acos(d1/R1) - d1 * math.sqrt( R1**2 - d1**2 )
            A_intersect += R2**2 * math.acos(d2/R2) - d2 * math.sqrt( R2**2 - d2**2 )

        return A_intersect



    def compute_unconfined_area(self):

        area_unconf = (self.globvar.Bx * self.globvar.By) - self.compute_confined_area()
        return area_unconf


    def set_project_name(self, s):

        try:
            name = str(s)
        except ValueError:
            warnings.warn("unsupported format of project name")
            name = "test"

        self.globvar.project_name = name
        logger.info(f"Setting project name to  = " + self.globvar.project_name)


    def select_oofem_folder(self):
        folder = str(QFileDialog.getExistingDirectory(self, "Select OOFEM Directory"))
        self.set_oofem_folder(folder)


        #self.globvar.oofem_folder =


    def set_oofem_folder(self, s):

        try:
            folder = str(s)
        except ValueError:
            warnings.warn("Unsupported format of oofem folder.")
            return

        self.globvar.oofem_folder = folder
        logger.info(f"Setting OOFEM folder to  = " + self.globvar.oofem_folder)


    def select_loading_file(self):
        #loading_file = "/home/pedro/Programs/oofem_official/python/"
        # temporary comment
        dialog = QFileDialog(self)
        dialog.setNameFilter(str("*.csv, tab-separated values [M,N] (*.csv)"))
        dialog.setViewMode(QFileDialog.Detail)
        dialog.setFileMode(QFileDialog.ExistingFile)
        if dialog.exec_():
            loading_file = dialog.selectedFiles()
        #loading_file = str(QFileDialog.getExistingDirectory(self, "Select *.csv file with [M,N] loading combinations"))
        #self.set_loading_file(loading_file[0])
        self.globvar.loading_file = loading_file[0]


    def set_loading_file(self, s):

        try:
            loading_file = str(s)
        except ValueError:
            warnings.warn("Unsupported format of loading file.")
            return

        if Path(loading_file).is_file():

            self.globvar.flag_loading_selected = True
            self.globvar.loading_file = loading_file
            logger.info(f"Setting loading file to  = " + self.globvar.loading_file)
            self.globvar.loading = pd.read_csv(loading_file, sep='\t', header=0, na_values=['nan'])
            self.checkbox_show_loading.setCheckState(Qt.Checked)
        else:

            self.globvar.flag_loading_selected = False
            self.globvar.loading_file = None
            self.globvar.loading = None

            self.checkbox_show_loading.setCheckState(Qt.Unchecked)

    def set_cpu_nr(self, s):

        try:
            cpu_nr = int(s)
        except ValueError:
            warnings.warn("unsupported format of cpu number")
            cpu_nr = 4

        self.globvar.cpu_nr = cpu_nr

    def update_ecc_selection(self):



        for checkbox in self.fem_checkboxes:
            index = self.fem_checkboxes.index(checkbox)

            checkbox_state = checkbox.checkState()
            current_task_status = self.globvar.tasks[index].status

            if (checkbox_state == Qt.Unchecked):
                if (current_task_status == Task_status.SELECTED):
                    self.globvar.tasks[index].status = Task_status.UNSELECTED
                else:
                    continue

            # Qt.Checked
            else:
                # was not selected before -> no results
                if (current_task_status == Task_status.UNSELECTED):
                    self.globvar.tasks[index].status = Task_status.SELECTED

        values = [
            task.eccentricity_normalized[1]
            for task in self.globvar.tasks
            if task.status == Task_status.SELECTED
        ]

    def report_progress(self, step, load, task):
        self.console.appendPlainText(f"Finished time step {(step)} at load level = {(load):.2f} MN")
        self.progressBar.setValue(int((load / self.globvar.estimate) * 100 ))


    def fail_dialog(self, message):
        self.warning_dialog_ok(message)
        self.console.setPlainText(f"-----------------------------------\n{message}\n-----------------------------------")

    def set_FEM_results(self, task):

        self.progressBar.setValue(0.)

        for load in task.load_level:
            task.N.append(load)
            e = ( 0.5 - task.eccentricity_normalized[1] / 2 ) * self.globvar.Bx
            task.M.append( load * e )
            task.d.append(task.eccentricity_actual[1])
            task.fb.append( load / (task.eccentricity_actual[1] * self.globvar.Bp))

        task.max_MN = [task.max_load * ( self.globvar.Bx / 2 - task.eccentricity_actual[1] / 2) , task.max_load]
        task.max_fb_d = [task.max_load / (task.eccentricity_actual[1] * self.globvar.Bp) ,task.eccentricity_actual[1] ]
        self.update_diagram_plot()


    def execute_next_task(self,tasks):



        searching = True

        while (searching):
            task = next(tasks,"end")
            if (task == "end"):
                self.write_console_report()
                self.progressBar.setValue(0)
                logger.info("Finished all selected tasks")
                return
            if (task.status == Task_status.SELECTED):
                searching = False

        logger.info(f"Starting FEM analysis for loading plate depth = {(task.eccentricity_normalized[1]):.3f} [-]")
        self.console.setPlainText(f"   Starting FEM analysis for loading plate depth = {(task.eccentricity_normalized[1]):.3f} [-]")

        if self.globvar.topology_layout == 1:
            N_max, M, fb, d = AS_TIE(task.eccentricity_actual[1], self.globvar)

        elif self.globvar.topology_layout == 0:
            N_max, M, fb, d = AS_MSR(task.eccentricity_actual[1], self.globvar)

        else:
            raise ValueError

        estimate = N_max
        self.console.appendPlainText(f"   estimated load =  {(estimate):.3f} MN")
        self.console.appendPlainText("---------------------------------------------------------------------------")
        self.globvar.estimate = estimate

        self.thread = QThread()
        self.worker = Worker(task,self.globvar.oofem_folder,self.globvar.cpu_nr, self.globvar.project_name)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.worker.progress.connect(self.report_progress)
        self.worker.fail.connect(self.fail_dialog)

        self.worker.result.connect( self.set_FEM_results )
        self.thread.finished.connect( lambda: self.execute_next_task(tasks) )

        self.thread.start()

    def run_oofem_problems(self):
        # multiprocessing - works but window freezes
        #https://superfastpython.com/run-function-in-new-process/
        #maybe in the future
        #https://stackoverflow.com/questions/26833093/how-to-terminate-qthread-in-python
        #https://realpython.com/python-pyqt-qthread/#using-qthread-to-prevent-freezing-guis

        logger.info("Started FEM analysis")

        if ( len(self.globvar.oofem_folder) ):
            oofem = Path(self.globvar.oofem_folder)

            if ( oofem.exists() ):
                self.globvar.flag_oofem_selected = True
            else:
                self.globvar.flag_oofem_selected = False

        if ( not self.globvar.flag_oofem_selected ):
            message = "Invalid or not selected OOFEM folder."
            warnings.warn(message)
            self.warning_dialog_ok(message)
            return

        elif (not self.globvar.flag_output_generated):
            message = "FEM input files have not been generated."
            warnings.warn(message)
            self.warning_dialog_ok(message)
            return

        else:

            self.thread = QThread()
            tasks = iter(self.globvar.tasks)
            self.execute_next_task(tasks)
            self.globvar.flag_analyses_run = True


    def change_status_to_progress(self,task):
        task.status = Task_status.PROGRESS

    def write_console_report(self):

        self.console.setPlainText("     SUMMARY OF FEM RESULTS:   ")
        self.console.appendPlainText("-----------------------------------------------------------------------")

        self.console.appendPlainText("d* [-]\tM [MNm]\tN [MN]\t% estimate")


        for task in self.globvar.tasks:
            if ( task.status == Task_status.COMPLETED ):

                if self.globvar.topology_layout == 1:
                    N_max, M, fb, d = AS_TIE(task.eccentricity_actual[1], self.globvar)

                elif self.globvar.topology_layout == 0:
                    N_max, M, fb, d = AS_MSR(task.eccentricity_actual[1], self.globvar)

                estimate = N_max
                self.console.appendPlainText(f"{(task.eccentricity_normalized[1]):.3f}\t{(task.max_MN[0]):.3f}\t{(task.max_MN[1]):.3f}\t{(100.*(task.max_load/estimate)):.1f}")


        self.console.appendPlainText("-----------------------------------------------------------------------")


if __name__ == "__main__":
    # LOGGER INITIALIZATION
    # logging warnings
    #  https://code-maven.com/python-warnings



    logger = logging.getLogger('py.warnings')
    logger.setLevel(logging.DEBUG)
    logging.captureWarnings(True)

    # console output
    sh = logging.StreamHandler()
    sh.setLevel(logging.DEBUG)
    sh.setFormatter(logging.Formatter('%(asctime)s %(levelname)-8s: %(message)s',
                                  "%Y-%m-%d %S:%M:%S") )
    # file output
    fh = logging.FileHandler(filename='multibear.log')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(asctime)s %(levelname)-8s: %(message)s',
                                  "%Y-%m-%d %S:%M:%S") )
    # add both outputs
    logger.addHandler(fh)
    logger.addHandler(sh)

    logger.info("Program started")



    app = QApplication(sys.argv)
    start = time()

    pixmap = QPixmap("multibear_logo_full_v2.png")

    message_top = "(c) Petr Janas, Petr Havlásek 2026"
    message_bottom = ('Acknowledgment: This software was created with the state '
                      'support of the Technology Agency of the Czech Republic under the '
                      'DELTA 2 Programme. TAČR project nr. TM04000013\n')

    if not pixmap.isNull():
        painter = QPainter(pixmap)
        painter.setPen(Qt.black)

        # Draw top text
        painter.drawText(pixmap.rect(), Qt.AlignLeft | Qt.AlignTop, message_top)

        # Make the font bold
        font = painter.font()
        font.setBold(True)
        #font.setPointSize(12)
        painter.setFont(font)


        # Draw bottom text
        painter.drawText(pixmap.rect(), Qt.AlignCenter | Qt.AlignBottom | Qt.TextWordWrap, message_bottom)

        painter.end()

    splash = QSplashScreen(pixmap, Qt.WindowStaysOnTopHint)
    splash.resize(1280,758)
    splash.setGeometry(
        QStyle.alignedRect(
        Qt.LeftToRight,
        Qt.AlignCenter,
        splash.size(),
        QGuiApplication.primaryScreen().availableGeometry(),
        ),
    )
    splash.show()

    while time() - start < 2.:
        sleep(0.001)
        app.processEvents()

    w = MainWindow()
    w.setFixedSize(1280, 720)
    w.setGeometry(
        QStyle.alignedRect(
        Qt.LeftToRight,
        Qt.AlignCenter,
        w.size(),
        QGuiApplication.primaryScreen().availableGeometry(),
        ),
    )
    w.show()

    splash.finish(w)

    app.setWindowIcon(QIcon('multibear_logo.png'))

    #splash.raise_()
    sys.exit(app.exec_())

    #logger.info("Program finished")
    #app.closeAllWindows()
    #app.exec_()


