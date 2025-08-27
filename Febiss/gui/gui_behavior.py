#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

import sys

import os
from .gui import Ui_Dialog
from PyQt6.QtWidgets import QDialog, QFileDialog
from Febiss.cpptraj_interface.gist import GistAnalyser
from Febiss.plotting.display import Plot
from Febiss.solvents import SOLVENT_DICT

class Window(QDialog):

    def __init__(self, mode, parent=None):

        QDialog.__init__(self, parent)

        # GUI mode of febiss_settings or febiss. 0 ... febiss_settings, 1 ... febiss
        self.mode = mode

        # Default solvent
        self.solvent_name = "Acetonitrile"

        # Instantiate GistAnalyser and Plot object
        self.analyser = GistAnalyser(**{
            'top': "Choose File ...",
            'trajectory_file': "Choose File ...",
            'solv_abb': SOLVENT_DICT[self.solvent_name][0],
            'com' : False,
            'solv_file': os.path.abspath(os.path.join(__file__, "../../solvents/{0}".format(SOLVENT_DICT[self.solvent_name][1]))),
            'refdens': SOLVENT_DICT[self.solvent_name][3],
            'ref_evv': SOLVENT_DICT[self.solvent_name][4],
            'rigid_atom_0': SOLVENT_DICT[self.solvent_name][2][1],
            'rigid_atom_1': SOLVENT_DICT[self.solvent_name][2][0],
            'rigid_atom_2': SOLVENT_DICT[self.solvent_name][2][2]})

        self.display = Plot()

        # Create an instance of the GUI
        self.ui = Ui_Dialog()

        # Run the .setupUi() method to show the GUI
        self.ui.setupUi(self)
        if self.mode == 0:
            self.ui.Run.setText('Save')

        # Set defaults
        self.set_simulation_values()
        self.set_solvent_values(self.solvent_name)
        self.set_cpptraj_values()
        self.set_plotting_values()

        # behavior
        self.ui.Browse_trajectory.clicked.connect(self.getTrajFile)
        self.ui.Browse_top.clicked.connect(self.getTopFile)
        self.ui.Browse_solvent.clicked.connect(self.getSolventFile)
        self.ui.Dropdown_solvent.currentIndexChanged.connect(self.changeValues)
        self.ui.com.stateChanged.connect(self.set_and_freeze_center)
        self.ui.Cancel.clicked.connect(self.on_cancel)
        self.ui.Run.clicked.connect(self.on_run)

        self.show()

    def set_simulation_values(self):
        self.ui.top.setText(self.analyser.top)
        self.ui.trajectory_file.setText(self.analyser.trajectory_file)
        self.ui.frame_selection.setPlaceholderText("e.g.: 1 1000 5")

    def set_solvent_values(self, solvent_name):
        self.ui.solv_abb.setText(SOLVENT_DICT[solvent_name][0])
        if not solvent_name == 'Custom':
            self.ui.solv_file.setText(os.path.abspath(os.path.join(__file__, "../../solvents/{0}".format(SOLVENT_DICT[solvent_name][1]))))
        else:
            self.ui.solv_file.setText("Choose File ...")
        self.ui.refdens.setText(str(SOLVENT_DICT[solvent_name][3]))
        self.ui.ref_evv.setText(str(SOLVENT_DICT[solvent_name][4]))
        if not self.ui.com.isChecked(): # Value shall only change upon choosing a different solvent
            self.ui.rigid_atom_0.setProperty('value', SOLVENT_DICT[solvent_name][2][1])
        self.ui.rigid_atom_1.setProperty('value', SOLVENT_DICT[solvent_name][2][0])
        self.ui.rigid_atom_2.setProperty('value', SOLVENT_DICT[solvent_name][2][2])

    def set_cpptraj_values(self):
        # Grid
        self.ui.grid_center.setPlaceholderText("e.g.: (10.0, 5.0, 0.0)")
        self.ui.grid_lenghts.setText(str(self.analyser.grid_lengths))
        self.ui.grid_spacing.setText(str(self.analyser.grid_spacing))

        #RDF
        self.ui.rdf_elements.setPlaceholderText("e.g.: C,N,O,center")
        self.ui.rdf_maximum.setText(str(self.analyser.rdf_maximum))
        self.ui.rdf_spacing.setText(str(self.analyser.rdf_spacing))

        # Solute residues
        self.ui.solute_residues.setText(str(self.analyser.solute_residues))

        # I/O
        self.ui.gist_cpptraj_command_file.setText(str(self.analyser.gist_cpptraj_command_file))
        self.ui.gist_grid_file.setText(str(self.analyser.gist_grid_file))
        self.ui.gist_out_file.setText(str(self.analyser.gist_out_file))

    def set_plotting_values(self):
        # Solvent Selection
        self.ui.cutoff1.setText(str(self.display.cutoff1))
        self.ui.cutoff2.setText(str(self.display.cutoff2))
        self.ui.displayed_solvents.setText(str(self.display.displayed_solvents))
        self.ui.solvent_selection.setPlaceholderText("e.g.: 1,3,4,5")
        self.ui.within.setText(str(self.display.within))
        self.ui.between.setText(str(self.display.between))
        self.ui.outside.setText(str(self.display.outside))
        self.ui.selected.setText(str(self.display.selected))
        self.ui.mark.setText(str(self.display.mark))

        #Bar Plot Settings
        self.ui.width.setText(str(self.display.width))
        self.ui.height.setText(str(self.display.height))
        self.ui.dpi.setText(str(self.display.dpi))
        self.ui.xlabel.setText(str(self.display.xlabel))
        self.ui.ylabel.setText(str(self.display.ylabel))
        self.ui.plotname.setText(str(self.display.plotname))
        self.ui.selected_plotname.setText((str(self.display.selected_plotname)))
        self.ui.fontsize.setText(str(self.display.fontsize))
        self.ui.number_xtics.setText(str(self.display.number_xtics))
        self.ui.y_numbers.setText(str(self.display.y_numbers))
        self.ui.febiss_file.setText(str(self.display.febiss_file))
        #self.ui.marks.setPlaceholderText("e.g.: [3,10]")

    def getTrajFile(self):
        response = QFileDialog.getOpenFileName(self, caption='Select trajectory file',
                                               filter="Trajectory (*.nc);;All Files (*)")
        self.ui.trajectory_file.setText(str(response[0]))

    def getTopFile(self):
        response = QFileDialog.getOpenFileName(self, caption='Select topology file',
                                               filter="Topology (*.prmtop);;All Files (*)")
        self.ui.top.setText(str(response[0]))

    def getSolventFile(self):
        response = QFileDialog.getOpenFileName(self, caption='Select solvent file',
                                               directory=os.path.abspath(os.path.join(__file__, "../../../solvents/")),
                                               filter="XYZ (*xyz);MOL2 (*mol2);;All Files (*)")
        self.ui.solv_file.setText(str(response[0]))

    def changeValues(self):
        self.solvent_name = self.ui.Dropdown_solvent.currentText()
        self.set_solvent_values(self.solvent_name)

    def set_and_freeze_center(self):
        if self.ui.com.isChecked():
            self.ui.rigid_atom_0.setMinimum(-1)
            self.ui.rigid_atom_0.setProperty('value', -1)
            self.ui.rigid_atom_0.setReadOnly(True)
        else:
            self.ui.rigid_atom_0.setMinimum(0)
            self.ui.rigid_atom_0.setProperty('value', SOLVENT_DICT[self.solvent_name][2][1])
            self.ui.rigid_atom_0.setReadOnly(False)

    def on_cancel(self):
        sys.exit()


    def on_run(self):
        # Write values of window back to GistAnalyser and Plot objects
        # analyser
        self.analyser.signal = True
        self.analyser.top = self.ui.top.text()
        self.analyser.trajectory_file = self.ui.trajectory_file.text()
        self.analyser.com = self.ui.com.isChecked()
        self.analyser.solv_file = self.ui.solv_file.text()
        self.analyser.solv_name =  self.ui.Dropdown_solvent.currentText()
        self.analyser.rigid_atom_0 = self.ui.rigid_atom_0.value()
        self.analyser.rigid_atom_1 = self.ui.rigid_atom_1.value()
        self.analyser.rigid_atom_2 = self.ui.rigid_atom_2.value()
        self.analyser.ref_evv = self.ui.ref_evv.text()
        self.analyser.solv_abb = self.ui.solv_abb.text()
        self.analyser.frame_selection = None if self.ui.frame_selection.text() == '' else self.ui.frame_selection.text()
        self.analyser.grid_center = None if self.ui.grid_center.text() == '' else self.ui.grid_center.text()
        self.analyser.grid_spacing = self.ui.grid_spacing.text()
        self.analyser.grid_lengths = self.ui.grid_lenghts.text()
        self.analyser.refdens = self.ui.refdens.text()
        self.analyser.solute_residues = self.ui.solute_residues.text()
        self.analyser.gist_cpptraj_command_file = self.ui.gist_cpptraj_command_file.text()
        self.analyser.gist_out_file = self.ui.gist_out_file.text()
        self.analyser.gist_grid_file = self.ui.gist_grid_file.text()
        self.analyser.rdf_elements = None if self.ui.rdf_elements.text() == '' else self.ui.rdf_elements.text()
        self.analyser.rdf_maximum = self.ui.rdf_maximum.text()
        self.analyser.rdf_spacing = self.ui.rdf_spacing.text()

        # display
        self.display.cutoff1 = self.ui.cutoff1.text() #distance between solute and solvent
        self.display.cutoff2 = self.ui.cutoff2.text() #distance between solute and solvent
        self.display.displayed_solvents = self.ui.displayed_solvents.text()
        self.display.within = self.ui.within.text()
        self.display.between = self.ui.between.text()
        self.display.outside = self.ui.outside.text()
        self.display.selected = self.ui.selected.text()
        self.display.mark = self.ui.mark.text()
        self.display.width = self.ui.width.text()
        self.display.height = self.ui.height.text()
        self.display.dpi = self.ui.dpi.text()
        self.display.xlabel = self.ui.xlabel.text()
        self.display.ylabel = self.ui.ylabel.text()
        self.display.plotname = self.ui.plotname.text()
        self.display.selected_plotname = self.ui.selected_plotname.text()
        #self.display.orientation = self.ui.orientation.currentText()
        self.display.fontsize = self.ui.fontsize.text()
        self.display.number_xtics = self.ui.number_xtics.text()
        self.display.y_numbers = self.ui.y_numbers.text()
        #self.display.marks = []
        self.display.febiss_file = self.ui.febiss_file.text()
        self.display.transparent = self.ui.transparent.isChecked()
        self.display.solvent_selection = None if self.ui.solvent_selection.text() == '' else self.ui.solvent_selection.text()

        # Check values
        self.analyser._sanity_check()
        self.display.sanity_checks(True)

        # Set Run available only after successful evaluation of all settings.
        if self.analyser.err_string == '' and self.display.err_string == '':
            self.close()

        else:
            print('The following settings did not pass the checks:'+self.analyser.err_string+self.display.err_string)
            self.analyser.err_string = ''
            self.display.err_string = ''