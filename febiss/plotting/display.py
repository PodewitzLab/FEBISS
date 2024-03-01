#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

from collections import OrderedDict
from matplotlib.ticker import FormatStrFormatter
from matplotlib.widgets import Button
from matplotlib.widgets import TextBox #new LM20240229
from typing import List
from warnings import warn
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import sys

from ..utilities.structures import Solute, Solvent, Reference
from ..utilities.distance_functions import distance_squared
from ..utilities.io_handling import write_pdb
from .rdf import ButtonActions


class Plot:
    def __init__(self, **kwargs):
        self._set_defaults()
        self.allowed_keys = {'cutoff1', 'cutoff2', 'displayed_solvents', 'colors', 'width', 'height', 'dpi', 'xlabel',
                             'ylabel', 'selected_plotname', 'plotname', 'orientation', 'fontsize', 'number_xtics',
                             'y_numbers', 'marks', 'transparent', 'display_once', 'testing', 'febiss_file', 'rdf_names',
                             'solvent_selection'}
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in self.allowed_keys)
        for key in kwargs.keys():
            if key not in self.allowed_keys:
                warn('WARNING: Did not recognize key: ' + str(key))

    def gui(self, abb, solute: Solute, solvent: Solvent, reference: Reference) -> str:
        if not self.testing and self.solvent_selection is None:
            self._determine_hetero_elements(solute)
            barcolors = self._determine_colors(solute, solvent)
            self._create_plot(barcolors, solvent, False)
        elif not self.testing:
            self._input_selection(self.solvent_selection)
        else:
            self.selected_solvents = [0, 1, 2, 3, 4]

        if self.display_once:
            sys.exit()
        # avoid bug of selecting out of range solvent
        if -1 in self.selected_solvents:
            self.selected_solvents.remove(-1)
        self._interactive_reselection(solute, solvent)
        filename = self._save_selection(abb, solute, solvent, reference)
        return filename

    def _set_defaults(self):
        # default values for bar chart
        self.selected_solvents = [] #get written in class ClickableBar #LM20231123
        self.cutoff1 = 3.0 #distance between solute and solvent
        self.cutoff2 = 6.0 #distance between solute and solvent
        self.displayed_solvents = 50
        self.colors = {'within': '#fdb462', 'between': '#80b1d3', 'outside': '#de2d26', 'selected': '#b3de69',
                       'mark': '#bc80bd'}
        self.width = 16
        self.height = 8
        self.dpi = 200
        self.xlabel = 'ID of solvent molecule (click bar or enter ID to select)' #changed LM20240301
        self.ylabel = '$-$ Free Energy / kcal mol$^{-1}$'
        self.selected_plotname = 'febiss-plot-selected.png'
        self.plotname = 'febiss-plot.png'
        self.orientation = 'landscape'
        self.fontsize = 18
        self.number_xtics = 10
        self.y_numbers = 0.25
        self.marks = []
        self.febiss_file = 'febiss.dat' #changed from febiss-solvents.pdb LM20231123
        self.rdf_names = {'center2': 'center of solute', 'C': 'carbon', 'O': 'oxygen', 'N': 'nitrogen', 'P': 'phosphor'}
        self.transparent = True
        self.display_once = False
        self.testing = False
        self.solvent_selection = None #new LM20240301
        self.drs = [] #not in allowed keys. replaces drs used for bar interaction below #LM20240229

    def _determine_hetero_elements(self, solute: Solute):
        self.existing_elements = []
        for symbol in self.rdf_names.keys():
            if symbol in solute.elements:
                self.existing_elements.append(True)
            elif symbol == 'center2':
                self.existing_elements.append(True)
            else:
                self.existing_elements.append(False)

    def _determine_colors(self, solute: Solute, solvent: Solvent) -> List[str]:
        """ set cutoffs """
        squared_cutoff1 = self.cutoff1 ** 2
        squared_cutoff2 = self.cutoff2 ** 2
        # if distance between two solvent atoms below this,
        # they belong to same solvent molecule, but they have to be in pdb file within 2 rows)
        # same_solvent_cutoff = 1.6 #not needed since per solvent only the com is given LM20231123
        #squared_same_solvent_cutoff = same_solvent_cutoff ** 2 #not needed since per solvent only the com is given LM20231123
        #if len(solute.polars) == 0: #not needed since the distance between solvent and solute will be measured between solvent COM and the nearest solute atom regardless of polarity
        #    solute.determine_polar_hydrogen_and_non_hydrogen() #gives all non-H atoms and polar hydrogen. Attention - unclear whether using the polar cut-off for H as done here is valid LM20231027
        within_cutoff = []
        between_cutoffs = []
        outside_cutoff2 = []
        """ cycle over solvents and assign to list for each polar solute atom """
        #skip_next = False #not needed LM20231123
        #skip_next_next = False #not needed LM20231123
        for count, sol in enumerate(solvent.coords): #this iterates through all atoms of the febiss_solvents.pdb which are marked with "HETATM" therefore either H or O information. LM20231027
            #if skip_next: #this and next 6 lines not needed LM20231123
            #    if skip_next_next:
            #        skip_next_next = False
            #        continue
            #    else:
            #        skip_next = False
            #        continue
            for atom in solute.coords: #changed from: "for polar in solute.polars:" LM20231123
                if distance_squared(atom, sol) < squared_cutoff1: #changed from polar to atom LM20231123
                    within_cutoff.append(count)
                    # relies on atoms of same solvent molecule to be right after each other in pdb!
                    #for i in range(-2, 3): #this and next 8 lines not needed because there is just one entry (the COM) per solvent and not several atoms LM20231123
                    #    # bound check and then distance
                    #    if i != 0 and len(solvent.elements) > count + i > 0 \
                    #            and distance_squared(sol, solvent.atoms[count + i]) < squared_same_solvent_cutoff:
                    #        within_cutoff.append(count + i)
                    #        if i == 1:
                    #            skip_next = True
                    #        elif i == 2:
                    #            skip_next_next = True
                    break  # close enough solute atom was found for solvent within cutoff -> break loop over solute atoms
                elif distance_squared(atom, sol) > squared_cutoff2: #changed from polar to atom LM20231123
                    outside_cutoff2.append(count)
                else:
                    between_cutoffs.append(count)



        """ delete multiple entries """
        within_cutoff = list(OrderedDict.fromkeys(within_cutoff))
        between_cutoffs = list(OrderedDict.fromkeys(between_cutoffs))
        outside_cutoff2 = list(OrderedDict.fromkeys(outside_cutoff2))
        """ remove duplicate between different lists """
        # within overrules between
        for within in within_cutoff:
            for between in between_cutoffs:
                if within == between:
                    between_cutoffs.remove(between)
        # within overrules outside (possible because distance of all solute atoms and solvents are calculated)
        for within in within_cutoff:
            for outside in outside_cutoff2:
                if within == outside:
                    outside_cutoff2.remove(outside)
        # between overrules outside
        for outside in outside_cutoff2:
            for between in between_cutoffs:
                if outside == between:
                    outside_cutoff2.remove(outside)
        """ set colors """
        barcolors = [self.colors['outside']] * len(solvent.coords) #changed from solvent.elements which contained H,H,O for one water molecule and therefore 3*nsolvent entries
        for within in within_cutoff:
            #if solvent.elements[within] == "O": #not needed since just one entry per solvent LM20231123
            #    index = int(round(within / 3)) #not needed since just one entry per solvent LM20231123
                barcolors[within] = self.colors['within']
        for outside in between_cutoffs:
            #if solvent.elements[outside] == "O": #same as above LM20231123
            #    index = int(round(outside / 3)) #same as above LM20231123
                barcolors[outside] = self.colors['between']

        return barcolors

    def _round_to_percentage(self, number: float, percentage: float) -> float:
        return round(number * percentage) / percentage

    def _axis_font(self):
        return {'size': str(self.fontsize)}

    def _create_legend(self, ax: plt.Axes, save_selected: bool):
        within = mpatches.Patch(color=self.colors['within'], label=u'd(solute-solvent) < %.1f Å' % self.cutoff1)
        between = mpatches.Patch(color=self.colors['between'],
                                 label=u'%.1f ≤ d(solute-solvent) ≤ %.1f Å' % (self.cutoff1, self.cutoff2))
        outside = mpatches.Patch(color=self.colors['outside'], label=u'd(solute-solvent) > %.1f Å' % self.cutoff2)
        if save_selected:
            # add green bars in legend
            selected = mpatches.Patch(color=self.colors['selected'], label='selected solvent molecules')
            ax.legend(handles=[within, between, outside, selected], prop=self._axis_font())
        else:
            ax.legend(handles=[within, between, outside], prop=self._axis_font())

    # creates interactive bar plot
    def _create_plot(self, barcolors: List[str], solvent: Solvent, save_selected: bool):
        plt.ioff()
        indices = np.arange(1, len(solvent.values) + 1) # x-values
        fig = plt.figure(figsize=(self.width, self.height))
        ax = fig.add_subplot(111)
        rects = ax.bar(indices, solvent.values, color=barcolors, picker=True)  # create bars

        # make bars interactive
        #drs = []  # necessary for interaction
        for rect, color in zip(rects, barcolors):
            dr = ClickableBar(rect, color, self)
            self.drs.append(dr)
        self.original_drs = self.drs

        # font specifications
        matplotlib.rcParams.update({'font.size': self.fontsize})
        if self.y_numbers == 0.5:
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        elif self.y_numbers == 1.0:
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        else:
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        # Additional plotting option for vertical lines to show cutoffs
        for x in self.marks:
            plt.axvline(x=x, color=self.colors['mark'], linestyle='--', lw=4)

        # label specifications
        plt.xlabel(self.xlabel, **self._axis_font())
        plt.ylabel(self.ylabel, **self._axis_font())

        # read user input for x-axis and determine range of y-axis
        ymax = np.max(solvent.values)
        fill_up = self._round_to_percentage(ymax, 1 / self.y_numbers)
        if fill_up < ymax:
            fill_up = self.y_numbers - (ymax - fill_up)
            ymax += fill_up
        else:
            ymax += (fill_up - ymax)

        if self.displayed_solvents == 'all':  # if user wants all solvents displayed
            xmax = len(solvent.values) + 0.5  # 0.5 ensures slight offset, so that last bar is seen completely
            ymin = np.min(solvent.values) - 0.5
        else:
            xmax = self.displayed_solvents + 0.5
            ymin = solvent.values[self.displayed_solvents]
            fill_up = self._round_to_percentage(ymin, 1 / self.y_numbers)
            if fill_up > ymin:
                fill_up = self.y_numbers - (fill_up - ymin)
                ymin -= fill_up
            else:
                ymin -= (ymin - fill_up)

        # ensures first bar is complete and touching the y-axis
        xmin = 0.5
        # determine stepsize on y-axis
        ystep = self.y_numbers
        if ystep > (ymax - ymin):
            ystep = (ymax - ymin) / 10

        # stepsize for x-axis
        xstep = int(round((xmax - xmin) / self.number_xtics))
        if xstep == 0:
            xstep = 1

        # set tics and limits
        plt.xticks(np.arange(0, xmax + xstep, step=xstep))
        if ymin == ymax == 0.0: #TODO: Change to try-except to catch if ymax == ymin. LM20231127
            ymax = 1
            ystep = 0.1
        plt.yticks(np.arange(ymin, ymax + ystep, step=ystep))
        ax.set_ylim([ymin, ymax])
        ax.set_xlim([xmin, xmax])
        ax.tick_params(axis='both', which='major', labelsize=self.fontsize)  # https://stackoverflow.com/questions/6390393/how-to-change-tick-label-font-size/11386056#11386056 (accessed 14 December 2023)

        # creates legend according with or without green bar
        self._create_legend(ax, save_selected)

        # saves figure without title and button
        if save_selected:
            plt.savefig(self.selected_plotname, dpi=self.dpi, transparent=self.transparent,
                        orientation=self.orientation)
        else:
            plt.savefig(self.plotname, dpi=self.dpi, transparent=self.transparent, orientation=self.orientation)

            # legend with green bars and title specifications for "GUI"
            self._create_legend(ax, save_selected=True)
            #plt.title('Click bars or enter solvent ID. Then close window.\n') #edited and eventually commented out LM20240301
            # specify button
            b = plt.axes((0.8, 0.9, 0.1, 0.075))  # set position
            button = Button(b, "Display RDF", color="0.85", hovercolor="0.95")  # set text and color
            callback_b = ButtonActions()  # init class with necessary functions and determine click actions
            button.on_clicked(lambda x: ButtonActions.plot_rdf(callback_b, self))

            # specify textbox
            t = plt.axes((0.125, 0.9, 0.1, 0.075)) #left, bottom, width, height
            global txt_box #needed to set txt_box to '' in self._input_selection. inspired by this: https://coderslegacy.com/python/matplotlib-textbox-widget/ (accessed 2024-03-01) LM20240301
            txt_box = TextBox(t, "IDs:")
            txt_box.on_submit(self._input_selection)

            # specify clear all
            c = plt.axes((0.5, 0.9, 0.1, 0.075))
            button_c = Button(c, "Clear all", color="0.85", hovercolor="0.95")  # set text and color
            button_c.on_clicked(self._deselect_all)

            plt.show()

    def _interactive_reselection(self, solute: Solute, solvent: Solvent):
        # save old values to tell user later
        old_cutoff1 = self.cutoff1
        old_cutoff2 = self.cutoff2
        old_displayed_solvents = self.displayed_solvents

        # while loop ensures selection or denial of selection by user
        # allows for changing distance bounds and displayed solvent number
        selected = False
        while not selected:
            if not self.selected_solvents:
                print("No solvents selected")
                inp = input("Do you want to change/set the cutoff value? [y/n] ")
                if inp in ['y', 'Y', 'yes', 'Yes', ' y']:
                    if old_displayed_solvents == "all":
                        print("Please enter the new values (previous distance cutoffs were {0:.2f} and {1:.2f} and "
                              "all solvent molecules were shown)".format(old_cutoff1, old_cutoff2))
                    else:
                        print("Please enter the new values (previous distance cutoffs were {0:.2f} and {1:.2f} and "
                              "{2:d} solvent molecules were shown)".format(old_cutoff1, old_cutoff2,
                                                                         old_displayed_solvents))

                    # ensures that user enters correct data type
                    selected_cutoff = False
                    while not selected_cutoff:
                        try:
                            self.cutoff1 = float(input("Choose first distance cutoff: "))
                            old_cutoff1 = self.cutoff1
                            selected_cutoff = True
                        except ValueError:
                            self.cutoff1 = float(input("Please enter a number: "))
                    # ensures that user enters correct data type
                    selected_cutoff = False
                    while not selected_cutoff:
                        try:
                            self.cutoff2 = float(input("Choose second distance cutoff: "))
                            old_cutoff2 = self.cutoff2
                            selected_cutoff = True
                        except ValueError:
                            self.cutoff2 = float(input("Please enter a number: "))

                    inp = input("Choose maximum number of solvent molecules displayed: ")
                    # user can enter 'all' or number -> first check all before ensuring datatype int
                    selected_solvent_number = False
                    while not selected_solvent_number:
                        if inp in ["all", "All", "ALL"]:
                            self.displayed_solvents = 'all'
                            old_displayed_solvents = 'all'
                            selected_solvent_number = True
                        else:
                            try:
                                self.displayed_solvents = int(inp)
                                selected_solvent_number = True
                            except ValueError:
                                inp = input("Please enter a number or 'all': ")

                    # user wants crazy number or is unaware of 'all' function --> print all
                    if self.displayed_solvents > len(solvent.values):
                        print('Chosen solvent number is greater than number of all placed solvent molecules.')
                        print('x-axis will stop at last solvent.')
                        self.displayed_solvents = 'all'
                        old_displayed_solvents = 'all'
                    # weird input and only one solvent is plotted
                    elif self.displayed_solvents < 1:
                        print('Chosen solvent number must be at least 1. Plotting 1 solvent molecule')
                        self.displayed_solvents = 1
                        old_displayed_solvents = 1
                    else:
                        old_displayed_solvents = self.displayed_solvents

                    # redetermine colors because of new distance bounds
                    barcolors = self._determine_colors(solute, solvent)
                    # display GUI again
                    self._create_plot(barcolors, solvent, False)
                elif inp in ['n', 'N', 'no', 'No', ' n']:
                    print("No solvents selected")
                    sys.exit()
                else:
                    print("Wrong input, just 'y' or 'n'.")
            else:  # solvent were selected, reselection of parameters are not necessary
                selected = True

    def _input_selection(self,text : str): #ids separated with ",". allows ranges with "-"
        if len(text) == 0: #this in connection with set_val('') prevents the double submission of the textbox when clicking on a bar after input selection. fyi: on_submit gets triggered with enter and with leaving the textbox
            return
        split = text.split(",")
        rm_list = []
        for i in range(len(split)):
            split[i] = split[i].split("-")
            if len(split[i]) == 1:
                try:
                    split[i] = int(split[i][0])
                except ValueError:
                    print("A non-integer was given. Please re-enter your selection!")
                    return
            else:
                rm_list.append(split[i])
                for j in range(len(split[i])):
                    try:
                        split[i][j] = int(split[i][j])
                    except ValueError:
                        print("A range containing a non-integer was given. Please re-enter your selection!")
                        return
                a = sorted(split[i])
                for k in range(a[0], a[-1] + 1):
                    split.append(k)

        for r in rm_list:
            split.remove(r)
        split = [x for x in split if x in range(1,self.displayed_solvents+1)]
        out = ','.join((str(x) for x in sorted(set(split))))
        print("Selected solvents: {0}".format(out))

        self.selected_solvents = sorted(set([x-1 for x in split]))  # overrules all previously clicked solvents

        if len(self.drs) != 0: #that is the case when there is no barplot due to predefined solvent selection
            for dr in self.drs:
                try:
                    if int(dr.rect.xy[0]) in self.selected_solvents:
                        dr.color = self.colors['selected']
                        canvas = dr.rect.figure.canvas
                        axes = dr.rect.axes
                        #dr.rect.set_animated(True)
                        canvas.draw()
                        dr.background = canvas.copy_from_bbox(dr.rect.axes.bbox)
                        dr.rect.set_color(self.colors['selected'])
                        axes.draw_artist(dr.rect)
                        canvas.blit(axes.bbox)
                except (ValueError,TypeError):
                    print("Something was wrong with re-coloring the bars upon input selection!")
                    return
            txt_box.set_val('')

    def _deselect_all(self,event):
        print("Resetting plot...")
        for dr in self.drs:
            dr.color = dr.original_color
            canvas = dr.rect.figure.canvas
            axes = dr.rect.axes
            dr.rect.set_color(dr.original_color)
            axes.draw_artist(dr.rect)
            canvas.blit(axes.bbox)
        self.selected_solvents = []
        print("Plot resetted!")

    def _save_selection(self, abb, solute: Solute, solvent: Solvent, reference: Reference) -> str:
        print('Number of solvents chosen: ' + str(len(self.selected_solvents)))
        filename = 'solvated_structure-' + str(len(self.selected_solvents)) + '.pdb'
        print('Your microsolvated structure is written to: ' + filename)
        # does not open GUI, but saves plot of selected bars
        barcolors = self._determine_colors(solute, solvent)
        for select in self.selected_solvents:
            barcolors[select] = self.colors['selected']
        self._create_plot(barcolors, solvent, True)

        # writes latest solvated structure to file to open with pymol
        with open("latest-solvation.log", "a") as latest:
            latest.write(filename + '\n')
        # writes file with solute and selected solvents
        write_pdb(filename, solute, abb, solute=True) #writes only the solute into the pdb-file LM20231123
        selected_solvent = Solvent() #info: this object finally contains all element labels, coords and values of all selected solvents. LM20231130

        for select in self.selected_solvents: #TODO: Check if selected_solvents order coincides with order of solvent.coord entries, i.e. are solvent.coord entries sorted wrt their energy
            voxel = int(solvent.data[select][0]) #new LM20231123. LM20231130: type conversion from str to int. TODO: Type conversion prone to ValueError
            quats = solvent.quats[voxel] #new LM20231123
            com = (float(solvent.data[select][1]), float(solvent.data[select][2]), float(solvent.data[select][3])) #new LM20231123. LM20231130: conversion from str to int. TODO: Prone to ValueError.
            elements, coords = reference._find_avg_solvent(voxel, quats, com, verbose=False) #new LM20231123. This finally determines the solvent to be placed.
            values = float(solvent.data[select][-1]) #new LM20231124. #LM20231130 conversion from str to float
            selected_solvent.elements.extend(elements) #changed from append which does not work since elements is a list itself. LM20231130
            selected_solvent.coords.extend(coords) #changed from select * 3 in brackets LM20231123. #changed from append which does not work since elements is a list itself. LM20231130
            #selected_solvent.atoms.append(solvent.atoms[select * 3 + 1]) #not needed since only the com coordinate is considered LM20231123
            #selected_solvent.atoms.append(solvent.atoms[select * 3 + 2])
            selected_solvent.values.extend([values]*len(elements)) #changed from select * 3 in brackets LM20231123 #needs as many entries as there are elements in the solvent. Therefore changed from append(values) to extend([values]*len(elements))
            #selected_solvent.values.append(solvent.all_values[select * 3 + 1]) #not needed since only the com coordinate is considered LM20231123
            #selected_solvent.values.append(solvent.all_values[select * 3 + 2]) #not needed since only the com coordinate is considered LM20231123
            #selected_solvent.elements.append('O') #not needed LM20231123
            #selected_solvent.elements.append('H') #not needed LM20231123
            #selected_solvent.elements.append('H') #not needed LM20231123

        write_pdb(filename, selected_solvent, abb, solute=False)
        return filename


class ClickableBar:
    lock = None  # only one can be clicked at a time

    def __init__(self, rect : matplotlib.patches.Rectangle, color, plot):
        self.rect = rect
        self.press = None
        self.background = None
        self.original_color = color
        self.color = color
        self.plot = plot
        self.cidpress = self.rect.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease = self.rect.figure.canvas.mpl_connect('button_release_event', self.on_release)

    def on_press(self, event):
        """on button press we will see if the mouse is over us and store some data"""
        if event.inaxes != self.rect.axes or ClickableBar.lock is not None:
            return
        contains, attrd = self.rect.contains(event)
        if not contains:
            return
        index = int(np.round(event.xdata)) - 1
        x0, y0 = self.rect.xy
        self.press = x0, y0, event.xdata, event.ydata
        ClickableBar.lock = self
        # draw everything but the selected rectangle and store the pixel buffer
        canvas = self.rect.figure.canvas
        axes = self.rect.axes
        self.rect.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.rect.axes.bbox)

        # now redraw just the rectangle
        # self.color is used to determine whether solvent gets selected or deselected
        # self.original_color is used to get previous color if solvent is deselected
        if self.color == self.plot.colors['selected']:
            print('solvent deselected:', int(np.round(event.xdata)))
            self.plot.selected_solvents.remove(index)
            self.rect.set_color(self.original_color)
            self.color = self.original_color
        else:
            print('solvent selected:', int(np.round(event.xdata)))
            self.plot.selected_solvents.append(index)
            self.rect.set_color(self.plot.colors['selected'])
            self.color = self.plot.colors['selected']
        axes.draw_artist(self.rect)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        """on release we reset the press data"""
        if ClickableBar.lock is not self:
            return

        self.press = None
        ClickableBar.lock = None

        # turn off the rect animation property and reset the background
        self.rect.set_animated(False)
        self.background = None

        # redraw the full figure
        self.rect.figure.canvas.draw()

    def disconnect(self):
        """disconnect all the stored connection ids"""
        self.rect.figure.canvas.mpl_disconnect(self.cidpress)
        self.rect.figure.canvas.mpl_disconnect(self.cidrelease)
