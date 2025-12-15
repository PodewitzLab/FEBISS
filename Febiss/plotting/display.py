#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

from collections import OrderedDict
from matplotlib.ticker import FormatStrFormatter
from matplotlib.widgets import Button
from matplotlib.widgets import TextBox
from typing import List
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import sys

from Febiss.cpptraj_interface.gist import GistAnalyser
from ..structures.solvent import Solvent
from ..structures.reference import Reference
from ..structures.solute import Solute
from Febiss.structures.write_pdb import write_pdb
from Febiss.utilities.check_settings import Checker
from Febiss.utilities.input import Input
from .rdf import ButtonActions


class Plot(Checker):
    def __init__(self, **kwargs):
        Checker.__init__(self)
        self._set_defaults()
        self.allowed_keys = OrderedDict({
            'cutoff1': ('type', 'float'),
            'cutoff2': ('type', 'float'),
            'displayed_solvents': ('format', '[0-9]+|(?i)all'),
            'within': ('format', '#[a-fA-F0-9]{6}'),
            'between': ('format', '#[a-fA-F0-9]{6}'),
            'outside': ('format', '#[a-fA-F0-9]{6}'),
            'selected': ('format', '#[a-fA-F0-9]{6}'),
            'mark': ('format', '#[a-fA-F0-9]{6}'),
            'width': ('type', 'int'),
            'height': ('type', 'int'),
            'dpi': ('type', 'int'),
            'xlabel': ('type', 'str'),
            'ylabel': ('type', 'str'),
            'selected_plotname': ('format', '[a-zA-Z0-9\._\-]+.png'),
            'plotname': ('format', '[a-zA-Z0-9\._\-]+.png'),
            #'orientation':('format', '(?i)landscape|portrait'),
            'fontsize': ('type', 'int'),
            'number_xtics': ('type', 'int'),
            'y_numbers': ('type', 'float'),
            #'marks': ('format', '\[[0-9]+.[0-9]+\]|\['),
            'transparent': ('type', 'bool'),
            'solvent_selection': ('format', '[0-9,\- ]*|(?i)None|^$'),
        })

        # Update attributes
        for k, v in kwargs.items():
            if k in self.allowed_keys.keys():
                if k == 'displayed_solvents':
                    if type(v) == str: # if 'all' is given
                        self.__dict__[k] = str(v).lower() # To account for misspelled "all"
                    else:
                        self.__dict__[k] = v
                elif k in ['within', 'between', 'outside', 'selected', 'mark']:
                    self.__dict__[k] = '#' + str(v) # places hashtag in front of color hex codes. yaml can't read '#'.
                                                    # conversion into str is necessary, since colors can appear as ints
                elif k == 'solvent_selection':
                    self.__dict__[k] = str(v) # converts the input into string
                else:
                    self.__dict__[k] = v
            else:
                print('WARNING: Did not recognize key: ' + str(k))

    def gui(self, analyser: GistAnalyser, solute: Solute, solvent: Solvent, reference: Reference) -> str:
        self._rdf_names = analyser._rdf_names
        if str(self.solvent_selection).lower() in ['', 'none']:
            self._determine_hetero_elements(solute)
            barcolors = self._determine_colors(solute, solvent)
            self._create_plot(barcolors, solvent, False)
        else:
            self._input_selection(self.solvent_selection)


        # avoid bug of selecting out of range solvent
        if -1 in self._selected_solvents:
            self._selected_solvents.remove(-1)
        self._interactive_reselection(solute, solvent)
        filename = self._save_selection(analyser.solv_abb, solute, solvent, reference)
        return filename

    def _distance_squared(self, a_array, b_array) -> float:
        return sum(((a - b) ** 2 for a, b in zip(a_array, b_array)))

    def _set_defaults(self):
        # default values for bar chart

        # allowed
        self.cutoff1 = 3.0 #distance between solute and solvent
        self.cutoff2 = 6.0 #distance between solute and solvent
        self.displayed_solvents = 30
        self.within = '000d98'
        self.between = '4f61ff'
        self.outside = 'd1d1ff'
        self.selected = 'b3de69'
        self.mark = 'bc80bd'
        self.width = 16
        self.height = 8
        self.dpi = 300
        self.xlabel = 'Solvent ID'
        self.ylabel = '$-$ Free Energy / kcal mol$^{-1}$'
        self.selected_plotname = 'febiss-plot-selected.png'
        self.plotname = 'febiss-plot.png'
        #self.orientation = 'landscape'
        self.fontsize = 18
        self.number_xtics = 6
        self.y_numbers = 0.5
        self.febiss_file = 'febiss.dat'
        #self.marks = str([])
        self.transparent = True
        self.solvent_selection = None

        #hidden
        self._selected_solvents = []
        self._drs = []


    def sanity_checks(self, gui = False):

        if gui:
            for k in ['within', 'between', 'outside', 'selected', 'mark']:
                self.allowed_keys[k] = ('format', '[a-fA-F0-9]{6}')

        for key, val in self.allowed_keys.items():
            if val[0] == 'path':
                if self.check_path(self.__dict__[key]):
                    continue

                else:
                    self.err_string += '\n\t-{0}. Path not found'.format(key)

            elif val[0] == 'type':
                if self.check_type(self.__dict__[key], eval(val[1])):
                    continue

                else:
                    self.err_string += '\n\t-{0}. Required type: {1}'.format(key, val[1])

            elif val[0] == 'format':
                if self.check_format(str(self.__dict__[key]), val[1]): # Conversion of __dict__[key] to string as required by check_format
                    continue

                else:
                    self.err_string += '\n\t-{0}. Required format: {1}'.format(key, val[1])

            else:
                print("You've just discovered a bug (required formats). Please reach out to us!")
                sys.exit()


    def _determine_hetero_elements(self, solute: Solute):
        self.existing_elements = []
        for symbol in self._rdf_names.keys():
            if symbol.lower() in [e.lower() for e in solute.elements]:
                self.existing_elements.append(True)
            elif symbol == 'center':
                self.existing_elements.append(True)
            else:
                self.existing_elements.append(False)

    def _determine_colors(self, solute: Solute, solvent: Solvent) -> List[str]:
        """ set cutoffs """
        squared_cutoff1 = self.cutoff1 ** 2
        squared_cutoff2 = self.cutoff2 ** 2
        within_cutoff = []
        between_cutoffs = []
        outside_cutoff2 = []
        """ cycle over solvents and assign to list for each polar solute atom """

        for count, sol in enumerate(solvent.coords):

            for atom in solute.coords:
                if self._distance_squared(atom, sol) < squared_cutoff1:
                    within_cutoff.append(count)
                    break  #close enough solute atom was found for solvent within cutoff -> break loop over solute atoms
                elif self._distance_squared(atom, sol) > squared_cutoff2:
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
        barcolors = [self.outside] * len(solvent.coords)
        for within in within_cutoff:
                barcolors[within] = self.within
        for outside in between_cutoffs:
                barcolors[outside] = self.between

        return barcolors

    def _round_to_percentage(self, number: float, percentage: float) -> float:
        return round(number * percentage) / percentage

    def _axis_font(self):
        return {'size': str(self.fontsize)}

    def _create_legend(self, ax: plt.Axes, save_selected: bool):
        within = mpatches.Patch(color=self.within, label=u'd(solute-solvent) < %.1f Å' % self.cutoff1)
        between = mpatches.Patch(color=self.between,
                                 label=u'%.1f ≤ d(solute-solvent) ≤ %.1f Å' % (self.cutoff1, self.cutoff2))
        outside = mpatches.Patch(color=self.outside, label=u'd(solute-solvent) > %.1f Å' % self.cutoff2)
        if save_selected:
            # add green bars in legend
            selected = mpatches.Patch(color=self.selected, label='selected solvent molecules')
            ax.legend(handles=[within, between, outside, selected], prop=self._axis_font())
        else:
            ax.legend(handles=[within, between, outside], prop=self._axis_font())

    # creates interactive bar plot
    def _create_plot(self, barcolors: List[str], solvent: Solvent, save_selected: bool):
        if len(solvent.values) == 0:
            sys.exit(print("No solvent data found. Please check your settings. Maybe the solvent abbreviation is wrong?"))
        plt.ioff()
        indices = np.arange(1, len(solvent.values) + 1) # x-values
        fig = plt.figure(figsize=(self.width, self.height))
        ax = fig.add_subplot(111)
        rects = ax.bar(indices, solvent.values, color=barcolors, picker=True)  # create bars

        # make bars interactive
        for rect, color in zip(rects, barcolors):
            dr = ClickableBar(rect, color, self)
            self._drs.append(dr)
        self.original_drs = self._drs

        # font specifications
        matplotlib.rcParams.update({'font.size': self.fontsize})
        if self.y_numbers == 0.5:
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        elif self.y_numbers == 1.0:
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        else:
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        # Additional plotting option for vertical lines to show cutoffs
        #for x in self.marks:
        #    plt.axvline(x=x, color=self.mark, linestyle='--', lw=4)

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
        if abs(ymax-ymin) < 1:
            ymax = ymin + 1
            ystep = 0.1
        plt.yticks(np.arange(ymin, ymax + ystep, step=ystep))
        ax.set_ylim([ymin, ymax])
        ax.set_xlim([xmin, xmax])
        ax.tick_params(axis='both', which='major', labelsize=self.fontsize)

        # creates legend according with or without green bar
        self._create_legend(ax, save_selected)

        # saves figure without title and button
        if save_selected:
            plt.savefig(self.selected_plotname, dpi=self.dpi, transparent=self.transparent,
                        #orientation=self.orientation
                        )
        else:
            plt.savefig(self.plotname, dpi=self.dpi, transparent=self.transparent,
            #orientation=self.orientation
            )

            # legend with green bars and title specifications for "GUI"
            self._create_legend(ax, save_selected=True)

            # specify button
            b = plt.axes((0.8, 0.9, 0.1, 0.075))  # set position
            button = Button(b, "Display RDF", color="0.85", hovercolor="0.95")  # set text and color
            callback_b = ButtonActions()  # init class with necessary functions and determine click actions
            button.on_clicked(lambda x: ButtonActions.plot_rdf(callback_b, self))

            # specify textbox
            t = plt.axes((0.125, 0.9, 0.1, 0.075)) #left, bottom, width, height
            global txt_box
            txt_box = TextBox(t, "IDs:")
            txt_box.on_submit(self._input_selection)

            # specify clear all
            c = plt.axes((0.5, 0.9, 0.1, 0.075))
            button_c = Button(c, "Clear all", color="0.85", hovercolor="0.95")  # set text and color
            button_c.on_clicked(lambda _: self._deselect_all())

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
            if not self._selected_solvents:
                print("No solvents selected")
                if Input("Do you want to change/set the cutoff value? [y/n] ").yn():
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
                        self.cutoff1 = Input("Choose first distance cutoff: ", type=float).input
                        old_cutoff1 = self.cutoff1
                        selected_cutoff = True

                    # ensures that user enters correct data type
                    selected_cutoff = False
                    while not selected_cutoff:
                        self.cutoff2 = Input("Choose second distance cutoff: ", type=float).input
                        old_cutoff2 = self.cutoff2
                        selected_cutoff = True

                    inp = Input("Choose maximum number of solvent molecules displayed: ")
                    # user can enter 'all' or number -> first check all before ensuring datatype int
                    selected_solvent_number = False
                    while not selected_solvent_number:
                        if inp.input.lower() == 'all':
                            self.displayed_solvents = 'all'
                            old_displayed_solvents = 'all'
                            selected_solvent_number = True
                        else:
                            try:
                                self.displayed_solvents = int(inp.input)
                                selected_solvent_number = True
                            except ValueError:
                                inp = Input("Please enter a number or 'all': ")

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

                    # redetermine _colors because of new distance bounds
                    barcolors = self._determine_colors(solute, solvent)
                    # display GUI again
                    self._create_plot(barcolors, solvent, False)
                else:
                    print("No solvents selected")
                    sys.exit()
            else:  # solvent were selected, reselection of parameters are not necessary
                selected = True

    def _input_selection(self, text: str = None):  # ids separated with ",". allows ranges with "-"
        # this in connection with set_val('') prevents the double submission of the textbox when clicking on a bar after
        # input selection. fyi: on_submit gets triggered with enter and with leaving the textbox
        if len(text) == 0:  # no input given
            return

        split_1 = text.split(",")
        sel_solv = [] # contains selected solvents

        for i in range(len(split_1)):
            split_2 = split_1[i].split("-")
            if len(split_2) == 1:  # split[i] is not a range
                try:  # try str-int conversion
                    sel_solv.append(int(split_2[0]))
                except ValueError:
                    print("Invalid input. Please re-enter your selection!")
                    return
            else:  # range was given. also takes care of weird ranges like 1-3-5 and 5-1, which are treated as 1-5
                range_list = []
                for j in range(len(split_2)):
                    try:  # try str-int conversion
                        range_list.append(int(split_2[j]))
                    except ValueError:
                        print("A range containing an invalid input was given. Please re-enter your selection!")
                        return
                a = sorted(range_list)

                for k in range(a[0], a[-1] + 1):
                    sel_solv.append(k)

        sel_solv = [x for x in sel_solv if (1 <= x <= self.displayed_solvents)]
        out = ','.join((str(x) for x in sorted(set(sel_solv))))
        print("Solvent selected: {0}".format(out))

        self._selected_solvents.extend(sorted(set([x - 1 for x in sel_solv])))

        if len(self._drs) != 0: #0 is only the case when there is no barplot due to predefined solvent selection
            for dr in self._drs:
                try:
                    if int(dr.rect.xy[0]) in self._selected_solvents:
                        dr.color = self.selected
                        canvas = dr.rect.figure.canvas
                        axes = dr.rect.axes
                        canvas.draw()
                        dr.background = canvas.copy_from_bbox(dr.rect.axes.bbox)
                        dr.rect.set_color(self.selected)
                        axes.draw_artist(dr.rect)
                        canvas.blit(axes.bbox)
                except (ValueError, TypeError):
                    print("Something was wrong with re-coloring the bars upon input selection!")
                    return
            txt_box.set_val('')

    def _deselect_all(self):
        print("Resetting plot...")
        for dr in self._drs:
            dr.color = dr.original_color
            canvas = dr.rect.figure.canvas
            axes = dr.rect.axes
            dr.rect.set_color(dr.original_color)
            axes.draw_artist(dr.rect)
            canvas.blit(axes.bbox)
        self._selected_solvents = []
        print("Plot resetted!")

    def _save_selection(self, abb, solute: Solute, solvent: Solvent, reference: Reference) -> str:
        self._selected_solvents = list(set(self._selected_solvents))
        print('\nSolvents chosen: ' + str([solv+1 for solv in sorted(self._selected_solvents)]) +
              ' (Total: {0})'.format(len(self._selected_solvents)))
        filename = 'solvated_structure-' + str(len(self._selected_solvents)) + '.pdb'

        # does not open GUI, but saves plot of selected bars
        barcolors = self._determine_colors(solute, solvent)
        for select in self._selected_solvents:
            barcolors[select] = self.selected
        self._create_plot(barcolors, solvent, True)

        # writes latest solvated structure to file to open with pymol
        with open("latest-solvation.log", "a") as latest:
            latest.write(filename + '\n')

        # writes file with solute and selected solvents
        write_pdb(filename, solute, abb, solute=True) #writes only the solute into the pdb-file

        # info: this object finally contains all element labels, coords and values of all selected solvents.
        selected_solvent = Solvent()

        # TODO: Check if _selected_solvents order coincides with order of solvent.coord entries,
        #  i.e. are solvent.coord entries sorted wrt their energy
        for select in self._selected_solvents:
            voxel = int(solvent.data[select][0])
            quats = solvent.quats[voxel]
            grid_point = (float(solvent.data[select][1]),
                          float(solvent.data[select][2]),
                          float(solvent.data[select][3]))

            #This finally determines the solvent to be placed.
            elements, coords = reference._find_avg_solvent(voxel, quats, grid_point)
            values = float(solvent.data[select][-1])
            selected_solvent.elements.extend(elements)
            selected_solvent.coords.extend(coords)
            selected_solvent.values.extend([values]*len(elements))

        write_pdb(filename, selected_solvent, abb, solute=False)
        print('Your microsolvated structure was written to: ' + filename)
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
        if self.color == self.plot.selected:
            print('Solvent deselected:', int(np.round(event.xdata)))
            self.plot._selected_solvents.remove(index)
            self.rect.set_color(self.original_color)
            self.color = self.original_color

        else:
            print('Solvent selected:', int(np.round(event.xdata)))
            self.plot._selected_solvents.append(index)
            self.rect.set_color(self.plot.selected)
            self.color = self.plot.selected
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
