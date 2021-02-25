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
from typing import List
from warnings import warn
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import sys

from ..utilities.structures import Solute, Water
from ..utilities.distance_functions import distance_squared
from ..utilities.io_handling import write_pdb
from .rdf import ButtonActions


class Plot:
    def __init__(self, **kwargs):
        self._set_defaults()
        self.allowed_keys = {'cutoff1', 'cutoff2', 'displayed_waters', 'colors', 'width', 'height', 'dpi', 'xlabel',
                             'ylabel', 'selected_plotname', 'plotname', 'orientation', 'fontsize', 'number_xtics',
                             'y_numbers', 'marks', 'transparent', 'display_once', 'testing', 'febiss_file', 'rdf_names'}
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in self.allowed_keys)
        for key in kwargs.keys():
            if key not in self.allowed_keys:
                warn('WARNING: Did not recognize key: ' + str(key))

    def gui(self, solute: Solute, water: Water) -> str:
        if not self.testing:
            self._determine_hetero_elements(solute)
            barcolors = self._determine_colors(solute, water)
            self._create_plot(barcolors, water, False)
        else:
            self.selected_waters = [0, 1, 2, 3, 4]
        if self.display_once:
            sys.exit()
        # avoid bug of selecting out of range water
        if -1 in self.selected_waters:
            self.selected_waters.remove(-1)
        self._interactive_reselection(solute, water)
        filename = self._save_selection(solute, water)
        return filename

    def _set_defaults(self):
        # default values for bar chart
        self.selected_waters = []
        self.cutoff1 = 3.0
        self.cutoff2 = 6.0
        self.displayed_waters = 50
        self.colors = {'within': '#fdb462', 'between': '#80b1d3', 'outside': '#de2d26', 'selected': '#b3de69',
                       'mark': '#bc80bd'}
        self.width = 16
        self.height = 8
        self.dpi = 200
        self.xlabel = 'Water molecules sorted by their Free Energy values'
        self.ylabel = '$-$ Free Energy / kcal mol$^{-1}$'
        self.selected_plotname = 'febiss-plot-selected.png'
        self.plotname = 'febiss-plot.png'
        self.orientation = 'landscape'
        self.fontsize = 18
        self.number_xtics = 10
        self.y_numbers = 0.25
        self.marks = []
        self.febiss_file = 'febiss-waters.pdb'
        self.rdf_names = {'center2': 'center of solute', 'C': 'carbon', 'O': 'oxygen', 'N': 'nitrogen', 'P': 'phosphor'}
        self.transparent = True
        self.display_once = False
        self.testing = False

    def _determine_hetero_elements(self, solute: Solute):
        self.existing_elements = []
        for symbol in self.rdf_names.keys():
            if symbol in solute.elements:
                self.existing_elements.append(True)
            elif symbol == 'center2':
                self.existing_elements.append(True)
            else:
                self.existing_elements.append(False)

    def _determine_colors(self, solute: Solute, water: Water) -> List[str]:
        """ set cutoffs """
        squared_cutoff1 = self.cutoff1 ** 2
        squared_cutoff2 = self.cutoff2 ** 2
        # if distance between two water atoms below this,
        # they belong to same water molecule, but they have to be in pdb file within 2 rows)
        same_water_cutoff = 1.6
        squared_same_water_cutoff = same_water_cutoff ** 2
        if len(solute.polars) == 0:
            solute.determine_polar_hydrogen_and_non_hydrogen()
        within_cutoff = []
        between_cutoffs = []
        outside_cutoff2 = []
        """ cycle over waters and assign to list for each polar solute atom """
        skip_next = False
        skip_next_next = False
        for count, wat in enumerate(water.atoms):
            if skip_next:
                if skip_next_next:
                    skip_next_next = False
                    continue
                else:
                    skip_next = False
                    continue
            for polar in solute.polars:
                if distance_squared(polar, wat) < squared_cutoff1:
                    within_cutoff.append(count)
                    # relies on atoms of same water molecule to be right after each other in pdb!
                    for i in range(-2, 3):
                        # bound check and then distance
                        if i != 0 and len(water.elements) > count + i > 0 \
                                and distance_squared(wat, water.atoms[count + i]) < squared_same_water_cutoff:
                            within_cutoff.append(count + i)
                            if i == 1:
                                skip_next = True
                            elif i == 2:
                                skip_next_next = True
                    break  # close enough solute atom was found for water within cutoff -> break loop over solute atoms
                elif distance_squared(polar, wat) > squared_cutoff2:
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
        # within overrules outside (possible because distance of all solute atoms and waters are calculated)
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
        barcolors = [self.colors['outside']] * len(water.elements)
        for within in within_cutoff:
            if water.elements[within] == "O":
                index = int(round(within / 3))
                barcolors[index] = self.colors['within']
        for outside in between_cutoffs:
            if water.elements[outside] == "O":
                index = int(round(outside / 3))
                barcolors[index] = self.colors['between']

        return barcolors

    def _round_to_percentage(self, number: float, percentage: float) -> float:
        return round(number * percentage) / percentage

    def _axis_font(self):
        return {'size': str(self.fontsize)}

    def _create_legend(self, ax: plt.Axes, save_selected: bool):
        within = mpatches.Patch(color=self.colors['within'], label=u'd(solute-water) < %.1f Å' % self.cutoff1)
        between = mpatches.Patch(color=self.colors['between'],
                                 label=u'%.1f ≤ d(solute-water) ≤ %.1f Å' % (self.cutoff1, self.cutoff2))
        outside = mpatches.Patch(color=self.colors['outside'], label=u'd(solute-water) > %.1f Å' % self.cutoff2)
        if save_selected:
            # add green bars in legend
            selected = mpatches.Patch(color=self.colors['selected'], label='selected water molecules')
            ax.legend(handles=[within, between, outside, selected], prop=self._axis_font())
        else:
            ax.legend(handles=[within, between, outside], prop=self._axis_font())

    # creates interactive bar plot
    def _create_plot(self, barcolors: List[str], water: Water, save_selected: bool):
        plt.ioff()
        indices = np.arange(1, len(water.values) + 1)  # x-values
        fig = plt.figure(figsize=(self.width, self.height))
        ax = fig.add_subplot(111)
        rects = ax.bar(indices, water.values, color=barcolors, picker=True)  # create bars

        # make bars interactive
        drs = []  # necessary for interaction
        for rect, color in zip(rects, barcolors):
            dr = ClickableBar(rect, color, self)
            drs.append(dr)

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
        ymax = np.max(water.values)
        fill_up = self._round_to_percentage(ymax, 1 / self.y_numbers)
        if fill_up < ymax:
            fill_up = self.y_numbers - (ymax - fill_up)
            ymax += fill_up
        else:
            ymax += (fill_up - ymax)

        if self.displayed_waters == 'all':  # if user wants all waters displayed
            xmax = len(water.values) + 0.5  # 0.5 ensures slight offset, so that last bar is seen completely
            ymin = np.min(water.values) - 0.5
        else:
            xmax = self.displayed_waters + 0.5
            ymin = water.values[self.displayed_waters]
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
        plt.yticks(np.arange(ymin, ymax + ystep, step=ystep))
        ax.set_ylim([ymin, ymax])
        ax.set_xlim([xmin, xmax])
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(self.fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(self.fontsize)

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
            plt.title('Choose water molecules for solvation by clicking bars and close window\n')
            # specify button
            b = plt.axes([0.85, 0.9, 0.1, 0.075])  # set position
            button = Button(b, "Display RDF", color="0.85", hovercolor="0.95")  # set text and color
            callback = ButtonActions()  # init class with necessary functions and determine click actions
            button.on_clicked(lambda x: ButtonActions.plot_rdf(callback, self))
            plt.show()

    def _interactive_reselection(self, solute: Solute, water: Water):
        # save old values to tell user later
        old_cutoff1 = self.cutoff1
        old_cutoff2 = self.cutoff2
        old_displayed_waters = self.displayed_waters

        # while loop ensures selection or denial of selection by user
        # allows for changing distance bounds and displayed water number
        selected = False
        while not selected:
            if not self.selected_waters:
                print("No waters selected")
                inp = input("Do you want to change/set the cutoff value? [y/n] ")
                if inp in ['y', 'Y', 'yes', 'Yes', ' y']:
                    if old_displayed_waters == "all":
                        print("Please enter the new values (previous distance cutoffs were {0:.2f} and {1:.2f} and "
                              "all water molecules were shown)".format(old_cutoff1, old_cutoff2))
                    else:
                        print("Please enter the new values (previous distance cutoffs were {0:.2f} and {1:.2f} and "
                              "{2:d} water molecules were shown)".format(old_cutoff1, old_cutoff2,
                                                                         old_displayed_waters))

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

                    inp = input("Choose maximum number of water molecules displayed: ")
                    # user can enter 'all' or number -> first check all before ensuring datatype int
                    selected_water_number = False
                    while not selected_water_number:
                        if inp in ["all", "All", "ALL"]:
                            self.displayed_waters = 'all'
                            old_displayed_waters = 'all'
                            selected_water_number = True
                        else:
                            try:
                                self.displayed_waters = int(inp)
                                selected_water_number = True
                            except ValueError:
                                inp = input("Please enter a number or 'all': ")

                    # user wants crazy number or is unaware of 'all' function --> print all
                    if self.displayed_waters > len(water.values):
                        print('Chosen water number is greater than number of all placed water molecules.')
                        print('x-axis will stop at last water.')
                        self.displayed_waters = 'all'
                        old_displayed_waters = 'all'
                    # weird input and only one water is plotted
                    elif self.displayed_waters < 1:
                        print('Chosen water number must be at least 1. Plotting 1 water molecule')
                        self.displayed_waters = 1
                        old_displayed_waters = 1
                    else:
                        old_displayed_waters = self.displayed_waters

                    # redetermine colors because of new distance bounds
                    barcolors = self._determine_colors(solute, water)
                    # display GUI again
                    self._create_plot(barcolors, water, False)
                elif inp in ['n', 'N', 'no', 'No', ' n']:
                    print("No waters selected")
                    sys.exit()
                else:
                    print("Wrong input, just 'y' or 'n'.")
            else:  # water were selected, reselection of parameters are not necessary
                selected = True

    def _save_selection(self, solute: Solute, water: Water) -> str:
        print('Number of waters chosen: ' + str(len(self.selected_waters)))
        filename = 'solvated_structure-' + str(len(self.selected_waters)) + '.pdb'
        print('Your microsolvated structure is written to: ' + filename)
        # does not open GUI, but saves plot of selected bars
        barcolors = self._determine_colors(solute, water)
        for select in self.selected_waters:
            barcolors[select] = self.colors['selected']
        self._create_plot(barcolors, water, True)

        # writes latest solvated structure to file to open with pymol
        with open("latest-solvation.log", "a") as latest:
            latest.write(filename + '\n')
        # writes file with solute and selected waters
        write_pdb(filename, solute, solute=True)
        selected_water = Water()
        for select in self.selected_waters:
            selected_water.atoms.append(water.atoms[select * 3])
            selected_water.atoms.append(water.atoms[select * 3 + 1])
            selected_water.atoms.append(water.atoms[select * 3 + 2])
            selected_water.values.append(water.all_values[select * 3])
            selected_water.values.append(water.all_values[select * 3 + 1])
            selected_water.values.append(water.all_values[select * 3 + 2])
            selected_water.elements.append('O')
            selected_water.elements.append('H')
            selected_water.elements.append('H')
        write_pdb(filename, selected_water)
        return filename


class ClickableBar:
    lock = None  # only one can be clicked at a time

    def __init__(self, rect, color, plot):
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
        # self.color is used to determine whether water gets selected or deselected
        # self.original_color is used to get previous color if water is deselected
        if self.color == self.plot.colors['selected']:
            print('water deselected:', int(np.round(event.xdata)))
            self.plot.selected_waters.remove(index)
            self.rect.set_color(self.original_color)
            self.color = self.original_color
        else:
            print('water selected:', int(np.round(event.xdata)))
            self.plot.selected_waters.append(index)
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
