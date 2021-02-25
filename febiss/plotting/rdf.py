#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

from scipy.signal import argrelextrema
from typing import List, Tuple
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os


class ButtonActions(object):
    def __init__(self):
        self.axs = []
        self.integrals = []
        self.scs = []
        self.annots = []

    def plot_rdf(self, display):
        matplotlib.rcParams.update({'font.size': 10})
        self.fig = plt.figure(figsize=(display.width, display.height))
        self.display = display

        rows, cols = self._get_rows_and_cols(display)

        count = 0  # only count existing -> not enumerate
        for existing, (symbol, name) in zip(display.existing_elements, display.rdf_names.items()):
            if existing:
                count += 1
                if os.path.exists('rdf-' + str(name) + '.dat'):
                    arr = np.loadtxt("rdf-" + str(name) + ".dat")
                else:
                    print("ERROR: RDF analysis for " + str(name) + " was not performed in this directory!")
                    ax = self.fig.add_subplot(rows, cols, count)
                    txt = ax.text(0.1, 0.5, '', transform=ax.transAxes)
                    txt.set_text("ERROR: RDF analysis for " + str(name) + "\nwas not performed in this directory!")
                    plt.plot()
                    continue

                x = arr[:, 0]
                y = arr[:, 1]
                ax = self.fig.add_subplot(rows, cols, count)
                self.axs.append(ax)

                # determine integrals
                sc_x, sc_y, integrals = self._find_local_minima_and_maxima(x, y, name)
                sc = plt.scatter(sc_x, sc_y, s=10, c=display.colors['mark'])
                self.integrals.append(integrals)
                self.scs.append(sc)
                annot = ax.annotate("", xy=(0, 0), xytext=(20, 20), textcoords="offset points",
                                    bbox=dict(boxstyle="round", fc="w"), arrowprops=dict(arrowstyle="->"))
                annot.set_visible(False)
                self.annots.append(annot)

                # title and label specifications
                plt.xlabel("Distance of " + str(name) + ' to oxygen atoms in water / \u00c5')
                plt.ylabel('RDF')
                plt.xticks(np.arange(0, np.max(x) + 0.5, step=0.5))
                ax.set_xlim([0, np.max(x)])
                ax.axhline(y=1, ls='--', color=display.colors['mark'])
                plt.plot(x, y, linestyle="-", color='#80b1d3')

        plt.ion()  # avoids 'The event loop is already running' error message
        self.fig.canvas.mpl_connect('motion_notify_event', lambda event: self._hover(event))
        plt.show()

    def _get_rows_and_cols(self, display) -> Tuple[int, int]:
        true_count = sum(display.existing_elements)
        if true_count % 2 == 0:
            rows = int(round(true_count / 2))
            cols = int(round(true_count / 2))
            if true_count == 2:
                rows = 2
        else:
            rows = int(round(true_count / 2 + 0.5))
            cols = int(round(true_count / 2 + 0.5))
            if true_count == 5:
                cols = 2
        return rows, cols

    def _find_local_minima_and_maxima(self, distances: np.array, values: np.array, name: str) -> Tuple[List[float],
                                                                                                       List[float],
                                                                                                       List[float]]:
        n_local = 5
        maxima = argrelextrema(values, np.greater, order=n_local)[0]
        minima = argrelextrema(values, np.less, order=n_local)[0]
        extrema = np.asarray(list(maxima) + list(minima))
        ext_distances = [distances[x] for x in extrema]
        ext_values = [values[x] for x in extrema]
        integrals = self._get_integrals(extrema, name)
        return ext_distances, ext_values, integrals

    def _get_integrals(self, indices: np.array, name: str) -> List[float]:
        arr = np.loadtxt("int-rdf-" + str(name) + ".dat")
        return [arr[:, 1][i] for i in indices]

    def _update_annot(self, ind, subplot_number: int):
        index = ind['ind'][0]
        integral = self.integrals[subplot_number][index]
        text = "{0:.2f} waters".format(integral)
        annot = self.annots[subplot_number]
        annot.xy = self.scs[subplot_number].get_offsets()[index]
        annot.set_text(text)
        annot.get_bbox_patch().set_facecolor(self.display.colors['mark'])
        annot.get_bbox_patch().set_alpha(0.4)

    def _hover(self, event):
        for i, a in enumerate(self.axs):
            if event.inaxes == a:
                contains, ind = self.scs[i].contains(event)
                annot = self.annots[i]
                visible = annot.get_visible()
                if contains:
                    self._update_annot(ind, i)
                    annot.set_visible(True)
                    self.fig.canvas.draw_idle()
                else:
                    if visible:
                        annot.set_visible(False)
                        self.fig.canvas.draw_idle()
