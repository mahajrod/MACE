#!/usr/bin/env python
__author__ = "tomarovsky"

import os
import re
import string
from copy import deepcopy
from functools import partial
from pathlib import Path

import distinctipy
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Rectangle
from RouToolPa.Collections.General import SynDict
from RouToolPa.Parsers.BED import CollectionBED
from RouToolPa.Parsers.BLAST import CollectionBLAST
from RouToolPa.Parsers.GFF import CollectionGFF
from RouToolPa.Parsers.STR import CollectionSTR
from RouToolPa.Parsers.VCF import CollectionVCF

from MACE.Routines import StatsVCF, Visualization

# TODO
# upgrade PSMC plot
# add legend props


class Plotter:
    def __init__(self):
        pass

    def set_figure_fontsize(self, font_scale):
        """
        plt.rcParams.update({'font.size': font_scale})
        """
        plt.rcParams.update({"font.size": font_scale})

    def set_paperticks_style(self, font_scale):
        """
        Configure a "ticks" style and "paper" context.
        """
        custom_params = {
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.grid": True,
            "axes.axisbelow": True,
            "grid.color": "#dfdfdf",
            "grid.linestyle": "--",
        }
        sns.set_theme(style="ticks", rc=custom_params)
        sns.set_context("paper", font_scale=font_scale)

    def annotate_subplot(self, ax, annotation, offset=(-0.1, 1.1), fontsize=12, fontweight="bold", color="black"):
        """
        Add an annotation to a specific subplot (Axes) in a figure.

        Parameters:
        - ax : matplotlib.axes.Axes
            The specific subplot (Axes) to label.
        - annotation : str
            The annotation to add (e.g., "A", "B", etc.).
        - offset : tuple of float, optional
            The x and y offsets for the label relative to the subplot in Axes coordinates.
        - fontsize : int, optional
            Font size for the label.
        - fontweight : str or int, optional
            Font weight for the label (e.g., "bold", "normal", or a numeric value).
        - color : str, optional
            The color used for the sign. Defaults to "black".
        """
        ax.text(
            offset[0],
            offset[1],
            annotation,
            transform=ax.transAxes,
            fontsize=fontsize,
            fontweight=fontweight,
            color=color,
            va="center",
            ha="center",
        )

    def annotate_all_subplots(self, axes, offset=(-0.1, 1.1), fontsize=12, fontweight="bold", color="black"):
        """
        Add an annotation to a specific subplot (Axes) in a figure.

        Parameters:
        - ax : matplotlib.axes.Axes
            All mtplotlib axes.
        - annotation : str
            The annotation to add (e.g., "A", "B", etc.).
        - offset : tuple of float, optional
            The x and y offsets for the label relative to the subplot in Axes coordinates.
        - fontsize : int, optional
            Font size for the label.
        - fontweight : str or int, optional
            Font weight for the label (e.g., "bold", "normal", or a numeric value).
        - color : str, optional
            The color used for the sign. Defaults to "black".
        """
        letters = string.ascii_uppercase  # 'A', 'B', 'C', ..., 'Z'
        letter_index = 0

        for row in axes:
            for axi in row:
                self.annotate_subplot(
                    axi,
                    letters[letter_index],
                    offset=(offset[0], offset[1]),
                    fontsize=fontsize,
                    fontweight=fontweight,
                    color=color,
                )
                letter_index += 1

    def scaled_histogram_with_extended_bins(self, df, bins, scale=0.45):
        hist, bin_edges = np.histogram(df, bins=bins, density=True)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_centers = np.insert(bin_centers, 0, 0)
        hist = np.insert(hist, 0, hist[0])
        bin_centers = np.insert(bin_centers, 0, 0)
        hist = np.insert(hist, 0, 0)
        bin_centers = np.append(bin_centers, bin_centers[-1])
        hist = np.append(hist, 0)
        scaling_factor = scale / max(hist, default=1)
        hist *= scaling_factor
        return hist, bin_centers

    def add_boxplot_classic(self, ax, df, position, width=0.15):
        ax.boxplot(
            df,
            vert=True,
            positions=[position],
            widths=width,
            showfliers=True,
            showmeans=True,
            patch_artist=False,
            meanprops=dict(marker="o", markerfacecolor="white", markeredgecolor="white", markersize=0.5),
            flierprops=dict(marker="o", markerfacecolor="#575757", markeredgecolor="#575757", markersize=0.5),
            medianprops=dict(color="#575757", solid_capstyle="butt"),
            boxprops=dict(color="#575757", linewidth=0.5),
            whiskerprops=dict(color="#575757", linewidth=0.5),
            capprops=dict(color="#575757", linewidth=0.5),
        )

    def add_boxplot_empty(self, ax, df, position, width=0.15):
        ax.boxplot(
            df,
            vert=True,
            positions=[position],
            widths=width,
            whis=1,
            showcaps=True,
            showmeans=True,
            showfliers=True,
            patch_artist=True,
            boxprops=dict(alpha=0),
            capprops=dict(color="black", linewidth=0.8),
            whiskerprops=dict(color="black", linewidth=0.8),
            meanprops=dict(marker="o", markerfacecolor="w", markeredgecolor="w", markersize=2),
            flierprops=dict(marker="o", markerfacecolor="black", markeredgecolor="none", markersize=.8),
            medianprops=dict(color="w", linewidth=0.8)
        )
        q25, q75 = df.quantile([0.25, 0.75])
        ax.hlines(y=q25, xmin=position-0.1, xmax=position+0.1, colors="black", linewidth=0.8)
        ax.hlines(y=q75, xmin=position-0.1, xmax=position+0.1, colors="black", linewidth=0.8)

    def add_double_boxplot(self, ax, df_1, df_2, position, width=0.04):
        for d, pos in zip([df_1, df_2], [position - 0.02, position + 0.02]):
            ax.boxplot(
                d,
                vert=True,
                positions=[pos],
                widths=width,
                whis=0,
                showfliers=False,
                showmeans=False,
                patch_artist=True,
                boxprops=dict(facecolor="#575757", color="#575757", linewidth=0),
                medianprops=dict(color="white", solid_capstyle="butt", linewidth=2),
                whiskerprops=dict(color="#575757", linewidth=0),
                capprops=dict(color="#575757", linewidth=0),
            )

    def process_variant_counts(self, file_paths, removed_chrX=None, reference=None, window_size=1, multiplicator=1):
        data_list = []
        for count in file_paths:
            df = pd.read_csv(count, sep="\t")
            file_path = Path(count)
            file_name = file_path.stem
            id = file_name.split(".")[0]

            if removed_chrX:
                # If removed_chrX is a list, remove the specified chromosomes; otherwise, remove a single string
                if isinstance(removed_chrX, list):
                    for chrX in removed_chrX:
                        if chrX is not None:
                            df = df[~df["CHROM"].str.contains(chrX)]
                else:
                    df = df[~df["CHROM"].str.contains(removed_chrX)]

            df.set_index("CHROM", inplace=True)
            density = df.iloc[:, -1] / window_size * multiplicator

            if reference:
                data_list.extend([{"density": d, "id": id, "Reference": reference} for d in density])
            else:
                data_list.extend([{"density": d, "id": id} for d in density])

        return data_list

    def draw_stripped_histograms(
        self,
        ax,
        data,
        ymin,
        ymax,
        yticklist,
        window_size=1000000,
        multiplicator=1000,
        removed_chrX="",
        title="",
        ylabel="Heterozygous SNPs/kbp",
        bin_width=0.1,
        boxplot_type="empty",
        boxplot_width=0.15,
        palette="turbo",
        rotation=45,
        statuses=dict(),
        sample_ids=[],  # same sorting as in data
        horizontal_lines=[],
        figure_grid=True,
        font_size=None,
        sort_by=None,  # 'mean' or 'median'
    ):
        """
        Draws stripped histograms and boxplots to visualize heterozygous SNP density for multiple samples.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the histograms.

        data : list of str
            A list of file paths to tab-separated files containing SNP density data.

        ymin : float
            The minimum value for the y-axis.

        ymax : float
            The maximum value for the y-axis.

        yticklist : list of float
            List of y-axis tick positions.

        window_size : int, optional
            The size of the window used for calculating density. Defaults to 1,000,000.

        multiplicator : float, optional
            A multiplier applied to the density values for scaling. Defaults to 1000.

        removed_chrX : str, optional
            Chromosome name to be removed from the analysis (e.g., "chrX"). If not provided, no chromosomes are removed.
            Defaults to an empty string.

        title : str, optional
            The title of the plot. Defaults to an empty string.

        ylabel : str, optional
            The label for the y-axis. Defaults to "Heterozygous SNPs/kbp".

        bin_width : float, optional
            The width of the bins for the histograms. Defaults to 0.1.

        boxplot_type : str, optional
            Type of boxplot to draw. Options are 'classic' or 'empty'. Defaults to 'empty'.

        boxplot_width : float, optional
            The width of the boxplot elements. Defaults to 0.15.

        palette : str, optional
            The color palette used for the histograms. Defaults to "turbo".

        rotation : int or float, optional
            The rotation angle for x-axis labels. Defaults to 45.

        statuses : dict, optional
            A dictionary mapping sample IDs to status labels. Status labels are displayed above the histograms.
            Defaults to an empty dictionary.

        sample_ids : list of str, optional
            A list of sample IDs for labeling the x-axis. If not provided, the unique IDs from the data are used.
            Defaults to an empty list.

        horizontal_lines : list of float, optional
            List of y-coordinates at which to draw horizontal lines across the plot. Defaults to an empty list.

        figure_grid : bool, optional
            Whether to display a grid on the y-axis. Defaults to True.

        font_size : int or float, optional
            Font size for the plot text. If not provided, the default font size is used.

        sort_by : str, optional
            Specifies how to sort the sample IDs. Options are 'mean' or 'median'. If not provided, no sorting is applied.
            Defaults to None.

        Notes:
        ------
        - The function reads density data from each file, optionally filters out a specified chromosome, and calculates
          SNP density. It then creates stripped histograms and optional boxplots for each sample.
        - Status labels are displayed above the histograms if provided in the `statuses` dictionary.
        """

        # ax.spines[["right", "top"]].set_visible(False)
        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)

        ax.set_axisbelow(True)

        if font_size is not None:
            plt.rcParams.update({"font.size": font_size})

        data_list = self.process_variant_counts(data, removed_chrX=removed_chrX, window_size=window_size, multiplicator=multiplicator)
        data = pd.DataFrame(data_list)

        # Calculate mean or median values for sorting if needed
        if sort_by in ["mean", "median"]:
            sort_values = data.groupby("id")["density"].agg(sort_by).sort_values()
            unique_ids = sort_values.index
        else:
            unique_ids = data["id"].unique()

        try:
            colors = sns.color_palette(palette, len(unique_ids))
        except ValueError:
            colors = palette

        for i, unique_id in enumerate(unique_ids):
            df = data[data["id"] == unique_id]["density"]

            # bins
            bins = np.arange(df.min(), df.max() + bin_width, bin_width)

            # draw stripped histograms
            hist, bin_centers = self.scaled_histogram_with_extended_bins(df, bins)
            ax.fill_betweenx(bin_centers, i, i - hist, color=colors[i], edgecolor=colors[i])
            ax.fill_betweenx(bin_centers, i, i + hist, color=colors[i], edgecolor=colors[i])

            # boxplot
            if boxplot_type == "classic":
                self.add_boxplot_classic(ax, df, i, width=boxplot_width)
            elif boxplot_type == "empty":
                self.add_boxplot_empty(ax, df, i, width=boxplot_width)

            # add statuses
            if statuses and unique_id in statuses:
                last_bin_center = bin_centers[-1]
                ax.text(i, last_bin_center + 0.05, statuses[unique_id], ha="center", va="bottom", color=colors[i])

        ax.set_xticks(range(len(unique_ids)))
        if sample_ids:
            ax.set_xticklabels(sample_ids, ha="right", rotation=rotation)
        else:
            ax.set_xticklabels(unique_ids, ha="right", rotation=rotation)
        ax.set_yticks(yticklist)
        if horizontal_lines:
            for ycoord in horizontal_lines:
                ax.axhline(y=ycoord, color="red", linestyle="--", linewidth=1)
        ax.set_ylim(ymin=ymin, ymax=ymax)
        ax.set_xlim(xmin=-0.5, xmax=len(unique_ids) - 0.5)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if figure_grid:
            ax.grid(axis="y", linestyle="--", alpha=0.5)

    def draw_double_stripped_histograms(
        self,
        ax,
        left_data,
        right_data,
        ymin,
        ymax,
        yticklist,
        window_size=1000000,
        multiplicator=1000,
        removed_chrX=[],
        sort_by=None,
        title="",
        ylabel="Heterozygous SNPs/kbp",
        bin_width=0.1,
        boxplot_width=0.04,
        references=["Species 1", "Species 2"],
        colors=["#E04B4B", "#6094C3"],
        rotation=45,
        horizontal_lines=[],
        figure_grid=True,
        font_size=None,
        show_legend=True,
        legend_title="Reference",
        legend_loc="upper left",
        legend_ncol=2,
    ):
        """
        Draw double stripped histograms comparing the density of data between two species across various IDs.

        This function plots histograms and boxplots comparing heterozygous SNP density from two different data sets,
        represented as `left_data` and `right_data`. The histograms are stripped and filled to visually compare
        the density distributions for each ID. The function also includes the ability to sort IDs by mean or
        median density values and customize various plot elements.

        Parameters
        ----------
        ax : matplotlib.axes._axes.Axes
            The axes object where the histograms and boxplots will be drawn.
        left_data : list of str
            List of file paths to data sets representing the left side (first species).
        right_data : list of str
            List of file paths to data sets representing the right side (second species).
        ymin : float
            Minimum value for the y-axis.
        ymax : float
            Maximum value for the y-axis.
        yticklist : list of float
            List of y-tick values to use on the y-axis.
        window_size : int, optional, default=1000000
            The size of the window used to normalize density values.
        multiplicator : int, optional, default=1000
            A multiplier applied to the density values for scaling purposes.
        removed_chrX : str, optional, default=""
            String pattern used to remove certain chromosomes from the data if needed.
        title : str, optional, default=""
            Title of the plot.
        ylabel : str, optional, default="Heterozygous SNPs/kbp"
            Label for the y-axis.
        bin_width : float, optional, default=0.1
            Width of the histogram bins.
        boxplot_width : float, optional, default=0.04
            Width of the boxplots.
        references : list of str, optional, default=["Species 1", "Species 2"]
            List of labels for the two data sets being compared.
        colors : list of str, optional, default=["#E04B4B", "#6094C3"]
            List of colors to use for the histograms of the two data sets.
        rotation : int or float, optional, default=45
            Angle of rotation for the x-axis labels.
        horizontal_lines : list of float, optional, default=[]
            List of y-coordinates for horizontal lines to be drawn across the plot.
        figure_grid : bool, optional, default=True
            Whether to display a grid on the y-axis.
        font_size : int, optional, default=None
            Font size for the plot text. If None, the default font size is used.
        sort_by : str, optional, default="mean"
            Method to sort the IDs, either 'mean' or 'median'.
            - 'mean': Sort IDs by the mean density value of the left data.
            - 'median': Sort IDs by the median density value of the left data.
        show_legend : bool, optional
            Whether to display a legend on the plot. Default is True.
        legend_title : str, optional
            Title of the legend. Default is "Reference".
        legend_loc : str, optional
            Location of the legend on the plot. Default is "upper left".
        legend_ncol : int, optional
            Number of columns in the legend. Default is 2.

        Notes
        -----
        - The function reads data from the provided file paths, calculates density values, and removes any data
          matching `removed_chrX` if specified.
        - It uses the `scaled_histogram_with_extended_bins` method to generate histograms and the
          `add_double_boxplot` method to add boxplots to the plot.
        - The x-axis labels are the unique IDs, and the y-axis limits, grid, and additional horizontal lines
          can be customized.
        """

        # ax.spines[["right", "top"]].set_visible(False)
        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)

        ax.set_axisbelow(True)

        if font_size is not None:
            plt.rcParams.update({"font.size": font_size})

        left_data_list = self.process_variant_counts(
            left_data,
            removed_chrX=[removed_chrX[0]] if removed_chrX else None,
            reference=references[0],
            window_size=window_size,
            multiplicator=multiplicator,
        )
        right_data_list = self.process_variant_counts(
            right_data,
            removed_chrX=[removed_chrX[1]] if removed_chrX else None,
            reference=references[1],
            window_size=window_size,
            multiplicator=multiplicator,
        )

        data = pd.DataFrame(left_data_list + right_data_list)

        # sorting
        if sort_by == "mean":
            sorting_values = data[data["Reference"] == references[0]].groupby("id")["density"].mean()
            unique_ids = sorting_values.sort_values().index
        elif sort_by == "median":
            sorting_values = data[data["Reference"] == references[0]].groupby("id")["density"].median()
            unique_ids = sorting_values.sort_values().index
        else:
            unique_ids = data["id"].unique()

        for i, unique_id in enumerate(unique_ids):
            df_1 = data[(data["Reference"] == references[0]) & (data["id"] == unique_id)]["density"]
            df_2 = data[(data["Reference"] == references[1]) & (data["id"] == unique_id)]["density"]

            # bins
            bins = np.arange(min([df_1.min(), df_2.min()]), max([df_1.max(), df_2.max()]) + bin_width, bin_width)

            # stripped histograms
            hist_1, bin_centers_1 = self.scaled_histogram_with_extended_bins(df_1, bins)
            hist_2, bin_centers_2 = self.scaled_histogram_with_extended_bins(df_2, bins)
            ax.fill_betweenx(
                bin_centers_1,
                i,
                i - hist_1,
                color=colors[0],
                edgecolor=colors[0],
                linewidth=0.5,
                label=r"$\mathit{" + references[0] + "}$",
            )
            ax.fill_betweenx(
                bin_centers_2,
                i,
                i + hist_2,
                color=colors[1],
                edgecolor=colors[1],
                linewidth=0.5,
                label=r"$\mathit{" + references[1] + "}$",
            )
            ax.plot([i, i], [0, max([len(hist_1) * bin_width, len(hist_2) * bin_width])], linewidth=0.6, color="#575757")

            # boxplot
            self.add_double_boxplot(ax, df_1, df_2, i, width=boxplot_width)

            if show_legend:
                ax.legend(title=legend_title, loc=legend_loc, ncol=legend_ncol, handlelength=0.7, frameon=True) if i == 0 else None

        ax.set_xticks(range(len(unique_ids)))
        ax.set_xticklabels(unique_ids, ha="right", rotation=rotation)
        ax.set_yticks(yticklist)
        if horizontal_lines:
            for ycoord in horizontal_lines:
                ax.axhline(y=ycoord, color="red", linestyle="--", linewidth=1)
        ax.set_ylim(ymax=ymax, ymin=ymin)
        ax.set_xlim(-0.5, len(unique_ids) - 0.5)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if figure_grid:
            ax.grid(axis="y", linestyle="--", alpha=0.5)

    def draw_single_admixture_barplot(
        self,
        ax,
        df,
        colors=None,
        rotation=45,
        yticks=[0, 25, 50, 75, 100],
        xlabel="",
        font_size=None,
        show_legend=True,
        legend_loc="upper right",
        legend_ncol=4,
    ):
        """
        Draw a single stacked barplot to visualize ADMIXTURE proportions for multiple samples.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the plot.

        df : pandas.DataFrame
            A DataFrame containing ADMIXTURE data. The first column must be `Sample_ID`,
            followed by columns with ADMIXTURE proportions.

        admixture_columns : list of str
            Column names in `df` corresponding to ADMIXTURE proportions.

        colors : list of str, optional
            Colors for the ADMIXTURE bars. Defaults to ["#e02828", "#4286c3"].

        rotation : int, optional
            Rotation angle for the x-axis labels. Defaults to 45 degrees.

        yticks : list of int, optional
            Y-axis tick values. Defaults to [0, 25, 50, 75, 100].

        xlabel : str, optional
            Label for the x-axis. Defaults to an empty string.

        font_size : int, optional
            Font size for the plot text. If `None`, the default font size is used.

        show_legend : bool, optional
            Whether to display a legend. Default is True.

        legend_loc : tuple of float, optional
            Position of the legend box within the plot. Defaults to (0.1, 0.95).

        legend_ncol : int, optional
            Number of columns in the legend. Defaults to 4.

        Notes:
        ------
        - The first column of `df` is used as the index (sample identifiers) for the plot.
        - A single stacked bar plot is created to visualize ADMIXTURE proportions.
        - If `show_legend` is True, a legend is added to the plot.
        """
        df_plot = df.copy()
        index_column = df_plot.columns[0]
        df_plot.set_index(index_column, inplace=True)

        admixture_columns = df_plot.columns.tolist()
        k_value = len(admixture_columns)

        # If colors is None or len(colors) < k_value:
        if colors is None or len(colors) < k_value:
            cmap = plt.get_cmap("tab10")
            colors = [cmap(i) for i in range(k_value)]

        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)
        ax.set_axisbelow(True)

        if font_size is not None:
            plt.rcParams.update({"font.size": font_size})

        # 4. Отрисовка
        df_plot.plot(
            kind="bar",
            stacked=True,
            width=0.9,
            color=colors,
            ax=ax,
            legend=False
        )

        if show_legend:
            ncol = legend_ncol if legend_ncol else min(k_value, 4)
            ax.legend(
                admixture_columns,
                loc=legend_loc,
                ncol=ncol,
                # bbox_to_anchor=(1.0, 1.0)
            )

        ax.set_xlim(-0.5, len(df_plot) - 0.5)
        ax.set_xticklabels(df_plot.index, ha="right", rotation=rotation)
        ax.set_xlabel(xlabel)
        ax.set_yticks(yticks)
        ax.set_ylim(0, 1)


    def legend_half_patch(self, color1, color2):
        class SplitPatchHandler(HandlerPatch):
            def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
                p1 = Rectangle((xdescent, ydescent), width / 2, height, facecolor=color1, transform=trans)
                p2 = Rectangle((xdescent + width / 2, ydescent), width / 2, height, facecolor=color2, transform=trans)
                return [p1, p2]

        return Rectangle((0, 0), 1, 1), SplitPatchHandler()

    def draw_double_admixture_barplot(
        self,
        ax,
        df,
        global_columns,
        local_columns,
        references,
        local_colors=["#e02828", "#4286c3"],
        global_colors=["#e06c6c", "#7da2c3"],
        rotation=45,
        yticks=[0, 25, 50, 75, 100],
        xlabel="",
        font_size=None,
        show_legend=True,
        legend_loc=(0.1, 0.95),
        legend_ncol=4,
    ):
        """
        Draw a double barplot to visualize global and local ADMIXTURE proportions for multiple samples.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the plot.

        df : pandas.DataFrame
            A DataFrame containing ADMIXTURE data. The first column must be `Sample_ID`,
            followed by columns with global and local ADMIXTURE proportions.

        global_columns : list of str
            Column names in `df` corresponding to global ADMIXTURE proportions.

        local_columns : list of str
            Column names in `df` corresponding to local ADMIXTURE proportions.

        references : list of str
            Names of the reference populations used in ADMIXTURE analysis. This should be a list of two items.

        local_colors : list of str, optional
            Colors for the local ADMIXTURE bars. Defaults to ["#e02828", "#4286c3"].

        global_colors : list of str, optional
            Colors for the global ADMIXTURE bars. Defaults to ["#e06c6c", "#7da2c3"].

        rotation : int, optional
            Rotation angle for the x-axis labels. Defaults to 45 degrees.

        yticks : list of int, optional
            Y-axis tick values. Defaults to [0, 25, 50, 75, 100].

        xlabel : str, optional
            Label for the x-axis. Defaults to an empty string.

        font_size : int, optional
            Font size for the plot text. If `None`, the default font size is used.

        show_legend : bool, optional
            Whether to display a legend. Default is True.

        legend_loc : tuple of float, optional
            Position of the legend box within the plot. Defaults to (0.1, 0.95).

        legend_ncol : int, optional
            Number of columns in the legend. Defaults to 2.

        Notes:
        ------
        - The first column of `df` is used as the index (sample identifiers) for the plot.
        - Two stacked bar plots are created side-by-side: one for global ADMIXTURE proportions
          and the other for local ADMIXTURE proportions.
        - A custom legend is generated to distinguish between global and local ADMIXTURE for
          each reference population.
        """

        # ax.spines[["right", "top"]].set_visible(False)
        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)
        ax.set_axisbelow(True)

        if font_size is not None:
            plt.rcParams.update({"font.size": font_size})

        index_column = df.columns[0]  # first column as index ('Sample_ID')

        df.set_index(index_column, inplace=True)

        stacked_data, stacked_data2 = df[global_columns], df[local_columns]

        stacked_data2.plot(kind="bar", stacked=True, width=0.46, color=local_colors, ax=ax, position=0, legend=False)
        stacked_data.plot(kind="bar", stacked=True, width=0.46, color=global_colors, ax=ax, position=1, legend=False)

        # Custom legend
        if show_legend:
            species_1 = self.legend_half_patch(global_colors[0], local_colors[0])
            species_2 = self.legend_half_patch(global_colors[1], local_colors[1])
            ax.legend(
                [species_1[0], species_2[0]],
                [
                    # f"Global and local admixture from $\\mathit{{{references[0]}}}$",
                    f"{references[0]} global and local ancestry",
                    # f"Глобальный и локальный ADMIXTURE от $\\mathit{{{references[0]}}}$",
                    # f"Global and local admixture from $\\mathit{{{references[1]}}}$",
                    f"{references[1]} global and local ancestry",
                    # f"Глобальный и локальный ADMIXTURE от $\\mathit{{{references[1]}}}$",
                ],
                handler_map={species_1[0]: species_1[1], species_2[0]: species_2[1]},
                loc=legend_loc,
                ncol=legend_ncol,
                frameon=False,
            )

        ax.set_xlim(right=len(stacked_data) - 0.5)
        ax.set_xticklabels(df.index, ha="right", rotation=rotation)
        ax.set_xlabel(xlabel)
        ax.set_yticks(yticks)

    def draw_roh_cumsum_plot(
        self,
        ax,
        data,
        genome_length,
        colors=None,
        xlim=(10e4, 150e6),
        ylim=(0, 0.35),
        xlabel="Homozygous tract length in base-pairs",
        ylabel="Cumulative genome fraction in ROHs",
        legend_title="",
        legend_loc="upper left",
        legend_ncol=2,
        linewidth=1.5,
        show_legend=True,
        figure_grid=True,
    ):
        """
        Draws a cumulative sum plot of runs of homozygosity (ROH) for multiple samples.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the plot.

        data : list of str
            A list of file paths to tab-separated files containing ROH data. Each file should have a column "length"
            indicating the length of homozygous tracts.

        genome_length : float
            The total length of the genome, used to normalize cumulative ROH values.

        colors : list of str or str, optional
            The colors used for the plot lines. If a list of colors is provided, it should match the number of samples.
            If a string is provided, it specifies a color palette name from seaborn. If None, distinct colors are
            automatically generated. Defaults to None.

        xlim : tuple of float, optional
            Limits for the x-axis (log scale), representing the range of tract lengths to display. Defaults to (10e4, 150e6).

        ylim : tuple of float, optional
            Limits for the y-axis, representing the range of cumulative genome fractions. Defaults to (0, 0.35).

        xlabel : str, optional
            The label for the x-axis. Defaults to "Homozygous tract length in base-pairs".

        ylabel : str, optional
            The label for the y-axis. Defaults to "Cumulative genome fraction in ROHs".

        legend_title : str, optional
            The title for the legend. Defaults to an empty string.

        legend_loc : str, optional
            The location of the legend box within the plot. Defaults to "upper left".

        legend_ncol : int, optional
            Number of columns in the legend. Defaults to 2.

        linewidth : float, optional
            The width of the plot lines. Defaults to 1.5.

        show_legend : bool, optional
            Whether to display a legend on the plot. Defaults to True.

        figure_grid : bool, optional
            Whether to display a grid on the plot. Defaults to True.

        Notes:
        ------
        - The function reads ROH data from each file, sorts it by tract length, and calculates the cumulative sum of
          tract lengths normalized by the genome length.
        - The plot uses a logarithmic x-axis to display a wide range of homozygous tract lengths.
        - A step plot is created for each sample, and a legend is added to differentiate between samples.
        """
        # ax.spines[["right", "top"]].set_visible(False)
        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)
        ax.set_axisbelow(True)

        if colors is None:
            colors = distinctipy.get_colors(len(data))
        else:
            if type(colors) == str:
                colors = sns.color_palette(colors, len(data))

        for count, f in enumerate(data):
            df = pd.read_csv(f, sep="\t", header=None, names=["scaffold", "start", "end", "length"])
            df.sort_values(by=["length"], inplace=True, ignore_index=True)

            df["cumsum"] = np.cumsum(df["length"]) / genome_length

            label = f.split("/")[-1].split(".")[0]
            print(f"{label}\t{len(df.index)}\t{df['length'].sum()}\t{df['length'].sum() / genome_length * 100}")

            df.loc[len(df.index)] = [None, None, None, xlim[1], df.iloc[-1]["cumsum"]]

            ax.step(df["length"], df["cumsum"], label=label, linewidth=linewidth, c=colors[count], where="post")

        ax.set_xscale("log")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        if show_legend:
            if legend_title is not None:
                # ax.legend(title=r"$\mathit{" + legend_title.replace(" ", r"\,") + "}$", loc=legend_loc, ncol=legend_ncol, fontsize=7.7)
                ax.legend(title=f"{legend_title}", loc=legend_loc, ncol=legend_ncol, fontsize=7.7)
            else:
                ax.legend(loc=legend_loc, ncol=legend_ncol, fontsize=7.7)

        if figure_grid:
            ax.grid(True, linestyle="--", alpha=0.5)

    def draw_busco_summary_plot(
        self,
        ax,
        data,
        colors=["#23b4e8", "#008dbf", "#fbbc04", "#ea4335"],
        xticks=[70, 80, 90, 100],
        xlim=(70, 100),
        bold_species_indices=None,
        vline_x_coord=77,
        show_legend=True,
        legend_labels=[
            "Complete and single-copy BUSCOs (S)",
            "Complete and duplicated BUSCOs (D)",
            "Fragmented BUSCOs (F)",
            "Missing BUSCOs (M)",
        ],
        legend_loc=(0.47, 0.95),
        legend_ncol=2,
    ):
        """
        Draws a horizontal bar plot to visualize BUSCO (Benchmarking Universal Single-Copy Orthologs) results
        for multiple species.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the plot.

        data : str or list of str
            Path to a tab-separated file with columns [file_path, species_name, species_id]
            or a list of paths to BUSCO short summary text files (short_summary_{species}.txt).

        colors : list of str, optional
            A list of four color hex codes used for the bar segments:
            - colors[0]: Color for "Complete and single-copy BUSCOs (S)"
            - colors[1]: Color for "Complete and duplicated BUSCOs (D)"
            - colors[2]: Color for "Fragmented BUSCOs (F)"
            - colors[3]: Color for "Missing BUSCOs (M)"
            Default is ["#23b4e8", "#008dbf", "#fbbc04", "#ea4335"].

        xticks : list of int, optional
            X-axis tick positions, representing percentages. Defaults to [70, 80, 90, 100].

        xlim : tuple of float, optional
            Limits for the x-axis (log scale), representing the range of tract lengths to display. Defaults to (70, 100).

        bold_species_indices : list of int, optional
            A list of indices (0-based) of species that should be labeled in bold. If specified,
            these species will have their labels printed in bold font. Indices are counted from the
            top to bottom, but can also be provided in reverse order (from bottom to top).
            For example, if `bold_species_indices=[1, 2, 3]`, the species at positions 1, 2, and 3 from the bottom
            will be labeled in bold.

        vline_x_coord : int or float, optional
            Vertical line to split barplots. Defaults to 77.

        show_legend : bool, optional
            Whether to display a legend. Default is True.

        legend_loc : tuple of float, optional
            Position of the legend box within the plot. Defaults to (-0.005, 0.97).

        legend_ncol : int, optional
            Number of columns in the legend. Defaults to 4.

        """
        # Check if data is a TSV file or a list of file paths
        if isinstance(data, str):
            df_input = pd.read_csv(data, sep="\t", header=None, names=["file_path", "Species", "ID"])
            df_input["Species"] = df_input["Species"].str.replace("_", " ")
            file_paths = df_input["file_path"].tolist()
        else:
            file_paths = data
            df_input = pd.DataFrame(
                {
                    "file_path": file_paths,
                    "Species": [file_path.split("/")[-1][14:-4].replace("_", " ") for file_path in file_paths],
                    "ID": [""] * len(file_paths),  # Placeholder if IDs are not provided
                }
            )

        # Extract data from each BUSCO summary file
        data_list = []
        for file_path, species, species_id in zip(df_input["file_path"], df_input["Species"], df_input["ID"]):
            with open(file_path, "r") as file:
                lines = file.readlines()
                target_line = lines[8]
                matches = re.findall(r"\d+\.\d+|\d+", target_line)
                numbers = [float(match) if "." in match else int(match) for match in matches]
                line = [species, species_id, numbers[1], numbers[2], numbers[3], numbers[4], numbers[5]]
                data_list.append(line)

        # Create DataFrame for plotting
        df = pd.DataFrame(data_list, columns=["Species", "ID", "S", "D", "F", "M", "N"]).iloc[::-1].reset_index(drop=True)

        # Customize plot
        # ax.spines[["left", "right", "top"]].set_visible(False)
        for spine in ["left", "right", "top"]:
            ax.spines[spine].set_visible(False)

        ax.set_axisbelow(True)
        position = range(len(file_paths))

        # Plot bars
        ax.barh(position, df["S"], height=0.9, label=legend_labels[0], color=colors[0])
        ax.barh(position, df["D"], height=0.9, left=df["S"], label=legend_labels[1], color=colors[1])
        ax.barh(position, df["F"], height=0.9, left=df["S"] + df["D"], label=legend_labels[2], color=colors[2])
        ax.barh(position, df["M"], height=0.9, left=df["S"] + df["D"] + df["F"], label=legend_labels[3], color=colors[3])

        # Add legend
        if show_legend:
            ax.legend(ncol=legend_ncol, loc=legend_loc, handlelength=0.8, frameon=False)

        # Remove default Y-axis ticks
        ax.set_yticks([])
        ax.set_yticklabels([])

        ax.set_xlim(xlim)
        ax.set_xticks(xticks)
        ax.set_xticklabels(["0%"] + [f"{i}%" for i in xticks[1:]])

        if vline_x_coord is not None:
            ax.vlines(x=vline_x_coord, ymin=-0.5, ymax=len(df.index) - 0.5, color="white", linewidth=10)
            ax.vlines(x=vline_x_coord, ymin=-0.75, ymax=len(df.index) - 0.25, color="gray", linewidth=1, linestyle="--")

        # Invert bold_species_indices
        if bold_species_indices:
            bold_species_indices = [len(df) - 1 - idx for idx in bold_species_indices]

        # Annotate species, IDs, and BUSCO scores
        ax.text(xlim[0] - 0.5, len(df.index), "Species", va="center", ha="right", fontweight="semibold")
        ax.text(xlim[0] + 0.5, len(df.index), "ID", va="center", ha="left", fontweight="semibold") if not df_input["ID"].eq("").all() else None
        for i, row in df.iterrows():
            fontweight = "bold" if bold_species_indices and i in bold_species_indices else "medium"
            ax.text(xlim[0] - 0.5, i, f"{row['Species']}", va="center", ha="right", fontweight=fontweight, style="italic")
            ax.text(xlim[0] + 0.5, i, f"{row['ID']}", va="center", ha="left", fontweight="bold", color="white")
            ax.text(
                df["S"].min() - 1,
                i,
                f"S: {row['S']}%   |   D: {row['D']}%   |   F: {row['F']}%   |   M: {row['M']}%",
                va="center",
                ha="right",
                fontweight="bold",
                color="white",
            )

    def draw_histogram_with_stats(self, ax, data, bin_width=None, show_stats=True, show_legend=True, legend_loc="upper right", legend_ncol=1):
        """
        Draws a histogram with optional statistics.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the histogram.

        data : array-like
            The dataset to plot as a histogram.

        bin_width : float, optional
            Width of each histogram bin. If None, matplotlib decides automatically.

        show_stats : bool, optional
            If True, shows mean, median, percentiles, and ±σ shaded areas.
            If False, plots only the histogram.

        show_legend : bool, optional
            Whether to display a legend. Default is True.

        legend_loc : str, optional
            Location of the legend on the plot. Default is "upper right".

        legend_ncol : int, optional
            Number of columns in the legend. Default is 1.
        """
        import numpy as np

        # Hide unnecessary spines
        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)

        # Determine bins based on bin_width
        if bin_width is not None:
            data_min, data_max = np.min(data), np.max(data)
            bins = np.arange(data_min, data_max + bin_width, bin_width)
        else:
            bins = "auto"

        # Draw histogram
        ax.hist(data, bins=bins, edgecolor="black")

        if show_stats:
            # Calculate statistics
            mean = np.mean(data)
            median = np.median(data)
            std_dev = np.std(data)
            percentile_5 = np.percentile(data, 5)
            percentile_95 = np.percentile(data, 95)

            # Plot key statistics
            ax.axvline(mean, color="red", linestyle="--", label=f"Mean: {mean:.2f}")
            ax.axvline(median, color="blue", linestyle="--", label=f"Median: {median:.2f}")
            ax.axvline(percentile_5, color="purple", linestyle="-.", label=f"5th percentile: {percentile_5:.2f}")
            ax.axvline(percentile_95, color="purple", linestyle="-.", label=f"95th percentile: {percentile_95:.2f}")

            # Shaded regions for ±1σ, ±2σ, ±3σ
            # for i, color in zip(range(1, 4), ["orange", "yellow", "green"]):
                # ax.axvspan(mean - i * std_dev, mean + i * std_dev, color=color, alpha=0.1)
                # left_value = mean - i * std_dev
                # right_value = mean + i * std_dev
                # ax.axvline(left_value, color=color, linestyle="--", label=f"Mean - {i}σ: {left_value:.2f}")
                # ax.axvline(right_value, color=color, linestyle="--", label=f"Mean + {i}σ: {right_value:.2f}")

        # Legend
        if show_legend and show_stats:
            ax.legend(ncol=legend_ncol, loc=legend_loc)

    def classify_roh(self, length):
        if length < 1_000_000:
            return "S"
        elif length >= 10_000_000:
            return "UL"
        else:
            return "L"

    def draw_roh_barplot(
        self,
        ax,
        data,
        genome_length,
        colors={"N": "#23b4e8", "S": "#008dbf", "L": "#fbbc04", "UL": "#ea4335"},
        xticks=[50, 60, 70, 80, 90, 100],
        xlim=(45, 100),
        vline_x_coord=65,
        sorting=True,
        groups=None,
        show_annotation=False,
        show_legend=True,
        legend_loc=(0.25, 0.97),
        legend_ncol=4,
    ):
        """
        Visualizes the distribution of ROHs (Runs of Homozygosity) across genome categories
        using a horizontal stacked bar plot.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the bar plot.

        data : list of str
            A list of file paths to BED files, where each file contains ROH data for a single sample.
            Each file should have tab-delimited columns: `scaffold`, `start`, `end`, `length`.

        genome_length : int or dict
            The total length of the genome, used to calculate percentages.
            If groups, dictionary with groups (keys) and genome sizes (values).

        colors : dict, optional
            A dictionary mapping ROH categories to their respective colors.
            Default is `{"N": "#23b4e8", "S": "#008dbf", "L": "#fbbc04", "UL": "#ea4335"}`.

        xticks : list of int, optional
            X-axis tick positions, representing percentages. Defaults to [0, 25, 50, 75, 100].

        xlim : tuple of float, optional
            Limits for the x-axis (log scale), representing the range of tract lengths to display. Defaults to (45, 100).

        sorting : bool, optional
            If True, bars are sorted by the percentage of Ultra Long ROHs (UL) in descending order.
            If False, bars follow the reverse order of `data`. Default is False.

        groups : dict, optional
            A dictionary where keys are group names and values are lists of sample names.
            This allows grouping samples on the y-axis with the group name followed by sample names.
            If None, no grouping is applied. Default is None.

        show_annotation : bool, optional
            Whether to display an annotation. Default is False.

        show_legend : bool, optional
            Whether to display a legend. Default is True.

        legend_loc : str, optional
            Location of the legend on the plot. Default is (0.25, 0.97).

        legend_ncol : int, optional
            Number of columns in the legend. Default is 4.

        Additional Details:
        --------------------
        - Each BED file represents a sample and contains ROH data with the following format:
            * `scaffold`: Chromosome or scaffold identifier.
            * `start`: Start position of the ROH.
            * `end`: End position of the ROH.
            * `length`: Length of the ROH (in base pairs).
        - Categories:
            * `"N"`: Non-ROHs (portion of the genome not covered by ROHs).
            * `"S"`: Short ROHs (<1,000,000 bp).
            * `"L"`: Long ROHs (1,000,000–10,000,000 bp).
            * `"UL"`: Ultra Long ROHs (≥10,000,000 bp).
        - The function calculates the percentage of the genome occupied by each ROH category for each sample.
        - The stacked bar plot ensures each bar represents 100% of the genome for a given sample.
        """
        # --- Customize axes
        for spine in ["left", "right", "top"]:
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_zorder(0)

        # --- Load data
        sample_data = {}
        sample_names = []
        for file_path in data:
            sample = os.path.basename(file_path).split(".")[0]
            sample_names.append(sample)
            df = pd.read_csv(file_path, sep="\t", header=None, names=["scaffold", "start", "end", "length"])
            df["classification"] = df["length"].apply(self.classify_roh)
            total = df.groupby("classification")["length"].sum().to_dict()

            group = next((g for g, s in (groups or {}).items() if sample in s), None)
            glen = genome_length[group] if groups else genome_length

            for cat in ["S", "L", "UL"]:
                total[cat] = total.get(cat, 0) / glen * 100
            total["N"] = max(0, 100 - sum(total.get(c, 0) for c in ["S", "L", "UL"]))

            sample_data[sample] = total

        # --- Build dataframe
        df = pd.DataFrame.from_dict(sample_data, orient="index", columns=["N", "S", "L", "UL"]).fillna(0)

        if sorting:
            df = df.sort_values(by=["UL", "L", "S", "N"], ascending=False)

        if not groups:
            groups = {"All": df.index.tolist()}

        # --- Order samples by groups
        ordered = []
        for g, s in groups.items():
            s = [x for x in s if x in df.index]
            if sorting:
                s = df.loc[s].sort_values(by=["UL", "L", "S", "N"], ascending=False).index.tolist()
            ordered.extend(s)
        df = df.loc[ordered]

        # --- Plot stacked bars
        labels = {
            "N": "Non-RoHs (N)",
            "S": "Short RoHs (S)",
            "L": "Long RoHs (L)",
            "UL": "Ultra Long RoHs (UL)",
        }

        left = pd.Series(0, index=df.index)
        for cat, color in colors.items():
            ax.barh(
                df.index,
                df[cat],
                left=left,
                color=color,
                label=labels[cat],
                height=0.9,
            )
            left += df[cat]

        # --- Initial shaded area (white span before first tick)
        if xticks[0] != 0:
            ax.axvspan(xlim[0], xticks[0], color="white")

        # --- Configure x/y axes
        ax.set_xlim(xlim)
        ax.set_xticks(xticks)
        ax.set_ylim(bottom=-1.25)
        ax.spines["bottom"].set_bounds(xticks[0], xlim[1])
        ax.yaxis.set_ticks_position("none")
        ax.set_yticks([])
        ax.set_yticklabels([])

        for i, sample in enumerate(df.index):
            ax.text(xticks[0] - 0.5, i, sample, va="center", ha="right")

        if vline_x_coord is not None:
            ax.set_xticklabels(["0%"] + [f"{i}%" for i in xticks[1:]])
            ax.vlines(x=vline_x_coord, ymin=-1.5, ymax=len(df.index) - 0.5, color="white", linewidth=10, clip_on=False)
            ax.vlines(x=vline_x_coord, ymin=-1.5, ymax=len(df.index) - 0.25, color="gray", linewidth=1, linestyle="--", clip_on=False)
        else:
            ax.set_xticklabels([f"{i}%" for i in xticks])

        # --- Legend
        if show_legend:
            ax.legend(loc=legend_loc, ncol=legend_ncol, handlelength=0.8, frameon=False)

        # --- Draw group separators
        if groups and "All" not in groups:
            for g, s in groups.items():
                s = [x for x in df.index if x in s]
                if not s:
                    continue
                start, end = df.index.get_loc(s[0]), df.index.get_loc(s[-1])
                ax.vlines(x=xlim[0], ymin=start - 0.25, ymax=end + 0.25, colors="black", linewidth=1)
                ax.text(xlim[0] - 1, (start + end) / 2, g, va="center", ha="right", fontstyle="italic")

        # --- Optional annotation per sample
        if show_annotation:
            for i, sample in enumerate(df.index):
                vals = df.loc[sample]
                text = f'  N: {vals["N"]:.1f}%   |   S: {vals["S"]:.1f}%   |   L: {vals["L"]:.1f}%   |   UL: {vals["UL"]:.1f}%'
                ax.text(
                    xticks[1],
                    i,
                    text,
                    va="center",
                    ha="left",
                    fontsize=8,
                    color="white",
                    fontweight="bold",
                )

    def read_series(self, s):
        return pd.read_csv(s, header=None).squeeze("columns") if os.path.exists(s) else pd.Series(s.split(","))

    def rgb_tuple_to_hex(self, rgb_tuple):
        color_code = "#"
        for i in [0, 1, 2]:
            color_code += "{:02X}".format(int(255 * rgb_tuple[i]))
        return color_code

    def draw_variant_window_densities(
        self,
        ax,
        input_file,
        input_type="vcf",
        output_prefix=None,
        output_formats=[],
        title=False,
        title_fontsize=10,
        window_size=1000000,
        window_step=None,
        density_multiplier=1000,
        scaffold_white_list=pd.Series(dtype=str),
        scaffold_ordered_list=pd.Series(dtype=str),
        scaffold_length_file=[],
        scaffold_syn_file=None,
        syn_file_key_column=0,
        syn_file_value_column=1,
        figure_width=10,
        figure_height_per_scaffold=0.5,
        figure_header_height=0,
        masking_track=None,
        masking_color="grey",
        masking_threshold=0.1,
        colormap="jet",
        custom_color_list=None,
        density_thresholds=(0.0, 0.1, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5),
        density_thresholds_expression_type="left_open",
        skip_top_interval=False,
        skip_bottom_interval=False,
        test_colormaps=False,
        hide_track_label=True,
        subplots_adjust_left=None,
        subplots_adjust_top=None,
        subplots_adjust_right=None,
        subplots_adjust_bottom=None,
        only_count=False,
        x_tick_fontsize=None,
        stranded=False,
        rounded=True,
        middle_break=False,
        stranded_end=False,
        feature_name="SNPs",
        centromere_bed=None,
        highlight_file=None,
        show_legend=True,
    ):
        # coverage=None,
        # scaffold_column_name="scaffold",
        # window_column_name="window",
        # coverage_column_name="median",
        # mean_coverage=None,
        # max_coverage_threshold=2.5,
        # min_coverage_threshold=0.5,
        # scaffold_black_list=pd.Series(dtype=str),
        # sort_scaffolds=False,
        # masking_gff_list=None,

        # if not output_prefix
        #     raise ValueError("Output prefix is required.")

        scaffold_white_list = self.read_series(scaffold_white_list)
        scaffold_ordered_list = self.read_series(scaffold_ordered_list)

        scaffold_ordered_list = scaffold_ordered_list[::-1]

        if isinstance(scaffold_ordered_list, list):
            if not scaffold_ordered_list:
                scaffold_ordered_list = scaffold_white_list
        else:
            if scaffold_ordered_list.empty:
                scaffold_ordered_list = scaffold_white_list


        if masking_track is not None:
            masking_track = pd.read_csv(masking_track, header=None, index_col=None, sep="\t")

            #detect type of masking track and assign colors
            if len(masking_track.columns) == 3:
                masking_track.columns = pd.Index(["scaffold", "start", "end"])
                masking_track["color"] = masking_color if masking_color is not None else masking_color
            elif len(masking_track.columns) == 4:
                masking_track.columns = pd.Index(["scaffold", "start", "end", "color"])
                if pd.api.types.is_integer_dtype(masking_track["color"]):
                    # assume that the fourth column contains counts of masked bases
                    threshold = int(window_size * masking_threshold)
                    masking_track = masking_track[masking_track["color"] >= threshold]
                    masking_track["color"] = masking_color
                elif pd.api.types.is_float_dtype(masking_track["color"]):
                    # assume that the fourth column contains counts of masked bases
                    masking_track = masking_track[masking_track["color"] >= masking_threshold]
                    masking_track["color"] = masking_color
                else:
                    # assume that the forth column is string like
                    raise ValueError("ERROR!!! Preset for string-like values in the fourth column is not implemented yet.")
                    #masking_track[masking_track["color"] == 'masked', "color"] = default_masking_color
            else:
                raise ValueError(f"ERROR! Masking track contains {len(masking_track.columns)} columns instead of 3 or 4!")

            masking_track["color"] = masking_track["color"].apply(mpl.colors.to_hex)
            masking_track = masking_track.set_index("scaffold")
            #print(masking_track)

        if custom_color_list is not None:
            if len(custom_color_list) != len(density_thresholds):
                raise ValueError(
                    "ERROR!!! Custom color list is set, but the number of colors ({0}) in the list is not equal to the number of the thresholds (1)!".format(
                        len(custom_color_list),
                    )
                )

        variants = CollectionVCF(input_file, parsing_mode="only_coordinates")

        chr_len_df = (
            pd.read_csv(scaffold_length_file, sep="\t", header=None, index_col=0) if scaffold_length_file else deepcopy(variants.scaffold_length)
        )
        chr_len_df.index = pd.Index(map(str, chr_len_df.index))
        chr_len_df.index.name = "scaffold"
        chr_len_df.columns = ["length"]

        chr_syn_dict = SynDict(filename=scaffold_syn_file, key_index=syn_file_key_column, value_index=syn_file_value_column)

        if centromere_bed:
            centromere_df = pd.read_csv(centromere_bed, usecols=(0, 1, 2), index_col=0, header=None, sep="\t", names=["scaffold_id", "start", "end"])
            centromere_df.rename(index=chr_syn_dict, inplace=True)
        else:
            centromere_df = None

        if input_type == "vcf":
            count_df = StatsVCF.count_variants_in_windows(
                variants,
                window_size,
                window_step,
                reference_scaffold_lengths=chr_len_df,
                ignore_scaffolds_shorter_than_window=True,
                # output_prefix=output_prefix,
                skip_empty_windows=False,
                expression=None,
                per_sample_output=False,
                scaffold_white_list=scaffold_white_list,
                scaffold_syn_dict=chr_syn_dict,
            )
            #print("Counts")
            #print(count_df)
            # instead of actual end of the window in track df start + window_step is used to avoid possible visualization artefacts.
            # Masking track must be modified to take it in account
            # It applies only for track_df calculated from vcf
            feature_df, track_df = StatsVCF.convert_variant_count_to_feature_df(count_df,
                                                                                window_size,
                                                                                window_step)
            if masking_track is not None:
                masking_track["end"] = masking_track["start"] + window_step
            #print("Feature df")
            #print(feature_df)
            #print("track_df")
            #print(track_df)

            # feature_df.to_csv("{}.features.counts".format(output_prefix), sep="\t", header=True, index=True)
            feature_df[feature_df.columns[-1]] = feature_df[feature_df.columns[-1]] * float(density_multiplier) / float(window_size)

            # feature_df.to_csv("{}.features.bed".format(output_prefix), sep="\t", header=True, index=True)

        elif input_type == "bedgraph":
            track_df = pd.read_csv(
                input_file,
                sep="\t",
                names=["scaffold", "start", "end", "value"],
                header=None,
                index_col=0,
                na_values=".",
                dtype={"scaffold": str, "start": int, "end": int, "value": float},
            )
            if scaffold_syn_file:
                track_df.rename(index=chr_syn_dict, inplace=True)
            track_df["value"] = track_df["value"].astype(float)

        if scaffold_syn_file:
            chr_len_df.rename(index=chr_syn_dict, inplace=True)
            if masking_track is not None:
                masking_track.rename(index=chr_syn_dict, inplace=True)
                #print(masking_track)
        #print(track_df)
        # scale counts
        track_df[track_df.columns[-1]] = track_df[track_df.columns[-1]] * float(density_multiplier) / float(window_size)

        if track_df.index.nlevels > 1:
            # drop second level of index if it was added by groupby
            track_df = (
                track_df.groupby("scaffold").apply(lambda df: df[df["end"] <= chr_len_df.loc[df.index[0], "length"]]).reset_index(level=1, drop=True)
            )

        if not only_count:
            if custom_color_list is not None:
                cmap_list = ["custom_list"]
            else:
                cmap_list = Visualization.colormap_list if test_colormaps else [colormap]

            for colormap in cmap_list:
                if colormap == "custom_list":
                    colors = custom_color_list
                else:
                    cmap = plt.get_cmap(colormap, len(density_thresholds))
                    colors = [self.rgb_tuple_to_hex(cmap(i)[:3]) for i in range(0, len(density_thresholds))]

                color_expression = partial(
                    Visualization.color_threshold_expression,
                    thresholds=density_thresholds,
                    colors=colors,
                    background="white",
                    interval_type=density_thresholds_expression_type,
                    skip_top_interval=skip_top_interval,
                    skip_bottom_interval=skip_bottom_interval,
                )

                track_with_colors_df = Visualization.add_color_to_track_df(
                    track_df,
                    color_expression,
                    value_column_index=-1,  # TODO fix it, add support for multiple tracks in the file
                )
                #apply track mask if exists
                if masking_track is not None:
                    #masking_track
                    track_with_colors_df = Visualization.apply_masking_to_track_df(track_with_colors_df, masking_track, masking_color=masking_color)
                    track_with_colors_df.to_csv("{}.{}.masked.track.bed".format(output_prefix,
                                                                                colormap), sep="\t", header=True, index=True)
                # track_with_colors_df.to_csv("{}.{}.track.bed".format(output_prefix, colormap), sep="\t", header=True, index=True)
                # print(feature_with_colors_df)
                # print(scaffold_ordered_list)
                Visualization.draw_features(
                    {"TR": track_with_colors_df},
                    chr_len_df,
                    scaffold_ordered_list,
                    output_prefix,
                    legend=(
                        Visualization.density_legend(
                            colors,
                            density_thresholds,
                            feature_name=feature_name,
                            interval_type=density_thresholds_expression_type,
                            skip_top_interval=skip_top_interval,
                            skip_bottom_interval=skip_bottom_interval,
                            masking_color=masking_color
                        )
                        if show_legend
                        else None
                    ),
                    centromere_df=centromere_df,
                    highlight_df=highlight_file,
                    figure_width=figure_width,
                    figure_height_per_scaffold=figure_height_per_scaffold,
                    figure_header_height=figure_header_height,
                    dpi=300,
                    default_color="red",
                    title=title,
                    extensions=output_formats,
                    feature_start_column_id="start",
                    feature_end_column_id="end",
                    feature_color_column_id="color",
                    feature_length_column_id="length",
                    subplots_adjust_left=subplots_adjust_left,
                    subplots_adjust_bottom=subplots_adjust_bottom,
                    subplots_adjust_right=subplots_adjust_right,
                    subplots_adjust_top=subplots_adjust_top,
                    show_track_label=not hide_track_label,
                    show_trackgroup_label=True,
                    close_figure=False,
                    subplot_scale=False,
                    track_group_scale=False,
                    # track_group_distance=2,
                    # xmax_multiplier=1.3,
                    # ymax_multiplier=1.00,
                    stranded_tracks=stranded,
                    rounded_tracks=rounded,
                    middle_break=middle_break,
                    stranded_end_tracks=stranded_end,
                    xtick_fontsize=x_tick_fontsize,
                    subplot_title_fontsize=title_fontsize,
                    subplot_title_fontweight="bold",
                    axes=ax,
                )

    def draw_features(
        self,
        ax,
        input_file,
        input_type="str",
        header=None,
        legend=None,
        output_prefix=None,
        output_formats=[],
        title=False,
        start_column_name=None,
        end_column_name=None,
        color_column_name=None,
        default_color="tab:blue",
        scaffold_white_list=pd.Series(dtype=str),
        scaffold_ordered_list=pd.Series(dtype=str),
        scaffold_length_file=None,
        scaffold_syn_file=None,
        syn_file_key_column=0,
        syn_file_value_column=1,
        colormap="jet",
        hide_track_label=True,
        x_tick_type="nucleotide",
        feature_shape="rectangle",
        subplots_adjust_left=None,
        subplots_adjust_top=None,
        subplots_adjust_right=None,
        subplots_adjust_bottom=None,
        figure_width=10,
        figure_height_per_scaffold=0.5,
        figure_header_height=0,
        verbose=False,
        subplot_scale=False,
        track_group_scale=False,
        x_tick_fontsize=None,
        stranded=False,
        rounded=True,
        stranded_end=False,
        fill_empty_tracks=False,
        empty_color="lightgrey",
        centromere_bed=None,
        highlight_file=None,
        title_fontsize=20,
    ):
        # scaffold_column_name=None,
        # scaffold_black_list=pd.Series(dtype=str),
        # sort_scaffolds=False,
        # figure_height_per_scaffold=0.5,
        # print(scaffold_ordered_list)
        scaffold_white_list = self.read_series(scaffold_white_list)
        scaffold_ordered_list = self.read_series(scaffold_ordered_list)
        scaffold_ordered_list = scaffold_ordered_list[::-1]

        chr_syn_dict = SynDict(filename=scaffold_syn_file, key_index=syn_file_key_column, value_index=syn_file_value_column)

        # print(scaffold_ordered_list)
        if isinstance(scaffold_ordered_list, list):
            if not scaffold_ordered_list:
                scaffold_ordered_list = deepcopy(scaffold_white_list)
                scaffold_ordered_list.replace(chr_syn_dict, inplace=True)
        else:
            # print("AAAA")
            if scaffold_ordered_list.empty:
                scaffold_ordered_list = deepcopy(scaffold_white_list)
                scaffold_ordered_list.replace(chr_syn_dict, inplace=True)

        if centromere_bed:
            centromere_df = pd.read_csv(centromere_bed, usecols=(0, 1, 2), index_col=0, header=None, sep="\t", names=["scaffold_id", "start", "end"])
            centromere_df.rename(index=chr_syn_dict, inplace=True)
        else:
            centromere_df = None
        try:
            if input_type == "str":
                feature_df = CollectionSTR(
                    in_file=input_file,
                    records=None,
                    format="filtered_str",
                    parsing_mode="all",
                    black_list=(),
                    white_list=(),
                    syn_dict=chr_syn_dict,
                )

                feature_df.records.set_index("scaffold_id", inplace=True)

                feature_start_column_id = start_column_name if start_column_name else "start"
                feature_end_column_id = end_column_name if end_column_name else "end"

            elif input_type == "gff":
                feature_df = CollectionGFF(in_file=input_file, parsing_mode="only_coordinates")

                feature_start_column_id = start_column_name if start_column_name else "start"
                feature_end_column_id = end_column_name if end_column_name else "end"

            elif input_type == "bed":
                feature_df = CollectionBED(in_file=input_file, parsing_mode="coordinates_only", format="bed")

                feature_start_column_id = "start"
                feature_end_column_id = "end"

            elif input_type == "bed_table":
                feature_df = CollectionBED(in_file=input_file, parsing_mode="complete", format="bed", header_in_file=False)
                feature_start_column_id = "start"
                feature_end_column_id = "end"

            elif input_type == "bedgraph":
                feature_df = CollectionBED(in_file=input_file, parsing_mode="all", format="bed")
                feature_df.records.columns = pd.Index(["start", "end", "value"])
                feature_df.records.index.name = "scaffold"
                feature_start_column_id = "start"
                feature_end_column_id = "end"
            elif input_type == "bed_track":
                feature_df = CollectionBED(
                    in_file=input_file,
                    parsing_mode="all",
                    format="bed_track",
                )
                print(feature_df.records)
                feature_df.records.columns = pd.Index(["start", "end", "value", "color"])
                feature_df.records.index.name = "scaffold"
                feature_start_column_id = "start"
                feature_end_column_id = "end"
                color_column_name = "color"
            elif input_type == "tab6":
                feature_df = CollectionBLAST(in_file=input_file, parsing_mode="complete")
                feature_df.records.reset_index(level="query_id", inplace=True)
                feature_start_column_id = start_column_name if start_column_name else "target_start"
                feature_end_column_id = end_column_name if end_column_name else "target_end"
            elif input_type == "tab6_colored":
                feature_df = CollectionBLAST(in_file=input_file, parsing_mode="complete", format="tab6_colored", header=header)
                feature_df.records.reset_index(level="query_id", inplace=True)
                feature_start_column_id = start_column_name if start_column_name else "target_start"
                feature_end_column_id = end_column_name if end_column_name else "target_end"
            else:
                raise ValueError("ERROR!!! Unrecognized input type ({}). ".format(input_type))
        except pd.errors.EmptyDataError:
            print(
                "Empty input file. Silent exit."
            )  # try-except added to handle case when input file is empty without raising exception. For use in snakemake
            exit(0)

        legend_df = pd.read_csv(legend, header=None, index_col=0, sep="\t") if legend else None

        # print(scaffold_white_list)
        # print(feature_df.records)
        scaffold_to_keep = StatsVCF.get_filtered_entry_list(
            feature_df.records.index.get_level_values(level=0).unique().to_list(), entry_white_list=scaffold_white_list
        )
        # print(scaffold_to_keep)
        # print(scaffold_to_keep)
        # remove redundant scaffolds
        # print(scaffold_white_list)
        # print(feature_df.records)
        # print(scaffold_to_keep)
        feature_df.records = feature_df.records[feature_df.records.index.isin(scaffold_to_keep)]
        # print("BBBBBBbb")
        # print(scaffold_white_list)
        # print(scaffold_ordered_list)
        if not scaffold_white_list.empty:
            scaffold_ordered_list = scaffold_ordered_list[scaffold_ordered_list.isin(pd.Series(scaffold_white_list).replace(chr_syn_dict))]
        # print("CCCC")
        # print(scaffold_ordered_list)
        # print(scaffold_to_keep)
        # print(pd.Series(scaffold_to_keep).replace(chr_syn_dict))
        chr_len_df = pd.read_csv(scaffold_length_file, sep="\t", header=None, names=("scaffold", "length"), index_col=0)
        chr_len_df.index = pd.Index(list(map(str, chr_len_df.index)))

        # print(chr_len_df)

        # print(feature_df.records)
        if scaffold_syn_file:
            chr_len_df.rename(index=chr_syn_dict, inplace=True)
            feature_df.records.rename(index=chr_syn_dict, inplace=True)
        if verbose:
            print(chr_syn_dict)
            print(feature_df.records)

        # print(chr_len_df)
        # print({"features": feature_df})
        # print(feature_df.records.columns)
        # print(feature_df.records)
        # print(chr_len_df)
        # print(scaffold_ordered_list)
        Visualization.draw_features(
            {"features": feature_df},
            chr_len_df,
            scaffold_ordered_list,
            output_prefix,
            legend=Visualization.feature_legend(legend_df, colormap=colormap),
            # legend_df=legend_df,
            centromere_df=centromere_df,
            highlight_df=highlight_file,
            figure_width=figure_width,
            figure_height_per_scaffold=figure_height_per_scaffold,
            figure_header_height=figure_header_height,
            dpi=300,
            # colormap=None, thresholds=None, colors=None, background=None,
            default_color=default_color,
            title=title,
            extensions=output_formats,
            feature_shape=feature_shape,
            feature_start_column_id=feature_start_column_id,
            feature_end_column_id=feature_end_column_id,
            feature_color_column_id=color_column_name,
            feature_length_column_id="length",
            subplots_adjust_left=subplots_adjust_left,
            subplots_adjust_bottom=subplots_adjust_bottom,
            subplots_adjust_right=subplots_adjust_right,
            subplots_adjust_top=subplots_adjust_top,
            show_track_label=not hide_track_label,
            show_trackgroup_label=True,
            subplot_scale=subplot_scale,
            track_group_scale=track_group_scale,
            stranded_tracks=stranded,
            rounded_tracks=rounded,
            stranded_end_tracks=stranded_end,
            fill_empty_tracks=fill_empty_tracks,
            empty_color=empty_color,
            xtick_fontsize=x_tick_fontsize,
            subplot_title_fontsize=title_fontsize,
            subplot_title_fontweight="bold",
            x_tick_type=x_tick_type,
            # xmax_multiplier=2,
            axes=ax,
        )

    def draw_psmc_plot(
        self,
        ax,
        diploid_data,
        round_data=None,
        ids=None,
        n=100,
        colorlist=None,
        xlim=(50000, 7000000),
        ylim=(0, 200),
        xticks=[100000, 200000, 500000, 1000000, 2000000, 3000000, 5000000],
        xtick_labels=["100k", "200k", "500k", "1M", "2M", "3M", "5M"],
        xlabel=r"Years Ago ($\mu=4.64e-9$, g=5)",
        ylabel=r"Effective population size, $10^4$",
        title=None,
        figure_grid=True,
        show_legend=True,
        legend_loc="upper right",
        legend_ncol=1,
    ):
        """
        Draw a PSMC plot using diploid and round data.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes object where the plot will be drawn.

        diploid_data : list of str
            A list of file paths containing the diploid PSMC data (each file with columns: 'x', 'y', 'z', 'w', 'h').
            The 'x' column represents time (in years ago), and 'y' represents the effective population size.

        round_data : list of str, optional
            A list of file paths containing the round PSMC data (default is None). Each file is expected to follow the
            same format as the diploid data. The round data will be visualized with a reduced alpha transparency.

        ids : list of str
            A list of sample IDs per diploid PSMC data.

        n : int, optional, default: 100
            Number of round data.

        colorlist : list of str or None, optional
            A list of colors for the plot. If None, colors are generated automatically. If a string is provided, it is
            interpreted as a seaborn color palette.

        xticks : list of int, optional, default: [50000, 70000, 100000, ..., 20000000]
            The ticks to display on the x-axis (years ago), with corresponding labels adjusted according to the `scale`.

        ylim : tuple of (int, int), optional, default: (0, 200)
            The y-axis limits, representing the effective population size.

        mu : float, optional, default: 2.2e-9
            The mutation rate used for PSMC analysis (used in axis label formatting).

        g : int, optional, default: 5
            The generation time used for PSMC analysis (used in axis label formatting).

        figure_grid : bool, optional
            Whether to show grid lines on the plot. Default is True.

        show_legend : bool, optional, default: True
            Whether to display the legend on the plot.

        legend_loc : str, optional, default: "upper right"
            The location of the legend on the plot.

        legend_ncol : int, optional, default: 1
            The number of columns in the legend.
        """
        if colorlist is None:
            colorlist = distinctipy.get_colors(len(diploid_data))
        else:
            if type(colorlist) is str:
                colorlist = sns.color_palette(colorlist, len(diploid_data))
            elif type(colorlist) is list or type(colorlist) is dict:
                colorlist = colorlist

        for i, diploid in enumerate(diploid_data):
            sample_name = diploid.split("/")[-1].split(".")[0] if ids is None else ids[i]
            data = pd.read_csv(diploid, names=["x", "y", "z", "w", "h"], sep="\t")
            data = data[data["x"] >= 0]
            if type(colorlist) is dict:
                color = colorlist[sample_name]
            else:
                color = colorlist[i]
            ax.step(data["x"], data["y"], where="post", color=color, linewidth=2, label=sample_name)

        if round_data is not None:
            for round, color in zip(round_data, colorlist):
                data = pd.read_csv(round, names=["x", "y", "z", "w", "h"], sep="\t")
                data = data[data["x"] >= 0]
                for i in range(0, len(data), len(data) // n):
                    ax.step(
                        data["x"].iloc[i : i + len(data) // n],
                        data["y"].iloc[i : i + len(data) // n],
                        where="post",
                        color=color,
                        linewidth=1,
                        alpha=0.1,
                    )

        # Customize plot
        # ax.spines[["right", "top"]].set_visible(False)
        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)

        ax.set_xscale("log")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if title:
            ax.set_title(title)

        if xticks:
            ax.set_xticks(xticks)
            ax.set_xticklabels(xtick_labels)

        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)

        if show_legend:
            ax.legend(loc=legend_loc, ncol=legend_ncol, frameon=False)

        if figure_grid:
            ax.grid(linestyle="--", alpha=0.2)

    def draw_pca_plot(
        self,
        ax,
        eigenvec_file,
        eigenval_file,
        colors=None,
        dot_size=70,
        highlight_samples=None,
        highlight_samples_fontsize=None,
        figure_grid=True,
        show_legend=True,
        legend_loc="upper left",
        legend_ncol=4,
    ):
        """
        Visualizes the results of Principal Component Analysis (PCA) from PLINK using a scatter plot.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the PCA plot.

        eigenvec_file : str
            Path to the eigenvectors file generated by PLINK. This file contains the PCA results for the samples, with one row per sample and one column per principal component.

        eigenval_file : str
            Path to the eigenvalues file generated by PLINK. This file contains the eigenvalues corresponding to each principal component.

        colors : list of str, optional
            A list of colors to be used for the dots representing each sample. If None, a default color map will be used. Default is None.

        dot_size : int, optional
            Size of the dots in the scatter plot. Default is 70.

        highlight_samples : list of str or bool, optional
            - If None or False: No samples are highlighted with text labels.
            - If True: All samples are labeled.
            - If a list of str: Only the samples listed are labeled.
            The labels' positions will be adjusted using adjust_text to minimize overlap. Default is None.

        highlight_samples_fontsize : int, optional
            Font size for the labels of the highlighted samples. If None, the default font size is used. Default is None.

        figure_grid : bool, optional
            Whether to show grid lines on the plot. Default is True.

        show_legend : bool, optional
            Whether to display a legend on the plot. Default is True.

        legend_loc : str, optional
            Location of the legend on the plot. Default is "upper left".

        legend_ncol : int, optional
            Number of columns in the legend. Default is 4.
        """
        pca_df = pd.read_csv(eigenvec_file, delim_whitespace=True, header=None)
        pca_df = pca_df.iloc[:, 1:].copy()
        pca_df.columns = ["Sample"] + [f"PC{i}" for i in range(1, pca_df.shape[1])]

        eigenvals = pd.read_csv(eigenval_file, header=None).squeeze("columns")
        explained_variance = eigenvals / eigenvals.sum() * 100

        if colors is None:
            colors = distinctipy.get_colors(len(pca_df))
        elif isinstance(colors, str):
            colors = sns.color_palette(colors, len(pca_df))[::-1]

        for i in range(pca_df.shape[0]):
            ax.scatter(
                x=pca_df.loc[i, "PC1"],
                y=pca_df.loc[i, "PC2"],
                color=colors[i] if len(colors) > i else colors[0],
                label=pca_df.loc[i, "Sample"],
                s=dot_size
            )

        # highlight samples
        texts = []
        samples_to_label = []

        if highlight_samples is True:
            samples_to_label = pca_df["Sample"].tolist()
        elif isinstance(highlight_samples, list):
            samples_to_label = highlight_samples

        if samples_to_label:
            for i, row in pca_df.iterrows():
                if row["Sample"] in samples_to_label:
                    texts.append(
                        ax.text(
                            row["PC1"],
                            row["PC2"],
                            row["Sample"],
                            fontsize=highlight_samples_fontsize,
                        )
                    )

            if texts:
                adjust_text(
                    texts,
                    ax=ax,
                    only_move={'points': 'xy', 'text': 'xy'},
                    arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
                    # force_text=1.0
                )

        # 5. Кастомизация графика
        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)

        ax.set_xlabel(f"Principal Component 1 ({explained_variance[0]:.2f}%)")
        ax.set_ylabel(f"Principal Component 2 ({explained_variance[1]:.2f}%)")

        if show_legend:
            ax.legend(loc=legend_loc, ncol=legend_ncol, handlelength=0.8, frameon=True, fontsize=8)

        if figure_grid:
            ax.set_axisbelow(True)
            ax.grid(linestyle="--", alpha=0.5)

    def draw_kmers_log_scale(
        self,
        ax,
        data,
        colors,
        linewidth=2,
        xlim=(1, 1e5),
        ylabel="Number of distinct 23-mers",
        xlabel=None,
        title=None,
        figure_grid=True,
        show_legend=True,
        legend_loc="upper right",
        legend_ncol=2,
    ):
        """
        Visualizes the k-mer distribution on a log scale on the specified subplot.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the k-mer distribution plot.

        data : list of str
            List of file paths to histogram files, each containing k-mer counts and their frequencies.

        colors : list of str or str, optional
            The colors used for the plot lines. If a list of colors is provided, it should match the number of samples.
            If a string is provided, it specifies a color palette name from seaborn. If None, distinct colors are
            automatically generated. Defaults to None.

        linewidth : float, optional
            The width of the plot lines. Default is 2.

        xlim : tuple of float, optional
            Limits for the x-axis (log scale), representing the range of tract lengths to display. Defaults to (10e4, 150e6).

        title : str, optional
            The title of the plot. If None, no title is displayed. Default is None.

        figure_grid : bool, optional
            Whether to show grid lines on the plot. Default is True.

        show_legend : bool, optional
            Whether to display a legend on the plot. Default is True.

        legend_loc : str, optional
            Location of the legend on the plot. Default is "upper right".

        legend_ncol : int, optional
            Number of columns in the legend. Default is 2.
        """
        # ax.spines[["right", "top"]].set_visible(False)
        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)
        ax.set_axisbelow(True)

        if colors is None:
            colors = distinctipy.get_colors(len(data))
        else:
            if type(colors) == str:
                colors = sns.color_palette(colors, len(data))

        for i, histo_file in enumerate(data):
            kmer_counts, frequencies = [], []
            with open(histo_file, "r") as f:
                for line in f:
                    kmer_count, frequency = map(int, line.strip().split())
                    kmer_counts.append(kmer_count)
                    frequencies.append(frequency)
            sample = histo_file.split("/")[-1].split(".")[0]
            ax.plot(kmer_counts, frequencies, "-", label=sample, color=colors[i], linewidth=linewidth)

        ax.set_xlim(xlim)
        ax.set_ylim(min(frequencies), max(frequencies))

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if title:
            ax.set_title(r"$\mathit{" + title.replace(" ", r"\,") + "}$")
        if show_legend:
            ax.legend(loc=legend_loc, ncol=legend_ncol, handlelength=0.8, frameon=True)
        if figure_grid:
            ax.set_axisbelow(True)
            ax.grid(linestyle="--", alpha=0.5)

    def draw_kmers_linear_scale(
        self,
        ax,
        data,
        colors=None,
        linewidth=2,
        xlim=(5, 100),
        ylim=(0, 1.5e8),
        ylabel="Number of distinct 23-mers",
        xlabel="23-mer coverage multiplicity",
        figure_grid=True,
        show_legend=False,
        legend_loc="upper right",
        legend_ncol=2,
    ):
        """
        Visualizes the k-mer distribution on a linear scale on the specified subplot.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes on which to draw the k-mer distribution plot.

        histo_files : list of str
            List of file paths to histogram files, each containing k-mer counts and their frequencies.

        colors : list of str or str, optional
            The colors used for the plot lines. If a list of colors is provided, it should match the number of samples.
            If a string is provided, it specifies a color palette name from seaborn. If None, distinct colors are
            automatically generated. Defaults to None.

        linewidth : float, optional
            The width of the plot lines. Default is 2.

        xlim : tuple of float, optional
            Limits for the x-axis (log scale), representing the range of tract lengths to display. Defaults to (10e4, 150e6).

        ylim : tuple of float, optional
            Limits for the y-axis, representing the range of cumulative genome fractions. Defaults to (0, 0.35).

        figure_grid : bool, optional
            Whether to show grid lines on the plot. Default is True.

        show_legend : bool, optional
            Whether to display a legend on the plot. Default is False.

        legend_loc : str, optional
            Location of the legend on the plot. Default is "upper right".

        legend_ncol : int, optional
            Number of columns in the legend. Default is 2.
        """
        # ax.spines[["right", "top"]].set_visible(False)
        for spine in ["right", "top"]:
            ax.spines[spine].set_visible(False)
        ax.set_axisbelow(True)

        if colors is None:
            colors = distinctipy.get_colors(len(data))
        else:
            if type(colors) == str:
                colors = sns.color_palette(colors, len(data))

        for i, histo_file in enumerate(data):
            kmer_counts, frequencies = [], []
            with open(histo_file, "r") as f:
                for line in f:
                    kmer_count, frequency = map(int, line.strip().split())
                    kmer_counts.append(kmer_count)
                    frequencies.append(frequency)

            sample = histo_file.split("/")[-1].split(".")[0]
            ax.plot(kmer_counts, frequencies, "-", label=sample, color=colors[i], linewidth=linewidth)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if show_legend:
            ax.legend(loc=legend_loc, ncol=legend_ncol, handlelength=0.8, frameon=True)
        if figure_grid:
            ax.set_axisbelow(True)
            ax.grid(linestyle="--", alpha=0.5)
