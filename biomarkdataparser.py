#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import sys
import argparse
import math
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as lsc
from matplotlib.colors import ListedColormap
import statsmodels.formula.api as smf
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

pd.options.mode.chained_assignment = None

def commandline():
    parser = argparse.ArgumentParser(
        prog="BioMarkqPCRdatavisualizer",
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="Write a summary from BioMark qPCR data",
        allow_abbrev=False)

    parser.add_argument(
        "-i", "--inputfile", type=str, help="qPCR data file")

    parser.add_argument(
        "-o", "--outdir", type=str, help="output directory",
        default=os.getcwd())

    parser.add_argument(
        "-v", "--validation", action="store_true", default=False,
        help="Validation format (species and strain) for samples")

    parser.add_argument(
        "-r", "--replicates", type=str, default="assay",
        choices=["assay", "sample", None],
        help="Where replicates were made (required for raw data display)")

    parser.add_argument(
        "-n", "--number", type=int, default=1,
        help="Number of replicates")

    parser.add_argument(
        "-t", "--typefilter", type=str, nargs="*", default=["all"],
        help="Filter for Sample Types, e.g. NTC, Unknown, Standard")

    parser.add_argument(
        "-s", "--samplefilter", type=str, nargs="*", default=["all"],
        help="Filter for Sample names")

    parser.add_argument(
        "-x", "--removesamples", type=str, nargs="*", default=[],
        help="Remove Samples that have input string(s) in Sample_Name")

    parser.add_argument(
        "-c", "--standardcurves", action="store_true", default=False,
        help="Draw standard curves")

    parser.add_argument(
        "-f", "--gformat", type=str, nargs="*", choices=["eps", "png"],
        default=["eps"],
        help="Picture format")

    parser.add_argument(
        "-a", "--assay_species", type=str,
        default=False,
        help="Csv file with first col Assay (Primer), second col species name")

    parser.add_argument(
        "-l", "--labelfile", type=str, default=False,
        help="A csv file with X and Y labels as columns")

    parser.add_argument(
        "--title", type=str, default=None,
        help="Title for final plots")

    parser.add_argument(
        "--transpose", action="store_true",
        help="change x and y axis, attention change columns in labelfile manually")

    parser.add_argument(
        "--italic", type=str, nargs="*", choices=["x", "y"], default=False,
        help="Select axis with species names (italic label style)")

    parser.add_argument(
        "--figsize", type=str, nargs="*", default=None,
        help="Figure size in format: Width Height")

    parser.add_argument(
        "--legendposition", type=str, nargs="*", default=False,
        help="Set legend position, format: X Y")

    parser.add_argument(
        "--deltact", type=str, nargs="*", default=False,
        help="First element preamp id, second element no preamp id")

    parser.add_argument(
        "--preamp", type=str, nargs="*", default=False,
        help="Include qualitative preamp data in results"
        "First element preamp id, second element no preamp id")

    parser.add_argument(
        "--datatoplot", type=str, nargs="*", choices=["cq", "copy", "None"],
        default=["copy"], help="Select the data to display in the plots")


    parser.add_argument(
        "--color", action="store_true", default=False,
        help="Use color in plots")

    parser.add_argument(
        "--rawdata", action="store_true", default=False,
        help="Create rawdata plots")

    parser.add_argument(
        "--reverse", action="store_true", default=False,
        help="Sort sample names in reverse order")

    parser.add_argument(
        "--outfile", type=str, help="Name of the Figure")


    return parser


class config:
    def __init__(
            self, outdir, validation, replicates,
            number, typefilter, samplefilter,
            removesamples, standardcurves, gformat,
            labelfile, title, assay_species,
            transpose, deltact, italic, figsize,
            legendposition, preamp, datatoplot,
            color, rawdata, reverse,
            outfile):
        self.outdir = outdir
        self.validation = validation
        self.replicates = replicates
        self.number = number
        self.typefilter = typefilter
        self.samplefilter = samplefilter
        self.removesamples = removesamples
        self.copythreshold = 800
        self.standardcurves = standardcurves
        self.format = gformat
        self.labelfile = labelfile
        self.title = title
        self.assay_species = assay_species
        self.transpose = transpose
        self.deltact = deltact
        self.italic = italic
        self.figsize = figsize
        self.legendposition = legendposition
        self.preamp = preamp
        self.datatoplot = datatoplot
        self.color = color
        self.rawdata = rawdata
        self.reverse = reverse
        self.outfile = outfile


def get_log(x):
    if x > 0:
        x = math.log10(x)
    return x

# read from fluidigm export file
def read_from_export(inputfile):
    def split_assayID(x):
        return str(x).split("-")[1]

    def split_sampleID(x):
        return str(x).split("-")[0]

    # Import data from Fluidigm export csv file
    # last header row 9, start of data row 10
    raw_df = pd.read_csv(inputfile, header=9)

    raw_df.columns = [
        'Chamber_ID', 'Sample_Name', 'Sample_Type', 'Sample_rConc', 'Assay_Name',
        'Assay_Type', 'Ct_value', 'Ct_Calibrated_rConc', 'Ct_Quality', 'Ct_Call',
        'Ct_Threshold', 'Tm_In_Range', 'Tm_Out_Range', 'Tm_Peak_Ratio']

    # Remove 'Not relevant' tagged assays from dataset
    raw_df = raw_df[raw_df.Assay_Name != "Not relevant"]
    # New Subset
    df = raw_df.copy()

    # Add ID to Assay name to make them unique
    df["assay_id"] = df["Chamber_ID"].apply(split_assayID)
    df["unique_assays"] = df[["Assay_Name", "assay_id"]].apply(lambda x: ' '.join(x), axis=1)

    # Add ID to Sample name to make them unique
    df["sample_id"] = df["Chamber_ID"].apply(split_sampleID)
    df["unique_samples"] = df[["Sample_Name", "sample_id"]].apply(lambda x: ' '.join(x), axis=1)

    # calculate log concentration
    df["log_Ct_Calibrated_rConc"] = df["Ct_Calibrated_rConc"].apply(get_log)
    df["log_Sample_rConc"] = df["Sample_rConc"].apply(get_log)

    return df

# Species Names extraction from Sample Names, Validation only
def get_standard_data(df):
    standarddf = df[df['Sample_Type'].isin(["Standard"])]
    subdf = standarddf[["Assay_Name", "Ct_value", "log_Sample_rConc", "Ct_Call"]]
    return subdf


def simple_reg_model(sta_data, conf):
    # create data for plotting
    datadict = {}
    assays = np.unique(sta_data[["Assay_Name"]])
    for i, assay in enumerate(assays):
        i += 1
        df = sta_data[sta_data.Assay_Name == assay]
        df.columns = ['Assay_Name', 'Cq', 'log_copy', "Ct_call"]
        excluded = df.copy()
        included = df.copy()
        included['Cq'] = np.where(included['Ct_call'] == "Flag", 999.1, included['Cq'])
        included = included[included.Cq < 999]

        excluded['Cq'] = np.where(excluded['Ct_call'] == "Pass", np.nan, excluded['Cq'])
        excluded['Cq'] = np.where(excluded['Cq'] >= 999, np.nan, excluded['Cq'])

        exvarx = excluded['log_copy'].values
        exvary = excluded['Cq'].values

        model = smf.ols('Cq ~ log_copy', data=included)
        model = model.fit()
        intercept = model.params[0]
        slope = model.params[1]
        equation = 'y = ' + str(round(slope, 2)) + 'x' + ' + ' + str(round(intercept, 2))
        rsquare = "R² = " + str(round(model.rsquared, 3))
        cq_pred = model.predict()

        datadict.update({
            assay: {
                "equation": equation, "rsquare": rsquare, "cq_pred":cq_pred,
                "xval": included['log_copy'], "yval": included['Cq'],
                "slope": slope, "intercept": intercept,
                "exvarx":  exvarx, "exvary": exvary}})

    return datadict


def plot_standard_curvedata(datadict, conf):
    ncols=6
    assays = list(datadict.keys())
    nrows = len(assays)//ncols
    assays.sort()
    if conf.title is None:
        title = "Standard calibration curves"
    else:
        title = conf.title

    # share axes
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', sharey='row', figsize=(ncols*2.5, nrows*2))

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.ylabel("Quantification cycle (Cq)", fontsize=16, labelpad=10)
    plt.xlabel("Gene copies (" + "$log_{10}$" + ")", fontsize=16, labelpad=10)
    fig.suptitle(title, x=0.5, y=1.02, fontsize=20)
    ind = 0
    # find xlim values
    x = datadict[assays[ind]]["xval"]
    exvarx = datadict[assays[ind]]["exvarx"]
    ex_min_x = min(exvarx)
    ex_max_x = max(exvarx)
    if ex_min_x < min(x):
        xmin = ex_min_x - 1
    else:
        xmin = min(x) - 1
    if ex_max_x > max(x):
        xmax = ex_max_x + 1
    else:
        xmax = max(x) + 1

    for r in range(0, nrows):
        for c in range(0, ncols):
            x = datadict[assays[ind]]["xval"]
            y = datadict[assays[ind]]["yval"]
            exvarx = datadict[assays[ind]]["exvarx"]
            exvary = datadict[assays[ind]]["exvary"]
            cq_pred = datadict[assays[ind]]["cq_pred"]
            equation = datadict[assays[ind]]["equation"]
            intercept = datadict[assays[ind]]["intercept"]
            rsquare = datadict[assays[ind]]["rsquare"]
            slope = datadict[assays[ind]]["slope"]
            ampeff = (-1 + math.pow(10, (-1/slope))) * 100
            eff = "Eff = " + str(ampeff).split(".")[0] + "%"
            ax[r, c].plot(x, y, '.', c="black", markersize=6)
            ax[r, c].plot(x, cq_pred, c="black", linewidth=0.8)
            ax[r, c].plot(exvarx, exvary, '.', c="red", markersize=6)
            ax[r, c].set_xlim(xmin, xmax)
            ax[r, c].set_ylim(0, 30)
            ax[r, c].set_title(assays[ind], fontsize=12)
            ax[r, c].annotate(equation, xy=(xmin + 0.25, 5), fontsize=10)
            ax[r, c].annotate(rsquare, xy=(xmin + 0.25, 1), fontsize=10)
            ax[r, c].annotate(eff, xy=(xmin + 0.25, 9), fontsize=10)

            # visualize cq cut-off
            cq_cut = slope * math.log10(800) + intercept
            ax[r, c].plot([xmin, xmax], [cq_cut, cq_cut], c="blue", linewidth=0.8, linestyle=":")
            ax[r, c].annotate(round(cq_cut, 1), xy=(xmax - 1.25, cq_cut - 4), fontsize=10, color="blue")

            for tick in ax[r, c].xaxis.get_major_ticks():
                tick.label.set_fontsize(10)
            for tick in ax[r, c].yaxis.get_major_ticks():
                tick.label.set_fontsize(10)

            ind += 1

    plt.subplots_adjust(
            left=None, bottom=None, right=None, top=None,
            wspace=0.05, hspace=0.2)

    if conf.outfile:
        filepath = conf.outfile
    else:
        filename = "standard_calibration_curves"
        filepath = os.path.join(conf.outdir, filename)

    for gformat in conf.format:
        graphpath = filepath + "." + gformat
        fig.savefig(graphpath, format=gformat, dpi=300, bbox_inches='tight')

    plt.show()


def summary(df, conf):
    index = "Sample_Name"
    cq_data_path = os.path.join(conf.outdir, "cq_heatmap_data.csv")
    copy_data_path = os.path.join(conf.outdir, "copy_heatmap_data.csv")

    nonzero_reps = pd.pivot_table(
                                    df, index=index, columns="Assay_Name",
                                    values="cq", dropna=False,
                                    aggfunc=np.count_nonzero)

    cqmean = pd.pivot_table(
                                    df, index=index, columns="Assay_Name",
                                    values="cq_nan", dropna=False,
                                    aggfunc=np.nanmean)

    cqstd = pd.pivot_table(
                                    df, index=index, columns="Assay_Name",
                                    values="cq_nan", dropna=False,
                                    aggfunc=np.nanstd)

    copymean = pd.pivot_table(
                                    df, index=index, columns="Assay_Name",
                                    values="log_calc_copies_nan", dropna=False,
                                    aggfunc=np.nanmean)

    copystd = pd.pivot_table(
                                    df, index=index, columns="Assay_Name",
                                    values="log_calc_copies_nan", dropna=False,
                                    aggfunc=np.nanstd)

    nonzero_reps.to_csv(os.path.join(conf.outdir, "cq_nonzero_reps.csv"))

    cqmean.loc[:, :] = np.where(nonzero_reps < conf.number - 1, np.nan, cqmean.loc[:, :])
    cqstd.loc[:, :] = np.where(nonzero_reps < conf.number - 1, np.nan, cqstd.loc[:, :])

    cqsummary = pd.concat([cqmean, cqstd], axis=1, sort=True)

    copymean.loc[:, :] = np.where(nonzero_reps < conf.number - 1, np.nan, copymean.loc[:, :])
    copystd.loc[:, :] = np.where(nonzero_reps < conf.number - 1, np.nan, copystd.loc[:, :])

    copysummary = pd.concat([copymean, copystd], axis=1, sort=True)

    annot = nonzero_reps
    annot.replace(to_replace = 0, value = "", inplace=True)
    annot.replace(to_replace = conf.number, value = "", inplace=True)
    for i in range(1, conf.number):
        annot = annot.replace(
                to_replace=i,
                value = "(" + str(i) + "/" + str(conf.number) + ")")

    if conf.reverse == True:
        cqmean.sort_index(inplace=True, ascending=False)
        cqstd.sort_index(inplace=True, ascending=False)
        cqsummary.sort_index(inplace=True, ascending=False)
        copymean.sort_index(inplace=True, ascending=False)
        copystd.sort_index(inplace=True, ascending=False)
        copysummary.sort_index(inplace=True, ascending=False)
        annot.sort_index(inplace=True, ascending=False)


    if conf.transpose == True:
        cqmean = cqmean.T
        cqstd = cqstd.T
        cqsummary = cqsummary.T
        copymean = copymean.T
        copystd = copystd.T
        copysummary = copysummary.T
        annot = annot.T

    # write results
    annot.to_csv(os.path.join(conf.outdir, "sample_n_annotation.csv"))

    cqsummary.to_csv(os.path.join(conf.outdir, "cq_data_summary.csv"))
    cqmean.to_csv(cq_data_path)

    copysummary.to_csv(os.path.join(conf.outdir, "copy_data_summary.csv"))
    copymean.to_csv(copy_data_path)

    assays = list(cqmean.columns)

    cqannot = (
            cqmean[assays].round(decimals=2).astype(str)
            + "\n" + u"\u00B1" + cqstd[assays].round(decimals=2).astype(str) + "\n" +
            annot[assays])
    cqannot = cqannot.replace(to_replace="nan", value=np.nan, regex=True)
    bool_df = cqannot.isnull() == cqmean[assays].isnull()
    cqannot = cqannot.mask(~bool_df, annot[assays])
    cqannot.to_csv(os.path.join(conf.outdir, "cq_annotation.csv"))

    copyannot = (
            copymean[assays].round(decimals=2).astype(str)
            + "\n" + u"\u00B1" + copystd[assays].round(decimals=2).astype(str) + "\n" +
            annot[assays])
    copyannot.to_csv(os.path.join(conf.outdir, "raw_copy_annotation.csv"))
    copyannot = copyannot.replace(to_replace="nan", value=np.nan, regex=True)
    copyannot = copyannot.mask(~bool_df, annot[assays])
    copyannot.to_csv(os.path.join(conf.outdir, "copy_annotation.csv"))

    return cq_data_path, copy_data_path, cqannot, copyannot, annot

def create_table(dataframe, conf, filename, mode="Ct_value"):
    if conf.replicates == "assay":
        index = "Sample_Name"
        columns = "unique_assays"
    elif conf.replicates == "sample":
        index = "unique_samples"
        columns = "Assay_Name"
    else:
        index = "Sample_Name"
        columns = "Assay_Name"

    newdf = dataframe[dataframe.Sample_Type != "Standard"]
    newdf[mode].replace(to_replace=0, value=np.nan, inplace = True)

    htdata = pd.pivot_table(newdf,index=index, columns=columns,
               values=mode, dropna=False)

    exportfile = os.path.join(conf.outdir, filename)
    htdata.to_csv(exportfile)

    return exportfile


def create_row_colors(df, conf):
    cols = df.columns.tolist()
    cols = cols[-2::] + cols[:-2]
    df = df[cols]
    df = df.drop('Sample_Name', 1)
    species = list(df.loc[:, "Species"])
    df = df.set_index("Species")
    once = np.unique(species)
    spec_pal = sns.color_palette("tab20b", once.size)

    # add better labels x-axis
    xtiklabel = list(df.columns)[1::]
    newxlabels = []
    for i, tik in enumerate(xtiklabel):
        if i % conf.number == 0:
            label = tik.split(" ")[0]
        else:
            label = ""
        newxlabels.append(label)

    spec_lut = dict(zip(map(str, once), spec_pal))

    newdf = df.set_index("Strains")
    colors = pd.Series(species).map(spec_lut)

    strains = list(newdf.index)
    for key in colors.keys():
        colors[strains[key]] = colors[key]
        del colors[key]

    return newdf, colors, newxlabels, once, spec_lut


def extract_labels(labelfile):
    df = pd.read_csv(labelfile, sep=",")
    xlabel = df.iloc[:, 0].dropna()
    ylabel = df.iloc[:, 1].dropna()
    return xlabel, ylabel


def change_assaynames_to_species(df, conf):
    specieslist = []
    adict = species_names_from_assaynames(conf.assay_species)
    if conf.transpose == True:
        assays = df.index
        for assay in assays:
            specieslist.append(adict[assay])
        return True, specieslist, "normal", "italic"
    else:
        assays = df.columns
        for assay in assays:
            specieslist.append(adict[assay])
        return specieslist, True, "italic", "normal"


def draw_heatmap(tablefile, conf, filename, final_plot=False, annotation=False, mode="cq", color=False):
    modedict = {
        "copy": {
                "vmin": 3, "vmax": 8, "cbarlabel": "Log copies / \u03BCl", "extend": "min"},
        "cq": {
                "vmin": 5, "vmax": 35, "cbarlabel": "Cq values", "extend": "neither"},
        "preamp": {
                "vmin": 3, "vmax": 8, "cbarlabel": "Log copies / \u03BCl", "extend": "min"},
        "deltact": {
                "vmin": 0, "vmax": 16, "cbarlabel": "\u0394 Cq", "extend": "neither"}}

    vmin = modedict[mode]["vmin"]
    vmax = modedict[mode]["vmax"]
    vdiff = vmax - vmin
    extend = modedict[mode]["extend"]
    cbar_label = modedict[mode]["cbarlabel"]

    if color:
        if mode == "cq":
            cmap = plt.get_cmap("RdYlBu", vdiff)
        else:
#            cmap = plt.get_cmap("RdYlBu_r", vdiff)
            cmap = ListedColormap(sns.color_palette("Blues", vdiff))
        cmap.set_under("lightgrey")
        kwscolor="black"
    else:
        cmap = lsc.from_list("", ["LightGrey", "Black"], vdiff)
        cmap.set_under("red")
        kwscolor = "snow"

    # overview all Ct's (filtered)
    sns.set_context("paper", font_scale=0.7)
    df = pd.read_csv(tablefile, sep=",")
    df = df.set_index(df.columns[0])
    if conf.figsize:
        figsize = (float(conf.figsize[0]), float(conf.figsize[1]))
    else:
        figsize = (8, 11)

    xstyle, ystyle = "normal", "normal"
    if df.shape[0] > 100 or df.shape[1] > 100:
        labelticksize = 3
        kwssize = 2
    if mode == "preamp":
        labelticksize = 8
        kwssize = 4
    else:
        labelticksize = 8
        kwssize = 6

    f, ax = plt.subplots(ncols=1, nrows=1, figsize=figsize)

    if conf.labelfile and final_plot:
        xlabels, ylabels = extract_labels(conf.labelfile)
        if conf.italic == ["x"]:
            xstyle = "italic"
        if conf.italic == ["y"]:
            ystyle = "italic"
        if conf.italic == ["x", "y"]:
            xstyle, ystyle = "italic", "italic"

    elif conf.assay_species and final_plot:
        xlabels, ylabels, xstyle, ystyle = change_assaynames_to_species(df, conf)
    else:
        xlabels = True
        ylabels = True

    if isinstance(annotation, str) == True:
        annotation = df.round(decimals=2).astype(str)

    x , y, w, h = -0.1, 0.77, 0.04, 0.2
    if conf.legendposition:
        x = float(conf.legendposition[0])
        y = float(conf.legendposition[1])
        if len(conf.legendposition) > 2:
            w = float(conf.legendposition[2])
            h = float(conf.legendposition[3])

    cax = f.add_axes([x, y, w, h])

    sns.heatmap(
            df, vmin=vmin, vmax=vmax, cmap=cmap, annot=annotation,
            annot_kws={"size": kwssize, "color": kwscolor}, fmt = '',
            cbar_kws={
                    'orientation': 'vertical', "label": cbar_label,
                    'extend': extend},
            ax=ax, cbar_ax=cax, xticklabels=xlabels, yticklabels=ylabels,
            linewidths=1, linecolor="black")

    ax.set_xticklabels(ax.get_xmajorticklabels(), rotation=90, fontsize=labelticksize, style=xstyle)
    ax.set_yticklabels(ax.get_ymajorticklabels(), rotation=0, fontsize=labelticksize, style=ystyle)
    ax.tick_params(axis='both', length=0.0, width=0.0)
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right')

    if conf.replicates == "assay":
        xlabel = "Assays n=" + str(conf.number)
        ylabel = "Samples"
    elif conf.replicates == "sample":
        xlabel = "Assays"
        ylabel = "Samples n=" + str(conf.number)
    else:
        xlabel = "Assays"
        ylabel = "Samples"

    labels = [xlabel, ylabel]

    if conf.transpose == True:
        ax.set_xlabel(labels[1], fontsize=10)
        ax.set_ylabel(labels[0], fontsize=10)
    else:
        ax.set_xlabel(labels[0], fontsize=10)
        ax.set_ylabel(labels[1], fontsize=10)

    ax.set_title(conf.title, fontsize=12)

    if conf.outfile:
        filepath = conf.outfile
    else:
        filepath = os.path.join(conf.outdir, "figures", filename)

    for gformat in conf.format:
        graphpath = filepath + "." + gformat
        f.savefig(graphpath, format=gformat, dpi=600, bbox_inches='tight')
    plt.show()


def filterdata(df, conf, datadict):
    def copy_number_threshold(conc):
        if int(conc) > conf.copythreshold or int(conc) == -1:
            return "Pass"
        else:
            return "Fail"

    if conf.samplefilter != ["all"] and conf.typefilter != ["all"]:
        print("Include Sample_Type and Sample_Names containing ")
        print(conf.typefilter, conf.samplefilter)
        df1 = df[df['Sample_Type'].isin(conf.typefilter)]
        fildf = df1[df1.Sample_Name.str.contains('|'.join(conf.samplefilter))]
    elif conf.samplefilter != ["all"]:
        print("Include Sample_Names containing ")
        print(conf.samplefilter)
        fildf = df[df.Sample_Name.str.contains('|'.join(conf.samplefilter))]
    elif conf.typefilter != ["all"]:
        print("Include only Sample_Type ")
        print(conf.typefilter)
        fildf = df[df['Sample_Type'].isin(conf.typefilter)]
    else:
        print("No filter")
        fildf = df

    if conf.removesamples != []:
        print("Remove Sample_Names containing ")
        print(conf.removesamples)
        fildf = fildf[~fildf.Sample_Name.str.contains('|'.join(conf.removesamples))]

    # Flag filter from Fluidigm (Quality Threshold or Multiple Peaks -> Threshold 0.8, Peak sensitivity 3)
    fildf.loc[:, 'Ct_value'] = np.where(fildf.loc[:, 'Ct_Call'] == "Flag", 999.1, fildf.loc[:, 'Ct_value'])
    # Calculate copy numbers for all Cq values below 999
    fildf = copynumber_calc(fildf, datadict)
    # Copy cutoff, adds column with "Fail" if estimated copy number below 800
    fildf["copy_cutoff"] = fildf["calc_copies"].apply(copy_number_threshold)
    fildf['calc_copies'] = np.where(fildf['copy_cutoff'] == "Fail", 0, fildf['calc_copies'])
    fildf['calc_copies_nan'] = fildf['calc_copies'].replace(0, np.nan)
    fildf['log_calc_copies'] = np.where(fildf['copy_cutoff'] == "Fail", 0, fildf['log_calc_copies'])
    fildf['log_calc_copies_nan'] = fildf['log_calc_copies'].replace(0, np.nan)
    # Filter results with lower concentrations than theoretically possible
    fildf['Ct_value'] = np.where(fildf['copy_cutoff'] == "Fail", 999.2, fildf['Ct_value'])
    # Replace non-valid ct values
    fildf["cq"] = fildf["Ct_value"].replace([999.0, 999.1, 999.2], [0, 0, 0,])
    fildf["cq_nan"] = fildf["Ct_value"].replace([999.0, 999.1, 999.2], [np.nan, np.nan, np.nan])

    return fildf


def calc_copies_per_ul(x, assaydict):
    assay = x.Assay_Name
    ct = x.Ct_value
    slope = assaydict[assay]['slope']
    intercept = assaydict[assay]['intercept']
    if ct < 999:
        copies = math.pow(10, (ct - intercept)/ slope)
    else:
        copies = 0
    return copies


def copynumber_calc(df, datadict):
    df["calc_copies"] = df.apply(calc_copies_per_ul, args=(datadict,), axis=1)
    df["log_calc_copies"] =  df["calc_copies"].apply(get_log)
    return df


def create_directory(path_to_dir):
    if not os.path.isdir(path_to_dir):
        try:
            os.makedirs(path_to_dir)
        except OSError:
            if not os.path.isdir(path_to_dir):
                raise


def species_names_from_assaynames(inputfile):
    adict = {}
    with open(inputfile) as f:
        reader = csv.reader(f)
        for row in reader:
            adict.update({row[0]: row[1]})
    return adict


def plot_cq_diff_preamp(cq_file, conf):
    deltactpath = os.path.join(conf.outdir, "preamp_deltact_data.csv")
    df = pd.read_csv(cq_file)
    df = df.set_index(df.columns[0])
    preamp = []
    no_preamp = []
    unknown = []
    if conf.transpose == True:
        samplelist = list(df.columns.values)
    else:
        samplelist = list(df.index.values)

    for item in samplelist:
        if conf.deltact[0] in item:
            preamp.append(item)
        elif conf.deltact[1] in item:
            no_preamp.append(item)
        else:
            unknown.append(item)

    if conf.transpose == True:
        preampdf = df[preamp]
        noampdf = df[no_preamp]
    else:
        preampdf = df.loc[preamp, :]
        noampdf = df.loc[no_preamp, :]

    if conf.transpose == True:
        namelist = list(noampdf.columns)
    else:
        namelist = list(noampdf.index)
    newnamelist = []
    for item in namelist:
        name = item.split(conf.deltact[1])[0]
        newnamelist.append(name)

    if conf.transpose == True:
        noampdf.columns = newnamelist
        preampdf.columns = newnamelist
    else:
        noampdf.index = newnamelist
        preampdf.index = newnamelist

    unknowndf = df[unknown]

    if conf.transpose == True:
        diffdf = noampdf.sub(preampdf, axis='columns')
    else:
        diffdf = noampdf.sub(preampdf, axis='rows')

    diffdf.append(unknowndf, sort=True)

    diffdf.to_csv(deltactpath)

    return deltactpath


def include_qualitative_preampdata(copy_file, conf, annot):
    df = pd.read_csv(copy_file)
    df = df.set_index(df.columns[0])
    preamp = []
    no_preamp = []
    unknown = []
    preamp_id = conf.preamp[0]

    if conf.transpose == True:
        samplelist = list(df.columns.values)
    else:
        samplelist = list(df.index.values)

    if len(conf.preamp) == 2:
        for item in samplelist:
            if preamp_id in item:
                preamp.append(item)
            elif conf.preamp[1] in item:
                no_preamp.append(item)
            else:
                unknown.append(item)

    else:
        for item in samplelist:
            if preamp_id in item:
                preamp.append(item)
            else:
                no_preamp.append(item)

    if conf.transpose == True:
        preampdf = df[preamp]
        noampdf = df[no_preamp]
        noamp_annot = annot[no_preamp]
        preamp_annot = annot[preamp]
    else:
        preampdf = df.loc[preamp, :]
        noampdf = df.loc[no_preamp, :]
        noamp_annot = annot.loc[no_preamp, :]
        preamp_annot = annot.loc[preamp, :]

    if conf.transpose == True:
        namelist = list(noampdf.columns)
    else:
        namelist = list(noampdf.index)

    newnamelist = []
    for item in namelist:
        if len(conf.preamp) == 2:
            name = item.split(conf.preamp[1])[0]
        else:
            name = item
        newnamelist.append(name)

    if conf.transpose == True:
        noampdf.columns = newnamelist
        preampdf.columns = newnamelist
        noamp_annot.columns = newnamelist
        preamp_annot.columns = newnamelist
    else:
        noampdf.index = newnamelist
        preampdf.index = newnamelist
        noamp_annot.index = newnamelist
        preamp_annot.index = newnamelist

    bool_df = noampdf.isnull() == preampdf.isnull()
    quali_mask = ~bool_df
    combined = noampdf.mask(quali_mask, 1.0)
    combi_annot = noamp_annot.mask(quali_mask, "preamp\n" + preamp_annot)
    preampcopypath = os.path.join(conf.outdir, "preamp_data.csv")
    combined.to_csv(preampcopypath)
    combi_annot.to_csv(os.path.join(conf.outdir, "preamp_annotation.csv"))

    return preampcopypath, combi_annot


def main():
    parser = commandline()
    args = parser.parse_args()

    if os.path.abspath(args.inputfile):
        inputfile = args.inputfile
    else:
        inputfile = os.path.join(os.getcwd(), args.inputfile)

    conf = config(
            args.outdir, args.validation, args.replicates, args.number,
            args.typefilter, args.samplefilter, args.removesamples,
            args.standardcurves, args.gformat, args.labelfile, args.title,
            args.assay_species, args.transpose, args.deltact, args.italic,
            args.figsize, args.legendposition, args.preamp, args.datatoplot,
            args.color, args.rawdata, args.reverse, args.outfile)

    create_directory(conf.outdir)
    df = read_from_export(inputfile)

    sta_data = get_standard_data(df)
    datadict = simple_reg_model(sta_data, conf)
    if conf.standardcurves:
        plot_standard_curvedata(datadict, conf)
    else:

        tmpfile = os.path.join(conf.outdir, "rawdata.csv")
        df.to_csv(tmpfile)

        filterdf = filterdata(df, conf, datadict)
        filterfile = os.path.join(conf.outdir, "filterdata.csv")
        filterdf.to_csv(filterfile)

        cq_data_path, copy_data_path, cqannot, copyannot, annot = summary(filterdf, conf)

        raw_cqs = create_table(filterdf, conf, "raw_cq_heatmap_data.csv")
        raw_copy = create_table(filterdf, conf, "raw_copy_heatmap_data.csv", "log_calc_copies")

        if conf.validation == False:
            create_directory(os.path.join(conf.outdir, "figures"))

            if "None" in conf.datatoplot:
                sys.exit()
                print("Exit")

            elif conf.deltact:
                deltactpath = plot_cq_diff_preamp(cq_data_path, conf)
                draw_heatmap(deltactpath, conf, "preamp_deltacts_heatmap", True, "custom", "deltact", conf.color)
                return

            elif conf.preamp:
                preampdatapath, newannot = include_qualitative_preampdata(copy_data_path, conf, copyannot)
                draw_heatmap(preampdatapath, conf, "qPCR_copies_preamp", True, newannot, "preamp", conf.color)
                return

            else:
                # draw heatmaps
                for datatype in conf.datatoplot:
                    if datatype == "cq":
                        annotation = cqannot
                        datapath = cq_data_path
                        rawpath = raw_cqs
                    else:
                        annotation = copyannot
                        datapath = copy_data_path
                        rawpath = raw_copy

                    if conf.rawdata == True:
                        draw_heatmap(rawpath, conf, "raw_qPCR_heatmap_" + datatype, False, "custom", datatype, conf.color)

                    draw_heatmap(datapath, conf, "qPCR_heatmap_" + datatype, True, annotation, datatype, conf.color)


if __name__ == "__main__":
    main()
