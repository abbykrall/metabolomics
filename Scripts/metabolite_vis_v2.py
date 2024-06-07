import math
import datetime
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mp
import json
import matplotlib.colors as mc
import matplotlib.cm as cm
from scipy import stats
import matplotlib.pyplot as plt
import Tkinter as tk
import tkFont as font
from tkFileDialog import askopenfilename
from pdfrw import PdfReader, PdfDict
from pdfrw.buildxobj import pagexobj
from pdfrw.toreportlab import makerl

from matplotlib.colors import Normalize
from matplotlib.collections import LineCollection

import io

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.enums import *
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, Table, Flowable, KeepTogether, TableStyle as TS, HRFlowable
from reportlab.lib.styles import ParagraphStyle as PS
from reportlab.lib.units import inch

import os
from matplotlib import font_manager as fm, rcParams

# --------------------------------------------------------------------------------------------------------

def validate_filename(*args):
    if ".pdf" in pdf_filename.get():
        tk.messagebox.showerror("Invalid Filename", "The filename should not include .pdf")
        pdf_filename.set(pdf_filename.get().replace(".pdf", ""))

# --------------------------------------------------------------------------------------------------------

pyenvdir = "C:\\Users\\AKrall-admin\\Desktop\\Py_Env\\christofk-lab"
save_data_filepath = os.path.join(pyenvdir, 'previous_data.json')

desired_width = 320
pd.set_option('display.width', desired_width)
pd.set_option('display.max_columns', 10)

open_sans_font = fm.FontProperties(fname="C:\\Users\\AKrall-admin\\Desktop\\Py_Env\\christofk-lab\\python_report\\fonts\\Open_Sans\\OpenSans-Regular.ttf")

def color_ramp(colors, n):
    colmap = mc.LinearSegmentedColormap.from_list("col", colors=colors, N=n)
    color_arr = []
    for i in range(0, n):
        color_arr.append(colmap((i / n)))
    return color_arr

def getCompoundDataFromData(data, group_num, cmpd_str):
    sample_cols = np.where(data[data['Compound'] == "group"].values == group_num)[1]
    sub_data = data[data['Compound'] == cmpd_str]
    sub_data.sort_values(by=sub_data.columns[2])
    sub_data = sub_data.iloc[:, sample_cols]
    return sub_data

def getNumTotalGroups(data):
    group_data = data[data['Compound'] == "group"]
    if group_data.empty:
        raise ValueError("No group data found.")
    max_group_number = group_data.iloc[:, 1:].max(axis=1)
    try:
        num_groups = int(max_group_number.values[0])
    except (ValueError, IndexError) as e:
        raise ValueError("Could not determine the number of groups from the data. Please check your input data.")
    return num_groups

def getGroupedAverageDataFromData(data, num_groups, cmpd_str, group_names={}):
    gp_data = pd.DataFrame()
    gp_std = pd.DataFrame()
    for group_id in range(1, num_groups + 1):
        group_iso_dist = getCompoundDataFromData(data, group_id, cmpd_str)
        group_mean = group_iso_dist.mean(axis=1).multiply(100)
        gp_data = gp_data.append(group_mean, ignore_index=True)
        group_stdev = group_iso_dist.std(axis=1).multiply(100)
        gp_std = gp_std.append(group_stdev, ignore_index=True)
    default_group_names = {}
    for i in range(1, num_groups + 1):
        default_group_names[i - 1] = "group " + repr(i)
    for group_i, group_name in group_names.items():
        default_group_names[group_i - 1] = group_name
    gp_data = gp_data.rename(default_group_names)
    gp_std = gp_std.rename(default_group_names)
    for i in range(0, len(gp_data.columns)):
        default_group_names[gp_data.columns[i]] = "M" + repr(i)
    gp_data = gp_data.rename(default_group_names, axis="columns")
    gp_std = gp_std.rename(default_group_names, axis="columns")
    return (gp_data, gp_std)

def getFractionalContributionDataFromIsoData(data, num_groups, cmpd_str, group_names={}):
    gp_data = pd.DataFrame()
    def frac(val):
        return 1 - val
    for group_id in range(1, num_groups + 1):
        group_iso_dist = getCompoundDataFromData(data, group_id, cmpd_str)
        if group_iso_dist.empty: continue
        new_row = group_iso_dist.iloc[0, :].apply(frac).multiply(100)
        mappings = {}
        for i in range(0, new_row.size):
            mappings[new_row.index[i]] = i;
        gp_data = gp_data.append(new_row.rename(mappings), ignore_index=True)
        gp_data = gp_data.replace(np.nan, 0)
    default_group_names = {}
    for i in range(1, num_groups + 1):
        default_group_names[i - 1] = "group " + repr(i)
    for group_i, group_name in group_names.items():
        default_group_names[group_i - 1] = group_name
    return gp_data.rename(default_group_names).transpose()

def getPoolSizeFromData(data, num_groups, cmpd_str, group_names={}):
    gp_data = pd.DataFrame()
    for group_id in range(1, num_groups + 1):
        new_row = getCompoundDataFromData(data, group_id, cmpd_str)
        if new_row.empty: continue
        new_row = new_row.iloc[0]
        mappings = {}
        for i in range(0, new_row.size):
            mappings[new_row.index[i]] = i;
        gp_data = gp_data.append(new_row.rename(mappings), ignore_index=True)
    default_group_names = {}
    for i in range(1, num_groups + 1):
        default_group_names[i - 1] = "group " + repr(i)
    for group_i, group_name in group_names.items():
        default_group_names[group_i - 1] = group_name
    return gp_data.rename(default_group_names).transpose()

def heatMapDatafromPoolData(data, num_groups, class_name):
    mets = ANALYTE_META_DATA[ANALYTE_META_DATA['pathway_class'] == class_name]['analyte_name']
    sample_sizes = []
    for group_id in range(1, num_groups + 1):
        sample_sizes.append(len(np.where(data[data['Compound'] == "group"].values == group_id)[1]))
    gp_data = pd.DataFrame(np.zeros((mets.size, sum(sample_sizes))))
    col_index = 0
    for i in range(0, len(sample_sizes)):
        for x in range(1, sample_sizes[i] + 1):
            gp_data = gp_data.rename(columns={col_index: str(i + 1) + "_" + str(x)})
            col_index += 1
    index_count = -1
    for met in mets:
        if ("<div>" in met):
            met = "[" + met.replace("<div>", "/") + "]"
        index_count += 1
        gp_data = gp_data.rename({index_count: met})
        met_data = data[data['Compound'] == met]
        if met_data.empty: continue
        met_data = met_data.iloc[:, 1:].astype(float)
        sample_cols = np.where(data[data['Compound'] == "group"].values == 1)[1] - 1
        if len(np.where(met_data.iloc[:, sample_cols] == 0)[0]) > 0: continue
        ctrl_grp_average = np.mean(met_data.iloc[:, sample_cols].values)
        if ctrl_grp_average == 0: continue
        met_data = met_data.divide(ctrl_grp_average)
        met_data = np.log2(met_data.replace(0, np.nan))
        insert_index = 0
        for group_id in range(1, num_groups + 1):
            sample_cols = np.where(data[data['Compound'] == "group"].values == group_id)[1] - 1
            group_data = met_data.iloc[0, sample_cols].values
            for val in group_data:
                gp_data.iloc[index_count, insert_index] = val
                insert_index += 1
    return gp_data.replace(0, np.nan)

HEATMAP_TILE_W = 0.12
HEATMAP_TILE_H = 0.12

def heatMapTableFromData(data, group_id, cmap, vmin=-2, vmax=2):
    norm = Normalize(vmin=vmin, vmax=vmax)
    table_data = []
    table_style = [('GRID', (0, 0), (-1, -1), 1, "white"), ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'), ('ALIGN', (0, 0), (-1, -1), 'CENTER')]
    for r in range(0, len(data.index)):
        table_row = []
        new_col_index = 0
        for c in range(0, len(data.columns)):
            if data.columns[c].split("_")[0] != str(group_id): continue
            datum = data.iloc
            datum = data.iloc[r, c]
            if np.isnan(datum):
                table_row.append('-')
                cell_color = colors.Color(red=1, green=1, blue=1)
            else:
                table_row.append('')
                col = cmap(norm(datum))
                cell_color = colors.Color(red=col[0], green=col[1], blue=col[2])
            table_style.append(('BACKGROUND', (new_col_index, r), (new_col_index, r), cell_color))
            new_col_index += 1
        table_data.append(table_row)
    table = Table(table_data, colWidths=HEATMAP_TILE_W * inch, rowHeights=HEATMAP_TILE_H * inch)
    table.setStyle(TS(table_style))
    return table


class PdfImage(Flowable):
    def __init__(self, fig, width=200, height=200, bgcolor="white"):
        self.img_width = width
        self.img_height = height
        idata = io.BytesIO()
        fig.savefig(idata, format='pdf', bbox_inches='tight', facecolor=bgcolor)
        idata.seek(0)
        self.img_data = self.form_xo_reader(idata)

    def form_xo_reader(self, imgdata):
        page, = PdfReader(imgdata).pages
        return pagexobj(page)

    def wrap(self, width, height):
        return self.img_width, self.img_height

    def drawOn(self, canv, x, y, _sW=0):
        if _sW > 0 and hasattr(self, 'hAlign'):
            a = self.hAlign
            if a in ('CENTER', 'CENTRE', TA_CENTER):
                x += 0.5 * _sW
            elif a in ('RIGHT', TA_RIGHT):
                x += _sW
            elif a not in ('LEFT', TA_LEFT):
                raise ValueError("Bad hAlign value " + str(a))
        canv.saveState()
        img = self.img_data
        if isinstance(img, PdfDict):
            xscale = self.img_width / img.BBox[2]
            yscale = self.img_height / img.BBox[3]
            canv.translate(x, y)
            canv.scale(xscale, yscale)
            canv.doForm(makerl(canv, img))
        else:
            canv.drawImage(img, x, y, self.img_width, self.img_height)
        canv.restoreState()

# Usual sheet of metabolites for the report
ANALYTE_META_DATA = pd.read_csv("C:\\Users\\AKrall-admin\\Desktop\\Py_Env\\christofk-lab\\metabolite_classes.csv", dtype={0: 'str', 1: 'str', 2: 'str'})

STYLE_TITLE = PS(name='Heading1', fontSize=28, leading=16, spaceAfter=0.2 * inch)
STYLE_TITLE_MIN = PS(name='Heading2', fontSize=22, leading=16, spaceAfter=0.2 * inch)
STYLE_H1 = PS(name='Heading1', fontSize=18, leading=16)
STYLE_H4 = PS(name='Heading4', fontSize=10, leading=16)
STYLE_MET_NAME = PS(name="Metname", fontSize=10, leading=10, spaceAfter=0.1 * inch)
STYLE_MET_RATIO_NAME = PS(name="Metname", fontSize=10, leading=10, spaceAfter=0.1 * inch)

HR_LINE = HRFlowable(width="100%", thickness=2, color="black", spaceBefore=1, spaceAfter=1, hAlign='CENTER', vAlign='CENTER', dash=None)
HR_LINE_MIN = HRFlowable(width="100%", thickness=0.8, color="black", spaceBefore=1, spaceAfter=1, hAlign='CENTER', vAlign='CENTER', dash=None)
VGAP_04 = Spacer(0, 0.4 * inch)
VGAP_02 = Spacer(0, 0.2 * inch)

GLOBAL_ABBRV_MAP = {"phosphogylcerate": "PG", "N-acetyl-glucosamine": "NAG", "phosphate": "P", "phosphoethanolamine": "PE"}

def ABBRV(word):
    for old, new in GLOBAL_ABBRV_MAP.items():
        word = word.replace(old, new)
    return word

class MetaboliteReport:
    def __init__(self, title, iso_df, ps_df, num_groups=-1, group_names={}, show_strct=True, show_iso=True, show_ps=True):
        self.title = title
        self.source_data_iso = iso_df
        self.source_data_pool = ps_df
        self.num_groups = num_groups
        if num_groups == -1: self.num_groups = getNumTotalGroups(ps_df)
        self.show_strct = show_strct
        self.show_iso = show_iso
        self.show_ps = show_ps
        self.title_params = {}
        self.group_names = group_names
        self.sample_sizes = []
        for group_id in range(1, self.num_groups + 1):
            self.sample_sizes.append(len(np.where(ps_df[ps_df['Compound'] == "group"].values == group_id)[1]))
        title = datetime.datetime.now().strftime("%Y_%m_%d") + "_" + title
        self.doc = SimpleDocTemplate(title + ".pdf", pagesize=letter, rightMargin=30, leftMargin=30, topMargin=30, bottomMargin=25)
        self.page_elements = []

    def add_element(self, element):
        if element is not None:
            self.page_elements.append(element)

    def configureTitlePage(self, title, project, exp_title, date, leader, exp_summary, normalization=[]):
        self.title_params["title"] = title
        self.title_params["date"] = date
        self.title_params["exp_title"] = exp_title
        self.title_params["project"] = project
        self.title_params["leader"] = leader
        self.title_params["exp_summary"] = exp_summary
        self.title_params["normalization"] = normalization

    def isotopomer_plot_image_for_data(self, data, error):
        fig = plt.figure(frameon=False)
        ax = fig.add_subplot(111)
        data.plot(ax=ax, kind="bar", stacked=True, rot=0, width=0.9, edgecolor="white", linewidth=1, yerr=error, capsize=2,
                  color=(
                      "#d9d9d9",
                      "#a6d2de",
                      "#5a8e9c",
                      "#e38e59",
                      "#2d67a8",
                      "#e0d699",
                      "#d496a9",
                      "#a5c797",
                      "#ceb4db",
                      "#cc4b64",
                      "#6aa695",
                      "#bf4a26",
                      "#bf7526",
                      "#c49f23",
                      "#bdc42f",
                      "#96c43f",
                      "#579e47",
                      "#308c75"
                  ))
        ax.set_ylim((0, 100))
        ax.set_ylabel("Composition of pool (%)")
        ax.yaxis.labelpad = 10
        ax.set_xlabel("Groups")
        ax.xaxis.labelpad = 10
        ax.tick_params(axis='x', pad=7, bottom=False)
        ax.legend(loc='lower right', frameon=False, bbox_to_anchor=(1, 0.1), bbox_transform=plt.gcf().transFigure)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.set_title('Mass Distribution (mean)', fontproperties=open_sans_font, y=1.05)
        ax.title.set_fontsize(17)
        size = fig.get_size_inches()
        ratio = size[1] / size[0]
        width = 2
        image = PdfImage(fig, width * inch, width * ratio * inch)
        plt.close(fig)
        return image

    def box_plot_image_for_data(self, data, title, ymax=100, yname=""):
        fig = plt.figure(frameon=False)
        ax = fig.add_subplot(111)
        ymin = data.replace(np.nan, 0).values.min()
        if ymin > 0: ymin = 0
        data = data.replace(0, np.nan)
        if ymax < 0: ymax = 0.2 * abs(ymax)
        data.plot(ax=ax, kind="box", color=dict(boxes='black', whiskers='black', caps='black'),
                  medianprops=dict(linewidth=0), showfliers=False,
                  meanprops=dict(linestyle='-', linewidth=2, color='black'), showmeans=True, meanline=True)
        ax.set_ylim((1.5 * ymin, ymax))
        ax.yaxis.labelpad = 10
        ax.set_xlabel("Groups")
        ax.set_ylabel(yname)
        ax.xaxis.labelpad = 10
        ax.tick_params(axis='x', pad=7, bottom=False)
        l1 = [(0, 0), (20, 0)]
        lc = LineCollection([l1, ], color=["black", "red", "red"], lw=1)
        ax.add_collection(lc)
        data_sets = []
        colors = ["red", "green", "blue"]
        for i in range(0, self.num_groups):
            y = data.iloc[:, i]
            data_sets.append(y.values.tolist())
            num_total = self.sample_sizes[i]
            valid_values = len(y) - y.isnull().sum()
            if valid_values / num_total < 1:
                plt.text(i + 1, 0, "(" + str(valid_values) + "/" + str(num_total) + ")", horizontalalignment='center',
                         fontsize=8, color="red", bbox={'facecolor': 'white', 'pad': 2, 'edgecolor': 'none'})
            y = y.dropna()
            if len(y) > 0:
                x = np.random.normal(i + 1, 0.04, size=len(y))
                plt.plot(x, y, color=colors[i % len(colors)], marker='.', linestyle='')

        border_color = 'black'
        border_width = 1

        f, p_val = stats.f_oneway(*data_sets)
        if p_val < 0.001:
            title += " ***"
            border_color = "red"
            border_width = 2
        elif p_val < 0.01:
            title += " **"
            border_color = "red"
            border_width = 2
        elif p_val < 0.05:
            title += " *"
            border_color = "red"
            border_width = 2

        ax.spines['bottom'].set_color(border_color)
        ax.spines['bottom'].set_linewidth(border_width)
        ax.spines['top'].set_color(border_color)
        ax.spines['top'].set_linewidth(border_width)
        ax.spines['left'].set_color(border_color)
        ax.spines['left'].set_linewidth(border_width)
        ax.spines['right'].set_color(border_color)
        ax.spines['right'].set_linewidth(border_width)
        ax.set_title(title, fontproperties=open_sans_font, y=1.05)
        ax.title.set_fontsize(17)

        size = fig.get_size_inches()
        ratio = size[1] / size[0]
        width = 2
        image = PdfImage(fig, width * inch, width * ratio * inch)
        plt.close(fig)

        return image

    def heatmap_for_data(self, data, title):
        if data.empty or not np.any(~np.isnan(data.values)):
            return None

        cmap = mc.LinearSegmentedColormap.from_list("col", colors=("#1f3e82", "white", "#e22d45"), N=100)
        num_groups = int(list(data)[-1].split("_")[0])

        if not data.empty and np.any(~np.isnan(data.values)):
            min_value = math.floor(min(-2, np.nanmin(data.values)))
            max_value = math.ceil(max(2, np.nanmax(data.values)))

            if max_value > abs(min_value):
                min_value = -1 * max_value
            else:
                max_value = -1 * min_value
        else:
            min_value = -2
            max_value = 2

        outer_table = []
        row = [""]

        for i in range(1, num_groups + 1):
            if i in self.group_names:
                row.append(self.group_names[i])
            else:
                row.append("group " + str(i))

        row.append("log2 foldchange")
        outer_table.append(row)

        plt.rcParams["font.weight"] = "bold"
        fig = plt.figure(frameon=False, figsize=(2, 10))
        ax = fig.add_subplot(111)
        cb = mp.colorbar.ColorbarBase(ax=ax, cmap=cmap, norm=Normalize(vmin=min_value, vmax=max_value), orientation='vertical')
        cb.set_ticks([min_value, 0, max_value])
        cb.set_ticklabels([min_value, " 0", " " + str(max_value)])
        cb.outline.set_linewidth(5)
        ax.tick_params(axis='both', which='major', labelsize=50)
        ax.tick_params(axis='y', pad=60, bottom=False)
        size = fig.get_size_inches()
        ratio = size[1] / size[0]
        width = 0.2
        image = PdfImage(fig, 1.65 * width * inch, width * ratio * inch)
        plt.close(fig)

        met_names_table = []
        table_style = [('ALIGN', (0, 0), (-1, -1), 'RIGHT'),
                    ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                    ('FONTSIZE', (0, 0), (-1, -1), 8)]
        for r in range(0, len(data.index)):
            name = data.index[r]
            pool_data = getPoolSizeFromData(self.source_data_pool, self.num_groups, name, self.group_names)
            data_sets = []

            if not pool_data.empty:
                pool_data = pool_data.replace(0, np.nan)
                for i in range(0, self.num_groups):
                    y = pool_data.iloc[:, i]
                    data_sets.append(y.values.tolist())

                f, p_val = stats.f_oneway(*data_sets)
                if p_val < 0.001:
                    name += " ***"
                elif p_val < 0.01:
                    name += " ** "
                elif p_val < 0.05:
                    name += " *  "
                else:
                    name += "    "

                if p_val < 0.05:
                    table_style.append(('TEXTCOLOR', (0, r), (0, r), colors.red))

            else:
                name += "    "

            met_names_table.append([ABBRV(name)])

        # Skip table creation if met_names_table is empty
        if not met_names_table:
            print("met_names_table is empty. Skipping table creation.")
            return None

        met_table = Table(met_names_table, rowHeights=HEATMAP_TILE_H * inch)
        met_table.setStyle(TS(table_style))

        row = [met_table]

        for i in range(1, num_groups + 1):
            row.append(heatMapTableFromData(data, i, cmap, vmin=min_value, vmax=max_value))

        row.append(image)
        outer_table.append(row)

        table_style = [('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                    ('ALIGN', (0, 0), (0, 1), 'RIGHT'),
                    ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                    ('FONTSIZE', (0, 0), (-1, -1), 8)]

        space_remaining = 7.3 * inch
        col_widths = [(1) * inch]
        space_remaining -= col_widths[0]
        num_samples = 0
        old_group_id = ""

        for c in range(0, len(data.columns)):
            group_id = str(data.columns[c].split("_")[0])
            num_samples += 1

            if old_group_id != "" and group_id != old_group_id:
                # Heatmap column spacing
                w = (HEATMAP_TILE_W * num_samples + 0.15) * inch
                col_widths.append(w)
                space_remaining -= w
                num_samples = 1

            old_group_id = group_id

        w = (HEATMAP_TILE_W * num_samples + 0.05) * inch
        col_widths.append(w)
        space_remaining -= w

        col_widths.append(space_remaining)

        table = Table(outer_table, colWidths=col_widths, rowSplitRange=(0, 0))
        table.setStyle(TS(table_style))

        return table

    def table_of_plots_for_metabolite(self, met_name):
        table_stct = []
        table_row = []

        if self.show_iso:
            iso_avg, iso_error = getGroupedAverageDataFromData(self.source_data_iso, self.num_groups, met_name, self.group_names)
            if not iso_avg.empty:
                table_row.append(self.isotopomer_plot_image_for_data(iso_avg, iso_error))

            frac_data = getFractionalContributionDataFromIsoData(self.source_data_iso, self.num_groups, met_name, self.group_names)
            if not frac_data.empty:
                table_row.append(self.box_plot_image_for_data(frac_data, title="Fractional Contribution", ymax=105, yname="Fractional Contribution (%)"))

        pool_data = getPoolSizeFromData(self.source_data_pool, self.num_groups, met_name, self.group_names)
        if not pool_data.empty:
            table_row.append(self.box_plot_image_for_data(pool_data, title="Pool Size", ymax=(1.2 * (np.nanmax(pool_data.values))), yname="Normalized Intensity"))

        if table_row:
            table_stct.append(table_row)
            plot_row = Table(table_stct, 2.8 * inch)

            if self.show_iso:
                plot_row.setStyle(TS([('ALIGN', (0, 0), (-1, -1), 'CENTER')]))

            return plot_row

        return None

    def append_new_page_title(self, title):
        self.page_elements.append(Paragraph('<b>' + title + '</b>', STYLE_H1))
        self.page_elements.append(Spacer(0, 0.2 * inch))
        self.page_elements.append(HR_LINE)
        self.page_elements.append(VGAP_04)


    def build_with_metabolite_classes(self, met_classes):
        if "title" in self.title_params:
            self.add_element(Paragraph('<b>' + self.title_params["title"] + '</b>', STYLE_TITLE))
            self.add_element(VGAP_04)
            self.add_element(HR_LINE)
            self.add_element(Spacer(0, 0.7 * inch))

        if "project" in self.title_params:
            self.add_element(Paragraph('<b>Project:  ' + self.title_params["project"] + '</b>', STYLE_H1))
            self.add_element(Spacer(0, 1 * inch))

        if "leader" in self.title_params:
            self.add_element(Paragraph('<b>' + self.title_params["leader"] + '</b>', STYLE_H4))
            self.add_element(Spacer(0, 0.15 * inch))

        if "exp_title" in self.title_params:
            self.add_element(Paragraph('<b>' + self.title_params["exp_title"] + '</b>', STYLE_H4))
            self.add_element(Spacer(0, 0.15 * inch))

        if "date" in self.title_params:
            self.add_element(Paragraph('<b>Date:  ' + self.title_params["date"] + '</b>', STYLE_H4))
            self.add_element(VGAP_04)

        if "exp_summary" in self.title_params:
            self.add_element(Paragraph(self.title_params["exp_summary"], STYLE_H4))

        self.add_element(Spacer(0, 1 * inch))

        self.add_element(Paragraph(
            "READ ME: All significance tests were performed on normalized data with a one-way ANOVA. * p-value < 0.05, ** p-value < 0.01, *** p-value < 0.001.",
            STYLE_H4))

        if "normalization" in self.title_params:
            norms = self.title_params["normalization"]

            if len(norms) > 0:
                self.add_element(Spacer(0, 1.5 * inch))
                self.add_element(Paragraph('<b>Pool sizes were normalized to field(s):</b>', STYLE_H4))

                for n in norms:
                    self.add_element(Paragraph(n, STYLE_H4))

        print("___________________FINISHED_TITLE_PAGE_________________________")

        self.add_element(PageBreak())

        self.add_element(Paragraph('<b> Pool Sizes </b>', STYLE_TITLE_MIN))
        self.add_element(VGAP_04)

        for met_class in met_classes:
            print("_________________STARTING_HEATMAP_FOR > " + met_class)

            hmap = self.heatmap_for_data(heatMapDatafromPoolData(self.source_data_pool, self.num_groups, met_class),
                                        title=r'$\log_2(fold-change-poolsize)$')

            if hmap is not None:
                group = KeepTogether([Paragraph('<b>' + met_class + '</b>', STYLE_H4),
                                    HR_LINE_MIN,
                                    VGAP_02,
                                    hmap,
                                    VGAP_02])

                print("_________________BUILT_HEATMAP_FOR > " + met_class)
                self.add_element(group)

        self.add_element(PageBreak())

        for met_class in met_classes:
            mets = ANALYTE_META_DATA[ANALYTE_META_DATA['pathway_class'] == met_class]['analyte_name']

            page_created = False

            for met in mets:
                ratio_value = False
                if "<div>" in met:
                    ratio_value = True
                    met = "[" + met.replace("<div>", "/") + "]"

                print("_________________STARTING > " + met)

                table = self.table_of_plots_for_metabolite(met)
                if table is not None:
                    if not page_created:
                        self.append_new_page_title(met_class)
                        page_created = True

                    if not ratio_value:
                        header = Paragraph(met, STYLE_MET_NAME)
                    else:
                        header = Paragraph(met, STYLE_MET_RATIO_NAME)

                    self.add_element(KeepTogether([header, table]))
                    self.add_element(VGAP_04)
                    print("_________________PROCESSED > " + met)

            if page_created:
                self.add_element(PageBreak())


    def output_build_to_path(self, file_path):
        try:
            # Build the PDF document
            self.doc.build(self.page_elements)
            print("PDF built successfully.")
        except Exception as e:
            print("Error while building PDF: ", e)

        try:
            # Move the file to the specified path
            os.rename(self.doc.filename, file_path)
            print("PDF moved to {}".format(file_path))
        except Exception as e:
            print("Error while moving PDF: ", e)


def generateReport():
    try:
        filename = askopenfilename()
        iso_dist_data = pd.read_excel(filename, sheet_name="Normalized", dtype={0: 'str'})
        pool_size_data = pd.read_excel(filename, sheet_name="PoolAfterDF", dtype={0: 'str'})

        pool_size_data.iloc[:, 0] = pool_size_data.iloc[:, 0].str.strip()

        print("____________________READ_IN_DATA_________________________")

        normalizing_rows = []
        normalizing_names = []
        row_index = 0
        for row in pool_size_data['Compound']:
            if "*" in row:
                normalizing_rows.append(row_index)
                normalizing_names.append(row)
            row_index += 1

        if len(normalizing_rows) > 0:
            for col in range(1, len(pool_size_data.columns)):
                factor = 0.0001
                for row in normalizing_rows:
                    factor *= float(pool_size_data.iloc[row, col])
                    pool_size_data.iloc[1:, col] /= factor

        met_names = ANALYTE_META_DATA['analyte_name'].dropna()

        for row in met_names:
            if "<div>" in row:
                mets = row.split("<div>")
                row1 = pool_size_data[pool_size_data['Compound'] == mets[0]]
                row2 = pool_size_data[pool_size_data['Compound'] == mets[1]]
                if not row1.empty and not row2.empty:
                    ratio_vals = [str("[" + mets[0] + "/" + mets[1] + "]")]
                    for col in range(1, len(pool_size_data.columns)):
                        v1 = row1.iloc[0, col]
                        v2 = row2.iloc[0, col]
                        val = 0
                        if not np.isnan(v1) and not np.isnan(v2) and v2 > 0:
                            val = v1 / v2
                        ratio_vals.append(val)
                    pool_size_data.loc[row_index] = ratio_vals
                    row_index += 1

        print(pool_size_data)

        print("____________________NORMALIZED_POOL_SIZE_________________________")

        group_names = {}
        i = 1
        for name in groups_str.get().split(","):
            group_names[i] = name
            i += 1

        report = MetaboliteReport(pdf_filename.get(),
                                  iso_dist_data, pool_size_data,
                                  show_iso=labelling_present.get(),
                                  group_names=group_names)

        report.configureTitlePage(title=title_str.get(),
                                  project=project_str.get(),
                                  exp_title=exp_title.get(),
                                  date=datetime.datetime.now().strftime('%b %d, %Y'),
                                  leader=leader_str.get(),
                                  exp_summary=description_box.get(1.0, "end-1c"),
                                  normalization=normalizing_names)

        report.build_with_metabolite_classes(("glycolysis",
                                              "TCA",
                                              #"pyruvate ratios",
                                              "amino acid",
                                              "PPP",
                                              "polyol-pathway",
                                              "nucleotide",
                                              "lipid metabolism",
                                              "energy-carrier",
                                              "redox",
                                              "urea cycle",
                                              "SAM-cycle",
                                              "vitamin",
                                              "small-molecule",
                                              "succinilation",
                                              "purine intermediate",
                                              "tryptophan metabolism",
                                              "histidine metabolism",
                                              "other"
                                              ))

        # Get the current datetime and format it
        current_datetime = datetime.datetime.now().strftime("%Y%m%d")
        # Insert the datetime at the beginning of the filename
        desktop_path = os.path.join(os.path.expanduser("~"), "Desktop", "{}_{}.pdf".format(current_datetime, pdf_filename.get()))
        
        report.output_build_to_path(desktop_path)


    except Exception as e:
        print("Error generating report:", e)

def savePreviousEntry():
    data = {}
    data['previous entry'] = []
    data['previous entry'].append({
        'project_str': project_str.get(),
        'exp_title': exp_title.get(),
        'leader_str': leader_str.get(),
        'description_box': description_box.get("1.0", "end-1c"),
        'groups_str': groups_str.get(),
        'pdf_filename': pdf_filename.get(),
        'labelling_present': labelling_present.get()
    })
    with open(save_data_filepath, 'w') as json_file:
        json.dump(data, json_file, indent=4)

def loadPreviousEntry():
    try:
        with open(save_data_filepath, 'r') as json_file:
            data = json.load(json_file)
            project_str.set(data["previous entry"][0]["project_str"])
            exp_title.set(data["previous entry"][0]["exp_title"])
            leader_str.set(data["previous entry"][0]["leader_str"])
            description_box.insert("1.0", data["previous entry"][0]["description_box"])
            groups_str.set(data["previous entry"][0]["groups_str"])
            pdf_filename.set(data["previous entry"][0]["pdf_filename"])
            labelling_present.set(data["previous entry"][0]["labelling_present"])
    except Exception as e:
        print("Error loading previous entry:", e)

window = tk.Tk()
window.title("Metabolite Report Generator v1")
text_font = font.Font(family='Arial')
window.option_add("*Font", text_font)

menu_title_list = ["Metabolomics Analysis", "Isotopomer Analysis", "Pool Size Analysis"]

title_str = tk.StringVar()
title_str.set(menu_title_list[0])

tk.Label(window, text="Analysis").grid(sticky="w", padx=(5, 20))
option = tk.OptionMenu(window, title_str, *menu_title_list).grid(row=0, column=1, padx=5, pady=10)

project_str = tk.StringVar()
tk.Label(window, text="Project").grid(sticky="w", row=1, padx=(5, 20))
tk.Entry(window, textvariable=project_str, width=50).grid(sticky="w", row=1, column=1, pady=10, padx=5)

exp_title = tk.StringVar()
tk.Label(window, text="Experiment Title").grid(sticky="w", row=2, padx=(5, 20))
tk.Entry(window, textvariable=exp_title, width=50).grid(sticky="w", row=2, column=1, pady=10, padx=5)

leader_str = tk.StringVar()
tk.Label(window, text="Leader").grid(sticky="w", row=3, padx=(5, 20))
tk.Entry(window, textvariable=leader_str, width=50).grid(sticky="w", row=3, column=1, pady=10, padx=5)

tk.Label(window, text="Experimental Details").grid(sticky="w", row=4, padx=(10, 20))
description_box = tk.Text(window, height=6, width=50)
description_box.grid(sticky="w", row=4, column=1, padx=5, pady=10)

groups_str = tk.StringVar()
tk.Label(window, text="Group Names (comma-separated)").grid(sticky="w", row=5, padx=(5, 20))
tk.Entry(window, textvariable=groups_str, width=50).grid(sticky="w", row=5, column=1, pady=10, padx=5)

pdf_filename = tk.StringVar()
tk.Label(window, text="PDF Filename (without .pdf)").grid(sticky="w", row=6, padx=(5, 20))
tk.Entry(window, textvariable=pdf_filename, width=50).grid(sticky="w", row=6, column=1, pady=10, padx=5)

labelling_present = tk.BooleanVar()
tk.Label(window, text="Labelling Present").grid(sticky="w", row=7, padx=(5, 20))
tk.Checkbutton(window, variable=labelling_present).grid(sticky="w", row=7, column=1, pady=10, padx=5)

tk.Button(window, text="Generate Report", command=generateReport).grid(sticky="e", row=8, column=1, padx=10, pady=10)
tk.Button(window, text="Save Entry", command=savePreviousEntry).grid(sticky="e", row=9, column=1, padx=5, pady=5)
tk.Button(window, text="Load Entry", command=loadPreviousEntry).grid(sticky="e", row=10, column=1, padx=5, pady=5)

window.mainloop()

