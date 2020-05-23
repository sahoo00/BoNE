# IMPORT STATEMENTS
import cv2
import re
import numpy as np
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from mpl_toolkits.axes_grid1 import SubplotDivider, Size
from mpl_toolkits.axes_grid1.mpl_axes import Axes
import matplotlib.patches as patches
import matplotlib.colors as colors
from matplotlib.transforms import *
import PIL
import math
#get_ipython().magic(u'matplotlib inline')
import pandas as pd
import seaborn as sns
import json
from sklearn.metrics import *
from scipy.stats import *
from pprint import pprint
import os
import pickle
import sys
sys.path.append("Hegemon")
import StepMiner as smn
import HegemonUtil as hu
try:
    reload  # Python 2.7
except NameError:
    try:
        from importlib import reload  # Python 3.4+
    except ImportError:
        from imp import reload  # Python 3.0 - 3.3

import bone
reload(bone)
acolor = ["#00CC00", "#D8A03D","#EC008C",
        'cyan', "#B741DC", "#808285",
        'blue', 'black', 'green', 'red',
        'orange', 'brown', 'pink', 'purple',
        'salmon', 'greenyellow', 'skyblue', 'plum',
        'sienna', 'darkseagreen', 'teal', 'magenta']


plt.rc('text', usetex=True)

def getPDF(cfile):
    import bone
    reload(bone)
    from matplotlib.backends.backend_pdf import PdfPages

    pdf = PdfPages(cfile)
    return pdf

def closePDF(pdf):
    import datetime
    d = pdf.infodict()
    d['Title'] = 'Prediction M1/M2 state'
    d['Author'] = 'Debashis Sahoo'
    d['Subject'] = 'COVID-19'
    d['Keywords'] = 'disease training validation ROC'
    d['CreationDate'] = datetime.datetime(2020, 1, 16)
    d['ModDate'] = datetime.datetime.today()
    pdf.close()

def T1():
    pdf = getPDF('results/heatma-1.pdf')
    import HegemonUtil as hu
    reload(hu)
    import bone
    reload(bone)
    ng = [0, 7, 6, 5, 1, 16, 17]
    genes, wt1, l1 = bone.getGeneGroups([ng[j] for j in [1, 2, 3]], [-3, -2, -1], 1)
    ana = bone.IBDAnalysis()
    ana.getPetersDf()
    ana.orderDataDf(l1, wt1)
    ofile = "results/heatmap-test.pdf"
    params = {'dx': 100, 'dy': 10, 'spaceAnn': 30, 'tAnn': 1, 'widthAnn':3,
              'sy': 35, 'thr': 1, 'w': 6, 'h': 6,
              'genes': genes, 'atypes': ana.atypes,'cval': ana.cval,
              'tl': 6, 'tw': 0.25, 'ts': 10, 'tsi': -100}
    i1 = ana.i1
    f_ranks = ana.f_ranks
    ana.params = {'genes': genes,'atypes': ana.atypes,'cval': ana.cval}
    ana.params.update(params)
    ax, divider = bone.plotHeatmap(ofile, ana.data, ana.col_labels,
            ana.row_labels, ana.params)
    pdf.savefig()

    closePDF(pdf)

T1()
