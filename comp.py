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
sys.path.append("/booleanfs2/sahoo/Hegemon/")
import StepMiner as smn
import HegemonUtil as hu
import bone
reload(bone)

plt.rc('text', usetex=True)

import bone
reload(bone)
ng = [0, 7, 6, 5, 1, 16, 17]
genes, wt1, l1 = bone.getGeneGroups([ng[j] for j in [1, 2, 3]], [-3, -2, -1],
        0)
#genes, wt1, l1 = bone.getGeneGroups([ng[j] for j in [4, 5, 6]], [1, 1, 2], 0)

acolor = ["#00CC00", "yellow", "#EC008C",
        'cyan', "#F7941D", "#808285",
        'blue', 'black', 'green', 'red']

ana = bone.IBDAnalysis()
ana.getPeters(2)

ana.orderData(l1, wt1)
fig = plt.figure(figsize=(4,4), dpi=100)
ax = plt.subplot2grid((4, 1), (0, 0))
params = {'spaceAnn': len(ana.order)/len(ana.atypes),
        'tAnn': 1, 'widthAnn':1, 'acolor': acolor, 'ax': ax,
        'w': 5, 'h': 0.8, 'atypes': ana.atypes ,'cval':
        ana.cval}
ax = ana.printTitleBar(params)
res = ana.getMetrics(ana.cval[0])
ax.text(len(ana.cval[0]), 4, ",".join(res))
ax = plt.subplot2grid((4, 1), (1, 0), rowspan=3)
ax = ana.densityPlot(ax)
plt.tight_layout()

fig.savefig("comp-1.pdf", dpi=100)

