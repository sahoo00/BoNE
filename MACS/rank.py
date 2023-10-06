
# coding: utf-8

# In[2]:

# IMPORT STATEMENTS
import cv2
import re
import numpy as np
from scipy import stats
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
from pprint import pprint
import os
import pickle
import sys
sys.path.append("/booleanfs/sahoo/scripts/")
sys.path.append("/booleanfs2/sahoo/Hegemon/")
import StepMinerNew as smn
import HegemonUtil as hu
reload(hu)

def readGenes(cfile):
    genes = "PRKAB1 PPARG PPP1CA PPARGC1A SIRT6 OCLN PARD3 IL23A IL6 S1PR5 ACTA2 CCDC88A COL1A1 CXCL10 ELMO1 HIF1A IL10 IL33 ITGA4 ITGB1 ITGB7 JAK1 MMP14 MMP2 MMP9 MRC1 NOD2 PML PRKCQ RIPK2 S1PR1 SNAI2 SPHK1 TGFB1 TIMP2 TLR2 TLR4 VIM CLDN2 IL11 MMP1 MMP3 CEMIP KIAA1199 PRKAA2 IL8 CXCL8 LGR5"
    if not os.path.isfile(cfile):
        print "Can't open file {0} <br>".format(cfile);
        exit()
    fp = open(cfile, "r")
    nodelist = re.split("[\[\]()\s]", genes)
    for line in fp:
        line = line.strip();
        ll = re.split("[\[\]()\s]", line);
        nodelist += ll
    fp.close();
    return [i for i in hu.uniq(nodelist) if i != '']


def plotSizes(csizes):
    w,h, dpi = (6, 4, 100)
    fig = plt.figure(figsize=(w,h), dpi=dpi)
    ax = fig.add_axes([70.0/w/dpi, 54.0/h/dpi, 1-2*70.0/w/dpi, 1-2*54.0/h/dpi])
    ax.loglog(range(len(csizes)), csizes, "r-", clip_on=False);
    ax.grid(False)
    ax.set_axis_bgcolor("white")
    for child in ax.get_children():
        if isinstance(child, matplotlib.spines.Spine):
            child.set_color('black')
            child.set_linewidth(0.5)
    ax.tick_params(direction='out', length=4, width=1, colors='k', top=False,
            right=False)
    ax.tick_params(which="minor", direction='out', length=2, width=0.5, colors='k', top=False, right=False)
    ax.set_xlabel("Clusters ranked by size")
    ax.set_ylabel("Cluster sizes")
    fig.savefig("Supplementary/cluster-sizes-1.pdf", dpi=200)

def getSynDictHs():
    cfile = "/booleanfs2/sahoo/Data/Sarah/CellCycle/database/gene-info-hs.txt"
    fp = open(cfile, "r")
    ghash = {}
    head = fp.readline()
    for line in fp:
        line = line.strip();
        ll = line.split("\t");
        if ll[2] == "0":
            continue
        for i in ll[3:]:
            if i not in ghash:
                ghash[i] = []
            ghash[i] += [ll[0]]
    fp.close();
    return ghash

def getGPL570():
    cfile = "/booleanfs/sahoo/Data/Annotation/GPL570.txt"
    fp = open(cfile, "r")
    ghash = {}
    head = fp.readline()
    for line in fp:
        line = line.strip();
        ll = line.split("\t");
        ghash[ll[0]] = ll
    fp.close();
    return ghash

reload(hu)
db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
h1 = hu.Hegemon(db.getDataset("PLP33"))
h1.init()
h1.initPlatform()
h1.initSurv()
with open('Supplementary/path-8.json') as data_file:
    data_item = json.load(data_file)
#pprint(data_item[7]) # 7, 6, 5, 1
cfile = "Supplementary/colon-network-clusters.txt"
if not os.path.isfile(cfile):
    print "Can't open file {0} <br>".format(cfile);
    exit()
fp = open(cfile, "r")
nodelist = {}
nhash = {}
csizes = []
for line in fp:
    line = line.strip();
    ll = line.split("\t");
    nodelist[ll[0]] = ll[2:]
    csizes.append(len(ll[2:]))
    for i in ll[2:]:
        nhash[i] = ll[0];
fp.close();
csizes.sort(reverse=True)
gene_groups = []
order = [1, 4, 5];
order = [8, 7, 6, 1, 2, 3];
order = [2, 1, 4, 5, 6, 7];
order = [8, 7, 6, 1, 2];
for i in range(len(order)):
    gene_groups.append(set())
    for g in data_item[order[i]][2]:
        gene_groups[i].add(g[0])
        if g[0] in nodelist:
            for k in nodelist[g[0]]:
                gene_groups[i].add(k)
    for g in data_item[order[i]][3]:
        gene_groups[i].add(g)
        if g in nodelist:
            for k in nodelist[g]:
                gene_groups[i].add(k)
#print([len(s) for s in gene_groups])
db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
h = hu.Hegemon(db.getDataset("PLP62.3"))
h.init()
h.initPlatform()
h.initSurv()
id = "PRKAB1"
e = h.getExprData(id)
atype = h.getSurvName("c Subtype")
btype = h.getSurvName("c Type")
#print hu.uniq(atype)
atypes = ['CFP', 'CAP']
ahash = {}
for i in range(len(atypes)):
    ahash[atypes[i]] = i
aval = [ahash[i] if i in ahash else None for i in atype]
cfp = [ i for i in h.aRange() if atype[i] == "CFP"]
cap = [ i for i in h.aRange() if atype[i] == "CAP"]
expg = cfp + cap
order = cfp + cap
expr = []
col_labels = [h.headers[i] for i in order]
row_labels = []
ghash = getSynDictHs()
gpl570 = getGPL570()
for gn in gene_groups[3]:
  gn1 = gn
  if gn == "IL8":
    gn = "CXCL8"
  if gn == "KIAA1199":
    gn = "CEMIP"
  if gn == "LOC730101":
    gn = "ENSG00000216775"
  if gn == "LOC645166":
    gn = "ENSG00000283196"
  if gn == "LOC254057":
    gn = "ENSG00000244300"
  l1 = h.getIDs(gn)
  if len(l1) == 0:
      if gn in ghash:
          for k in ghash[gn]:
              l1.update(h.getIDs(k))
  if len(l1) == 0:
      id1 = h1.getSimpleName(gn)
      if id1 in gpl570:
          l2 = gpl570[id1]
          for k in l2[3:]:
              ll = re.split("[\s]", k)
              for m in ll:
                  l1.update(h.getIDs(m))
  if len(l1) == 0:
      id1 = h1.getSimpleName(gn)
      print gn, id1, "not found"
  for id in l1:
    e = h.getExprData(id);
    t = h.getThrData(id);
    v1 = np.array([float(e[i]) for i in cfp])
    v2 = np.array([float(e[i]) for i in cap])
    t, p = stats.ttest_ind(v1,v2)
    print id, gn, t, p, np.mean(v1)-np.mean(v2), h1.getSimpleName(gn1)
