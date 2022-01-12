
# coding: utf-8

# In[1]:

import cv2
import re
import numpy as np
import matplotlib
#matplotlib.use('agg')
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
from scipy.stats import fisher_exact, ttest_ind
import scipy.stats
from pprint import pprint
import os
import pickle
import sys
sys.path.append("/booleanfs2/sahoo/Hegemon/")
import StepMiner as smn
import HegemonUtil as hu
acolor = ["#00CC00", "#D8A03D","#EC008C",
        'cyan', "#B741DC", "#808285",
        'blue', 'black', 'green', 'red']

def getGhash(l1):
    ghash = {}
    for i in range(len(l1)):
        for name in l1[i]:
            ghash[name] = i
    return ghash

def getBooleanInfo(cfile, dbid, l1):
    ghash = getGhash(l1)
    db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
    h = hu.Hegemon(db.getDataset(dbid))
    h.init()
    h.initPlatform()
    h.initSurv()
    ifile = h.rdataset.getInfo()
    fp = open(ifile, "r")
    ihash = {}
    head = fp.readline().strip();
    ll= re.split("[\t]", head)
    hh1 = {}
    for i in range(len(ll)):
        hh1[ll[i]] = i
    for line in fp:
        line = line.strip();
        ll = re.split("[\t]", line);
        ihash[ll[0]] = ll
    fp.close();

    qd = [0, 1, 2, 3]
    fp = open(cfile, "r")
    res = {}
    head = fp.readline().strip();
    for line in fp:
        line = line.strip();
        ll = re.split("[\t]", line);
        m1 = float(ihash[ll[1]][hh1["min"]])
        m2 = float(ihash[ll[1]][hh1["max"]])
        sd = float(ihash[ll[1]][hh1["sd"]])
        dr = float("%.2f" % (m2 - m1))
        test = 1
        name = h.getSimpleName(ll[1])
        if (name not in ghash and ll[1] not in ghash):
            continue
        if name not in res:
            res[name] = [{ll[1]:[dr, sd]}]
            for i in qd:
                s = float(ll[2 + 4 * 1 + i])
                p = float(ll[2 + 4 * 2 + i])
                res[name] += [[s, p, ll[1]]]
        else:
            res[name][0][ll[1]] = [dr, sd]
            j = 1
            for i in qd:
                s = float(ll[2 + 4 * 1 + i])
                p = float(ll[2 + 4 * 2 + i])
                if (s > res[name][j][0] and p < res[name][j][1]):
                    res[name][j][0] = s
                    res[name][j][1] = p
                    res[name][j][2] = ll[1]
                j += 1
    fp.close();
    return res

def getGeneInfo(dbid, l1):
    ghash = getGhash(l1)
    db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
    h = hu.Hegemon(db.getDataset(dbid))
    h.init()
    h.initPlatform()
    h.initSurv()
    ifile = h.rdataset.getInfo()
    fp = open(ifile, "r")
    ihash = {}
    head = fp.readline().strip();
    ll= re.split("[\t]", head)
    hh1 = {}
    for i in range(len(ll)):
        hh1[ll[i]] = i
    for line in fp:
        line = line.strip();
        ll = re.split("[\t]", line);
        if ll[1] not in ghash:
            continue
        m1 = float(ll[hh1["min"]])
        m2 = float(ll[hh1["max"]])
        sd = float(ll[hh1["sd"]])
        thr = float(ll[hh1["thr"]])
        dr = float("%.2f" % (m2 - m1))
        if ll[1] not in ihash:
            ihash[ll[1]] = [dr, sd, thr, ll[0]]
        else:
            if dr > ihash[ll[1]][0]:
                ihash[ll[1]] = [dr, sd, thr, ll[0]]
    fp.close()
    return ihash

def getInfo(dbid):
    db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
    h = hu.Hegemon(db.getDataset(dbid))
    h.init()
    ifile = h.rdataset.getInfo()
    fp = open(ifile, "r")
    ihash = {}
    head = fp.readline().strip();
    ll= re.split("[\t]", head)
    hh1 = {}
    for i in range(len(ll)):
        hh1[ll[i]] = i
    for line in fp:
        line = line.strip();
        ll = re.split("[\t]", line);
        ihash[ll[0]] = ll
    fp.close()
    return ihash

def getCorrInfo(cfile, thr, l1):
    ghash = getGhash(l1)
    fp = open(cfile, "r")
    res = {}
    for line in fp:
        line = line.strip();
        ll = re.split("[\t]", line);
        l1 = re.split(": ", ll[3])
        l2 = re.split(" /// ", l1[0])
        test= 0
        name = l2[0]
        for k in l2:
            if k in ghash:
                test = 1
                name = k
        if test == 0:
            continue
        if name not in res:
            res[name] = [float(ll[0]), ll[2]]
        else:
            if abs(float(ll[0])) > abs(res[name][0]):
                res[name] = [float(ll[0]), ll[2]]
    fp.close();
    return res

def getCellInfo(l1, ahash):
    ghash = getGhash(l1)
    db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
    h = hu.Hegemon(db.getDataset("G19"))
    h.init()
    h.initPlatform()
    h.initSurv()
    atype = h.getSurvName("c cell type")
    order = [i for i in h.aRange() if atype[i] in ahash]
    ifile = h.rdataset.getInfo()
    fp = open(ifile, "r")
    ihash = {}
    head = fp.readline().strip();
    ll= re.split("[\t]", head)
    hh1 = {}
    for i in range(len(ll)):
        hh1[ll[i]] = i
    for line in fp:
        line = line.strip();
        ll = re.split("[\t]", line);
        if ll[1] not in ghash:
            continue
        m1 = float(ll[hh1["min"]])
        m2 = float(ll[hh1["max"]])
        sd = float(ll[hh1["sd"]])
        thr = float(ll[hh1["thr"]])
        dr = float("%.2f" % (m2 - m1))
        ex = h.getExprData(ll[0])
        v1 = [float(ex[i]) for i in order]
        vm1 = np.mean(v1)
        vdr = float("%.2f" % (vm1 - thr))
        if ll[1] not in ihash:
            ihash[ll[1]] = [vdr, sd, thr, ll[0]]
        else:
            if vdr > ihash[ll[1]][0]:
                ihash[ll[1]] = [vdr, sd, thr, ll[0]]
    fp.close()
    return ihash

def getTCell(l1):
    ahash = {'CD4+ Central Memory':1, 'CD4+ Effector Memory':1,
             'Naive CD4+ T-cell':1, 'CD8+ Central Memory':1,
             'Naive CD8+ T-cell':1, 'CD8+ Effector Memory RA':1,
             'CD8+ Effector Memory':1}
    return getCellInfo(l1, ahash)

def getBCell(l1):
    ahash = {'Mature B-cell class able to switch':1,
              'Na\xc3\xafve B-cells':1, 'Early B-cell':1,
              'Pro B-cell':1, 'Mature B-cell class switched':1,
              'Mature B-cells':1}
    return getCellInfo(l1, ahash)

def getCode(p):
    if p <= 0:
        return '0'
    if p <= 0.001:
        return '***'
    if p <= 0.01:
        return '**'
    if p <= 0.05:
        return '*'
    if p <= 0.1:
        return '.'
    return ''

def printOLS(fm, df1):
    import statsmodels.formula.api as smf
    lm1 = smf.ols(formula=fm, data=df1).fit()
    print(lm1.summary())
    idx = lm1.params.index
    ci = lm1.conf_int()
    ci_1 = [ ci[0][i] for i in range(len(idx))]
    ci_2 = [ ci[1][i] for i in range(len(idx))]
    c_1 = [ getCode(p) for p in lm1.pvalues]
    df = pd.DataFrame({'Name': idx, 
	'coeff' : lm1.params, 'lower 0.95' : ci_1,
        'upper 0.95' : ci_2, 'pvalues' : lm1.pvalues, 'codes': c_1},
        columns=['Name', 'coeff', 'lower 0.95',
	    'upper 0.95', 'pvalues', 'codes'])
    for i in range(len(idx)):
        print('%s\t%.2f\t(%0.2f - %0.2f)\t%0.3f' % \
                (idx[i], lm1.params[i], ci[0][i], ci[1][i], lm1.pvalues[i]))
    print(df.to_string(formatters={'coeff':'{:,.2f}'.format,
        'lower 0.95':'{:,.2f}'.format, 'upper 0.95':'{:,.2f}'.format,
        'pvalues': '{:,.3f}'.format}))
    return df

def printStats(cfile, thr):
    if not os.path.isfile(cfile):
        print("Can't open file {0} <br>".format(cfile));
        exit()
    fp = open(cfile, "r")
    numhigh = 0
    numlow = 0
    total = 0
    for line in fp:
        line = line.strip();
        ll = re.split("[\t]", line);
        if float(ll[2]) >= 0 and float(ll[3]) < thr:
            numhigh += 1
        if float(ll[2]) < 0  and float(ll[3]) < thr:
            numlow += 1
        total += 1
    fp.close();
    print(cfile, numhigh, numlow, total)

def getStats(cfile, thr, index):
    if not os.path.isfile(cfile):
        print("Can't open file {0} <br>".format(cfile));
        exit()
    fp = open(cfile, "r")
    high = set()
    low = set()
    for line in fp:
        line = line.strip();
        ll = re.split("[\t]", line);
        if float(ll[2]) >= 0 and float(ll[3]) < thr:
            high.add(ll[index])
        if float(ll[2]) < 0  and float(ll[3]) < thr:
            low.add(ll[index])
    fp.close();
    return high, low

def getEntries(cfile, index):
    if not os.path.isfile(cfile):
        print("Can't open file {0} <br>".format(cfile));
        exit()
    fp = open(cfile, "r")
    res = []
    for line in fp:
        line = line.strip();
        ll = re.split("[\t]", line);
        res += [ll[index]]
    fp.close();
    return res

def getPVal(cfile):
    return [float(i) for i in getEntries(cfile, 3)]

def getFdrStats(cfile, thr, index):
    pval = getPVal(cfile)
    ids = getEntries(cfile, index)
    stat = getEntries(cfile, 2)
    from statsmodels.stats import multitest
    mstat = multitest.multipletests(pval, thr, 'fdr_bh')
    high = set()
    low = set()
    for i in range(len(pval)):
        if mstat[0][i] and float(stat[i]) >= 0:
            high.add(ids[i])
        if mstat[0][i] and float(stat[i]) < 0:
            low.add(ids[i])
    return high, low

def readGenes(cfile):
    genes = "PRKAB1 PPARG PPP1CA PPARGC1A SIRT6 OCLN PARD3 IL23A IL6 S1PR5 ACTA2     CCDC88A COL1A1 CXCL10 ELMO1 HIF1A IL10 IL33 ITGA4 ITGB1 ITGB7 JAK1 MMP14 MMP2     MMP9 MRC1 NOD2 PML PRKCQ RIPK2 S1PR1 SNAI2 SPHK1 TGFB1 TIMP2 TLR2 TLR4 VIM CLDN2     IL11 MMP1 MMP3 CEMIP KIAA1199 PRKAA2 IL8 CXCL8 LGR5"
    genes = ""
    if not os.path.isfile(cfile):
        print("Can't open file {0} <br>".format(cfile));
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
    ax.tick_params(which="minor", direction='out', length=2, width=0.5, 
                   colors='k', top=False, right=False)
    ax.set_xlabel("Clusters ranked by size")
    ax.set_ylabel("Cluster sizes")
    fig.savefig("Supplementary/cluster-sizes-1.pdf", dpi=200)

def printReport(actual, predicted, score, target_names):
    print(classification_report(actual, predicted, target_names=target_names))
    fpr, tpr, _ = roc_curve(actual, score)
    roc_auc = auc(fpr, tpr)
    print('ROC AUC', roc_auc)
    print('ROC AUC', roc_auc_score(actual, score))
    print('Accuracy', accuracy_score(actual, predicted))
    wi,hi, dpi = (4, 4, 100)
    fig = plt.figure(figsize=(wi,hi), dpi=dpi)
    ax = fig.add_axes([70.0/wi/dpi, 54.0/hi/dpi, 1-2*70.0/wi/dpi, 1-2*54.0/hi/dpi])
    ax.plot(fpr, tpr, color='darkorange',
            lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
    ax.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('Receiver operating characteristic')
    ax.legend(loc="lower right")
    return ax

def convertScore(mylist):
    hs = dict()
    for x in mylist:
        if (x not in hs):
            hs[x] = 1
        else:
            hs[x] += 1
    keys = hs.keys()
    values = [0] + list(np.cumsum(hs.values()))
    for i in range(len(keys)):
        hs[keys[i]] = values[i]
    hh = dict()
    res = []
    for x in mylist:
        if (x not in hh):
            hh[x] = 0
        else:
            hh[x] += 1
        res += [hh[x] + hs[x]]
    return res

def mergeRanks(group, start, exp, weight):
    X = np.array([[e[k-start] for e in exp] for k in group])
    arr = np.dot(X, np.array(weight))
    return arr

def mergeRanks2(group, exp, weight):
    X = np.array([[e[k] for e in exp] for k in range(len(group))])
    arr = np.dot(X, np.array(weight))
    return arr

def getOrder(group, start, exp, weight):
    arr = mergeRanks(group, start, exp, weight)
    return [group[i] for i in np.argsort(arr)]

def getRanks(gene_groups, h):
    expr = []
    row_labels = []
    ranks = []
    for s in gene_groups:
        print(len(s), s)
        count = 0
        avgrank = [0 for i in h.aRange()]
        for gn in s:
          for id in h.getIDs(gn):
            e = h.getExprData(id);
            t = h.getThrData(id);
            if e[-1] == "":
                continue
            v = np.array([float(e[i]) for i in h.aRange()])
            te = []
            for i in h.aRange():
                v1 = (float(e[i]) - t[3]) / 3;
                if np.std(v) > 0:
                    v1 = v1 / np.std(v)
                avgrank[i-h.start] += v1
                te.append(v1)
            expr.append(te)
            #row_labels.append(h.getSimpleName(id))
            row_labels.append(gn)
            count += 1
            #if count > 100:
            #    break
        ranks.append(avgrank)
    return ranks, row_labels, expr

def getRanks2(gene_groups, h):
    expr = []
    row_labels = []
    row_ids = []
    row_numhi = []
    ranks = []
    g_ind = 0
    counts = []
    for s in gene_groups:
        count = 0
        avgrank = [0 for i in h.aRange()]
        for gn in s:
          for id in h.getIDs(gn):
            e = h.getExprData(id);
            t = h.getThrData(id);
            if e[-1] == "":
                continue
            v = np.array([float(e[i]) if e[i] != "" else 0 for i in h.aRange()])
            te = []
            sd = np.std(v)
            for i in h.aRange():
                if (e[i] != ""):
                    v1 = (float(e[i]) - t[3]) / 3;
                    if sd > 0:
                        v1 = v1 / sd
                else:
                    v1 = -t[3]/3/sd
                avgrank[i-h.start] += v1
                te.append(v1)
            expr.append(te)
            row_labels.append(h.getSimpleName(id))
            row_ids.append(id)
            v1 = [g_ind, sum(v > t[3])]
            if g_ind > 3:
                v1 = [g_ind, sum(v <= t[3])]
            else:
                v1 = [g_ind, sum(v > t[3])]
            row_numhi.append(v1)
            count += 1
            #if count > 200:
            #    break
        ranks.append(avgrank)
        g_ind += 1
        counts += [count]
    print(counts)
    return ranks, row_labels, row_ids, row_numhi, expr

def getRanks3(gene_groups, h, order):
    expr = []
    row_labels = []
    row_ids = []
    row_numhi = []
    ranks = []
    g_ind = 0
    for s in gene_groups:
        count = 0
        avgrank = [0 for i in order]
        for gn in s:
          for id in h.getIDs(gn):
            e = h.getExprData(id);
            if e[-1] == "":
                continue
            v = np.array([float(e[i]) for i in order])
            t = hu.getThrData(v)
            te = []
            for i in range(len(order)):
                v1 = (float(e[order[i]]) - t[3]) / 3;
                if np.std(v) > 0:
                    v1 = v1 / np.std(v)
                avgrank[i] += v1
                te.append(v1)
            expr.append(te)
            row_labels.append(h.getSimpleName(id))
            row_ids.append(id)
            v1 = [g_ind, sum(v > t[3])]
            if g_ind > 3:
                v1 = [g_ind, sum(v <= t[3])]
            else:
                v1 = [g_ind, sum(v > t[3])]
            row_numhi.append(v1)
            count += 1
            #if count > 200:
            #    break
        ranks.append(avgrank)
        g_ind += 1
    return ranks, row_labels, row_ids, row_numhi, expr

def saveList(ofile, l1):
    of = open(ofile, "w")
    for i in l1:
        of.write("\t".join([i]) + "\n")
    of.close()

def readList(cfile):
    if not os.path.isfile(cfile):
        print("Can't open file {0} <br>".format(cfile));
        exit()
    fp = open(cfile, "r")
    f_order = []
    for line in fp:
        line = line.strip();
        ll = re.split("[\s]", line);
        f_order += [ll[0]]
    fp.close();
    return f_order

def saveCData(ofile, h, i1, f_ranks):
    f_order = dict()
    for i in h.aRange():
        f_order[i] = ""
    for i in range(len(i1)):
        f_order[i1[i]] = str(i)
    of = open(ofile, "w")
    for i in h.aRange():
        id1 = h.headers[i]
        of.write("\t".join([id1, f_order[i], str(f_ranks[i - h.start])]) + "\n")
    of.close()

def saveHeatmapData(ofile, row_labels, row_numhi, row_ids, index, expr):
    ind_r = np.array(sorted(range(len(row_labels)), key=lambda x: (row_numhi[x][0], row_numhi[x][1])))
    of = open(ofile, "w")
    for i in ind_r:
        id1 = row_ids[i]
        of.write("\t".join([id1, row_labels[i], str(row_numhi[i][0]),                             str(row_numhi[i][1])] + [str(expr[i][j]) for j in index]) + "\n")
    of.close()

def readCData(cfile):
    if not os.path.isfile(cfile):
        print("Can't open file {0} <br>".format(cfile));
        exit()
    fp = open(cfile, "r")
    f_order = []
    f_ranks = []
    for line in fp:
        line = line.strip();
        ll = re.split("[\s]", line);
        f_order += [ll[1]]
        f_ranks += [float(ll[2])]
    fp.close();
    return f_order, f_ranks

def readHeatmapData(cfile):
    if not os.path.isfile(cfile):
        print("Can't open file {0} <br>".format(cfile));
        exit()
    fp = open(cfile, "r")
    row_labels, row_numhi, row_ids, expr = [], [], [], []
    for line in fp:
        line = line.strip();
        ll = re.split("\t", line);
        row_ids += [ll[0]]
        row_labels += [ll[1]]
        row_numhi += [ [int(ll[2]), int(ll[3])] ]
        expr += [[float(k) for k in ll[4:]]]
    fp.close();
    return row_labels, row_numhi, row_ids, expr

def barTop(tax, atypes, color_sch1, params):
    spaceAnn = 70
    widthAnn = 3
    tAnn = 1
    if 'spaceAnn' in params:
        spaceAnn = params['spaceAnn']
    if 'widthAnn' in params:
        widthAnn = params['widthAnn']
    if 'tAnn' in params:
        tAnn = params['tAnn']
    for i in range(len(atypes)):
        tax.add_patch(patches.Rectangle( (i *spaceAnn, 0), widthAnn, 3,
                                        facecolor=color_sch1[i], edgecolor="none", alpha=1.0))
        tax.text(i * spaceAnn + widthAnn + tAnn, 1, atypes[i], rotation='horizontal',
                 ha='left', va='center', fontsize=12)

def plotHeatmap(ofile, data, col_labels, row_labels, params):
    
    genes = []
    atypes = []
    cval = []
    dpi, tl, tw, ts, tsi = (100, 3, 0.25, 0.5, 0)
    if 'genes' in params:
        genes = params['genes']
    if 'atypes' in params:
        atypes = params['atypes']
    if 'cval' in params:
        cval = params['cval']
    if 'dpi' in params:
        dpi = params['dpi']
    if 'tl' in params:
        tl = params['tl']
    if 'tw' in params:
        tw = params['tw']
    if 'ts' in params:
        ts = params['ts']
    if 'tsi' in params:
        tsi = params['tsi']
    
    w,h = (12, 12)
    dx, dy = (10, 10)
    if 'dx' in params:
        dx = params['dx']
    if 'dy' in params:
        dy = params['dy']
    if 'w' in params:
        w = params['w']
    if 'h' in params:
        h = params['h']
        
    nAt, nGt = (len(col_labels), len(row_labels))
    fig = plt.figure(figsize=(w,h), dpi=dpi)
    ax = fig.add_axes([70.0/w/dpi, 54.0/h/dpi, 1-2*70.0/w/dpi, 1-2*54.0/h/dpi])
    extent = [0, nAt*dx, 0, nGt*dy]

    cvals  = [-1, -0.7, -0.4, 0, 0, 0.2, 1]
    clrs = ["#210B61","#0B614B","#04B45F", "#D8F781", "#F2F5A9", "red", "#DF0101"]
    norm=plt.Normalize(min(cvals),max(cvals))
    tuples = list(zip(map(norm,cvals), clrs))
    cmap = colors.LinearSegmentedColormap.from_list("BGYR1", tuples)
    plt.register_cmap(cmap=cmap)

    cvals  = [-1, -0.7, -0.4, 0, 0, 0.8, 1]
    clrs = ["#210B61","#0B614B","#04B45F", "#D8F781", "#F2F5A9", "red", "#DF0101"]
    norm=plt.Normalize(min(cvals),max(cvals))
    tuples = list(zip(map(norm,cvals), clrs))
    cmap = colors.LinearSegmentedColormap.from_list("BGYR2", tuples)
    plt.register_cmap(cmap=cmap)

    im = ax.imshow(data, cmap="bwr", interpolation='nearest', vmin=-2.0, vmax=2.0, extent = extent)

    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    ax.set_xticklabels([])
    
    yticks = []
    ylabels = []
    if 'rowlabels' in params:
        for i in range(len(row_labels)):
            g = row_labels[i]
            if g in genes:
                yticks += [-dy/2 + (len(row_labels) - i) * dy]
                ylabels += [ row_labels[i] ]
    else:
        for g in genes:
            if g in row_labels:
                i = row_labels.index(g)
                yticks += [-dy/2 + (len(row_labels) - i) * dy]
                ylabels += [ row_labels[i] ]
    si = np.argsort(np.array(yticks))
    yiticks = np.array(yticks)[si]
    yoticks = np.array(yticks)[si]
    ylabels = np.array(ylabels)[si]
    sy = 5
    if 'sy' in params:
        sy = params['sy']
    for i in range(1, len(yoticks)):
        diff = yoticks[i] - yoticks[i - 1]
        if diff < sy*dy:
            yoticks[i] = yoticks[i - 1] + sy*dy
    for i in range(len(yoticks)):
        yoticks[i] = yoticks[i] + tsi
    ax.set_yticks(yiticks)
    ax.set_yticklabels([])
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right", rotation_mode="anchor")
    for edge, spine in ax.spines.items():
        spine.set_visible(False)
    ax.grid(False)
    ax.tick_params(top=False, left=True, bottom=False, right=False,
            length=tl, width=tw)
    plt.xlim(xmin=0)
    trans = blended_transform_factory(ax.transData, ax.transData)
    fx, fy =  ax.transData.transform((1, 1)) - ax.transData.transform((0, 0))
    fx = dpi/fx/72
    fy = dpi/fy/72
    print(fx, fy)
    fx = max(fx, fy)
    fy = max(fx, fy)
    oo = 2
    for i in range(len(yoticks)):
        ax.annotate(str(ylabels[i]), xy=(0.0, yoticks[i]),
                xycoords=trans,
                xytext=(-(2*tl+ts), 0), textcoords='offset points', color="black",
                fontsize=8, ha="right")
        ax.plot((-(2*tl+ts)*fx, -(tl+ts+oo)*fx, -(tl+oo)*fx, -tl*fx),
                (yoticks[i]+4*fy, yoticks[i]+4*fy, yiticks[i], yiticks[i]),
                transform=trans,
                linestyle='solid', linewidth=tw, color='black', clip_on=False)
        oo += 0.5
        if (oo > 2):
            oo = 0

    # Create colorbar
    aspect = 20
    pad_fraction = 0.5
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)
    cbar = plt.colorbar(im, cax=cax)
    cbar.ax.set_ylabel("Expression", rotation=-90, va="bottom")

    color_sch1 = ["#3B449C", "#B2509E","#EA4824"]
    color_sch1 = ["#00CC00", "#EFF51A","#EC008C", "#F7941D", "#808285",
            'cyan', 'blue', 'black', 'green', 'red']
    if 'acolor' in params:
        color_sch1 = params['acolor']

    if len(cval) > 0:
        width = axes_size.AxesX(ax, aspect=1./aspect)
        pad = axes_size.Fraction(pad_fraction, width)
        tax = divider.append_axes("top", size=width, pad=pad)
        extent = [0, nAt, 0, 5]
        tax.axis(extent)
        cmap = colors.ListedColormap(color_sch1)
        boundaries = range(len(color_sch1) + 1)
        norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)
        tax.imshow(cval, interpolation='nearest', cmap=cmap, norm=norm, extent=extent, aspect="auto")
        tax.set_xticklabels([])
        tax.set_yticklabels([])
        tax.tick_params(top=False, left=False, bottom=False, right=False)
        if 'tline' in params and params['tline'] == 1:
            tax.set_xticks(np.arange(0, nAt, 1))
            tax.grid(which='major', alpha=1, linestyle='-', linewidth='1',
                    color='black')
        else:
            tax.grid(False)

    pad = axes_size.Fraction(0.2, width)
    lax = divider.append_axes("top", size=width, pad=pad, frame_on=False)
    lax.axison = False
    lax.axis(extent)
    lax.set_xticklabels([])
    lax.set_yticklabels([])
    lax.grid(False)
    lax.tick_params(top=False, left=False, bottom=False, right=False)
    barTop(lax, atypes, color_sch1, params)

    fig.savefig(ofile, dpi=dpi)
    return ax, divider

def plotTitleBar(cval, atypes, params):
    dpi = 100
    if 'dpi' in params:
        dpi = params['dpi']
    w,h = (5, 0.8)
    if 'w' in params:
        w = params['w']
    if 'h' in params:
        h = params['h']
    color_sch1 = ["#3B449C", "#B2509E","#EA4824"]
    color_sch1 = ["#00CC00", "#EFF51A","#EC008C", "#F7941D", "#808285",
            'cyan', 'blue', 'black', 'green', 'red']
    if 'acolor' in params:
        color_sch1 = params['acolor']
    if 'cval' in params:
        cval = params['cval']

    ax = None
    if 'ax' in params:
        ax = params['ax']
    if ax is None:
        fig = plt.figure(figsize=(w,h), dpi=dpi)
        ax = fig.add_subplot(1, 1, 1)
    nAt = len(cval[0])
    extent = [0, nAt, 0, 5]
    ax.axis(extent)
    cmap = colors.ListedColormap(color_sch1)
    boundaries = range(len(color_sch1) + 1)
    norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)
    ax.imshow(cval, interpolation='nearest', cmap=cmap, \
                      norm=norm, extent=extent, aspect="auto")
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.tick_params(top=False, left=False, bottom=False, right=False)
    ax.set_xticks(np.arange(0, nAt, 1))
    ax.grid(which='major', alpha=0.2, linestyle='-', linewidth=0.5,
            color='black')
    for edge, spine in ax.spines.items():
                spine.set_visible(False)
    divider = make_axes_locatable(ax)
    width = axes_size.AxesX(ax, aspect=1./20)
    spaceAnn = 70
    widthAnn = 3
    tAnn = 1
    if 'spaceAnn' in params:
        spaceAnn = params['spaceAnn']
    if 'widthAnn' in params:
        widthAnn = params['widthAnn']
    if 'tAnn' in params:
        tAnn = params['tAnn']
    pad = axes_size.Fraction(0.1, width)
    lax = divider.append_axes("top", size="100%", pad="20%", frame_on=False)
    lax.axison = False
    lax.axis(extent)
    lax.set_xticklabels([])
    lax.set_yticklabels([])
    lax.grid(False)
    lax.tick_params(top=False, left=False, bottom=False, right=False)
    if 'atypes' in params:
        atypes = params['atypes']
    barTop(lax, atypes, color_sch1, params)
    return ax

def plotDensity(x, atypes, ax = None, color = None):
    if color is None:
        color = acolor
    df = pd.Series(x)
    for i in range(len(atypes)):
        df1 = pd.DataFrame(pd.Series(df[df == i].index),
                columns=[atypes[i]])
        if len(df1[atypes[i]]) <= 1:
            continue
        if ax is None:
            ax = df1.plot.kde(bw_method=1.0, c=color[i])
        else:
            ax = df1.plot.kde(bw_method=1.0, ax = ax, c=color[i])

    ax.set_title("Density plot")
    ax.set_xlabel("Sample rank")
    return ax

def adj_light(color, l=1, s=1):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, l * c[1])), 
		 max(0, min(1, s * c[2])))

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, h

def plotScores(data, atypes, params):
    dpi = 100
    if 'dpi' in params:
        dpi = params['dpi']
    vert = 0
    if 'vert' in params:
        vert = params['vert']
    if vert == 0:
        w,h = (2, 0.25 * len(atypes))
    else:
        w,h = (0.75 * len(atypes), 2)
    if 'w' in params:
        w = params['w']
    if 'h' in params:
        h = params['h']
    color_sch1 = ["#3B449C", "#B2509E","#EA4824"]
    color_sch1 = ["#00CC00", "#EFF51A","#EC008C", "#F7941D", "#808285",
            'cyan', 'blue', 'black', 'green', 'red']
    if 'acolor' in params:
        color_sch1 = params['acolor']
    if 'cval' in params:
        cval = params['cval']

    ax = None
    if 'ax' in params:
        ax = params['ax']
    if ax is None:
        fig = plt.figure(figsize=(w,h), dpi=dpi)
        ax = fig.add_subplot(1, 1, 1)

    meanpointprops = dict(marker='D', markerfacecolor='none', markersize=5,
                              linestyle='none')
    meanlineprops = dict(linestyle='--', linewidth=1, color='purple')
    cols = [ color_sch1[i % len(color_sch1)] for i in range(len(data)) ]
    bp = ax.boxplot(data, notch=0, sym='+', vert=vert, whis=1.5,
            showmeans=1, meanline=0, meanprops=meanpointprops, widths=0.7)
    if vert == 0:
        ax.set_yticklabels(atypes)
        ax.set_xticklabels([])
        for edge, spine in ax.spines.items():
            spine.set_visible(False)
    else:
        ax.set_xticklabels(atypes)
    ax.grid(False)
    ax.tick_params(top=False, left=True, bottom=False, right=False)

    from matplotlib.patches import Polygon,Arrow
    for i in range(len(data)):
        col1 = adj_light(cols[i], 1, 0.5)
        col2 = 'blue'
        col3 = adj_light(cols[i], 0.5, 1)
        plt.setp(bp['medians'][i], color='red')
        plt.setp(bp['boxes'][i], color=col1)
        plt.setp(bp['means'][i], markeredgecolor=col2)
        plt.setp(bp['whiskers'][2*i], color=col1)
        plt.setp(bp['whiskers'][2*i+1], color=col1)
        plt.setp(bp['caps'][2*i], color=col1)
        plt.setp(bp['caps'][2*i+1], color=col1)
        plt.setp(bp['fliers'][i], markeredgecolor=cols[i], marker='o')
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        box_coords = np.column_stack([boxX, boxY])
        ax.add_patch(Polygon(box_coords, facecolor=cols[i], alpha=0.2))
        m, h = mean_confidence_interval(data[i])
        if vert == 0:
            ax.add_patch(Arrow(m, i + 1,
                h, 0, width=1, facecolor=col3, alpha=0.4))
            ax.add_patch(Arrow(m, i + 1,
                -h, 0, width=1, facecolor=col3, alpha=0.4))
        else:
            ax.add_patch(Arrow(i + 1, m,
                0, h, width=1, facecolor=col3, alpha=0.4))
            ax.add_patch(Arrow(i + 1, m,
                0, -h, width=1, facecolor=col3, alpha=0.4))

    return ax, bp

def getGroupsMm(gene_groups):
    cfile = "/booleanfs2/sahoo/Data/SeqData/genome/Homo_sapiens.GRCh38.95.chr_patch_hapl_scaff.len.txt"
    fp = open(cfile, "r")
    hsdict = {}
    for line in fp:
        line = line.strip();
        ll = re.split("\t", line);
        hsdict[ll[0]] = ll[1]
    fp.close();

    cfile = "/booleanfs2/sahoo/Data/Sarah/CellCycle/database/mart-export-hs-mm.txt"
    fp = open(cfile, "r")
    mmdict = {}
    for line in fp:
        line = line.strip();
        ll = re.split("\t", line);
        if len(ll) > 3 and ll[0] in hsdict:
            g = hsdict[ll[0]]
            if g not in mmdict:
                mmdict[g] = []
            mmdict[g] += [ll[3]]
    fp.close();

    gene_groups_mm = []
    for s in gene_groups:
        s1 = set()
        for g in s:
            if g in mmdict:
                for k in mmdict[g]:
                    s1.add(k)
        gene_groups_mm.append(s1)
    return gene_groups_mm


def getSimpleName(gene_groups, h):
    res = []
    for s in gene_groups:
        s1 = set()
        for g in s:
            for id1 in h.getIDs(g):
                s1.add(h.getSimpleName(id1))
        res.append(s1)
    return res

def getGeneGroups(order = None, weight = None, debug = 1):
    reload(hu)
    db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
    h = hu.Hegemon(db.getDataset("PLP7"))
    h.init()
    h.initPlatform()
    h.initSurv()
    data_item = []
    with open('Supplementary/path-1.json') as data_file:
        data_item += json.load(data_file)
    with open('Supplementary/path-2.json') as data_file:
        data_item += json.load(data_file)
    cfile = "Supplementary/ibd-network-clusters.txt"
    if not os.path.isfile(cfile):
        print("Can't open file {0} <br>".format(cfile));
        exit()
    fp = open(cfile, "r")
    nodelist = {}
    nhash = {}
    for line in fp:
        line = line.strip();
        ll = line.split("\t");
        nodelist[ll[0]] = ll[2:]
        for i in ll[2:]:
            nhash[i] = ll[0];
    fp.close();
    gene_groups = []
    for i in range(len(data_item)):
        gene_groups.append(set())
        gn = data_item[i][2][0][0]
        for g in data_item[i][2]:
            gene_groups[i].add(g[0])
            if g[0] in nodelist:
                for k in nodelist[g[0]]:
                    gene_groups[i].add(k)
        for g in data_item[i][3]:
            gene_groups[i].add(g)
            if g in nodelist:
                for k in nodelist[g]:
                    gene_groups[i].add(k)
        if debug == 1:
            print(i, gn, h.getSimpleName(gn), data_item[i][0], len(gene_groups[i]))
    print([len(s) for s in gene_groups])
    if order is None:
        order = [7, 6, 5, 1];
        order = [7, 6, 5];
    gene_groups = [gene_groups[i] for i in order]
    print([len(s) for s in gene_groups])
    gene_groups = getSimpleName(gene_groups, h)
    print([len(s) for s in gene_groups])
    if weight is None:
        weight = [-3, -2, -1]
    print(weight)
    genes = readGenes("cluster-names.txt")
    return genes, weight, gene_groups

def getGeneGroups2(order = None, weight = None, debug = 1):
    reload(hu)
    db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
    h = hu.Hegemon(db.getDataset("PLP11"))
    h.init()
    h.initPlatform()
    h.initSurv()
    data_item = []
    with open('CD/info2/path-1.json') as data_file:
        data_item += json.load(data_file)
    with open('CD/info2/path-2.json') as data_file:
        data_item += json.load(data_file)
    cfile = "CD/cd-network-2-cls.txt"
    if not os.path.isfile(cfile):
        print("Can't open file {0} <br>".format(cfile));
        exit()
    fp = open(cfile, "r")
    nodelist = {}
    nhash = {}
    for line in fp:
        line = line.strip();
        ll = line.split("\t");
        nodelist[ll[0]] = ll[2:]
        for i in ll[2:]:
            nhash[i] = ll[0];
    fp.close();
    gene_groups = []
    for i in range(len(data_item)):
        gene_groups.append(set())
        gn = data_item[i][2][0][0]
        for g in data_item[i][2]:
            gene_groups[i].add(g[0])
            if g[0] in nodelist:
                for k in nodelist[g[0]]:
                    gene_groups[i].add(k)
        for g in data_item[i][3]:
            gene_groups[i].add(g)
            if g in nodelist:
                for k in nodelist[g]:
                    gene_groups[i].add(k)
        if debug == 1:
            print(i, gn, h.getSimpleName(gn), data_item[i][0], len(gene_groups[i]))
    print([len(s) for s in gene_groups])
    if order is None:
        order = [7, 6, 5, 1];
        order = [7, 6, 5];
    gene_groups = [gene_groups[i] for i in order]
    print([len(s) for s in gene_groups])
    gene_groups = getSimpleName(gene_groups, h)
    print([len(s) for s in gene_groups])
    if weight is None:
        weight = [-3, -2, -1]
    print(weight)
    genes = readGenes("cluster-names.txt")
    return genes, weight, gene_groups

def getC4genes():
    cfile = "Supplementary/c4-ranks.txt"
    fp = open(cfile, "r")
    g4 = set()
    for line in fp:
        line = line.strip();
        ll = line.split(" ");
        if float(ll[2]) > 0:
            continue
        g4.add(ll[1])
    fp.close();
    return g4

def processGeneGroups(ana, l1, wt1, debug = 0, fthr = None):
    ana.orderData(l1, wt1); print("ROC-AUC", ana.getMetrics())
    actual = [1 if ana.aval[i] >= 1 else 0 for i in ana.i1]
    score = [ana.f_ranks[i - ana.h.start] for i in ana.i1]
    thr = hu.getThrData(score)
    nm = (np.max(ana.f_ranks) - np.min(ana.f_ranks))/16
    if fthr is None:
        fthr = thr[0]
    if fthr == "thr0":
        fthr = thr[0] - nm
    if fthr == "thr2":
        fthr = thr[0] + nm
    if fthr == "thr3":
        fthr = thr[0] + 3 * nm
    print(thr)
    print(nm, fthr)
    predicted = [1 if ana.f_ranks[i - ana.h.start] >= fthr else 0 for i in ana.i1]
    c_dict = {}
    for i in ana.order:
        c_dict[i] = ana.f_ranks[i - ana.h.start]
        c_dict[i] = 0
        if ana.f_ranks[i - ana.h.start] >= fthr:
            c_dict[i] = 1
    fpr, tpr, thrs = roc_curve(actual, score)
    roc_auc = auc(fpr, tpr)
    if debug == 1:
        print(actual)
        print(predicted)
        print(score)
        print(thr[0], thr, nm)
        print("ROC-AUC", roc_auc)
        print('Accuracy', accuracy_score(actual, predicted))
        print(classification_report(actual, predicted, target_names=ana.atypes))
    return c_dict, fpr, tpr, roc_auc

class IBDAnalysis:

    def __init__(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.normal = []
        self.uc = []
        self.cd = []
        self.otype = 0
        self.axes = []

    def addAxes(self, ax):
        self.axes += [ax]


    def printInfo(self):
        print(self.h.rdataset.getName() + " (n = " + str(self.h.getNum()) + ")")
        url = "http://hegemon.ucsd.edu/Tools/explore.php?key=polyps&id="
        if self.dbid.startswith("CRC"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=colon&id="
        print(len(self.order), len(self.normal), \
                len(self.uc), len(self.cd), self.h.getSource(), \
                url + self.dbid, self.dbid)

    def convertMm(self, gene_groups, genes):
        self.gene_groups = getGroupsMm(gene_groups)
        self.genes = getGroupsMm([genes])[0]

    def orderData(self, gene_groups, weight):
        self.col_labels = [self.h.headers[i] for i in self.order]
        #ranks, row_labels, expr = getRanks(gene_groups, h)
        ranks, row_labels, row_ids, row_numhi, expr = getRanks2(gene_groups,
                self.h)
        i1 = getOrder(self.order, self.h.start, ranks, weight)
        index = np.array([i - self.h.start for i in i1])
        self.cval = np.array([[self.aval[i] for i in i1]])
        #self.data = np.array(expr)[:,index]
        self.ind_r = np.array(sorted(range(len(row_labels)),
            key=lambda x: (row_numhi[x][0], row_numhi[x][1])))
        row_labels = [row_labels[i] for i in self.ind_r]
        row_ids = [row_ids[i] for i in self.ind_r]
        self.data = np.array([expr[i] for i in self.ind_r])
        if self.data.shape[0] > 0:
            self.data = self.data[:,index]
        self.f_ranks = mergeRanks(self.h.aRange(), self.h.start, ranks, weight)
        self.ranks = ranks
        self.row_labels = row_labels
        self.row_ids = row_ids
        self.row_numhi = row_numhi
        self.expr = expr
        self.i1 = i1
        self.index = index
        self.otype = 1

    def orderData2(self, gene_groups, weight):
        self.col_labels = [self.h.headers[i] for i in self.order]
        #ranks, row_labels, expr = getRanks(gene_groups, h)
        ranks, row_labels, row_ids, row_numhi, expr = getRanks3(gene_groups,
                self.h, self.order)
        self.f_ranks = mergeRanks2(self.order, ranks, weight)
        i1 = [self.order[i] for i in np.argsort(self.f_ranks)]
        index = np.argsort(self.f_ranks)
        self.cval = np.array([[self.aval[i] for i in i1]])
        #self.data = np.array(expr)[:,index]
        self.ind_r = np.array(sorted(range(len(row_labels)),
            key=lambda x: (row_numhi[x][0], row_numhi[x][1])))
        row_labels = [row_labels[i] for i in self.ind_r]
        row_ids = [row_ids[i] for i in self.ind_r]
        self.data = np.array([expr[i] for i in self.ind_r])
        if self.data.shape[0] > 0:
            self.data = self.data[:,index]
        self.ranks = ranks
        self.row_labels = row_labels
        self.row_ids = row_ids
        self.row_numhi = row_numhi
        self.expr = expr
        self.i1 = i1
        self.index = index
        self.otype = 2

    def saveData(self, ofile1, ofile2):
        saveHeatmapData(ofile1, self.row_labels, self.row_numhi,
                self.row_ids, self.index, self.expr)
        saveCData(ofile2, self.h, self.i1, self.f_ranks)

    def readData(self, file1, file2):
        row_labels, row_numhi, row_ids, data = readHeatmapData(file1)
        i1, f_ranks = readCData(file2)
        self.ranks = ranks
        self.row_labels = row_labels
        self.row_ids = row_ids
        self.row_numhi = row_numhi
        self.expr = expr
        self.i1 = i1
        self.f_ranks = f_ranks

    def printHeatmap(self, ofile, genes, params):
        i1 = self.i1
        f_ranks = self.f_ranks
        self.params = {'genes': genes,'atypes': self.atypes,'cval': self.cval}
        self.params.update(params)
        ax, divider = plotHeatmap(ofile, self.data, self.col_labels, 
                self.row_labels, self.params)
        actual = [1 if self.aval[i] >= 1 else 0 for i in i1]
        if "actual" in self.params:
            actual = self.params["actual"]
        thr = hu.getThrData(f_ranks)
        nm = (np.max(f_ranks) - np.min(f_ranks))/15
        predicted = [1 if f_ranks[i - self.h.start] >= thr[0] - nm else 0 for i in i1]
        if "thr" in self.params:
            if self.params["thr"] == 1:
                predicted = [1 if f_ranks[i - self.h.start] >= thr[0] else 0 for i in i1]  
            if self.params["thr"] == 2:
                predicted = [1 if f_ranks[i - self.h.start] >= thr[0] + nm else 0 for i in i1]  
        if "predicted" in self.params:
            predicted = self.params["predicted"]
        data_list = {'x' : predicted}
        df = pd.DataFrame(data_list)
        df['y'] = pd.Series(np.array(actual))
        target_names = ['Normal', 'IBD']
        score=convertScore(predicted)
        print(list(actual))
        print(list(predicted))
        ax1 = printReport(actual, predicted, score, target_names)
        tab = pd.crosstab(df.y > 0, df.x > 0)
        print(tab)
        print(fisher_exact(tab))
        print('Fisher Exact pvalue =', fisher_exact(tab)[1])
        return ax, ax1, divider

    def printHeatmap2(self, ofile, genes, params):
        i1 = self.i1
        f_ranks = self.f_ranks
        self.params = {'genes': genes,'atypes': self.atypes,'cval': self.cval}
        self.params.update(params)
        ax, divider = plotHeatmap(ofile, self.data, self.col_labels, 
                self.row_labels, self.params)
        actual = [1 if self.aval[i] >= 1 else 0 for i in i1]
        if "actual" in self.params:
            actual = self.params["actual"]
        thr = hu.getThrData([f_ranks[i - self.h.start] for i in i1])
        nm = (np.max(f_ranks) - np.min(f_ranks))/15
        predicted = [1 if f_ranks[i - self.h.start] >= thr[0] - nm else 0 for i in i1]
        if "thr" in self.params:
            if self.params["thr"] == 1:
                predicted = [1 if f_ranks[i - self.h.start] >= thr[0] else 0 for i in i1]  
            if self.params["thr"] == 2:
                predicted = [1 if f_ranks[i - self.h.start] >= thr[0] + nm else 0 for i in i1]  
        if "predicted" in self.params:
            predicted = self.params["predicted"]
        data_list = {'x' : predicted}
        df = pd.DataFrame(data_list)
        df['y'] = pd.Series(np.array(actual))
        target_names = ['Normal', 'IBD']
        score=convertScore(predicted)
        print(list(actual))
        print(list(predicted))
        ax1 = printReport(actual, predicted, score, target_names)
        tab = pd.crosstab(df.y > 0, df.x > 0)
        print(tab)
        print(fisher_exact(tab))
        print('Fisher Exact pvalue =', fisher_exact(tab)[1])
        return ax, ax1, divider

    def printHeatmap3(self, ofile, genes, params):
        i1 = self.i1
        f_ranks = self.f_ranks
        self.params = {'genes': genes,'atypes': self.atypes,'cval': self.cval}
        self.params.update(params)
        ax, divider = plotHeatmap(ofile, self.data, self.col_labels, 
                self.row_labels, self.params)
        actual = [1 if self.aval[i] >= 1 else 0 for i in i1]
        if "actual" in self.params:
            actual = self.params["actual"]
        thr = hu.getThrData(f_ranks)
        nm = (np.max(f_ranks) - np.min(f_ranks))/15
        i2 = np.argsort(f_ranks)
        predicted = [1 if f_ranks[i] >= thr[0] - nm else 0 for i in i2]
        if "thr" in self.params:
            if self.params["thr"] == 1:
                predicted = [1 if f_ranks[i] >= thr[0] else 0 for i in i2]  
            if self.params["thr"] == 2:
                predicted = [1 if f_ranks[i] >= thr[0] + nm else 0 for i in i2]  
        if "predicted" in self.params:
            predicted = self.params["predicted"]
        data_list = {'x' : predicted}
        df = pd.DataFrame(data_list)
        df['y'] = pd.Series(np.array(actual))
        target_names = ['Normal', 'IBD']
        score=convertScore(predicted)
        print(list(actual))
        print(list(predicted))
        ax1 = printReport(actual, predicted, score, target_names)
        tab = pd.crosstab(df.y > 0, df.x > 0)
        print(tab)
        print(fisher_exact(tab))
        print('Fisher Exact pvalue =', fisher_exact(tab)[1])
        return ax, ax1, divider

    def printTitleBar(self, params):
        self.params = {'atypes': self.atypes,'cval': self.cval}
        self.params.update(params)
        ax = plotTitleBar(self.params['cval'], \
                self.params['atypes'], self.params)
        return ax

    def getMetrics(ana, actual = None, ahash = None):
        if actual is None:
            actual = [ana.aval[i] for i in ana.i1]
        if ana.otype == 2:
            score = [ana.f_ranks[i] for i in ana.i1]
        else:
            if ahash is None:
                score = [ana.f_ranks[i - ana.h.start] for i in ana.i1]
            else:
                score = [ana.f_ranks[ana.i1[i] - ana.h.start] \
                        for i in range(len(ana.i1)) \
                        if ana.cval[0][i] in  ahash ]
        res = None
        fpr, tpr, thrs = roc_curve(actual, score, pos_label=1)
        roc_auc = auc(fpr, tpr)
        res = "%.2f" % roc_auc
        return res

    def getMetrics2(ana, actual, ahash = None, fthr = None):
        if ana.otype == 2:
            score = [ana.f_ranks[i] for i in ana.i1]
        else:
            if ahash is None:
                score = [ana.f_ranks[i - ana.h.start] for i in ana.i1]
            else:
                score = [ana.f_ranks[ana.i1[i] - ana.h.start] \
                        for i in range(len(ana.i1)) \
                        if ana.cval[0][i] in  ahash ]
        res = [None, None, None]
        for i in range(3):
            fpr, tpr, thrs = roc_curve(actual, score, pos_label=i)
            roc_auc = auc(fpr, tpr)
            if roc_auc < 0.5:
                roc_auc = 1 - roc_auc
            res[i] = "%.2f" % roc_auc
        thr = hu.getThrData(score[:-2])
        nm = (np.max(score) - np.min(score))/15
        if fthr is None:
            fthr = thr[0]
        if fthr == "thr0":
            fthr = thr[0] - nm
        if fthr == "thr2":
            fthr = thr[0] + nm
        if fthr == "thr3":
            fthr = thr[0] + 3 * nm
        if fthr < np.min(score) or fthr > np.max(score):
            fthr = thr[0]
        predicted = [1 if ana.f_ranks[i - ana.h.start] >= fthr else 0 for i in ana.i1]
        print(thr)
        print(nm, fthr)
        print(score)
        print(actual)
        print(predicted)
        x0 = [i for i in predicted if i == 0]
        x1 = [i for i in predicted if i == 1]
        if len(x0) == 0:
            predicted[0] = 0
        if len(x1) == 0:
            predicted[0] = 1
        res[1] = "%.2f" % accuracy_score(actual, predicted)
        data_list = {'x' : predicted}
        df = pd.DataFrame(data_list)
        df['y'] = pd.Series(np.array(actual))
        tab = pd.crosstab(df.y > 0, df.x > 0)
        res[2] = "%.3g" % fisher_exact(tab)[1]
        return res

    def densityPlot(self, ax=None, color = None):
        if (color is None):
            color = acolor
        ax = plotDensity(self.cval[0], self.atypes, ax, color)
        self.addAxes(ax)
        return ax

    def getScores(ana, ahash = None):
        lval = [[] for i in ana.atypes]
        cval = ana.cval[0]
        if ana.otype == 2:
            score = [ana.f_ranks[i] for i in ana.i1]
        else:
            if ahash is None:
                score = [ana.f_ranks[i - ana.h.start] for i in ana.i1]
            else:
                score = [ana.f_ranks[ana.i1[i] - ana.h.start] \
                        for i in range(len(ana.i1)) \
                        if ana.cval[0][i] in  ahash ]
                cval = [ana.cval[0][i] \
                        for i in range(len(ana.i1)) \
                        if ana.cval[0][i] in  ahash ]
            for i in range(len(cval)):
                lval[cval[i]] += [score[i]]
        return lval, score

    def printAllPvals(self, ahash = None, params = None):
        lval, score = self.getScores(ahash=ahash)
        cAllPvals(lval, self.atypes)

    def getPvals(self, label, ahash = None, params = None):
        lval, score = self.getScores(ahash=ahash)
        actual = [self.aval[i] for i in self.i1]
        if "actual" in self.params:
            actual = self.params["actual"]
        thr = hu.getThrData(score)
        nm = (np.max(score) - np.min(score))/15
        i2 = np.argsort(score)
        predicted = [1 if score[i] >= thr[0] - nm else 0 for i in i2]
        if "thr" in self.params:
            if self.params["thr"] == 1:
                predicted = [1 if score[i] >= thr[0] else 0 for i in i2]  
            if self.params["thr"] == 2:
                predicted = [1 if score[i] >= thr[0] + nm else 0 for i in i2]  
        if "predicted" in self.params:
            predicted = self.params["predicted"]
        pval_fe = getFisher(predicted, actual);
        t, p = ttest_ind(lval[0],lval[label], equal_var=False)
        fpr, tpr, thrs = roc_curve(actual, score, pos_label=label)
        roc_auc = auc(fpr, tpr)
        desc = "%.3g, %.3g, %.3g" % (pval_fe, p, roc_auc)
        return desc

    def printScores(self, ahash = None, params = None):
        self.params = {'atypes': self.atypes,'cval': self.cval}
        if params is not None:
            self.params.update(params)
        lval, score = self.getScores(ahash=ahash)
        ax,bp = plotScores(lval, self.params['atypes'], self.params)
        ax.text(ax.get_xlim()[1], ax.get_ylim()[1], self.h.getSource(),
                horizontalalignment='left', verticalalignment='center')
        if ('vert' not in self.params or  self.params['vert'] == 0):
            for i in range(1, len(lval)):
                desc = self.getPvals(i, ahash, params)
                ax.text(ax.get_xlim()[1], i + 1, desc,
                        horizontalalignment='left', verticalalignment='center')
        self.addAxes(ax)
        self.addAxes(bp)
        return ax,bp

    def printViolin(self, ahash = None, params = None):
        self.params = {'atypes': self.atypes,'cval': self.cval}
        if params is not None:
            self.params.update(params)
        lval, score = self.getScores(ahash=ahash)
        df = pd.DataFrame()
        df['score'] = score
        atypes = [str(self.atypes[i]) + "("+str(len(lval[i]))+")"
                for i in range(len(self.atypes))]
        df['category'] = [atypes[self.aval[i]] for i in self.i1]
        m1 = []
        pvals = []
        for i in range(1, len(lval)):
            if len(lval[i]) <= 0:
                m1 += [0]
                pvals += [""]
                continue
            m1 += [max(lval[i]) + (max(lval[i]) - min(lval[i])) * 0.1]
            t, p = ttest_ind(lval[0],lval[i], equal_var=False)
            if (p < 0.05):
                pvals += ["p=%.3g" % p]
            else:
                pvals += [""]
        params = self.params
        dpi = 100
        if 'dpi' in params:
            dpi = params['dpi']
        w,h = (1.5 * len(self.atypes), 4)
        if 'w' in params:
            w = params['w']
        if 'h' in params:
            h = params['h']
        color_sch1 = acolor
        if 'acolor' in params:
            color_sch1 = params['acolor']
        ax = None
        if 'ax' in params:
            ax = params['ax']
        if ax is None:
            fig,ax = plt.subplots(figsize=(w,h), dpi=dpi)
        sns.set()
        sns.set_style("white")
        sns.set_style({'xtick.color':'.5', 'ytick.color':'.5', 'axes.labelcolor': '.5'})
        sns.set_context("notebook")
        sns.set_palette([adj_light(c, 1.5, 1) for c in color_sch1])
        width = 1
        height = 1
        if 'width' in params:
            width = params['width']
        if 'vert' in params and params['vert'] == 1:
            ax = sns.violinplot(x="category", y="score", inner='quartile',
                    linewidth=0.5, width=width, ax = ax, data=df,
                    order = atypes)
            ax = sns.swarmplot(x="category", y="score", color = 'blue', alpha=0.2,
                    ax=ax, data=df, order = atypes)
            ax.set_xlabel("")
            pos = range(len(atypes))
            for tick,label in zip(pos[1:],ax.get_xticklabels()[1:]):
                ax.text(pos[tick], m1[tick - 1], pvals[tick - 1],
                        horizontalalignment='center', size=12,
                        color='0.3')
            ax.yaxis.grid(True, clip_on=False)
        else:
            ax = sns.violinplot(x="score", y="category", inner='quartile',
                    linewidth=0.5, width=width, ax = ax, data=df,
                    order = atypes)
            ax = sns.swarmplot(x="score", y="category", color = 'blue', alpha=0.2,
                    ax=ax, data=df, order = atypes)
            ax.set_ylabel("")
            pos = range(len(atypes))
            for tick,label in zip(pos[1:],ax.get_yticklabels()[1:]):
                ax.text(m1[tick - 1], pos[tick]-0.5, pvals[tick - 1],
                        horizontalalignment='center', size=12,
                        color='0.3')
            ax.xaxis.grid(True, clip_on=False)
        return ax

    def printGene(self, name, ahash = None, params = None):
        self.params = {'atypes': self.atypes, 'vert':1}
        if params is not None:
            self.params.update(params)
        id1 = self.h.getBestID(self.h.getIDs(name).keys())
        expr = self.h.getExprData(id1)
        if expr is None:
            print("Not Found")
            return None, None
        lval = [[] for i in self.params['atypes']]
        aval = self.aval
        if ahash is not None:
            aval = [ahash[i] if i in ahash else None for i in self.atype]
        for i in self.h.aRange():
            if aval[i] is None:
                continue
            lval[aval[i]] += [float(expr[i])]
        ax,bp = plotScores(lval, self.params['atypes'], self.params)
        if self.params['vert'] == 0:
            ax.text(ax.get_xlim()[1], 1, self.h.getSource(),
                    horizontalalignment='left', verticalalignment='center')
        else:
            title = self.h.rdataset.getName() + " (" + self.h.getSource() + "; n = " + str(self.h.getNum()) + ")"
            ax.set_title(title)
            ax.set_ylabel(self.h.getSimpleName(id1))
        self.addAxes(ax)
        self.addAxes(bp)
        cAllPvals(lval, self.params['atypes'])
        return ax,bp

    def getROCspecific(ana, m=0, n=1):
        actual = [ana.aval[i] for i in ana.i1
                if ana.aval[i] == m or ana.aval[i] == n]
        score = [ana.f_ranks[i - ana.h.start] for i in ana.i1
                if ana.aval[i] == m or ana.aval[i] == n]
        fpr, tpr, thrs = roc_curve(actual, score, pos_label=n)
        df = pd.DataFrame()
        df['fpr'] = fpr
        df['tpr'] = tpr
        df['thrs'] = thrs
        df['group'] = "-".join([str(m), str(n)])
        return df

    def getROCAUCspecific(ana, m=0, n=1):
        actual = [ana.aval[i] for i in ana.i1
                if ana.aval[i] == m or ana.aval[i] == n]
        score = [ana.f_ranks[i - ana.h.start] for i in ana.i1
                if ana.aval[i] == m or ana.aval[i] == n]
        fpr, tpr, thrs = roc_curve(actual, score, pos_label=n)
        roc_auc = auc(fpr, tpr)
        return "%.2f" % roc_auc

    def getROCAUC(ana):
        res = []
        for k in range(1, len(ana.atypes)):
            v = ana.getROCAUCspecific(0, k)
            res += [v]
        return ",".join(res)

    def getROC(ana):
        res = None
        for k in range(1, len(ana.atypes)):
            v = ana.getROCspecific(0, k)
            if (res is None):
                res = v
            else:
                res = pd.concat([res, v], ignore_index=True, sort=False)
        return res

    def getTPvals(ana):
        res = []
        lval, score = ana.getScores()
        for k in range(1, len(ana.atypes)):
            t, p = ttest_ind(lval[0],lval[k], equal_var=False)
            res += ["%.3g" % p]
        return ",".join(res)

    def getStats(self, l1, wt1, annotation=[]):
        src = re.sub(" .*", "", self.h.getSource())
        species = annotation[1]
        if species == 'Hs' or species == 'Rm' :
            self.orderData(l1, wt1)
        else:
            genes = []
            self.convertMm(l1, genes)
            self.orderData(self.gene_groups, wt1)
        roc = self.getROCAUC()
        p = self.getTPvals()
        lval, score = self.getScores()
        n1 = np.sum([len(lval[i])for i in range(1, len(lval))])
        return [src, roc, p, len(lval[0]), n1] + annotation
        
    def getSurvival(self, dbid = "CRC35.3"):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = dbid
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("status")
        atypes = ['Censor', 'Relapse']
        ahash = {"0": 0, "1":1}
        aval = [ahash[atype[i]] if atype[i] in ahash \
                        else None for i in range(self.h.getNum()+2)]
        st0 = [ i for i in self.h.aRange() if aval[i] == 0]
        st1 = [ i for i in self.h.aRange() if aval[i] == 1]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = st0 + st1
        self.printInfo()

    def printSurvival(self, fthr = None, pG = None, genex = "CDX2",
            ct = None, ax = None):
        f_ranks = self.f_ranks
        order = self.order
        thr = hu.getThrData(f_ranks)
        nm = (np.max(f_ranks) - np.min(f_ranks))/16
        if fthr is None:
            fthr = thr[0]
        if fthr == "thr0":
            fthr = thr[0] - nm
        if fthr == "thr2":
            fthr = thr[0] + nm
        if fthr == "thr3":
            fthr = thr[0] + 3 * nm
        print(thr)
        print(nm, fthr)
        g1 = [i for i in order if f_ranks[i - self.h.start] < fthr]
        g2 = [i for i in order if f_ranks[i - self.h.start] >= fthr]
        if pG is None:
            pG = [ ["Low", "green", g1], ["High", "red", g2]]
        obj = hu.getHegemonPatientData(self.dbid, 'time')
        time = obj[1]
        obj = hu.getHegemonPatientData(self.dbid, 'status')
        status = obj[1]
        if ct is not None:
            time, status = hu.censor(time, status, ct)
        sax = hu.survival(time, status, pG, ax)
        df = pd.DataFrame()
        df["f_ranks"] = pd.to_numeric(pd.Series(f_ranks))
        e = self.h.getExprData(genex)
        df[genex] = pd.to_numeric(pd.Series(e[2:]))
        ax = df.plot.scatter(x=genex, y='f_ranks')
        return sax, ax

    def printCoxSurvival(ana, fthr = None):
        score = [ana.f_ranks[i - ana.h.start] for i in ana.h.aRange()]
        thr = hu.getThrData(score)
        nm = (np.max(ana.f_ranks) - np.min(ana.f_ranks))/16
        if fthr is None:
            fthr = thr[0]
        if fthr == "thr0":
            fthr = thr[0] - nm
        if fthr == "thr2":
            fthr = thr[0] + nm
        if fthr == "thr3":
            fthr = thr[0] + 3 * nm
        print(thr)
        print(nm, fthr)
        c1 = ["", ""] + [1 if ana.f_ranks[i - ana.h.start] >= fthr else 0 for i in
                ana.h.aRange()]
        obj = hu.getHegemonPatientData(ana.dbid, 'time')
        time = obj[1]
        obj = hu.getHegemonPatientData(ana.dbid, 'status')
        status = obj[1]
        mdata = {"time": time, "status": status, "c1" : c1}
        df = pd.DataFrame(mdata)
        df.drop([0, 1], inplace=True)
        df.replace('', np.nan, inplace=True)
        df.dropna(inplace=True)

        import rpy2.robjects as ro
        ro.r('library(survival)')
        ro.r("time <- c(" + ",".join(df['time']) + ")")
        ro.r("status <- c(" + ",".join(df['status']) + ")")
        ro.r("c1 <- c(" + ",".join([str(i) for i in df['c1']]) + ")")
        ro.r('x <- coxph(Surv(time, status) ~ c1)')
        ro.r('s <- summary(x)')
        print(ro.r('s'))
        t = [float(i) for i in df['time']]
        s = [int(i) for i in df['status']]
        g = [int(i) for i in df['c1']]
        from lifelines.statistics import multivariate_logrank_test
        from matplotlib.legend import Legend
        res = multivariate_logrank_test(t, g, s)
        print("p = %.2g" % res.p_value)

    def getJSTOM(self):
        self.getSurvival("CRC35.3")

    def getBos(self):
        self.getSurvival("BC20")

    def getPeters(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP7"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c clinical condition")
        atypes = ['Normal', 'UC', 'CD']
        ahash = {"control": 0, "Ulcerative Colitis":1, "Crohn's disease":2}
        if (tn == 2):
            atypes = ['Normal', 'IBD']
            ahash = {"control": 0, "Ulcerative Colitis":1, "Crohn's disease":1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.normal = [ i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [ i for i in self.h.aRange() if aval[i] == 1]
        self.cd = [ i for i in self.h.aRange() if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getNoble(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP6"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c disease")
        atypes = ['Normal', 'UC']
        ahash = {}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.normal = [ i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [ i for i in self.h.aRange() if aval[i] == 1]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getArijs2018(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP10"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        rtype = self.h.getSurvName("c induction therapy_maintenance therapy")
        rtypes = ['CO', 'PLAC', 'VDZ', 'IFX']
        ahash = {'plac_plac':1, 'vdz_vdz4w':2, 'vdz_plac':2, 'vdz_vdz8w':2,
                'vdz4w':2}
        for i in range(len(rtypes)):
            ahash[rtypes[i]] = i
        rval = [ahash[i] if i in ahash else None for i in rtype]
        atype = self.h.getSurvName("c Response")
        atypes = ['Control', 'UC R', 'UC NR', 'Active UC', 'UC other']
        ahash = {}
        if (tn == 2):
            atypes = ['Normal', 'IBD']
            ahash = {'Control':0, 'UC R':1, 'UC NR':1, 'Active UC':1,
                    'UC other':1}
        if (tn == 3 or tn == 4):
            atypes = ['R', 'NR']
            ahash = {'UC R':0, 'UC NR':1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        tissue = self.h.aRange()
        if (tn == 4):
            tissue = [i for i in self.h.aRange() if rval[i] == 3]
        self.normal = [ i for i in tissue if aval[i] == 0]
        self.uc = [ i for i in tissue if aval[i] > 0]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.rval = rval
        self.rtype = rtype
        self.rtypes = rtypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getWu2007(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP12"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c type")
        atypes = ['N', 'UC-un', 'CD-un', 'IC-un', 'INF', 
                'UC-Aff', 'CD-Aff', 'IC-Aff']
        atypes = ['N', 'UC', 'CD']
        ahash = {'N': 0, 'UC-Aff': 1, 'CD-Aff': 2}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'N': 0, 'UC-Aff': 1, 'CD-Aff': 1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [ i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [i for i in self.h.aRange() if aval[i] == 1]
        self.cd = [i for i in self.h.aRange() if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getVancamelbeke(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP16"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c src1")
        atypes = ['N', 'UC', 'CD']
        ahash = {'Biopsy from inflamed colonic mucosa of active UC patient':1,
                'Biopsy from inflamed colonic mucosa of active CD patient':2,
                'Biopsy from normal colonic mucosa of control individual':0}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'Biopsy from inflamed colonic mucosa of active UC patient':1,
                    'Biopsy from inflamed colonic mucosa of active CD patient':1,
                    'Biopsy from normal colonic mucosa of control individual':0}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [ i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [i for i in self.h.aRange() if aval[i] == 1]
        self.cd = [i for i in self.h.aRange() if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getDePreter(self, t1 = 1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP24"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c src1")
        atypes = ['UC Rp', 'UC Rb', 'UC p', 'UC b']
        ahash = {'Colonic mucosal biopsy from UC patient in remission before probiotics intake':1,
                'Colonic mucosal biopsy from UC patient in remission before placebo intake':0,
                'Colonic mucosal biopsy from active UC patient before placebo intake':2,
                'Colonic mucosal biopsy from active UC patient before probiotics intake':3}
        if t1 == 2:
            atypes = ['UC R', 'UC']
            ahash = {'Colonic mucosal biopsy from UC patient in remission before placebo intake':0,
                    'Colonic mucosal biopsy from active UC patient before placebo intake':1}
        if t1 == 3:
            atypes = ['UC R', 'UC']
            ahash = {'Colonic mucosal biopsy from UC patient in remission before probiotics intake':0,
                    'Colonic mucosal biopsy from active UC patient before probiotics intake':1}
        if t1 == 4:
            atypes = ['UC R', 'UC']
            ahash = {'Colonic mucosal biopsy from UC patient in remission before probiotics intake':0,
                    'Colonic mucosal biopsy from UC patient in remission before placebo intake':0,
                    'Colonic mucosal biopsy from active UC patient before placebo intake':1,
                    'Colonic mucosal biopsy from active UC patient before probiotics intake':1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [ i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [i for i in self.h.aRange() if aval[i] == 1]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getArijs2009(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP27"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        tissue = self.h.getSurvName("c tissue")
        treatment = self.h.getSurvName("c before or after first infliximab treatment")
        response = self.h.getSurvName("c response to infliximab")
        atype = self.h.getSurvName("c disease")
        atypes = ['Control', 'UC', 'CD', "R", "NR", "N/A"]
        ahash = {"Yes": 3, "No": 4, "Not applicable": 5}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {"Control": 0, "UC": 1, "CD": 1}
        if (tn == 3):
            atypes = ['R', 'NR']
            ahash = {"Yes": 0, "No": 1}
            atype = response
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        colon = [i for i in self.h.aRange() if tissue[i] == "Colon" and \
                treatment[i] != "After first infliximab treatment"]
        aval = [ahash[i] if i in ahash else None for i in atype]
        rval = [ahash[i] if i in ahash else None for i in response]
        self.normal = [ i for i in colon if aval[i] == 0]
        self.uc = [ i for i in colon if aval[i] == 1]
        self.cd = [ i for i in colon if aval[i] == 2]
        self.aval = aval
        self.rval = rval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getHaberman2014(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP11"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        ulcer = self.h.getSurvName("c deep ulcer")
        atype = self.h.getSurvName("c diagnosis")
        atypes = ['Control', 'UC', 'CD', "Yes", "No", "NA"]
        ahash = {'Not IBD':0, 'not IBD':0, "Yes": 3, "No": 4, "no":4, "NA": 5}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        rval = [ahash[i] if i in ahash else None for i in ulcer]
        self.normal = [ i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [ i for i in self.h.aRange() if aval[i] == 1]
        self.cd = [ i for i in self.h.aRange() if aval[i] == 2]
        self.aval = aval
        self.rval = rval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getHaberman2018(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP14"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c diagnosis")
        atypes = ['Control', 'CD']
        ahash = {'Non-IBD':0}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.normal = [ i for i in self.h.aRange() if aval[i] == 0]
        self.uc = []
        self.cd = [ i for i in self.h.aRange() if aval[i] == 1]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getVanhove(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP23"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        activity = self.h.getSurvName("c disease activity")
        atype = self.h.getSurvName("c disease")
        atypes = ['Control', 'UC', 'CD', 'I', 'A', 'N']
        ahash = {'ulcerative colitis':1, "Crohn's disease":2, 'control':0,
                'active':4, 'inactive':3, 'normal': 5}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'ulcerative colitis':1, "Crohn's disease":1, 'control':0}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        rval = [ahash[i] if i in ahash else None for i in activity]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        if (tn == 3):
            atypes = ['I', 'A']
            ahash = {'active':1, 'inactive':0}
            aval = [ahash[i] if i in ahash else None for i in activity]
            ahash = {'ulcerative colitis':1, "Crohn's disease":2, 'control':0}
            rval = [ahash[i] if i in ahash else None for i in atype]
            expg = [i for i in self.h.aRange() if aval[i] is not None and
                    rval[i] == 1]
            self.normal = [i for i in self.h.aRange() if aval[i] == 0]
            self.uc = [i for i in self.h.aRange() if aval[i] == 1]
            self.cd = []
        else:
            self.normal = [i for i in self.h.aRange() if aval[i] == 0]
            self.uc = [i for i in self.h.aRange() if aval[i] == 1 and rval[i] == 4]
            self.cd = [i for i in self.h.aRange() if aval[i] == 2 and rval[i] == 4]
        self.aval = aval
        self.rval = rval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getVanderGoten(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP25"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        activity = self.h.getSurvName("c disease activity")
        atype = self.h.getSurvName("c disease")
        atypes = ['Control', 'UC', 'CD', 'I', 'A', 'NA']
        ahash = {'control':0,
                'active':4, 'inactive':3, 'not applicable': 5}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'control':0, 'UC':1, 'CD':1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        rval = [ahash[i] if i in ahash else None for i in activity]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        if (tn == 3):
            atypes = ['I', 'A']
            ahash = {'active':1, 'inactive':0}
            aval = [ahash[i] if i in ahash else None for i in activity]
            ahash = {'control':0, 'UC':1, 'CD':2}
            rval = [ahash[i] if i in ahash else None for i in atype]
            expg = [i for i in self.h.aRange() if aval[i] is not None and
                    rval[i] == 1]
            self.normal = [i for i in self.h.aRange() if aval[i] == 0]
            self.uc = [i for i in self.h.aRange() if aval[i] == 1]
            self.cd = []
        else:
            self.normal = [i for i in self.h.aRange() if aval[i] == 0]
            self.uc = [i for i in self.h.aRange() if aval[i] == 1 and rval[i] == 4]
            self.cd = [i for i in self.h.aRange() if aval[i] == 2 and rval[i] == 4]
        self.aval = aval
        self.rval = rval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getPekow(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP59"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c src1")
        atypes = ['C', 'qUC', 'nUC']
        ahash = {'normal control': 0, 
                'quiescent ulcerative colitis': 1,
                'ulcerative colitis with neoplasia': 2}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'normal control': 0, 
                    'ulcerative colitis with neoplasia': 1}
        if (tn == 3):
            atypes = ['qUC', 'nUC']
            ahash = {'quiescent ulcerative colitis': 0,
                    'ulcerative colitis with neoplasia': 1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.normal = [ i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [ i for i in self.h.aRange() \
                if aval[i] == 1 or aval[i] == 2]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getGao(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP34"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c group")
        atypes = ['C', 'D', 'A', 'A/D']
        ahash = {'control': 0, 'DSS': 1, 'AOM': 2, 'AOM/DSS': 3}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.normal = [ i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [ i for i in self.h.aRange() \
                if aval[i] == 1 or aval[i] == 2 or aval[i] == 3]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getWatanabe(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP57"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c desc")
        atypes = ['UC-NonCa', 'UC-Ca']
        ahash = {}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.normal = []
        self.uc = [i for i in self.h.aRange() if aval[i] is not None]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getEColi(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "CRC141"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        series = self.h.getSurvName("c Series")
        atype = self.h.getSurvName("c treatment")
        atypes = ['C', 'K12', 'O157']
        ahash = {'control, 60min':0, 'K-12, 60min':1,
                'O157:H7, 120min':2, 'control, 90min':0,
                'O157:H7, 60min':2, 'O157:H7, 90min':2,
                'control, 120min':0, 'K12, 120min':1, 'K-12, 90min':1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.normal = [ i for i in self.h.aRange() \
                if aval[i] == 0 or aval[i] == 1 and series[i] == "GSE50040"]
        self.uc = [ i for i in self.h.aRange() \
                if aval[i] == 2 and series[i] == "GSE50040"]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getMatsuki(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "CRC142"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c src1")
        atypes = ['C', 'Lc', 'Bb']
        ahash = {'Caco-2 cells cultured with B. breve':2,
                'Caco-2 cells alone':0,
                'Caco-2 cells cultured with L. casei': 1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.normal = [ i for i in self.h.aRange() \
                if aval[i] == 0]
        self.uc = [ i for i in self.h.aRange() \
                if aval[i] == 1 or aval[i] == 2]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getArbibe(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "CRC143"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c Title")
        atypes = ['C', 'M90T', 'OspF', 'OspF C']
        ahash = {'Caco2_OspF complementation infected_Rep1':3,
                'Caco2_Non infected_Rep2':0,
                'Caco2_OspF mutant infected_Rep1':2,
                'Caco2_OspF complementation infected_Rep2':3,
                'Caco2_M90T wild type infected_Rep1':1,
                'Caco2_Non infected_Rep1':0,
                'Caco2_M90T wild type infected_Rep2':1,
                'Caco2_OspF mutant infected_Rep2':2}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.normal = [ i for i in self.h.aRange() \
                if aval[i] == 0 or aval[i] == 1]
        self.uc = [ i for i in self.h.aRange() \
                if aval[i] == 2 or aval[i] == 3]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getPereiraCaro(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "CRC141"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        series = self.h.getSurvName("c Series")
        atype = self.h.getSurvName("c agent")
        atypes = ['C', 'HTy', 'EHTy']
        ahash = {'control':0}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.normal = [ i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [ i for i in self.h.aRange() if aval[i] == 1]
        self.cd = [i for i in self.h.aRange() if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getKonnikova(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP68"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        media = self.h.getSurvName("c collection media")
        atype = self.h.getSurvName("c status")
        atypes = ['uninflamed', 'inflamed', 'D', 'R']
        ahash = {'RNA Later':3, 'DMSO':2}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        rval = [ahash[i] if i in ahash else None for i in media]
        self.normal = [ i for i in self.h.aRange() \
                if aval[i] == 0]
        self.uc = [ i for i in self.h.aRange() \
                if aval[i] == 1]
        self.cd = []
        self.aval = aval
        self.rval = rval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getKarns2019(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP69"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c disease subtype")
        atypes = ['Control', 'UC', 'CD', 'iCD']
        ahash = {'ileal CD':3, 'not IBD':0, 'colonic CD':2}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [i for i in self.h.aRange() if aval[i] == 1]
        self.cd = [i for i in self.h.aRange() if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getPeck2015(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP70"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        inflamed = self.h.getSurvName("c inflamed")
        etype = self.h.getSurvName("c Type")
        atype = self.h.getSurvName("c disease_stage")
        atypes = ['NA', 'B1', 'B2', 'B3']
        ahash = {'B1/non-strictuing, non-penetrating':1,
                'B3/penetrating':3,
                'NA':0,
                'B2/stricturing':2}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if etype[i] == "RNA-Seq"]
        self.rval = [ahash[i] if i in ahash else None for i in inflamed]
        self.normal = [i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [i for i in self.h.aRange() if aval[i] == 1]
        self.cd = [i for i in self.h.aRange() if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getCorraliza2018(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP71"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        response = self.h.getSurvName("c hsct responder")
        atype = self.h.getSurvName("c disease")
        atypes = ['Control', 'CD', "YES", "NO", "C"]
        ahash = {'Healthy non-IBD':0, "Crohn's disease (CD)":1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.rval = [ahash[i] if i in ahash else None for i in response]
        self.normal = [i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [i for i in self.h.aRange() if aval[i] == 2]
        self.cd = [i for i in self.h.aRange() if aval[i] == 1]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getArze2019(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP72"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c disease status")
        atypes = ['Control', 'UC', "CD"]
        ahash = {'Non IBD':0, 'Ulcerative Colitis':1, "Crohn's Disease":2}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.normal = [i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [i for i in self.h.aRange() if aval[i] == 1]
        self.cd = [i for i in self.h.aRange() if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getVerstockt2019(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP73"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c clinical history")
        atypes = ['R', 'NR']
        ahash = {'responder':0, 'non-responder':1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.normal = [i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [i for i in self.h.aRange() if aval[i] == 1]
        self.cd = [i for i in self.h.aRange() if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getHasler(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP74"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        tissue = self.h.getSurvName("c tissue")
        inflammation = self.h.getSurvName("c inflammation")
        atype = self.h.getSurvName("c diagnosis")
        atypes = ['Control', 'UC', "CD", 'non inflamed', 'inflamed']
        ahash = {'disease control':0, 'healthy':0}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.rval = [ahash[i] if i in ahash else None for i in inflammation]
        self.normal = [i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [i for i in self.h.aRange() if aval[i] == 1]
        self.cd = [i for i in self.h.aRange() if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getZhao2019(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP75"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c disease state")
        atype = self.h.getSurvName("c src1")
        atypes = ['Normal', "CD u", "CD i"]
        ahash = {'control':0,
                'Crohn\xe2\x80\x99s disease uninvolved':1,
                'Crohn\xe2\x80\x99s disease involved':2}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [i for i in self.h.aRange() if aval[i] == 3]
        self.cd = [i for i in self.h.aRange() if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getKugathasan2008(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP76"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c Type")
        atypes = ['Normal', 'UC', "CD", "CD i"]
        ahash = {'Healthy control':0,
                'Colon-only CD':2,
                'Ileo-colonic CD':3,
                'Ulcerative colitis':1,
                'Internal control':0}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [i for i in self.h.aRange() if aval[i] == 1]
        self.cd = [i for i in self.h.aRange() if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getZhao2015(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP77"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c disease state")
        atypes = ['Normal', 'UC I', "UC A"]
        ahash = {'ulcerative colitis inactive':1,
                'healthy control':0,
                'ulcerative colitis active':2}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [i for i in self.h.aRange() if aval[i] == 2]
        self.cd = [i for i in self.h.aRange() if aval[i] == 3]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getTang2017(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP78"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        state = self.h.getSurvName("c inflammation")
        atype = self.h.getSurvName("c disease state")
        atypes = ['Normal', 'UC', "CD", 'Inactive', "Active"]
        ahash = {'non-IBD control':0}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        rval = [ahash[i] if i in ahash else None for i in state]
        for i in self.h.aRange():
            if rval[i] == 0:
                aval[i] = 0
        expg = [i for i in self.h.aRange() if rval[i] is not None]
        self.normal = [i for i in self.h.aRange() if rval[i] == 0]
        self.uc = [i for i in self.h.aRange() if aval[i] == 1 and rval[i] == 4]
        self.cd = [i for i in self.h.aRange() if aval[i] == 2 and rval[i] == 4]
        self.aval = aval
        self.rval = rval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getCarey2008(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP79"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c Type")
        atypes = ['Normal', 'UC', "CD", "CD t"]
        ahash = {'healthy control reference':0,
                'CD':2,
                'treated CD':3,
                'UC':1,
                'Internal Control':0}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in self.h.aRange() if aval[i] == 0]
        self.uc = [i for i in self.h.aRange() if aval[i] == 1]
        self.cd = [i for i in self.h.aRange() if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getDotti2017(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP80"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        culture = self.h.getSurvName("c organoid_culture")
        atype = self.h.getSurvName("c case_phenotype")
        atypes = ['Normal', 'UC'];
        ahash = {'ulcerative colitis (UC) patient':1, 'non-IBD control':0}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        epoc = [i for i in self.h.aRange() if culture[i] == "EpOC"]
        depoc = [i for i in self.h.aRange() if culture[i] == "d-EpOC"]
        select = depoc
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.normal = [i for i in select if aval[i] == 0]
        self.uc = [i for i in select if aval[i] == 1]
        self.cd = [i for i in select if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getDenson2018(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP86"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        response = self.h.getSurvName("c week 4 remission")
        atype = self.h.getSurvName("c diagnosis")
        atypes = ['Control', 'UC', 'Yes', 'No', 'NA'];
        ahash = {'Ulcerative Colitis':1}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'Control':0, 'Ulcerative Colitis':1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        self.rval = [ahash[i] if i in ahash else None for i in response]
        self.normal = [i for i in select if aval[i] == 0]
        self.uc = [i for i in select if aval[i] == 1]
        self.cd = [i for i in select if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = self.normal + self.uc + self.cd
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getBoyd2018(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP87"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        response = self.h.getSurvName("c condition")
        atype = self.h.getSurvName("c condition")
        atypes = ['Control', 'UC', 'CD', 'Inactive', 'Active', 'NA'];
        ahash = {'CD active':2, 'CD inactive':2, 'control':0,
                'UC active':1, 'UC inactive':1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        ahash = {'CD active':4, 'CD inactive':3, 'control':5,
                'UC active':4, 'UC inactive':3}
        rval = [ahash[i] if i in ahash else None for i in response]
        self.normal = [i for i in select if aval[i] == 0]
        self.uc = [i for i in select if aval[i] == 1 and rval[i] == 4]
        self.cd = [i for i in select if aval[i] == 2 and rval[i] == 4]
        self.aval = aval
        self.rval = rval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getBreynaert2013(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP38"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c colitis group")
        atypes = ['C', 'DA', 'D1', 'D2', 'D3', 'A'];
        ahash = {'2 cycles DSS with additional recovery period':1,
                'control':0,
                '1 cycle DSS':2,
                'acute colitis':5,
                '3 cycles DSS':4,
                '2 cycles DSS':3}
        if (tn == 2):
            atypes = ['N', 'C']
            ahash = {'2 cycles DSS with additional recovery period':1,
                    'control':0,
                    '1 cycle DSS':1,
                    'acute colitis':1,
                    '3 cycles DSS':1,
                    '2 cycles DSS':1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] != 0]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getGerstgrasser2017(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP36"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c Title")
        atypes = ['C', 'DSS'];
        ahash = {'colon wt DSS rep3':1,
                'colon wt DSS rep2':1,
                'colon wt DSS rep4':1,
                'colon wt DSS rep1':1,
                'colon wt water rep3':0,
                'colon wt water rep2':0,
                'colon wt water rep1':0,
                'colon wt water rep4':0}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] != 1]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getTang2012(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP40"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c disease state")
        atypes = ['N', 'I', 'LD', 'HD', 'C'];
        ahash = {'low grade dysplasia lesion':2,
                'inflamed colorectal mucosa':1,
                'high grade dysplasia':3,
                'normal':0,
                'colorectal adenocarcinoma':4}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] != 1]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getJensen2017(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP66"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c disease")
        atypes = ['normal', 'colitis']
        ahash = {}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] != 0]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getGkouskou2016(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP84"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        tissue = self.h.getSurvName("c tissue")
        atype = self.h.getSurvName("c src1")
        atypes = ['normal', 'AD2', 'AD4', 'proximal', 'distal']
        ahash = {'UNTREATED':0, 'AOM, 4 DSS CYCLES':2, 'AOM, 2 DSS CYCLES':1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        ahash = {'proximal colon':3, 'distal colon':4}
        rval = [ahash[i] if i in ahash else None for i in tissue]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] != 1]
        self.cd = []
        self.aval = aval
        self.rval = rval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getLopezDee2012(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP49"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        gtype = self.h.getSurvName("c genotype/variation")
        atype = self.h.getSurvName("c src1")
        atypes = ['WC', 'TC', 'WD', 'TD', 'WDS', 'WDST2', 'WDST3', 'WRFK',
                'WT', 'TN']
        ahash = {'Wt, water control':0,
                'TSP-null-water':1,
                'Wt, DSS treated':2,
                'DSS-treated TSP-null':3,
                'Wt, DSS-saline treated':4,
                'Wt, DSS-TSR2 treated':5,
                'Wt, DSS-3TSR treated':6,
                'Wt, DSS-TSR+RFK treated':7}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        ahash = {'Wild-type':8, 'TSP-null':9}
        rval = [ahash[i] if i in ahash else None for i in gtype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] != 1]
        self.cd = []
        self.aval = aval
        self.rval = rval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getSarvestani2018(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "ORG16"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c Pathology")
        atypes = ['N', 'UC', 'D', 'C']
        ahash = {'Normal':0,
                'Chronic active colitis':1,
                'Colitis':1,
                'Colitis with benign strictures':1,
                'Colitis with fibrosis and treatment effect':1,
                'Low-grade dysplasia':2,
                'T3N1a':3}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] == 1]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getPlanell2013(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP88"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c src1")
        atypes = ['C', 'NI', 'Re', 'I']
        ahash = {
                'Human colon biopsies from UC patient with active disease (involved mucosa)':3,
                'Human colon biopsies from non-inflammatory control':0,
                'Human colon biopsies from UC patient with active disease (non-involved mucosa)':1,
                'Human colon biopsies from UC patient in remission (involved mucosa)':2}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0 or aval[i] == 1]
        self.uc = [i for i in expg if aval[i] == 2 or aval[i] == 3]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getLyons2018(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP89"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c inflammation level")
        atypes = ['NI', 'M', 'S']
        ahash = {'severe':2, 'moderate':1, 'non-inflamed':0}
        if (tn ==2):
            atypes = ['NI', 'I']
            ahash = {'severe':1, 'moderate':1, 'non-inflamed':0}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] == 1 or aval[i] == 2]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getFang2012(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP90"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c time")
        atypes = ['W0', 'W2', 'W4', 'W6']
        ahash = {'0 week':0, '4 weeks':2, '6 weeks':3, '2 weeks':1}
        if (tn ==2):
            atypes = ['W0', 'W2-6']
            ahash = {'0 week':0, '4 weeks':1, '6 weeks':1, '2 weeks':1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] == 1 or aval[i] == 2 or  aval[i] == 3]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getSchiering2014(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP91"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        gtype = self.h.getSurvName("c genotype")
        atype = self.h.getSurvName("c cell type")
        atypes = ['TC', 'TP', 'TN', 'WT', 'I23', 'F3']
        ahash = {
                'TCR\xce\xb2+CD4+ T cells from colon':0,
                'TCR\xce\xb2+CD4+Foxp3+ from colon lamina propria (cLP)':1,
                'TCR\xce\xb2+CD4+Foxp3+ from mesenteric lymph node (MLN)':2}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        ahash = {'wild type':3, 'Il23r-/-':4, 'Foxp3gfp':5}
        rval = [ahash[i] if i in ahash else None for i in gtype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] != 1]
        self.cd = []
        self.aval = aval
        self.rval = rval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getKremer2012(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP92"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c disease status")
        atypes = ['N', 'UC']
        ahash = {'TNBS colitis':1, 'Healthy Control':0}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] == 1]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getHo2014(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP93"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        gtype = self.h.getSurvName("c tissue")
        atype = self.h.getSurvName("c src1")
        atypes = ['N', 'UC', 'UCt', 'Sp', 'CO']
        ahash = {
                'Mock':0,
                'EA treatment/TNBS-induced colitis':2,
                'TNBS-induced colitis':1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        ahash = {'Spleen':3, 'Colon':4}
        rval = [ahash[i] if i in ahash else None for i in gtype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] != 1]
        self.cd = []
        self.aval = aval
        self.rval = rval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getDohi2014(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP94"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        gtype = self.h.getSurvName("c treated with")
        atype = self.h.getSurvName("c injected with")
        atypes = ['N', 'UC', 'T0', 'T1', 'T2', 'T3', 'T4']
        atypes = ['N', 'UC', 'T0']
        ahash = {'trinitrobenzene sulfonic acid (TNBS)':1,
                'none (na\xc3\xafve control)':0}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        ahash = {
                'none (untreated control)':2,
                '10 mg/kg control IgG2a mAb (anti-human CD20)':3,
                '0.3 mg/kg TNFR-Ig':4,
                '10 mg/kg anti-TWEAK mP2D10':5,
                'combination of TNFR-Fc (0.3 mg/kg) and anti-TWEAK mP2D10 (10 mg/kg)':6}
        rval = [ahash[i] if i in ahash else None for i in gtype]
        expg = [i for i in self.h.aRange() if rval[i] == 2]
        if (tn == 2):
            expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] == 1]
        self.cd = []
        self.aval = aval
        self.rval = rval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getDeBuhr2006(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP95"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c Title")
        atypes = ['B6-WT', 'B6-IL10', 'C3-WT', 'C3-IL10']
        ahash = {'C57BL/6J, sample 1':0,
                'C57BL/6J, 4 week old, sample 2':0,
                'C57BL/6J-Il10tm1Cgn, sample 2':1,
                'C57BL/6J-Il10tm1Cgn, sample 1':1,
                'C3H/HeJBir, sample 2':2,
                'C3H/HeJBir, sample-1':2,
                'C3H/HeJBir-Il10tm1Cgn, sample 2':3,
                'C3H/HeJBir-Il10tm1Cgn, sample 1':3}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0 or aval[i] == 2]
        self.uc = [i for i in expg if aval[i] == 1 or aval[i] == 3]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getRuss2013(self, tval):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP96"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        tissue = self.h.getSurvName("c tissue")
        atype = self.h.getSurvName("c genotype/variation")
        atypes = ['WT', 'IL10']
        ahash = {'IL10-/-':1, 'wildtype':0,
                'colon epithelium':2, 'colon':3}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        rval = [ahash[i] if i in ahash else None for i in tissue]
        select = [i for i in self.h.aRange() if rval[i] == tval]
        expg = [i for i in select if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] == 1]
        self.cd = []
        self.aval = aval
        self.rval = rval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getPunit2015(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP97"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c genotype")
        atypes = ['WT', 'TNFR2']
        ahash = {'Wildtype':0, 'TNFR2-knockout':1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] == 1]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getTam2019(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP98"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c genotype")
        atypes = ['WT', 'Tnfr1']
        ahash = {'Tnfrsf1a-/-':1, 'WT':0}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] == 1]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getLamas2018(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP99"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c src1")
        atypes = ['D0', 'D4', 'D12', 'D22', 'WT', 'KO']
        ahash = {'Cecum_Flore WT_Day0':0,
                'Cecum_Flore KO_Day0':0,
                'Cecum_Flore WT_Day4':1,
                'Cecum_Flore KO_Day4':1,
                'Cecum_Flore WT_Day12':2,
                'Cecum_Flore KO_Day12':2,
                'Cecum_Flore WT_Day22':3,
                'Cecum_Flore KO_Day22':3}
        if (tn == 2):
            atypes = ['D0-4', 'D12-22']
            ahash = {'Cecum_Flore WT_Day0':0,
                    'Cecum_Flore WT_Day4':0,
                    'Cecum_Flore WT_Day12':1,
                    'Cecum_Flore WT_Day22':1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        select = self.h.aRange()
        aval = [ahash[i] if i in ahash else None for i in atype]
        ahash = {'Cecum_Flore WT_Day0':4,
                'Cecum_Flore KO_Day0':5,
                'Cecum_Flore WT_Day4':4,
                'Cecum_Flore KO_Day4':5,
                'Cecum_Flore WT_Day12':4,
                'Cecum_Flore KO_Day12':5,
                'Cecum_Flore WT_Day22':4,
                'Cecum_Flore KO_Day22':5}
        rval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [i for i in expg if aval[i] == 0]
        self.uc = [i for i in expg if aval[i] != 0]
        self.cd = []
        self.aval = aval
        self.rval = rval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getFuso1(self):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "CRC112.2"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c treatment")
        atypes = ['C', 'FN']
        ahash = {"non-infected": 0, "F. nucleatum":1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [ i for i in expg if aval[i] == 0]
        self.uc = [ i for i in expg if aval[i] == 1]
        self.cd = []
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getCampanioni2017(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "GS8"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c tissue")
        atypes = ['HC', 'CC', 'IC', 'CG', 'IG']
        ahash = {'Gastric mucosa of Healthy Control':0,
                'Gastric mucosa of CIM Control':1,
                'Gastric mucosa of CIM-GC':3,
                'Gastric mucosa of IIM-GC':4,
                'Gastric mucosa of IIM Control':2}
        if (tn == 2):
            atypes = ['HC', 'IIM-C', 'IIM-GC']
            ahash = {'Gastric mucosa of Healthy Control':0,
                    'Gastric mucosa of IIM-GC':2,
                    'Gastric mucosa of IIM Control':1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [ i for i in expg if aval[i] == 0]
        self.uc = [ i for i in expg if aval[i] == 1]
        self.cd = [ i for i in expg if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getPG2019(self, tn=1):
        #self.db = hu.Database("/Users/sataheri/public_html/Hegemon/explore.conf")
        #self.dbid = "ibd.9"
        self.db = hu.Database("/Users/mahdi/public_html/Hegemon/explore.conf")
        self.dbid = "IBD1"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        #atype = self.h.getSurvName("c treatment groups (ch1)")
        atype = self.h.getSurvName("c type")
        atypes = ['PBS', 'DSS', 'PAR']
        #ahash = {'PBS':0, 'DSS':1, 'DSS+PAR':2}
        ahash = {'DSS':1, 'DSSPAR':2, 'PBS':0}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [ i for i in expg if aval[i] == 0]
        self.uc = [ i for i in expg if aval[i] == 1]
        self.cd = [ i for i in expg if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getChoteau2016(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "ORG28"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c treatment")
        atypes = ['PIO', 'GED-1', 'GED-5', '5-ASA']
        ahash = {'pioglitazone (1 uM)':0,
                'GED (1mM)':1,
                'GED (30mM)':2,
                '5-ASA (30mM)':3}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [ i for i in expg if aval[i] == 0]
        self.uc = [ i for i in expg if aval[i] == 1]
        self.cd = [ i for i in expg if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getBunger2007(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP109"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c desc")
        gtype = [str(k).split(" ")[2] if len(str(k).split(" ")) > 2 else None
                for k in atype]
        ttype = [str(k).split(" ")[5] if len(str(k).split(" ")) > 5 else None
                for k in atype]
        atype = [ str(gtype[i]) + " " + str(ttype[i]) for i in
                range(len(atype))]
        atypes = ['WC', 'WT', 'KC', 'KT']
        ahash = {'PPARalpha control':2,
                'wild control':0,
                'PPARalpha treated':3,
                'wild treated':1}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [ i for i in expg if aval[i] == 0]
        self.uc = [ i for i in expg if aval[i] == 1]
        self.cd = [ i for i in expg if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

    def getdeVogel2008(self, tn=1):
        self.db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
        self.dbid = "PLP110"
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        atype = self.h.getSurvName("c src1")
        gtype = [str(k).split(" ")[3] if len(str(k).split(" ")) > 3 else None
                for k in atype]
        ttype = [str(k).split(" ")[5] if len(str(k).split(" ")) > 5 else None
                for k in atype]
        atype = [ str(gtype[i]) + " " + str(ttype[i]) for i in
                range(len(atype))]
        atypes = ['WTP', 'WTH', 'WTO', 'WWY', 'KTP', 'KTH', 'KTO', 'KWY']
        ahash = {'PPARalpha trieicosapentaenoin':4,
                'PPARalpha WY14643':7,
                'PPARalpha triolein':6,
                'PPARalpha tridocosahexaenoin':5,
                'wild trieicosapentaenoin':0,
                'wild tridocosahexaenoin':1,
                'wild triolein':2,
                'wild WY14643':3}
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.h.aRange() if aval[i] is not None]
        self.normal = [ i for i in expg if aval[i] == 0]
        self.uc = [ i for i in expg if aval[i] == 1]
        self.cd = [ i for i in expg if aval[i] == 2]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.ibd = self.uc + self.cd
        self.printInfo()

