
# coding: utf-8

# In[1]:

import cv2
import re
import numpy as np
import scipy
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
import array
#get_ipython().magic(u'matplotlib inline')
import pandas as pd
import seaborn as sns
import json
from sklearn.metrics import *
from scipy.stats import fisher_exact, ttest_ind
from pprint import pprint
import os
import pickle
import sys
sys.path.append("Hegemon")
import StepMiner as smn
import HegemonUtil as hu

acolor = ["#00CC00", "#D8A03D","#EC008C",
          'cyan', "#B741DC", "#808285",
          'blue', 'black', 'green', 'red',
          'orange', 'brown', 'pink', 'purple']

try:
    reload  # Python 2.7
except NameError:
    try:
        from importlib import reload  # Python 3.4+
    except ImportError:
        from imp import reload  # Python 3.0 - 3.3

def getRealpath(cfile):
    return os.path.realpath(os.path.join(os.getcwd(), 
        os.path.dirname(__file__), cfile))

def asciiNorm(ah):
    if sys.version_info[0] >= 3:
        keys = list(ah.keys())
        for k in keys:
            ah[bytes(k, encoding='latin-1').decode('utf-8')] = ah[k]
    return ah

def reactome(idlist):
    import requests
    reactomeURI = 'http://www.reactome.org/AnalysisService/identifiers/projection?pageSize=100&page=1';
    response = requests.post(reactomeURI, data = idlist, \
                headers = { "Content-Type": "text/plain",  "dataType" : "json" })
    obj = json.loads(response.text)
    df = pd.DataFrame()
    df['name'] = [p["name"] for  p in obj["pathways"]]
    df['pValue'] = [p["entities"]["pValue"] for  p in obj["pathways"]]
    df['fdr'] = [p["entities"]["fdr"] for  p in obj["pathways"]]
    return df

def getBoolean(cfile, sthr, pthr, code):
    res = []
    with open(cfile, "r") as bFile:
        for ln in bFile:
            ll = ln.strip().split("\t")
            bs = [ [int(ll[i]) for i in range(2, 6)] ]
            bs += [ [int(ll[i]) for i in range(2, 6)] ]
            bs += [ [float(ll[i]) for i in range(6, 10)] ]
            bs += [ [float(ll[i]) for i in range(10, 14)] ]
            rel, stats = hu.getBooleanRelationType(bs, sthr, pthr)
            if rel == code:
                res.append(ll[1])
    return res

def plotViolinBar(ana, desc=None):
    fig = plt.figure(figsize=(4,4), dpi=100)
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    ax1 = plt.subplot2grid((4, 1), (0, 0))
    ax2 = plt.subplot2grid((4, 1), (1, 0), rowspan=3)
    params = {'spaceAnn': len(ana.order)/len(ana.atypes), 'tAnn': 1, 'widthAnn':1,
              'genes': [], 'ax': ax1, 'acolor': acolor}
    ax = ana.printTitleBar(params)
    res = ana.getROCAUC()
    ax.text(len(ana.cval[0]), 4, res)
    if desc is not None:
        ax.text(-1, 2, desc, horizontalalignment='right',
                    verticalalignment='center')
    params = {'spaceAnn': len(ana.order)/len(ana.atypes), 'tAnn': 1, 'widthAnn':1,
            'genes': [], 'ax': ax2, 'acolor': acolor, 'vert': 0}
    ax = ana.printViolin(None, params)
    return fig

def plotDensityBar(ana, desc=None):
    fig = plt.figure(figsize=(4,4), dpi=100)
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    ax1 = plt.subplot2grid((4, 1), (0, 0))
    ax2 = plt.subplot2grid((4, 1), (1, 0), rowspan=3)
    params = {'spaceAnn': len(ana.order)/len(ana.atypes), 'tAnn': 1, 'widthAnn':1,
              'genes': [], 'ax': ax1, 'acolor': acolor}
    ax = ana.printTitleBar(params)
    res = ana.getMetrics(ana.cval[0])
    ax.text(len(ana.cval[0]), 4, ",".join(res))
    if desc is not None:
        ax.text(-1, 2, desc, horizontalalignment='right',
                    verticalalignment='center')
    ax = ana.densityPlot(ax2, acolor)
    return fig

def processData(ana, l1, wt1, desc=None, violin=1):
    ana.orderData(l1, wt1)
    if (violin == 1):
        return plotViolinBar(ana, desc)
    return plotDensityBar(ana, desc)

def rugplot(data, pos=0, height=.1, ax=None, **kwargs):
    from matplotlib.collections import LineCollection
    ax = ax or plt.gca()
    zero = np.zeros_like(data)
    kwargs.setdefault("linewidth", 1)
    segs = np.stack((np.c_[data, data],
                     np.c_[zero+pos*height, zero+(pos+1)*height]),
                    axis=-1)
    lc = LineCollection(segs, transform=ax.get_xaxis_transform(), **kwargs)
    ax.add_collection(lc)
    return

def plotSingle(expr, pG):
    fig = plt.figure(figsize=(6,4), dpi=100)
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    ax1 = plt.subplot2grid((4, 1), (0, 0), rowspan=3)
    ax2 = plt.subplot2grid((4, 1), (3, 0))
    ax2.axison = False
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax2.grid(False)
    ax2.tick_params(top=False, left=False, bottom=False, right=False)

    for i in range(len(pG)):
        name, col, order = pG[i]
        if len(order) <= 0:
            continue
        vals = [expr[j] for j in order]
        ax = sns.kdeplot(vals, bw = 1, cut = 2, color=col, label=name, ax=ax1)
        ax1.axvline(x=np.mean(vals), c=col)
        rugplot(vals, pos=i, color=col, ax=ax2, height=1/len(pG))

    lims = ax2.axis(ax1.axis())
    return fig

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
        print("Can't open file {0} <br>".format(cfile))
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
        print("Can't open file {0} <br>".format(cfile))
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
        print("Can't open file {0} <br>".format(cfile))
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
    genes = ""
    if not os.path.isfile(cfile):
        print("Can't open file {0} <br>".format(cfile))
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
            nm = h.getSimpleName(id)
            row_labels.append(nm)
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

def getSName(name):
    l1 = re.split(": ", name)
    l2 = re.split(" /// ", l1[0])
    return l2[0]

def getRanksDf(df_e, df_t):
    expr = []
    row_labels = []
    row_ids = []
    row_numhi = []
    ranks = []
    g_ind = 0
    counts = []
    for k in range(len(df_e)):
        count = 0
        order = range(2, df_e[k].shape[1])
        avgrank = [0 for i in order]
        for j in range(df_e[k].shape[0]):
            e = df_e[k].iloc[j,:]
            t = df_t[k]['thr2'][j]
            if e[-1] == "":
                continue
            v = np.array([float(e[i]) if e[i] != "" else 0 for i in order])
            te = []
            sd = np.std(v)
            for i in order:
                if (e[i] != ""):
                    v1 = (float(e[i]) - t) / 3;
                    if sd > 0:
                        v1 = v1 / sd
                else:
                    v1 = -t/3/sd
                avgrank[i-2] += v1
                te.append(v1)
            expr.append(te)
            nm = getSName(e[1])
            row_labels.append(nm)
            row_ids.append(e[0])
            v1 = [g_ind, sum(v > t)]
            if g_ind > 3:
                v1 = [g_ind, sum(v <= t)]
            else:
                v1 = [g_ind, sum(v > t)]
            row_numhi.append(v1)
            count += 1
            #if count > 200:
            #    break
        ranks.append(avgrank)
        g_ind += 1
        counts += [count]
    print(counts)
    return ranks, row_labels, row_ids, row_numhi, expr

def saveList(ofile, l1):
    of = open(ofile, "w")
    for i in l1:
        of.write("\t".join([i]) + "\n")
    of.close()

def readList(cfile):
    if not os.path.isfile(cfile):
        print("Can't open file {0} <br>".format(cfile))
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
        print("Can't open file {0} <br>".format(cfile))
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
        print("Can't open file {0} <br>".format(cfile))
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
        idx = df[df == i].index
        n = len(idx)
        l = str(atypes[i]) + "(" + str(n) + ")"
        df1 = pd.DataFrame(pd.Series(idx),
                columns=[l])
        if n <= 1:
            continue
        if ax is None:
            ax = df1.plot.kde(bw_method=1.0, c=color[i], label=l)
        else:
            ax = df1.plot.kde(bw_method=1.0, ax = ax, c=color[i], label=l)
    for i in range(len(atypes)):
        idx = df[df == i].index
        n = len(idx)
        l = str(atypes[i]) + "(" + str(n) + ")"
        df1 = pd.DataFrame(pd.Series(idx),
                columns=[l])
        if n != 1:
            continue
        df1['y'] = 1
        if ax is None:
            ax = df1.plot.line(x=l, y='y', c=color[i], label=l)
            ax.axvline(x=idx[0], c=color[i])
        else:
            ax = df1.plot.line(x=l, y='y', ax = ax, c=color[i], label=l)
            ax.axvline(x=idx[0], c=color[i])

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

def setPlotStyle(params=None):
    color_sch1 = acolor
    if params is not None and 'acolor' in params:
        color_sch1 = params['acolor']
    sns.set()
    sns.set_style("white")
    sns.set_style({'text.color': '.5', 
        'xtick.color':'.5', 'ytick.color':'.5', 'axes.labelcolor': '.5'})
    sns.set_context("notebook")
    sns.set_palette([adj_light(c, 1.5, 1) for c in color_sch1])

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, h

def cAllPvals(lval, atypes):
    for i in range(len(lval)):
        for j in range(i +1, len(lval)):
            if len(lval[i]) <= 0:
                continue
            if len(lval[j]) <= 0:
                continue
            #print(lval[i])
            #print(lval[j])
            t, p = ttest_ind(lval[i],lval[j], equal_var=False)
            desc = "%s vs %s %.3g, %.3g" % (atypes[i], atypes[j], t, p)
            print(desc)

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

def plotViolin(data, atypes, params):
    df = pd.DataFrame()
    df['score'] = [k for i in range(len(data)) for k in data[i]]
    df['category'] = [atypes[i] for i in range(len(data)) for k in data[i]]
    m1 = []
    pvals = []
    for i in range(1, len(data)):
        if len(data[i]) <= 0:
            m1 += [0]
            pvals += [""]
            continue
        m1 += [max(data[i]) + (max(data[i]) - min(data[i])) * 0.1]
        t, p = ttest_ind(data[0],data[i], equal_var=False)
        if (p < 0.05):
            pvals += ["p=%.3g" % p]
        else:
            pvals += [""]
    dpi = 100
    if 'dpi' in params:
        dpi = params['dpi']
    w,h = (1.5 * len(atypes), 4)
    if 'w' in params:
        w = params['w']
    if 'h' in params:
        h = params['h']
    color_sch1 = acolor
    if 'acolor' in params:
        color_sch1 = params['acolor']
    sns.set()
    sns.set_style("white")
    sns.set_style({'text.color': '.5', 
        'xtick.color':'.5', 'ytick.color':'.5', 'axes.labelcolor': '.5'})
    sns.set_context("notebook")
    sns.set_palette([adj_light(c, 1.5, 1) for c in color_sch1])
    ax = None
    if 'ax' in params:
        ax = params['ax']
    if ax is None:
        fig,ax = plt.subplots(figsize=(w,h), dpi=dpi)
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

def getGroupsMmv1(gene_groups):
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

def getGroupsMm(gene_groups):
    cfile = getRealpath("data/ensembl-GRCh38.p13-100-hs-mm.txt")
    fp = open(cfile, "r")
    mmdict = {}
    for line in fp:
        line = line.strip();
        ll = re.split("\t", line);
        if len(ll) > 3 and ll[2] != '' and ll[3] != '':
            g = ll[3]
            if g not in mmdict:
                mmdict[g] = []
            mmdict[g] += [ll[2]]
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

def getGroupsHsv1(gene_groups):
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
            if ll[3] not in mmdict:
                mmdict[ll[3]] = []
            mmdict[ll[3]] += [g]
    fp.close();

    gene_groups_hs = []
    for s in gene_groups:
        s1 = set()
        for g in s:
            if g in mmdict:
                for k in mmdict[g]:
                    s1.add(k)
        gene_groups_hs.append(s1)
    return gene_groups_hs

def getGroupsHs(gene_groups):
    cfile = getRealpath("data/ensembl-GRCm38.p6-100-mm-hs.txt")
    fp = open(cfile, "r")
    mmdict = {}
    for line in fp:
        line = line.strip();
        ll = re.split("\t", line);
        if len(ll) > 3 and ll[1] != '' and ll[2] != '':
            g = ll[1]
            if g not in mmdict:
                mmdict[g] = []
            mmdict[g] += [ll[2]]
    fp.close();

    gene_groups_hs = []
    for s in gene_groups:
        s1 = set()
        for g in s:
            if g in mmdict:
                for k in mmdict[g]:
                    s1.add(k)
        gene_groups_hs.append(s1)
    return gene_groups_hs


def getSimpleName(gene_groups, h):
    res = []
    for s in gene_groups:
        s1 = set()
        for g in s:
            for id1 in h.getIDs(g):
                s1.add(h.getSimpleName(id1))
        res.append(s1)
    return res

def PathView(file1,  debug = 1):
    db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
    h = hu.Hegemon(db.getDataset("PLP7"))
    h.init()
    h.initPlatform()
    h.initSurv()
    data_item = []
    with open(file1) as data_file:
        data_item += json.load(data_file)
    cfile = "data/ibd-network-g-eq-cls-4.txt"
    if not os.path.isfile(cfile):
        print("Can't open file {0} <br>".format(cfile))
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

def getNBGeneGroups(order = None, weight = None, debug = 1):
    db = hu.Database("explore.conf")
    h = hu.Hegemon(db.getDataset("NB4"))
    h.init()
    h.initPlatform()
    h.initSurv()
    data_item = []
    with open('path-1.json') as data_file:
        data_item += json.load(data_file)
    with open('path-2.json') as data_file:
        data_item += json.load(data_file)
    with open('path-3.json') as data_file:
        data_item += json.load(data_file)
    cfile = "nb-net-cls.txt"
    if not os.path.isfile(cfile):
        print("Can't open file {0} <br>".format(cfile))
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
    print([len(s) for s in gene_groups])
    if order is None:
        order = [1, 3, 4, 5];
        order = [35]
        order = [43, 44, 45];
        order = [8, 9, 10]
    gene_groups = [gene_groups[i] for i in order]
    print([len(s) for s in gene_groups])
    gene_groups = getSimpleName(gene_groups, h)
    print([len(s) for s in gene_groups])
    if weight is None:
        weight = [-1, 1, 2, 3]
        weight = [-1, -2, -3]
        weight = [-1]
        weight = [-1, -2, -3]
    print(weight)
    genes = []
    return genes, weight, gene_groups

def getGeneGroups(order = None, weight = None, debug = 1):
    data_item = []
    with open('data/path-1.json') as data_file:
        data_item += json.load(data_file)
    with open('data/path-2.json') as data_file:
        data_item += json.load(data_file)
    cfile = "data/ibd-network-g-eq-cls-4.txt"
    if not os.path.isfile(cfile):
        print("Can't open file {0} <br>".format(cfile))
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
    print([len(s) for s in gene_groups])
    if order is None:
        order = [7, 6, 5, 1];
        order = [7, 6, 5];
    gene_groups = [gene_groups[i] for i in order]
    print([len(s) for s in gene_groups])
    if weight is None:
        weight = [-3, -2, -1]
    print(weight)
    genes = readGenes("data/cluster-names.txt")
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
        print("Can't open file {0} <br>".format(cfile))
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
    score = [ana.f_ranks[i - ana.start] for i in ana.i1]
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
    predicted = [1 if ana.f_ranks[i - ana.start] >= fthr else 0 for i in ana.i1]
    c_dict = {}
    for i in ana.order:
        c_dict[i] = ana.f_ranks[i - ana.start]
        c_dict[i] = 0
        if ana.f_ranks[i - ana.start] >= fthr:
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

def processGeneGroupsDf(ana, l1, wt1, debug = 0, fthr = None):
    ana.orderDataDf(l1, wt1); print("ROC-AUC", ana.getMetrics())
    actual = [1 if ana.aval[i] >= 1 else 0 for i in ana.i1]
    score = [ana.f_ranks[i - ana.start] for i in ana.i1]
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
    predicted = [1 if ana.f_ranks[i - ana.start] >= fthr else 0 for i in ana.i1]
    c_dict = {}
    for i in ana.order:
        c_dict[i] = ana.f_ranks[i - ana.start]
        c_dict[i] = 0
        if ana.f_ranks[i - ana.start] >= fthr:
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
        self.state = []
        self.params = {}
        self.otype = 0
        self.start = 2
        self.end = 2
        self.axes = []

    def addAxes(self, ax):
        self.axes += [ax]

    def aRange(self):
        return range(self.start, self.end + 1)

    def getTitle(self):
        title = self.name + " (" + self.source + "; n = " + str(self.num) + ")"
        return title

    def printInfo(self):
        print(self.name + " (n = " + str(self.num) + ")")
        url = "http://hegemon.ucsd.edu/Tools/explore.php?key=polyps&id="
        if self.dbid.startswith("NB"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=nb&id="
        if self.dbid.startswith("PLP"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=polyps&id="
        if self.dbid.startswith("CRC"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=colon&id="
        if self.dbid.startswith("MAC"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=mac&id="
        if self.dbid.startswith("MACV"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=macv&id="
        if self.dbid.startswith("LIV"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=liver&id="
        if self.dbid.startswith("G16"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=gbm&id="
        if self.dbid.startswith("GL"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=global&id="
        if self.dbid.startswith("GS"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=gastric&id="
        if self.dbid.startswith("AD"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=ad&id="
        if self.dbid.startswith("COV"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=covid&id="
        if self.dbid.startswith("LU"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=lung&id="
        if self.dbid.startswith("HRT"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=heart&id="
        print(self.source + " " + url + self.dbid)
        print(len(self.order), [len(i) for i in self.state], \
                self.source, url + self.dbid, self.dbid)

    def prepareDataDf(self, dbid):
        self.dbid = dbid
        self.dataset = hu.getHegemonDataset(self.dbid)
        self.num = self.dataset[2]
        self.name = self.dataset[1]
        self.source = self.dataset[3]
        obj = hu.getHegemonPatientData(self.dbid, 'time')
        self.headers = obj[0]
        self.hhash = {}
        self.start = 2;
        self.end = len(self.headers) - 1
        for i in range(len(self.headers)):
            self.hhash[self.headers[i]] = i

    def prepareData(self, dbid, cfile =
            "/booleanfs2/sahoo/Hegemon/explore.conf"):
        self.db = hu.Database(cfile)
        self.dbid = dbid
        self.h = hu.Hegemon(self.db.getDataset(self.dbid))
        self.h.init()
        self.h.initPlatform()
        self.h.initSurv()
        self.num = self.h.getNum()
        self.start = self.h.getStart()
        self.end = self.h.getEnd()
        self.name = self.h.rdataset.getName()
        self.source = self.h.getSource()
        self.headers = self.h.headers
        self.hhash = {}
        for i in range(len(self.headers)):
            self.hhash[self.headers[i]] = i

    def initData(self, atype, atypes, ahash):
        for i in range(len(atypes)):
            ahash[atypes[i]] = i
        aval = [ahash[i] if i in ahash else None for i in atype]
        expg = [i for i in self.aRange() if aval[i] is not None]
        self.state = [[i for i in range(len(atype)) if aval[i] == k] 
                for k in range(len(atypes))]
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = expg
        self.printInfo()

    def getSurvName(self, name):
        return hu.getHegemonPatientData(self.dbid, name)[1]

    def convertMm(self, gene_groups, genes):
        self.gene_groups = getGroupsMm(gene_groups)
        self.genes = getGroupsMm([genes])[0]

    def orderData(self, gene_groups, weight):
        self.col_labels = [self.h.headers[i] for i in self.order]
        #ranks, row_labels, expr = getRanks(gene_groups, h)
        ranks, row_labels, row_ids, row_numhi, expr = getRanks2(gene_groups,
                self.h)
        self.f_ranks = mergeRanks(self.h.aRange(), self.h.start, ranks, weight)
        #i1 = getOrder(self.order, self.h.start, ranks, weight)
        arr = [self.f_ranks[i - self.h.start] for i in self.order]
        i1 = [self.order[i] for i in np.argsort(arr)]
        index = np.array([i - self.h.start for i in i1])
        self.cval = np.array([[self.aval[i] for i in i1]])
        #self.data = np.array(expr)[:,index]
        self.ind_r = np.array(sorted(range(len(row_labels)),
            key=lambda x: (row_numhi[x][0], row_numhi[x][1])))
        row_labels = [row_labels[i] for i in self.ind_r]
        row_ids = [row_ids[i] for i in self.ind_r]
        if len(self.ind_r) > 0:
            self.data = np.array([expr[i] for i in self.ind_r])[:,index]
        else:
            self.data = np.array([expr[i] for i in self.ind_r])
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
        self.data = np.array([expr[i] for i in self.ind_r])[:,index]
        self.ranks = ranks
        self.row_labels = row_labels
        self.row_ids = row_ids
        self.row_numhi = row_numhi
        self.expr = expr
        self.i1 = i1
        self.index = index
        self.otype = 2

    def orderDataDf(self, gene_groups, weight):
        data_e = []
        data_t = []
        for k in gene_groups:
            df_e = hu.getHegemonDataFrame(self.dbid, k, None)
            df_t = hu.getHegemonThrFrame(self.dbid, k)
            rhash = {}
            for i in range(df_t.shape[0]):
                rhash[df_t.iloc[i,0]] = i
            order = [rhash[df_e.iloc[i,0]] for i in range(df_e.shape[0])]
            df_t = df_t.reindex(order)
            df_t.reset_index(inplace=True)
            data_e.append(df_e)
            data_t.append(df_t)
        self.col_labels = self.headers[self.start:]
        if len(gene_groups) > 0:
            self.col_labels = data_e[0].columns[self.start:]
        self.chash = {}
        for i in range(len(self.col_labels)):
            self.chash[self.col_labels[i]] = i
        ranks, row_labels, row_ids, row_numhi, expr = getRanksDf(data_e, data_t)
        i1 = getOrder(self.order, self.start, ranks, weight)
        index = np.array([i - self.start for i in i1])
        self.cval = np.array([[self.aval[i] for i in i1]])
        #self.data = np.array(expr)[:,index]
        self.ind_r = np.array(sorted(range(len(row_labels)),
            key=lambda x: (row_numhi[x][0], row_numhi[x][1])))
        row_labels = [row_labels[i] for i in self.ind_r]
        row_ids = [row_ids[i] for i in self.ind_r]
        self.data = np.array([expr[i] for i in self.ind_r])[:,index]
        self.f_ranks = mergeRanks(range(self.start, len(self.headers)),
                self.start, ranks, weight)
        self.ranks = ranks
        self.row_labels = row_labels
        self.row_ids = row_ids
        self.row_numhi = row_numhi
        self.expr = expr
        self.i1 = i1
        self.index = index
        self.otype = 1

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
        predicted = [1 if f_ranks[i - self.start] >= thr[0] - nm else 0 for i in i1]
        if "thr" in self.params:
            if self.params["thr"] == 1:
                predicted = [1 if f_ranks[i - self.start] >= thr[0] else 0 for i in i1]  
            if self.params["thr"] == 2:
                predicted = [1 if f_ranks[i - self.start] >= thr[0] + nm else 0 for i in i1]  
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
        thr = hu.getThrData([f_ranks[i - self.start] for i in i1])
        nm = (np.max(f_ranks) - np.min(f_ranks))/15
        predicted = [1 if f_ranks[i - self.start] >= thr[0] - nm else 0 for i in i1]
        if "thr" in self.params:
            if self.params["thr"] == 1:
                predicted = [1 if f_ranks[i - self.start] >= thr[0] else 0 for i in i1]  
            if self.params["thr"] == 2:
                predicted = [1 if f_ranks[i - self.start] >= thr[0] + nm else 0 for i in i1]  
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
                score = [ana.f_ranks[i - ana.start] for i in ana.i1]
            else:
                score = [ana.f_ranks[ana.i1[i] - ana.start] \
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
                score = [ana.f_ranks[i - ana.start] for i in ana.i1]
            else:
                score = [ana.f_ranks[ana.i1[i] - ana.start] \
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
        #print(thr, nm)
        if fthr is None or fthr == "thr1":
            fthr = thr[0]
        if fthr == "thr0":
            fthr = thr[0] - nm
        if fthr == "thr2":
            fthr = thr[0] + nm
        if fthr == "thr3":
            fthr = thr[0] + 3 * nm
        predicted = [1 if ana.f_ranks[i - ana.start] >= fthr else 0 for i in ana.i1]
        if "predicted" in ana.params:
            predicted = ana.params["predicted"]
        print(list(actual))
        print(list(predicted))
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
                score = [ana.f_ranks[i - ana.start] for i in ana.i1]
            else:
                score = [ana.f_ranks[ana.i1[i] - ana.start] \
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
        atypes = self.params['atypes']
        atypes = [str(atypes[i]) + "("+str(len(lval[i]))+")"
                for i in range(len(atypes))]
        ax,bp = plotScores(lval, atypes, self.params)
        ax.text(ax.get_xlim()[1], ax.get_ylim()[1], self.source,
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
                pvals += [""]
                m1 += [0.1]
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
            ax.text(ax.get_xlim()[1], 1, self.source,
                    horizontalalignment='left', verticalalignment='center')
        else:
            title = self.getTitle()
            ax.set_title(title)
            ax.set_ylabel(self.h.getSimpleName(id1))
        self.addAxes(ax)
        self.addAxes(bp)
        cAllPvals(lval, self.params['atypes'])
        return ax,bp

    def getROCAUCspecific(ana, m=0, n=1):
        actual = [ana.aval[i] for i in ana.i1
                if ana.aval[i] == m or ana.aval[i] == n]
        score = [ana.f_ranks[i - ana.start] for i in ana.i1
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

    def getTPvals(ana):
        res = []
        lval, score = ana.getScores()
        for k in range(1, len(ana.atypes)):
            t, p = ttest_ind(lval[0],lval[k], equal_var=False)
            res += ["%.3g" % p]
        return ",".join(res)

    def getStats(self, l1, wt1, annotation=[]):
        src = re.sub(" .*", "", self.source)
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
        self.prepareData(dbid)
        atype = self.h.getSurvName("status")
        atypes = ['Censor', 'Relapse']
        ahash = {"0": 0, "1":1}
        self.initData(atype, atypes, ahash)

    def getSurvivalDf(self, dbid = "CRC35.3"):
        self.prepareDataDf(dbid)
        atype = self.getSurvName("status")
        atypes = ['Censor', 'Relapse']
        ahash = {"0": 0, "1":1}
        self.initData(atype, atypes, ahash)

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
        g1 = [i for i in order if f_ranks[i - self.start] < fthr]
        g2 = [i for i in order if f_ranks[i - self.start] >= fthr]
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
        score = [ana.f_ranks[i - ana.start] for i in ana.h.aRange()]
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
        c1 = ["", ""] + [1 if ana.f_ranks[i - ana.start] >= fthr else 0 for i in
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

    def printMeanAbsoluteDeviation(ana, ofile):
        from scipy.stats import median_absolute_deviation
        fp = ana.h.fp;
        fp.seek(0, 0);
        head = fp.readline();
        of = open(ofile, "w")
        of.write("\t".join(["ArrayID", "MAD"]) + "\n")
        index = 0
        for line in fp:
          line = re.sub("[\r\n]", "", line)
          ll = line.split("\t")
          if len([i for i in ana.order if ll[i] == '']) > 0:
              continue
          v1 = [float(ll[i]) for i in ana.order]
          mad = median_absolute_deviation(v1)
          of.write("\t".join([ll[0], str(mad)]) +"\n")
          index += 1
        of.close()

    def getJSTOM(self):
        self.getSurvival("CRC35.3")

    def getBos(self):
        self.getSurvival("BC20")

    def getPeters(self, tn=1):
        self.prepareData("PLP7")
        atype = self.h.getSurvName("c clinical condition")
        ahash = {'Ulcerative Colitis':1, 'control':0, "Crohn's disease":2}
        atypes = ['Normal', 'UC', 'CD']
        if (tn == 2):
            atypes = ['Normal', 'IBD']
            ahash = {"control": 0, "Ulcerative Colitis":1, "Crohn's disease":1}
        self.initData(atype, atypes, ahash)

    def getPetersDf(self, tn=1):
        self.prepareDataDf("PLP7")
        atype = hu.getHegemonPatientData(self.dbid, 'c clinical condition')[1]
        ahash = {'Ulcerative Colitis':1, 'control':0, "Crohn's disease":2}
        atypes = ['Normal', 'UC', 'CD']
        if (tn == 2):
            atypes = ['Normal', 'IBD']
            ahash = {"control": 0, "Ulcerative Colitis":1, "Crohn's disease":1}
        self.initData(atype, atypes, ahash)

    def getWu2007(self, tn=1):
        self.prepareData("PLP12")
        atype = self.h.getSurvName("c type")
        atypes = ['N', 'UC-un', 'CD-un', 'IC-un', 'INF', 
                'UC-Aff', 'CD-Aff', 'IC-Aff']
        atypes = ['N', 'UC', 'CD']
        ahash = {'N': 0, 'UC-Aff': 1, 'CD-Aff': 2}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'N': 0, 'UC-Aff': 1, 'CD-Aff': 1}
        self.initData(atype, atypes, ahash)

    def getWu2007Df(self, tn=1):
        self.prepareDataDf("PLP12")
        atype = hu.getHegemonPatientData(self.dbid, 'c type')[1]
        atypes = ['N', 'UC-un', 'CD-un', 'IC-un', 'INF', 
                'UC-Aff', 'CD-Aff', 'IC-Aff']
        atypes = ['N', 'UC', 'CD']
        ahash = {'N': 0, 'UC-Aff': 1, 'CD-Aff': 2}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'N': 0, 'UC-Aff': 1, 'CD-Aff': 1}
        self.initData(atype, atypes, ahash)

    def getVerstockt2019(self, tn=1):
        self.prepareData("PLP73")
        treatment = self.h.getSurvName("c treatment")
        status = self.h.getSurvName("c organism part")
        res = hu.uniq(status)
        res1 = hu.uniq(treatment)
        atype = self.h.getSurvName("c clinical history")
        atypes = ['R', 'NR']
        ahash = {'responder':0, 'non-responder':1}
        if (tn == 2):
            atype = [atype[i] if status[i] == res[3] and treatment[i] == res1[3]
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getVerstockt2019Df(self, tn=1):
        self.prepareDataDf("PLP73")
        treatment = self.getSurvName("c treatment")
        status = self.getSurvName("c organism part")
        res = hu.uniq(status)
        res1 = hu.uniq(treatment)
        atype = self.getSurvName("c clinical history")
        atypes = ['R', 'NR']
        ahash = {'responder':0, 'non-responder':1}
        if (tn == 2):
            atype = [atype[i] if status[i] == res[3] and treatment[i] == res1[3]
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getArijs2018(self, tn=1):
        self.prepareData("PLP10")
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
        if (tn == 4):
            atype = [atype[i] if rval[i] == 3
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getArijs2018Df(self, tn=1):
        self.prepareDataDf("PLP10")
        rtype = self.getSurvName("c induction therapy_maintenance therapy")
        rtypes = ['CO', 'PLAC', 'VDZ', 'IFX']
        ahash = {'plac_plac':1, 'vdz_vdz4w':2, 'vdz_plac':2, 'vdz_vdz8w':2,
                'vdz4w':2}
        for i in range(len(rtypes)):
            ahash[rtypes[i]] = i
        rval = [ahash[i] if i in ahash else None for i in rtype]
        atype = self.getSurvName("c Response")
        atypes = ['Control', 'UC R ', 'UC NR ', 'Active UC ', 'UC other ']
        ahash = {}
        if (tn == 2):
            atypes = ['Normal', 'IBD']
            ahash = {'Control':0, 'UC R ':1, 'UC NR ':1, 'Active UC ':1,
                    'UC other ':1}
        if (tn == 3 or tn == 4):
            atypes = ['R', 'NR']
            ahash = {'UC R ':0, 'UC NR ':1}
        if (tn == 4):
            atype = [atype[i] if rval[i] == 3
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getArijs2009(self, tn=1):
        self.prepareData("PLP27")
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
        atype = [atype[i] if tissue[i] == "Colon" and \
                treatment[i] != "After first infliximab treatment"
                else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getArijs2009Df(self, tn=1):
        self.prepareDataDf("PLP27")
        tissue = self.getSurvName("c tissue")
        treatment = self.getSurvName("c before or after first infliximab treatment")
        response = self.getSurvName("c response to infliximab")
        atype = self.getSurvName("c disease")
        atypes = ['Control', 'UC', 'CD', "R", "NR", "N/A"]
        ahash = {"Yes": 3, "No": 4, "Not applicable": 5}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {"Control": 0, "UC": 1, "CD": 1}
        if (tn == 3):
            atypes = ['R', 'NR']
            ahash = {"Yes": 0, "No": 1}
            atype = response
        atype = [atype[i] if tissue[i] == "Colon" and \
                treatment[i] != "After first infliximab treatment"
                else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getVanhove(self, tn=1):
        self.prepareData("PLP23")
        activity = self.h.getSurvName("c disease activity")
        atype = self.h.getSurvName("c disease")
        atypes = ['Control', 'UC', 'CD', 'I', 'A', 'N']
        ahash = {'ulcerative colitis':1, "Crohn's disease":2, 'control':0,
                'active':4, 'inactive':3, 'normal': 5}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'ulcerative colitis':1, "Crohn's disease":1, 'control':0}
        if (tn == 3):
            atype = activity
            atypes = ['I', 'A']
            ahash = {'active':1, 'inactive':0}
        self.initData(atype, atypes, ahash)

    def getVanhoveDf(self, tn=1):
        self.prepareDataDf("PLP23")
        activity = self.getSurvName("c disease activity")
        atype = self.getSurvName("c disease")
        atypes = ['Control', 'UC', 'CD', 'I', 'A', 'N']
        ahash = {'ulcerative colitis':1, "Crohn's disease":2, 'control':0,
                'active':4, 'inactive':3, 'normal': 5}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'ulcerative colitis':1, "Crohn's disease":1, 'control':0}
        if (tn == 3):
            atype = activity
            atypes = ['I', 'A']
            ahash = {'active':1, 'inactive':0}
        self.initData(atype, atypes, ahash)

    def getVanderGoten(self, tn=1):
        self.prepareData("PLP25")
        activity = self.h.getSurvName("c disease activity")
        atype = self.h.getSurvName("c disease")
        atypes = ['Control', 'UC', 'CD', 'I', 'A', 'NA']
        ahash = {'control':0,
                'active':4, 'inactive':3, 'not applicable': 5}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'control':0, 'UC':1, 'CD':1}
        if (tn == 3):
            atype = activity
            atypes = ['I', 'A']
            ahash = {'active':1, 'inactive':0}
        self.initData(atype, atypes, ahash)

    def getVanderGotenDf(self, tn=1):
        self.prepareDataDf("PLP25")
        activity = self.getSurvName("c disease activity")
        atype = self.getSurvName("c disease")
        atypes = ['Control', 'UC', 'CD', 'I', 'A', 'NA']
        ahash = {'control':0,
                'active':4, 'inactive':3, 'not applicable': 5}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'control':0, 'UC':1, 'CD':1}
        if (tn == 3):
            atype = activity
            atypes = ['I', 'A']
            ahash = {'active':1, 'inactive':0}
        self.initData(atype, atypes, ahash)

    def getDePreter(self, t1=1):
        self.prepareData("PLP24")
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
        self.initData(atype, atypes, ahash)

    def getDePreterDf(self, t1=1):
        self.prepareDataDf("PLP24")
        atype = self.getSurvName("c src1")
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
        self.initData(atype, atypes, ahash)

    def getPekow(self, tn=1):
        self.prepareData("PLP59")
        atype = self.h.getSurvName("c src1")
        atypes = ['C', 'qUC', 'nUC']
        ahash = {'normal control': 0, 
                'quiescent ulcerative colitis': 1,
                'ulcerative colitis with neoplasia': 2}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'normal control': 0, 
                    'quiescent ulcerative colitis': 1,
                    'ulcerative colitis with neoplasia': 1}
        if (tn == 3):
            atypes = ['qUC', 'nUC']
            ahash = {'quiescent ulcerative colitis': 0,
                    'ulcerative colitis with neoplasia': 1}
        self.initData(atype, atypes, ahash)

    def getPekowDf(self, tn=1):
        self.prepareDataDf("PLP59")
        atype = self.getSurvName("c src1")
        atypes = ['C', 'qUC', 'nUC']
        ahash = {'normal control': 0, 
                'quiescent ulcerative colitis': 1,
                'ulcerative colitis with neoplasia': 2}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'normal control': 0, 
                    'quiescent ulcerative colitis': 1,
                    'ulcerative colitis with neoplasia': 1}
        if (tn == 3):
            atypes = ['qUC', 'nUC']
            ahash = {'quiescent ulcerative colitis': 0,
                    'ulcerative colitis with neoplasia': 1}
        self.initData(atype, atypes, ahash)

    def getBreynaert2013(self, tn=1):
        self.prepareData("PLP38")
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
        self.initData(atype, atypes, ahash)

    def getBreynaert2013Df(self, tn=1):
        self.prepareDataDf("PLP38")
        atype = self.getSurvName("c colitis group")
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
        self.initData(atype, atypes, ahash)

    def getJensen2017(self):
        self.prepareData("PLP66")
        atype = self.h.getSurvName("c disease")
        atypes = ['normal', 'colitis']
        ahash = {}
        self.initData(atype, atypes, ahash)

    def getJensen2017Df(self):
        self.prepareDataDf("PLP66")
        atype = self.getSurvName("c disease")
        atypes = ['normal', 'colitis']
        ahash = {}
        self.initData(atype, atypes, ahash)

    def getDohi2014(self, tn=1):
        self.prepareData("PLP94")
        gtype = self.h.getSurvName("c treated with")
        ahash = {
                'none (untreated control)':2,
                '10 mg/kg control IgG2a mAb (anti-human CD20)':3,
                '0.3 mg/kg TNFR-Ig':4,
                '10 mg/kg anti-TWEAK mP2D10':5,
                'combination of TNFR-Fc (0.3 mg/kg) and anti-TWEAK mP2D10 (10 mg/kg)':6}
        rval = [ahash[i] if i in ahash else None for i in gtype]
        atype = self.h.getSurvName("c injected with")
        atypes = ['N', 'UC']
        ahash = {'trinitrobenzene sulfonic acid (TNBS)':1,
                'none (na\xc3\xafve control)':0}
        ahash = asciiNorm(ahash)
        if (tn == 2):
            atype = [atype[i] if rval[i] == 2
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getDohi2014Df(self, tn=1):
        self.prepareDataDf("PLP94")
        gtype = self.getSurvName("c treated with")
        ahash = {
                'none (untreated control)':2,
                '10 mg/kg control IgG2a mAb (anti-human CD20)':3,
                '0.3 mg/kg TNFR-Ig':4,
                '10 mg/kg anti-TWEAK mP2D10':5,
                'combination of TNFR-Fc (0.3 mg/kg) and anti-TWEAK mP2D10 (10 mg/kg)':6}
        rval = [ahash[i] if i in ahash else None for i in gtype]
        atype = self.getSurvName("c injected with")
        atypes = ['N', 'UC']
        ahash = {'trinitrobenzene sulfonic acid (TNBS)':1,
                'none (na\xc3\xafve control)':0}
        ahash = asciiNorm(ahash)
        if (tn == 2):
            atype = [atype[i] if rval[i] == 2
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getLamas2018(self, tn=1):
        self.prepareData("PLP99")
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
        self.initData(atype, atypes, ahash)

    def getLamas2018Df(self, tn=1):
        self.prepareDataDf("PLP99")
        atype = self.getSurvName("c src1")
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
        self.initData(atype, atypes, ahash)

    def getLyons2018(self, tn=1):
        self.prepareData("PLP89")
        atype = self.h.getSurvName("c inflammation level")
        atypes = ['NI', 'M', 'S']
        ahash = {'severe':2, 'moderate':1, 'non-inflamed':0}
        if (tn ==2):
            atypes = ['NI', 'I']
            ahash = {'severe':1, 'moderate':1, 'non-inflamed':0}
        self.initData(atype, atypes, ahash)

    def getLyons2018Df(self, tn=1):
        self.prepareDataDf("PLP89")
        atype = self.getSurvName("c inflammation level")
        atypes = ['NI', 'M', 'S']
        ahash = {'severe':2, 'moderate':1, 'non-inflamed':0}
        if (tn ==2):
            atypes = ['NI', 'I']
            ahash = {'severe':1, 'moderate':1, 'non-inflamed':0}
        self.initData(atype, atypes, ahash)

    def getFang2012(self, tn=1):
        self.prepareData("PLP90")
        atype = self.h.getSurvName("c time")
        atypes = ['W0', 'W2', 'W4', 'W6']
        ahash = {'0 week':0, '4 weeks':2, '6 weeks':3, '2 weeks':1}
        if (tn ==2):
            atypes = ['W0', 'W2-6']
            ahash = {'0 week':0, '4 weeks':1, '6 weeks':1, '2 weeks':1}
        self.initData(atype, atypes, ahash)

    def getFang2012Df(self, tn=1):
        self.prepareDataDf("PLP90")
        atype = self.getSurvName("c time")
        atypes = ['W0', 'W2', 'W4', 'W6']
        ahash = {'0 week':0, '4 weeks':2, '6 weeks':3, '2 weeks':1}
        if (tn ==2):
            atypes = ['W0', 'W2-6']
            ahash = {'0 week':0, '4 weeks':1, '6 weeks':1, '2 weeks':1}
        self.initData(atype, atypes, ahash)

    def getRuss2013(self, tn=1):
        self.prepareData("PLP96")
        tissue = self.h.getSurvName("c tissue")
        ahash = {'colon epithelium':2, 'colon':3}
        rval = [ahash[i] if i in ahash else None for i in tissue]
        atype = self.h.getSurvName("c genotype/variation")
        atypes = ['WT', 'IL10']
        ahash = {'IL10-/-':1, 'wildtype':0}
        if (tn == 2):
            atype = [atype[i] if rval[i] == 2
                    else None for i in range(len(atype))]
        if (tn == 3):
            atype = [atype[i] if rval[i] == 3
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getRuss2013Df(self, tn=1):
        self.prepareDataDf("PLP96")
        tissue = self.getSurvName("c tissue")
        ahash = {'colon epithelium':2, 'colon':3}
        rval = [ahash[i] if i in ahash else None for i in tissue]
        atype = self.getSurvName("c genotype/variation")
        atypes = ['WT', 'IL10']
        ahash = {'IL10-/-':1, 'wildtype':0}
        if (tn == 2):
            atype = [atype[i] if rval[i] == 2
                    else None for i in range(len(atype))]
        if (tn == 3):
            atype = [atype[i] if rval[i] == 3
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getTam2019(self):
        self.prepareData("PLP98")
        atype = self.h.getSurvName("c genotype")
        atypes = ['WT', 'Tnfr1']
        ahash = {'Tnfrsf1a-/-':1, 'WT':0}
        self.initData(atype, atypes, ahash)

    def getTam2019Df(self):
        self.prepareDataDf("PLP98")
        atype = self.getSurvName("c genotype")
        atypes = ['WT', 'Tnfr1']
        ahash = {'Tnfrsf1a-/-':1, 'WT':0}
        self.initData(atype, atypes, ahash)

    def getPunit2015(self):
        self.prepareData("PLP97")
        atype = self.h.getSurvName("c genotype")
        atypes = ['WT', 'TNFR2']
        ahash = {'Wildtype':0, 'TNFR2-knockout':1}
        self.initData(atype, atypes, ahash)

    def getPunit2015Df(self):
        self.prepareDataDf("PLP97")
        atype = self.getSurvName("c genotype")
        atypes = ['WT', 'TNFR2']
        ahash = {'Wildtype':0, 'TNFR2-knockout':1}
        self.initData(atype, atypes, ahash)

    def getVancamelbeke(self, tn=1):
        self.prepareData("PLP16")
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
        self.initData(atype, atypes, ahash)

    def getVancamelbekeDf(self, tn=1):
        self.prepareDataDf("PLP16")
        atype = self.getSurvName("c src1")
        atypes = ['N', 'UC', 'CD']
        ahash = {'Biopsy from inflamed colonic mucosa of active UC patient':1,
                'Biopsy from inflamed colonic mucosa of active CD patient':2,
                'Biopsy from normal colonic mucosa of control individual':0}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'Biopsy from inflamed colonic mucosa of active UC patient':1,
                    'Biopsy from inflamed colonic mucosa of active CD patient':1,
                    'Biopsy from normal colonic mucosa of control individual':0}
        self.initData(atype, atypes, ahash)

    def getDenson2018(self, tn=1):
        self.prepareData("PLP86")
        response = self.h.getSurvName("c week 4 remission")
        atype = self.h.getSurvName("c diagnosis")
        atypes = ['Control', 'UC', 'Yes', 'No', 'NA'];
        ahash = {'Ulcerative Colitis':1}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'Control':0, 'Ulcerative Colitis':1}
        self.initData(atype, atypes, ahash)

    def getDenson2018Df(self, tn=1):
        self.prepareDataDf("PLP86")
        response = self.getSurvName("c week 4 remission")
        atype = self.getSurvName("c diagnosis")
        atypes = ['Control', 'UC', 'Yes', 'No', 'NA'];
        ahash = {'Ulcerative Colitis':1}
        if (tn == 2):
            atypes = ['N', 'IBD']
            ahash = {'Control':0, 'Ulcerative Colitis':1}
        self.initData(atype, atypes, ahash)

class NBAnalysis(IBDAnalysis):

    def __init__(self):
        IBDAnalysis.__init__(self)

    def getZage2020(self, tn=1):
        self.prepareDataDf("NB14")
        atype = self.getSurvName('c Type')
        ahash = {'Kely':0, 'SKNSH':1, 'SKNBE2':2, 'P134':3, 'SKNAS':4,
                'IMR32':5, 'NGP':6, 'LAN1':7}
        rval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName('c Group')
        atypes = ['C', 'RA']
        ahash = {}
        if (tn >= 3):
            atype = [atype[i] if rval[i] == (tn - 3) else None 
                    for i in range(len(atype))]
        if (tn == 2):
            bhash = {1:1, 2:1, 4:1, 5:1}
            atype = [atype[i] if rval[i] in bhash else None 
                    for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getTARGET2013(self, tn=1):
        self.prepareDataDf("NB12")
        atype = self.getSurvName('c Grade')
        atypes = ['U', 'D']
        ahash = {'Undifferentiated or Poorly Differentiated':0,
                'Differentiating':1}
        self.initData(atype, atypes, ahash)

    def getJanoueixLerosey2008(self, tn=1):
        self.prepareDataDf("NB15")
        atype = self.getSurvName("c Title")
        atype = [re.sub(" .*", "", str(k)) for k in atype]
        atypes = ['NB', 'GGNB', 'GN']
        ahash = {}
        self.initData(atype, atypes, ahash)

    def getAlbino2008(self, tn=1):
        self.prepareDataDf("NB16")
        atype = self.getSurvName("c Phenotype")
        atype = [re.sub("Borderline g", "G", str(k)) for k in atype]
        atype = [re.sub("[, ].*", "", str(k)) for k in atype]
        atypes = ['NB', 'GGNB', 'GN']
        ahash = {'Neuroblastoma':0, 'Ganglioneuroma':2,
                'Ganglioneuroblastoma':1}
        self.initData(atype, atypes, ahash)

    def getNishida2008(self, tn=1, tb=1):
        self.prepareData("NB17")
        atype = self.h.getSurvName("c src1")
        atype = [re.sub(".*RA\), ", "RA ", str(k)) for k in atype]
        atype = [re.sub(".*NF\), ", "BDNF ", str(k)) for k in atype]
        ahash = {'RA LY294002, 1day':1, 'RA LY294002, 3day':3,
                'RA LY294002, 2day':2, 'RA LY294002, 0hour':0,
                'RA LY294002, 6hour':0.25, 'RA LY294002, 5day':5,
                'BDNF 6hour':0.25, 'BDNF 1day':1, 'BDNF 2day':2, 'BDNF 3day':3,
                'RA 1day':1, 'RA 0hour':0, 'RA 5day':5,
                'RA 3day':3, 'RA 6hour':0.25, 'RA 2day':2}
        tval = [ahash[i] if i in ahash else None for i in atype]
        atypes = ['C', 'RA', 'RAi', 'BDNF']
        ahash = {'RA LY294002, 1day':0, 'RA LY294002, 3day':2,
                'RA LY294002, 2day':0, 'RA LY294002, 0hour':0,
                'RA LY294002, 6hour':0, 'RA LY294002, 5day':2,
                'BDNF 6hour':0, 'BDNF 1day':0, 'BDNF 2day':0, 'BDNF 3day':3,
                'RA 1day':0, 'RA 0hour':0, 'RA 5day':1,
                'RA 3day':1, 'RA 6hour':0, 'RA 2day':0}
        if (tn == 2):
            atypes = ['C', 'RA']
            ahash = {'RA 1day':1, 'RA 0hour':0, 'RA 5day':1,
                    'RA 3day':1, 'RA 6hour':1, 'RA 2day':1}
            atype = [atype[i] if tval[i] == tb or tval[i] == 0
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getNishida2008Df(self, tn=1, tb=1):
        self.prepareDataDf("NB17")
        atype = self.getSurvName("c src1")
        atype = [re.sub(".*RA\), ", "RA ", str(k)) for k in atype]
        atype = [re.sub(".*NF\), ", "BDNF ", str(k)) for k in atype]
        ahash = {'RA LY294002, 1day':1, 'RA LY294002, 3day':3,
                'RA LY294002, 2day':2, 'RA LY294002, 0hour':0,
                'RA LY294002, 6hour':0.25, 'RA LY294002, 5day':5,
                'BDNF 6hour':0.25, 'BDNF 1day':1, 'BDNF 2day':2, 'BDNF 3day':3,
                'RA 1day':1, 'RA 0hour':0, 'RA 5day':5,
                'RA 3day':3, 'RA 6hour':0.25, 'RA 2day':2}
        tval = [ahash[i] if i in ahash else None for i in atype]
        atypes = ['C', 'RA', 'RAi', 'BDNF']
        ahash = {'RA LY294002, 1day':0, 'RA LY294002, 3day':2,
                'RA LY294002, 2day':0, 'RA LY294002, 0hour':0,
                'RA LY294002, 6hour':0, 'RA LY294002, 5day':2,
                'BDNF 6hour':0, 'BDNF 1day':0, 'BDNF 2day':0, 'BDNF 3day':3,
                'RA 1day':0, 'RA 0hour':0, 'RA 5day':1,
                'RA 3day':1, 'RA 6hour':0, 'RA 2day':0}
        if (tn == 2):
            atypes = ['C', 'RA']
            ahash = {'RA 1day':1, 'RA 0hour':0, 'RA 5day':1,
                    'RA 3day':1, 'RA 6hour':1, 'RA 2day':1}
            atype = [atype[i] if tval[i] == tb or tval[i] == 0
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getOhtaki2010(self, tn=1):
        self.prepareData("NB19")
        atype = self.h.getSurvName("c outcome of the patient")
        atypes = ['D', 'A']
        ahash = {'Alive':1, 'Died of disease':0}
        if (tn == 2):
            atype = [None if self.headers[i] == 'GSM408940'
                    else atype[i] for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getOhtaki2010Df(self, tn=1):
        self.prepareDataDf("NB19")
        atype = self.getSurvName("c outcome of the patient")
        atypes = ['D', 'A']
        ahash = {'Alive':1, 'Died of disease':0}
        if (tn == 2):
            atype = [None if self.headers[i] == 'GSM408940'
                    else atype[i] for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getAsgharzadeh2006(self, tn=1):
        self.prepareDataDf("NB11")
        atype = self.getSurvName("c src1")
        atype = [re.sub(".* ", "", str(k)) for k in atype]
        atypes = ['D', 'R']
        ahash = {'diagnosis':0, 'relapse':1}
        self.initData(atype, atypes, ahash)

    def getAckermann2018(self, tn=1):
        self.prepareDataDf("NB9")
        atype = self.getSurvName("c mycn status")
        ahash = {'not amplified':0, 'amplified':1, 'N/A':2}
        mval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c tert rearrangement")
        atypes = ['+', '-']
        ahash = {}
        if (tn == 2):
            atype = [atype[i] if (mval[i] == 0)
                    else None for i in range(len(atype))]
        if (tn == 3):
            atype = [atype[i] if (mval[i] == 1)
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getRifatbegovic2018(self, tn=1):
        self.prepareDataDf("NB21")
        atype = self.getSurvName("c Title")
        ctype = [re.sub(" .*", "", str(k)) for k in atype]
        ahash = {'DTCs':0, 'MNCs':1, 'Tumor':2}
        tval = [ahash[i] if i in ahash else None for i in ctype]
        atype = [re.sub(".* at ", "", str(k)) for k in atype]
        atype = [re.sub(". .*", "", str(k)) for k in atype]
        atypes = ['D', 'R']
        ahash = {'diagnosis':0, 'relapse':1}
        if (tn == 2):
            atype = [atype[i] if (tval[i] == 0)
                    else None for i in range(len(atype))]
        if (tn == 3):
            atype = [atype[i] if (tval[i] == 1)
                    else None for i in range(len(atype))]
        if (tn == 4):
            atype = tval
            atypes = ['DTCs', 'MNCs', 'Tumor']
            ahash = {0:0, 1:1, 2:2}
        self.initData(atype, atypes, ahash)

    def getClaeys2019I(self, tn=1):
        self.prepareDataDf("NB22")
        atype = self.getSurvName("c transduced with")
        atype = [re.sub(" .*", "", str(k)) for k in atype]
        atypes = ['none', 'HBP1']
        ahash = {}
        self.initData(atype, atypes, ahash)

    def getClaeys2019II(self, tn=1):
        self.prepareDataDf("NB22.2")
        atype = self.getSurvName("c Array ID")
        atype = [re.sub("_.$", "", str(k)) for k in atype]
        atypes = ['DMSO', 'BEZ', 'SAHA', 'BEZ_SAHA']
        ahash = {}
        if (tn == 2):
            atypes = ['DMSO', 'BEZ', 'SAHA']
            ahash = {}
        if (tn == 3):
            atypes = ['DMSO', 'BEZ']
            ahash = {}
        self.initData(atype, atypes, ahash)

    def getWesterlund2017II(self, tn=1):
        self.prepareDataDf("NB23.2")
        atype = self.getSurvName("c time point")
        ahash = {'EP':0, 'day10':1, 'day14':2}
        tval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c injected cells")
        ahash = {'CHP212':0, 'LAN-1':1}
        gval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c treatment")
        atypes = ['C', 'RA', 'AZA', 'RA+AZA']
        ahash = {'AZA':2, 'AZA + RA':3, 'DMSO (CTRL)':0, 'RA':1}
        if (tn == 2):
            atype = [atype[i] if (gval[i] == 0)
                    else None for i in range(len(atype))]
        if (tn == 3):
            atype = [atype[i] if gval[i] == 1
                    else None for i in range(len(atype))]
        if (tn == 4):
            atypes = ['C', 'RA+AZA']
            ahash = {'AZA + RA':1, 'DMSO (CTRL)':0}
            atype = [atype[i] if gval[i] == 1 and tval[i] == 2
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getWesterlund2017III(self, tn=1):
        self.prepareDataDf("NB23.3")
        atype = self.getSurvName("c time point")
        ahash = {'EP':0}
        tval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c injected cells")
        ahash = {'SK-N-AS':0}
        gval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c treatment")
        atypes = ['C', 'RA', 'AZA', 'RA+AZA']
        ahash = {'AZA':2, 'AZA + RA':3, 'DMSO (CTRL)':0, 'RA':1}
        if (tn == 2):
            atypes = ['C', 'RA+AZA']
            ahash = {'AZA + RA':1, 'DMSO (CTRL)':0}
        self.initData(atype, atypes, ahash)

    def getFrumm2013(self, tn=1):
        self.prepareDataDf("NB24")
        atype = self.getSurvName("c treatment duration")
        ahash = {'6':6, '72':72, '24':24}
        tval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c cell line")
        ahash = {'BE(2)-C':0}
        gval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c treatment")
        atypes = ['C', 'RA', 'VPA', 'RA+VPA']
        ahash = {'5uM ATRA':1, '1mM valproic acid (VPA)':2, 'DMSO':0,
                '5uM ATRA + 1mM valproic acid (VPA)':3}
        if (tn == 2):
            atypes = ['C', 'RA']
            ahash = {'5uM ATRA':1, 'DMSO':0}
            atype = [atype[i] if tval[i] == 24 or tval[i] == 72
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getPezzini2017(self, tn=1):
        self.prepareDataDf("NB25")
        atype = self.getSurvName("c cell line")
        ahash = {'SH-SY5Y':0}
        gval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c treatment")
        atypes = ['C', 'RA']
        ahash = {'RA-NBM':1, '5% FBS':0}
        self.initData(atype, atypes, ahash)

class BINetwork:

    def __init__(self, filename):
        if not os.path.isfile(filename):
            print("Can't open file {0} <br>".format(filename));
            exit()
        self.fp = open(filename, "rb")
        self.file = filename
        self.magic = -1;
        self.major = -1;
        self.minor = -1;
        self.num = -1;
        self.numBits = -1;
        self.numBytes = -1;
        self.startPtr = -1;
        self.name = None;
        self.low = []
        self.high = []
        self.balanced = []
        self.balancedhash = {}

    def getCounts(res):
        count = [0, 0, 0, 0, 0, 0, 0] #count each relationships
        for k in res:
            count[k] += 1
        return count

    def getCountsString(res):
        count = BINetwork.getCounts(res)
        return ", ".join([str(k) for k in count])

    def readVal(buffer, i, nBits, debug=0):
        index = int(i * nBits/8)
        v1 = 0
        v2 = 0
        if index < len(buffer): 
            v1 = buffer[index]
        if (index + 1) < len(buffer): 
            v2 = buffer[index + 1]
        shift = (i * nBits) % 8
        mask = (1 << nBits) - 1
        val = ((v1 | v2 << 8) >> shift) & mask
        if (debug == 1):
            print(v1, v2, shift, mask, val)
        return val

    #
    # Codes for Boolean Relationships:
    # 0 - No relation
    # 1 - X low -> i high
    # 2 - X low -> i low
    # 3 - X high -> i high
    # 4 - X high -> i low
    # 5 - Equivalent
    # 6 - Opposite
    #
    # X - buffer1 and buffer2 corresponds to the two lines for probe X
    # i - probe i
    def readCode(buffer1, buffer2, i, numBits, debug=0):
        index = 2 * i * numBits;
        v1 = BINetwork.readVal(buffer1, index, numBits);
        v2 = BINetwork.readVal(buffer1, index + numBits, numBits);
        v3 = BINetwork.readVal(buffer2, index, numBits);
        v4 = BINetwork.readVal(buffer2, index + numBits, numBits);
        if (debug == 1):
            print (index, index + numBits, numBits, v1, v2, v3, v4)
        total = v1 + v2 + v3 + v4;
        if (total == 1):
            if (v1 == 1):
                return 1
            if (v2 == 1):
                return 2
            if (v3 == 1):
                return 3
            if (v4 == 1):
                return 4
        if (total == 2):
            if (v2 == 1 and v3 == 1):
                return 5
            if (v1 == 1 and v4 == 1):
                return 6
        return 0;

    def getLow(self):
        return self.readList(0)
    def getHigh(self):
        return self.readList(1)
    def getBalanced(self):
        return self.readList(2)

    def readList(self, num):
        fh = self.fp
        fh.seek(3 + num * 4)
        buffer = fh.read(4)
        ptr = array.array("I", buffer)[0]
        fh.seek(ptr)
        buffer = fh.read(4);
        length = array.array("I", buffer)[0]
        buffer = fh.read(length)
        name = buffer.decode('utf-8')
        buffer = fh.read(4);
        size = array.array("I", buffer)[0]
        res = []
        for i in range(size):
            buffer = fh.read(4);
            length = array.array("I", buffer)[0]
            buffer = fh.read(length)
            res += [buffer.decode('utf-8')]
        return res
    
    def init(self):
        fh = self.fp
        fh.seek(0)
        self.magic = array.array("B", fh.read(1))[0]
        self.major = array.array("B", fh.read(1))[0]
        self.minor = array.array("B", fh.read(1))[0]
        self.low = self.getLow()
        self.high = self.getHigh();
        self.balanced = self.getBalanced();
        self.balancedhash = {}
        for i in range(len(self.balanced)):
            self.balancedhash[self.balanced[i]] = i
        fh.seek(3 + 3 * 4)
        buffer = fh.read(4)
        ptr = array.array("I", buffer)[0]
        buffer = fh.read(4)
        self.num = array.array("I", buffer)[0]
        buffer = fh.read(4)
        self.numBits = array.array("I", buffer)[0]
        self.numBytes = int(self.num * self.numBits/8) + 1
        fh.seek(ptr)
        self.startPtr = ptr
        return

    def readBlock(self):
        fh = self.fp
        buffer1 = fh.read(self.numBytes)
        buffer2 = fh.read(self.numBytes)
        res = []
        for i in range(int(self.num/2)):
            code = BINetwork.readCode(buffer1, buffer2, i, self.numBits);
            res += [code]
        return res

    def readBlockIndex(self, a):
        ptr = self.startPtr + 2 * a * self.numBytes;
        self.fp.seek(ptr)
        return self.readBlock()

    def readBlockID(self, id1):
        if id1 in self.balancedhash:
            a = self.balancedhash[id1]
            return self.readBlockIndex(a)
        return None

    def printDetails(self):
        print("Network: ", self.file)
        print(self.magic, self.major, self.minor)
        print(self.num, self.numBits, self.numBytes)
        print ("Low(", len(self.low), "):")
        #print (" ".join(self.low))
        print ("High(", len(self.high), "):")
        #print (" ".join(self.high))
        print ("Balanced(", len(self.balanced), "):")
        #print (" ".join(self.balanced))
        for i in range(10):
            count = [0, 0, 0, 0, 0, 0, 0] #count each relationships
            res = self.readBlock(); #get all relationships code for one probe
            for k in res:
                count[k] += 1
            #print the counts
            print (self.balanced[i], ", ".join([str(k) for k in count]))
    
class BIGraph:
    def readEqGraph(cfile):
        edges = {}
        nodes = {}
        count = 0
        with open(cfile, "r") as netFile:
            for ln in netFile:
                if (not ln.startswith("Found : 5")):
                    continue
                pln = ln.strip().split("\t")
                id1, id2 = pln[3], pln[4]
                nodes[id1] = 1
                nodes[id2] = 1
                if (id1 not in edges):
                    edges[id1] = {}
                count += 1
                edges[id1][id2] = 1
        print(str(count) + " edges Processed")
        return nodes, edges
    def gamma(u, edges, hash1):
        if (hash1 and u in hash1):
            return hash1[u]
        res = [u] + list(edges[u].keys())
        if hash1:
            hash1[u] = res
        return res
    def rho(u, v, edges, hash1, shash):
        if (shash and u in shash and v in shash[u]):
            return shash[u][v]
        gu = BIGraph.gamma(u, edges, hash1)
        gv = BIGraph.gamma(v, edges, hash1)
        g_union = set(gu).union(gv)
        g_int = set(gu).intersection(gv)
        res = 0
        if len(g_union) == 0:
            res = 0
        else:
            res = len(g_int) / len(g_union)
        if shash:
            shash[u][v] = res
        return res
    def pruneEqGraph(edges):
        from networkx.utils.union_find import UnionFind
        uf = UnionFind()
        hash1 = {}
        shash = {}
        num = len(edges)
        count = 0
        eqscores = []
        for u in edges:
            print(count, num, end='\r', flush=True)
            data = []
            scores = []
            for v in edges[u]:
                r = BIGraph.rho(u, v, edges, hash1, shash)
                data += [r]
                scores += [[r, v]]
                uf[u], uf[v]
            filter1 =  sorted(scores, reverse=True)
            for s in filter1:
                if (uf[u] != uf[s[1]]):
                    uf.union(u, s[1])
                    eqscores.append([u, s[0], s[1]])
                    #print(u, s[0], s[1])
                    break
            count += 1
        return pd.DataFrame(eqscores)

    def getClusters(df, thr=0.5):
        from networkx.utils.union_find import UnionFind
        uf = UnionFind()
        edges = {}
        for i in df.index:
            id1 = df[0][i]
            id2 = df[2][i]
            if id1 not in edges:
                edges[id1] = {}
            edges[id1][id1] = df[1][i] 
            if id2 not in edges:
                edges[id2] = {}
            edges[id2][id2] = df[1][i] 
            uf[id1], uf[id2]
            if (df[1][i] > thr):
                uf.union(id1, id2)
                edges[id1][id2] = df[1][i] 
                edges[id2][id1] = df[1][i]
        rank = {}
        for k in edges:
            rank[k] = len(edges[k])
        cls = {}
        for k in uf.to_sets():
            l = sorted(k, key=lambda x: rank[x], reverse=True)
            cls[l[0]] = [len(l), l]
        return cls

    def getNCodes(net, l):
        relations = {}
        for u in l:
            res = net.readBlockID(u)
            i = net.balancedhash[u]
            for j in range(len(net.balanced)):
                v = net.balanced[j]
                code = res[j]
                if code <= 0:
                    continue
                if code not in relations:
                    relations[code] = {}
                if v not in relations[code]:
                    relations[code][v] = 0
                relations[code][v] += 1
        return relations

    def getClustersGraph(net, cls):
        nodes = cls
        ids = [k for k in cls]
        eqgraph = []
        count = 0
        for u in ids:
            if (count % 100) == 0:
                print(count, end='\r', flush=True)
            nl = nodes[u][0]
            mid = int(nl / 2)
            l = [u]
            if (nl > 2):
                l = [u, nodes[u][1][1], nodes[u][1][mid]]
            if (nl > 10):
                l += [nodes[u][1][i] for i in [int(mid/4), int(mid/2), mid - 1]]
            ru = BIGraph.getNCodes(net, l)
            for c in range(1, 7):
                if c not in ru:
                    continue
                for v in ru[c]:
                    if v not in nodes:
                        continue
                    if (nl > 10 and ru[c][v] < 3):
                        continue
                    eqgraph.append([u, str(c), v, ru[c][v]])
            count += 1
        return pd.DataFrame(eqgraph)

    def readClustersGraph(cls, cg):
        edges = {}
        clusters = {}
        nodep = {}
        for u in cls:
            clusters[u] = cls[u][1]
            for k in cls[u][1]:
                nodep[k] = u
        print(str(len(clusters)) + " Processed")
        print("Processing edges...")
        count = 0
        for i in cg.index:
            if (count % 100) == 0:
                print(count, end='\r', flush=True)
            u, c, v, ru = cg.loc[i, :]
            if(u not in clusters or v not in clusters):
                continue
            if (u not in edges):
                edges[u] = {}
            if c not in edges[u]:
                edges[u][c] = {}
            count += 1
            edges[u][c][v] = 1
        print(str(count) + " edges processed")
        return edges, clusters, nodep

    def readClustersGraphFile(netprm):
        edges = {}
        clusters = {}
        nodep = {}
        with open(netprm + "-cls.txt", "r") as clusterFile:
            for ln in clusterFile:
                pln = ln.strip().split("\t")
                if(int(pln[1]) > 0): # min size threshold
                    clusters[pln[0]] = pln[2:]
                    for k in pln[2:]:
                        nodep[k] = pln[0]
        print(str(len(clusters)) + " Processed")
        print("Processing edges...")
        count = 0
        with open(netprm + "-eq-g.txt", "r") as edgeFile:
            for ln in edgeFile:
                if (count % 1000) == 0:
                    print(count, end='\r', flush=True)
                pln = ln.strip().split("\t")
                if(pln[0] not in clusters or
                   pln[2] not in clusters):
                    continue
                if (pln[0] not in edges):
                    edges[pln[0]] = {}
                if pln[1] not in edges[pln[0]]:
                    edges[pln[0]][pln[1]] = {}
                count += 1
                edges[pln[0]][pln[1]][pln[2]] = 1
        print(str(count) + " edges processed")
        return edges, clusters, nodep

    def getPath(G, s, p):
        visited = {}
        visited[s] = 1
        res = [s]
        l1 = list(G.neighbors(s))
        l1 = list(set(p).intersection(l1)) + list(set(l1).difference(p))
        while (len(l1) > 0):
            if l1[0] in visited:
                break
            visited[l1[0]] = 1
            res += [l1[0]]
            l1 = list(G.neighbors(l1[0]))
            l1 = list(set(p).intersection(l1)) + list(set(l1).difference(p))
        return res

    def getBINGraph(ana, edges, clusters, dr):
        import networkx as nx
        from networkx.drawing.nx_agraph import write_dot, graphviz_layout

        sclusters = {k:len(clusters[k]) for k in clusters}
        dclusters = {k:dr.loc[k] for k in clusters}
        def getW(k):
            return (-sclusters[k], -dclusters[k])
        def getS(k):
            return np.log(sclusters[k]+1)/np.log(1.5) + 2
        keys = sorted(sclusters, key=lambda k:getW(k))

        net = nx.DiGraph()
        for id1 in edges:
            if '2' in edges[id1]:
                l1 = sorted(edges[id1]['2'], key=lambda k:getW(k))
                for id2 in l1:
                    net.add_edge(id1, id2, rel='2')
            if 2 in edges[id1]:
                l1 = sorted(edges[id1][2], key=lambda k:getW(k))
                for id2 in l1:
                    net.add_edge(id1, id2, rel='2')

        G = nx.DiGraph()
        for id1 in keys[0:10]:
            l1 = []
            if id1 in edges and '4' in edges[id1]:
                l1 += list(edges[id1]['4'].keys())
            if id1 in edges and '6' in edges[id1]:
                l1 += list(edges[id1]['6'].keys())
            if id1 in edges and 4 in edges[id1]:
                l1 += list(edges[id1][4].keys())
            if id1 in edges and 6 in edges[id1]:
                l1 += list(edges[id1][6].keys())
            l1 = list(set(l1))
            l2 = sorted(l1, key=lambda k:getW(k))[0:2]
            print (sclusters[id1], dclusters[id1], l2)
            for id2 in l2:
                G.add_node(id1, label=ana.h.getSimpleName(id1), 
                           size=getS(id1), title='hilo', group=1)
                G.add_node(id2, label=ana.h.getSimpleName(id2),
                           size=getS(id2), title='hilo', group=1)
                G.add_edge(id1, id2, rel='4', color='red')

        l1 = list(G)
        for id1 in l1:
            l2 = BIGraph.getPath(net, id1, l1)
            print(l2)
            t1 = id1
            for id2 in l2[1:]:
                G.add_node(t1, label=ana.h.getSimpleName(t1), 
                           size=getS(t1), title='lolo', group=2)
                G.add_node(id2, label=ana.h.getSimpleName(id2),
                           size=getS(id2), title='lolo', group=2)
                G.add_edge(t1, id2, rel='2', color='blue',
                          arrows={'to':{'enabled':True, 'type':'arrow'}})
                t1 = id2
                
        return G

    def visualizeNetwork(G):
        from pyvis.network import Network
        nt = Network("500px", "500px", heading="BoNE", notebook=True)
        nt.from_nx(G)
        return nt.show("network/nx.html")

    def getBINGraphGML(ana, edges, clusters, dr):
        import networkx as nx
        from networkx.drawing.nx_agraph import write_dot, graphviz_layout

        sclusters = {k:len(clusters[k]) for k in clusters}
        dclusters = {k:dr.loc[k] for k in clusters}
        def getW(k):
            return (-sclusters[k], -dclusters[k])
        def getS(k):
            return np.log(sclusters[k]+1)/np.log(1.5) + 2
        keys = sorted(sclusters, key=lambda k:getW(k))

        net = nx.DiGraph()
        for id1 in edges:
            if '2' in edges[id1]:
                l1 = sorted(edges[id1]['2'], key=lambda k:getW(k))
                for id2 in l1:
                    net.add_edge(id1, id2, rel='2')
            if 2 in edges[id1]:
                l1 = sorted(edges[id1][2], key=lambda k:getW(k))
                for id2 in l1:
                    net.add_edge(id1, id2, rel='2')

        G = nx.DiGraph()
        for id1 in keys[0:10]:
            l1 = []
            if id1 in edges and '4' in edges[id1]:
                l1 += list(edges[id1]['4'].keys())
            if id1 in edges and '6' in edges[id1]:
                l1 += list(edges[id1]['6'].keys())
            if id1 in edges and 4 in edges[id1]:
                l1 += list(edges[id1][4].keys())
            if id1 in edges and 6 in edges[id1]:
                l1 += list(edges[id1][6].keys())
            l1 = list(set(l1))
            l2 = sorted(l1, key=lambda k:getW(k))[0:2]
            print (sclusters[id1], dclusters[id1], l2)
            for id2 in l2:
                G.add_node(id1, name=ana.h.getSimpleName(id1),
                           graphics = {'x': 0, 'y': 0, 'w': getS(id1), 'h': getS(id1),
                                       'type': 'ellipse', 'fill': '#889999', 'outline': '#666666',
                                       'outline_width': 1.0})
                G.add_node(id2, name=ana.h.getSimpleName(id2),
                           graphics = {'x': 0, 'y': 0, 'w': getS(id2), 'h': getS(id2),
                                       'type': 'ellipse', 'fill': '#889999', 'outline': '#666666',
                                       'outline_width': 1.0})
                G.add_edge(id1, id2, rel='4', 
                           graphics={'width': 1.0, 'fill': '#ff0000', 'type': 'line',
                                    'source_arrow': 0, 'target_arrow': 0})

        l1 = list(G)
        for id1 in l1:
            l2 = BIGraph.getPath(net, id1, l1)
            print(l2)
            t1 = id1
            for id2 in l2[1:]:
                G.add_node(t1, name=ana.h.getSimpleName(t1), 
                           graphics = {'x': 0, 'y': 0, 'w': getS(t1), 'h': getS(t1),
                                       'type': 'ellipse', 'fill': '#889999', 'outline': '#666666',
                                       'outline_width': 1.0})
                G.add_node(id2, name=ana.h.getSimpleName(id2),
                           graphics = {'x': 0, 'y': 0, 'w': getS(id2), 'h': getS(id2),
                                       'type': 'ellipse', 'fill': '#889999', 'outline': '#666666',
                                       'outline_width': 1.0})
                G.add_edge(t1, id2, rel='2', 
                          graphics={'width': 1.0, 'fill': '#0000ff', 'type': 'line',
                                    'source_arrow': 0, 'target_arrow': 1})
                t1 = id2
        return G

    def writeGML(G, ofile="network/nx.gml"):
        import networkx as nx
        nx.write_gml(G, ofile)

    def getDiff(ana, l1):
        res = []
        for k in l1:
            for u in ana.h.getIDs(k):
                expr = ana.h.getExprData(u)
                v1 = [float(expr[i]) for i in ana.state[0]]
                v2 = [float(expr[i]) for i in ana.state[1]]
                t, p = ttest_ind(v1,v2, equal_var=False)
                m1 = np.mean(v1)
                m2 = np.mean(v2)
                diff = m2 - m1
                res.append([k, u, m1, m2, diff, t, p])
        cl = ['Name', 'ID', 'm1', 'm2', 'diff', 't', 'p']
        return pd.DataFrame(res, columns=cl)

    def saveClusters(ofile, cls):
        fp = open(ofile, "w")
        for k in sorted(cls, key=lambda x: cls[x][0], reverse=True):
            l1 = [k, cls[k][0]] + cls[k][1]
            l1 = [str(k) for k in l1]
            fp.write("\t".join(l1) + "\n")
        fp.close()
        return
    def readClusters(cfile):
        cls = {}
        fp = open(cfile, "r")
        for line in fp:
            l1 = line.strip().split("\t")
            cls[l1[0]] = [int(l1[1]), l1[2:]]
        fp.close()
        return cls
    def getVolcano(ana, cfile, genes):
        h = ana.h
        list1 = []
        for g in genes:
            id1 = h.getBestID(h.getIDs(g).keys())
            if id1 is None:
                print (g)
                continue
            list1 += [id1]
        idlist = bone.getEntries(cfile, 0)
        idhash = {}
        for i in range(len(idlist)):
            idhash[idlist[i]] = i
        order = [idhash[k] for k in list1]
        pval = [float(i) for i in bone.getEntries(cfile, 3)]
        fc = [float(i) for i in bone.getEntries(cfile, 4)]
        c = ["black" for i in fc]
        df = pd.DataFrame()
        df["-log10(p)"] = -np.log10(pval)
        df["log(FC)"] = fc
        for i in range(len(pval)):
            if df["-log10(p)"][i] >= -np.log10(0.05) and abs(df["log(FC)"][i]) >= 1:
                c[i] = "green"
            if df["-log10(p)"][i] >= -np.log10(0.05) and abs(df["log(FC)"][i]) < 1:
                c[i] = "red"
            if df["-log10(p)"][i] < -np.log10(0.05) and abs(df["log(FC)"][i]) >= 1:
                c[i] = "orange"
        df["color"] = c
        df.dropna(inplace = True)
        ax = df.plot.scatter("log(FC)", "-log10(p)", c = df["color"],
                s=8, alpha=0.2, figsize=(6, 4))
        ax.set_xlim([-9, 9])
        ax.set_ylim([0, 9])
        for i in order:
            ax.text(fc[i], -np.log10(pval[i]),  h.getSimpleName(idlist[i]),
                        horizontalalignment='left', verticalalignment='top')
            ax.plot([fc[i]], [-np.log10(pval[i])], "bo", ms=4)
        return ax
