import requests
import sys
import argparse
import json
import numpy as np

import os
import os.path
import re
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import bone
import cStringIO
import base64

acolor = ["#00CC00", "#D8A03D","#EC008C",
        'cyan', "#B741DC", "#808285",
        'blue', 'black', 'green', 'red']
acolor = ["#00CC00", "#EFF51A","#EC008C",
        'cyan', "#B741DC", "#808285",
        'blue', 'black', 'green', 'red']

ap = argparse.ArgumentParser()
ap.add_argument('-c','--cmd', help = 'Command')
ap.add_argument('-i','--id', help = 'id')
ap.add_argument('-p','--pdf', dest='otype', action='store_const',
        const="pdf", default="jpg", help = 'pdf')
args = vars(ap.parse_args())

cmd = ""
id1 = ""
otype = "jpg"

if (args['cmd']):
    cmd = args['cmd'].strip()
if (args['id']):
    id1 = args['id']
if (args['otype']):
    otype = args['otype']

def printPeters(gene):
    ana = bone.IBDAnalysis()
    ana.getPeters()
    gn = []
    gi = []
    for g in ana.h.getIDs(gene):
        expr = ana.h.getExprData(g)
        m1 = np.mean([float(expr[i]) for i in ana.normal])
        m2 = np.mean([float(expr[i]) for i in ana.ibd])
        if m1 > m2:
            gn += [g]
        else:
            gi += [g]
    l1 = [gn, gi]
    wt1 = [-1, 1]
    c_dict, fpr, tpr, roc_auc = bone.processGeneGroups(ana, l1, wt1)
    if (otype == "jpg"):
        params = {'spaceAnn': len(ana.order)/len(ana.atypes),
                'tAnn': 1, 'widthAnn':6, 'acolor': acolor,
                'w': 5, 'h': 0.8, 'atypes': ana.atypes ,'cval': ana.cval}
        ax = ana.printTitleBar(params)
        res = ana.getMetrics2(ana.cval[0], fthr="thr2");
        my_stringIObytes = cStringIO.StringIO()
        ax.get_figure().savefig(my_stringIObytes, format='jpg')
        my_stringIObytes.seek(0)
        bar = base64.b64encode(my_stringIObytes.read())
        ax = ana.densityPlot();
        my_stringIObytes = cStringIO.StringIO()
        ax.get_figure().savefig(my_stringIObytes, format='jpg')
        my_stringIObytes.seek(0)
        density = base64.b64encode(my_stringIObytes.read())
        data = { "roc-auc": roc_auc, "accuracy": res[1],
                "fisher": res[2], "roc": res[0], "bar": bar,
                "density": density, "source": ana.h.getSource(),
                "title": ana.h.getTitle(), "dbid": ana.dbid,
                "key": "polyps", "otype": otype}
    if (otype == "pdf"):
        atypes = []
        for i in range(len(ana.atypes)):
            l = [k for k in ana.aval if k == i]
            atypes += [str(ana.atypes[i]) + "("+str(len(l))+")"]
        ana.atypes = atypes
        fig = plt.figure(figsize=(4,4), dpi=100)
        plt.subplots_adjust(hspace=0.5, wspace=0.5)
        ax1 = plt.subplot2grid((4, 1), (0, 0))
        ax2 = plt.subplot2grid((4, 1), (1, 0), rowspan=3)
        params = {'spaceAnn': len(ana.order)/len(ana.atypes),
                'tAnn': 1, 'widthAnn':6, 'acolor': acolor, 'ax':ax1,
                'w': 5, 'h': 0.8, 'atypes': atypes ,'cval': ana.cval}
        ax = ana.printTitleBar(params)
        res = ana.getMetrics2(ana.cval[0], fthr="thr2");
        #ax.text(len(ana.cval[0]), 4, ",".join(res))
        #ax = ana.densityPlot(ax2);
        params['ax'] = ax2
        params['vert'] = 0
        ax = ana.printViolin(None, params)
        ax.set_title(",".join(res))
        my_stringIObytes = cStringIO.StringIO()
        fig.savefig(my_stringIObytes, format='pdf')
        my_stringIObytes.seek(0)
        density = base64.b64encode(my_stringIObytes.read())
        data = { "roc-auc": roc_auc, "accuracy": res[1],
                "fisher": res[2], "roc": res[0], "otype": otype,
                "density": density, "source": ana.h.getSource(),
                "title": ana.h.getTitle(), "dbid": ana.dbid,
                "key": "polyps"}
    print "\n--data--"
    print json.dumps(data);

def getGeneSet(gene):
    ana = bone.IBDAnalysis()
    ana.getPeters()
    gn = []
    gi = []
    for g in ana.h.getIDs(gene):
        expr = ana.h.getExprData(g)
        m1 = np.mean([float(expr[i]) for i in ana.normal])
        m2 = np.mean([float(expr[i]) for i in ana.ibd])
        if m1 > m2:
            gn += [g]
        else:
            gi += [g]
    l1 = [gn, gi]
    wt1 = [-1, 1]

    return l1, wt1

def printOutcome(gene):
    l1, wt1 = getGeneSet(gene)
    data = { "roc-auc": [], "accuracy": [], "fisher": [], "roc": [], "bar": [],
            "source": [], "title": [], "dbid": [], "key": []}
    index = 0
    if otype == "pdf":
        fig = plt.figure(figsize=(12,6), dpi=100)
        n1 = 7
        axes = []
        for i in range(n1):
            ax = plt.subplot2grid((n1, 2), (i, 1))
            axes.extend([ax])
    def processHs(ana, data):
        idx = index
        c_dict, fpr, tpr, roc_auc = bone.processGeneGroups(ana, l1, wt1)
        params = {'spaceAnn': len(ana.order)/len(ana.atypes),
                'tAnn': 1, 'widthAnn':1, 'acolor': acolor,
                'w': 5, 'h': 0.8, 'atypes': ana.atypes ,'cval': ana.cval}
        if (otype == "pdf"):
            params['ax'] = axes[index]
        ax = ana.printTitleBar(params)
        res = ana.getMetrics2(ana.cval[0], fthr="thr2");
        if (otype == "pdf"):
            ax.text(-0.7, 0.5, ",".join(res), horizontalalignment='left',
                    verticalalignment='center', transform=ax.transAxes)
            desc = ana.h.getSource()
            ax.text(-0.1, 0.5, desc, horizontalalignment='right',
                    verticalalignment='center', transform=ax.transAxes)
            desc = ana.h.rdataset.getName() + " (n = " + str(ana.h.getNum()) + ")"
            ax.text(-1.3, 0.5, desc, horizontalalignment='left',
                    verticalalignment='center', transform=ax.transAxes)
        if (otype == "jpg"):
            my_stringIObytes = cStringIO.StringIO()
            ax.get_figure().savefig(my_stringIObytes, format='jpg')
            my_stringIObytes.seek(0)
            bar = base64.b64encode(my_stringIObytes.read())
            data["bar"] += [bar]
        data["roc-auc"] += [float("%.2f" % roc_auc)]
        data["accuracy"] += [res[1]]
        data["fisher"] += [res[2]]
        data["roc"] += [res[0]]
        data["source"] += [ana.h.getSource()]
        data["title"] += [ana.h.getTitle()]
        data["dbid"] += [ana.dbid]
        data["key"] += ["polyps"]

    ana = bone.IBDAnalysis()
    ana.getArijs2018(2)
    processHs(ana, data)
    index += 1

    ana = bone.IBDAnalysis()
    ana.getArijs2009(tn=2)
    processHs(ana, data)
    index += 1

    ana = bone.IBDAnalysis()
    ana.getVanhove(tn=2)
    processHs(ana, data)
    index += 1

    ana = bone.IBDAnalysis()
    ana.getVanderGoten(2)
    processHs(ana, data)
    index += 1

    ana = bone.IBDAnalysis()
    ana.getDePreter(4)
    processHs(ana, data)
    index += 1

    ana = bone.IBDAnalysis()
    ana.getPekow(3)
    processHs(ana, data)
    index += 1

    ana = bone.IBDAnalysis()
    ana.getVerstockt2019()
    processHs(ana, data)
    index += 1

    if otype == "pdf":
        my_stringIObytes = cStringIO.StringIO()
        fig.savefig(my_stringIObytes, format='pdf')
        my_stringIObytes.seek(0)
        bar = base64.b64encode(my_stringIObytes.read())
        data["outcome"] = bar
    data["otype"] =  otype

    print "\n--data--"
    print json.dumps(data);

def printGender(gene):
    ana = bone.IBDAnalysis()
    ana.getPeters()
    atype = ana.h.getSurvName("c gender")
    ahash = {'male':0, 'female': 1}
    aval = [ahash[i] if i in ahash else None for i in atype]
    male = [i for i in ana.h.aRange() if aval[i] == 0]
    female = [i for i in ana.h.aRange() if aval[i] == 1]
    lval = []
    atypes = []
    alist = []
    for g in ana.h.getIDs(gene):
        expr = ana.h.getExprData(g)
        lval += [ [ float(expr[i]) for i in male ] ]
        atypes += [ "male" ]
        lval += [ [ float(expr[i]) for i in female ] ]
        atypes += [ "female" ]
        alist += [g]
    params = {'vert': 1}
    ax,bp = bone.plotScores(lval, atypes, params)
    for i in range(len(alist)):
        ax.text(i* 2 + 1.5, ax.get_ylim()[0] + 1, alist[i],
                horizontalalignment='center', verticalalignment='center')
    if otype == "jpg":
        my_stringIObytes = cStringIO.StringIO()
        ax.get_figure().savefig(my_stringIObytes, format='jpg')
        my_stringIObytes.seek(0)
        bar = base64.b64encode(my_stringIObytes.read())
    if otype == "pdf":
        my_stringIObytes = cStringIO.StringIO()
        ax.get_figure().savefig(my_stringIObytes, format='pdf')
        my_stringIObytes.seek(0)
        bar = base64.b64encode(my_stringIObytes.read())
    data = {"bar": bar, "atypes":atypes, "alist":alist,
            "lval": lval, "source": ana.h.getSource(),
            "title": ana.h.getTitle(), "dbid": ana.dbid,
            "key": "polyps"}
    data["otype"] =  otype
    print "\n--data--"
    print json.dumps(data);

def printMouse(gene):
    l1, wt1 = getGeneSet(gene)
    data = { "roc-auc": [], "accuracy": [], "fisher": [], "mtype": [], "model": [],
            "source": [], "title": [], "dbid": [], "key": [], "bar": None}
    genes = []
        
    fig = plt.figure(figsize=(10,5), dpi=100)
    n1 = 11
    axlist = []
    for i in range(n1):
        ax = plt.subplot2grid((n1, 3), (i, 2))
        axlist.extend([ax])

    def processHeader(ax):
        nAt = 10
        extent = [0, nAt, 0, 5]
        ax.axis(extent)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.tick_params(top=False, left=False, bottom=False, right=False)
        for edge, spine in ax.spines.items():
                    spine.set_visible(False)
        ax.text(-2.9, 0.5, "Model Type",
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
        ax.text(-2.5, 0.5, "Mouse Model",
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
        ax.text(-1.9, 0.5, "Dataset",
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
        ax.text(-1.35, 0.5, "#Samples\n(N, Colitis)",
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
        ax.text(-0.9, 0.5, "ROC-\nAUC",
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
        ax.text(-0.65, 0.5, "Accu",
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
        ax.text(-0.4, 0.5, "Fisher p",
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
        
    def processMm(ana, genes, wt1, l1, ax, desc1, desc2,
            acolor = acolor, data = data, fthr = "thr0"):
        ana.convertMm(l1, genes)
        ana.orderData(ana.gene_groups, wt1)
        nAt = len(ana.cval[0])
        extent = [0, nAt, 0, 5]
        ax.axis(extent)
        cmap = colors.ListedColormap(acolor)
        boundaries = range(len(acolor) + 1)
        norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)
        ax.imshow(ana.cval, interpolation='nearest', cmap=cmap, \
                          norm=norm, extent=extent, aspect="auto")
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.tick_params(top=False, left=False, bottom=False, right=False)
        ax.set_xticks(np.arange(0, nAt, 1))
        ax.grid(which='major', alpha=0.2, linestyle='-', linewidth=0.5,
                color='black')
        for edge, spine in ax.spines.items():
                    spine.set_visible(False)
        res = ana.getMetrics2(ana.cval[0], fthr=fthr)
        ax.text(-2.9, 0.5, desc1,
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
        ax.text(-2.5, 0.5, desc2,
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
        ax.text(-1.9, 0.5, ana.h.getSource(),
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
        ax.text(-1.35, 0.5, "(" + str(len(ana.normal)) + ", " + str(len(ana.uc)) + ")",
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
        ax.text(-0.9, 0.5, res[0],
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
        ax.text(-0.65, 0.5, res[1],
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
        ax.text(-0.4, 0.5, res[2],
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
        data["roc-auc"] += [res[0]]
        data["accuracy"] += [res[1]]
        data["fisher"] += [res[2]]
        data["source"] += [ana.h.getSource()]
        data["title"] += [ana.h.getTitle()]
        data["dbid"] += [ana.dbid]
        data["key"] += ["polyps"]
        data["mtype"] += [desc1]
        data["model"] += [desc2]

    processHeader(axlist[0])
    ana = bone.IBDAnalysis()
    ana.getBreynaert2013(2)
    desc1 = "Chemical"; desc2 = "DSS (bulk)"
    processMm(ana, genes, wt1, l1, axlist[1], desc1, desc2)

    ana = bone.IBDAnalysis()
    ana.getJensen2017()
    desc1 = "Chemical"; desc2 = "DSS (epithelium)"
    processMm(ana, genes, wt1, l1, axlist[2], desc1, desc2)

    ana.getDohi2014(2)
    desc1 = "Chemical"; desc2 = "TNBS"
    processMm(ana, genes, wt1, l1, axlist[3], desc1, desc2)

    ana.getLamas2018(2)
    desc1 = "Infection"; desc2 = "Citrobacter"
    processMm(ana, genes, wt1, l1, axlist[4], desc1, desc2)

    ana.getLyons2018(2)
    desc1 = "T cell"; desc2 = "ACT"
    processMm(ana, genes, wt1, l1, axlist[5], desc1, desc2)

    ana.getFang2012(2)
    desc1 = "T cell"; desc2 = "ACT"
    processMm(ana, genes, wt1, l1, axlist[6], desc1, desc2)

    ana.getRuss2013(2)
    desc1 = "Genetic"; desc2 = "IL10 -/- (epithelium)"
    processMm(ana, genes, wt1, l1, axlist[7], desc1, desc2)

    ana.getRuss2013(3)
    desc1 = "Genetic"; desc2 = "IL10 -/- (bulk)"
    processMm(ana, genes, wt1, l1, axlist[8], desc1, desc2)

    ana.getTam2019()
    desc1 = "Genetic"; desc2 = "TNFR1 -/-"
    processMm(ana, genes, wt1, l1, axlist[9], desc1, desc2, fthr = "thr2")

    ana.getPunit2015()
    desc1 = "Genetic"; desc2 = "TNFR2 -/-"
    processMm(ana, genes, wt1, l1, axlist[10], desc1, desc2, fthr = "thr2")

    my_stringIObytes = cStringIO.StringIO()
    if otype == "jpg":
        fig.savefig(my_stringIObytes, format='jpg')
    if otype == "pdf":
        fig.savefig(my_stringIObytes, format='pdf')
    my_stringIObytes.seek(0)
    bar = base64.b64encode(my_stringIObytes.read())
    data["bar"] = bar
    data["otype"] =  otype
    print "\n--data--"
    print json.dumps(data);

def printColonTissue(gene):
    l1, wt1 = getGeneSet(gene)
    data = {}
    for i in range(len(l1)):
        for name in l1[i]:
            data[name] = { "top": 0, "bottom": 0, "macrophage": 0,
                          "lymphocyte": 0, "fibroblast": 0}
    k20info = bone.getBooleanInfo("k20-boolean.txt", "CRC115", l1)
    for i in range(len(l1)):
        for name in l1[i]:
            if name in k20info:
                if k20info[name][2][0] > 3 and k20info[name][2][1] < 0.16:
                    data[name]["top"] = 1
    cinfo = bone.getGeneInfo("CRC102", l1)
    rinfo = bone.getCorrInfo("jung-2015-sc-lgr5-corr.txt", 0.8, l1)
    for i in range(len(l1)):
        for name in l1[i]:
            if name in rinfo:
                if rinfo[name][0] > 0.8:
                    data[name]["bottom"] = 1
            if name in k20info and name in cinfo:
                if k20info[name][1][0] > 3 and k20info[name][1][1] < 0.16 \
                and cinfo[name][0] > 3:
                    data[name]["bottom"] = 1
    cinfo = bone.getInfo("GL1")
    rinfo = bone.getGeneInfo("MAC59", l1)
    for i in range(len(l1)):
        for name in l1[i]:
            if name in rinfo:
                id1 = rinfo[name][3]
                if rinfo[name][2] > float(cinfo[id1][2]) + 0.5:
                    data[name]["macrophage"] = 1
    rinfo = bone.getGeneInfo("PLP114", l1)
    for i in range(len(l1)):
        for name in l1[i]:
            if name in rinfo:
                id1 = rinfo[name][3]
                if rinfo[name][2] > float(cinfo[id1][2]) + 0.5:
                    data[name]["fibroblast"] = 1
    rinfo = bone.getTCell(l1)
    for i in range(len(l1)):
        for name in l1[i]:
            if name in rinfo:
                id1 = rinfo[name][3]
                if rinfo[name][0] > 0.5 and rinfo[name][1] > 0.5:
                    data[name]["lymphocyte"] = 1
    rinfo = bone.getBCell(l1)
    for i in range(len(l1)):
        for name in l1[i]:
            if name in rinfo:
                id1 = rinfo[name][3]
                if rinfo[name][0] > 0.5 and rinfo[name][1] > 0.5:
                    data[name]["lymphocyte"] = 1
    print "\n--data--"
    print json.dumps(data);
        
if (cmd == 'ibd-1'):
    printPeters(id1)
if (cmd == 'ibd-outcome'):
    printOutcome(id1)
if (cmd == 'ibd-gender'):
    printGender(id1)
if (cmd == 'ibd-mouse'):
    printMouse(id1)
if (cmd == 'colon-tissue'):
    printColonTissue(id1)
