import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import numpy as np
import re
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import sys
sys.path.append("../Hegemon/")
sys.path.append("../")
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


#plt.rc('text', usetex=True)

def getPDF(cfile):
    import bone
    reload(bone)
    from matplotlib.backends.backend_pdf import PdfPages

    pdf = PdfPages(cfile)
    return pdf

def closePDF(pdf):
    import datetime
    d = pdf.infodict()
    d['Title'] = 'Analysis of Merck compound'
    d['Author'] = 'Debashis Sahoo'
    d['Subject'] = 'COVID-19'
    d['Keywords'] = 'disease training validation ROC'
    d['CreationDate'] = datetime.datetime(2020, 1, 16)
    d['ModDate'] = datetime.datetime.today()
    pdf.close()

class Analysis(bone.IBDAnalysis):

    def __init__(self):
        bone.IBDAnalysis.__init__(self)

    def getPG2020LungHamNamir(self, tn=1):
        self.prepareData("COV323", "explore.conf")
        atype = self.h.getSurvName('c info')
        atypes = ['U', 'V', 'N', '3', '4']
        ahash = {'Merck Drug':2, 'Merck Drug Vehicle Control':1, 'UN':0}
        if (tn == 2):
            atypes = ['U', 'V', 'N']
            ahash = {'Merck Drug':2, 'Merck Drug Vehicle Control':1, 'UN':0}
            ah = {'3', '4'}
            atype = [None if k in ah else k for k in atype]
        self.initData(atype, atypes, ahash)


def T1():
    pdf = getPDF("namir.pdf")
    import bone
    reload(bone)
    ana = Analysis()
    ana.getPG2020LungHamNamir(2)
    ana.printMeanAbsoluteDeviation("namir-mad.txt")
    df = pd.read_csv("namir-mad.txt", sep="\t")
    thr = hu.getThrData(list(df['MAD']))
    print(len(df[df['MAD'] > thr[2]]['ArrayID']))
    thr = hu.getThrData(list(df[df['MAD'] > thr[2]]['MAD']))
    print(len(df[df['MAD'] > thr[2]]['ArrayID']))
    l1 = [list(df[df['MAD'] > thr[2]]['ArrayID'])]
    ranks, row_labels, row_ids, row_numhi, expr = bone.getRanks3(l1,ana.h,
            ana.order)
    row_labels = [re.sub("_", "\\_", k) for k in row_labels]
    df = pd.DataFrame(expr, columns=[ana.atypes[ana.aval[i]] for i in ana.order])
    df.index = row_labels
    cvals  = [np.min(np.min(df)), 0, np.max(np.max(df))]
    clrs = ["blue","#FFFFFF", "red"]
    norm=plt.Normalize(min(cvals),max(cvals))
    tuples = list(zip(map(norm,cvals), clrs))
    import matplotlib.colors as colors
    cmap = colors.LinearSegmentedColormap.from_list("BGYR1", tuples)
    plt.register_cmap(cmap=cmap)
    sns.set(font_scale=1.0)
    g = sns.clustermap(df, cmap=cmap)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    pdf.savefig(transparent=True)
    closePDF(pdf)

T1()
