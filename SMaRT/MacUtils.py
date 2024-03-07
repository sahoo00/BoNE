import os
import sys
import pandas as pd
import requests
import json
try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3

import numpy as np
import re
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
from sklearn.metrics import *
from scipy.stats import fisher_exact, ttest_ind
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import matplotlib.patches as patches
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

try:
    reload  # Python 2.7
except NameError:
    try:
        from importlib import reload  # Python 3.4+
    except ImportError:
        from imp import reload  # Python 3.0 - 3.3

acolor = ["#00CC00", "#D8A03D","#EC008C",
        'cyan', "#B741DC", "#808285",
        'blue', 'black', 'green', 'red',
        'orange', 'brown', 'pink', 'purple']

basedir = os.path.join(os.getcwd(), os.path.dirname(__file__)) + "/"
def getRealpath(cfile):
    return os.path.realpath(os.path.join(os.getcwd(), 
        os.path.dirname(__file__), cfile))

urlbase="http://hegemon.ucsd.edu/Tools/explore.php"
class HegemonUtil:
    @staticmethod
    def uniq(mylist):
      used = set()
      unique = [x for x in mylist if x not in used and (used.add(x) or True)]
      return unique
    @staticmethod
    def getHegemonDataset(dbid, urlbase=urlbase):
        url = urlbase + "?go=getdatasetjson&id=" + dbid
        response = requests.get(url)
        obj = json.loads(response.text)
        return  obj
    @staticmethod
    def getHegemonPatientInfo(dbid, urlbase=urlbase):
        url = urlbase + "?go=getpatientinfojson" + "&id=" + dbid
        response = requests.get(url)
        obj = json.loads(response.text)
        return  obj
    @staticmethod
    def getHegemonPatientData(dbid, name, urlbase=urlbase):
        hdr = HegemonUtil.getHegemonPatientInfo(dbid, urlbase)
        clinical = 0
        if name in hdr:
            clinical = hdr.index(name)
        url = urlbase + "?go=getpatientdatajson" + \
            "&id=" + dbid + "&clinical=" + str(clinical)
        response = requests.get(url)
        obj = json.loads(response.text)
        return  obj
    @staticmethod
    def getHegemonDataFrame(dbid, genelist=None, pGroups=None, urlbase=urlbase):
        genes =''
        if genelist is not None:
            genes = ' '.join(genelist)
        groups = ''
        if pGroups is not None:
            for i in range(len(pGroups)):
                str1 = "=".join([str(i), pGroups[i][0], ":".join(pGroups[i][2])])
                if i == 0:
                    groups += str1
                else:
                    groups = groups + ';' + str1
        url = urlbase
        opt = {'go': 'dataDownload', 'id': dbid, 'genes': genes, 'groups' : groups}
        response = requests.post(url, opt)
        data = StringIO(response.text)
        df = pd.read_csv(data, sep="\t")
        return df
    @staticmethod
    def getHegemonThrFrame(dbid, genelist=None, urlbase=urlbase):
        genes =''
        if genelist is not None:
            genes = ' '.join(genelist)
        url = urlbase
        opt = {'go': 'dataDownload', 'id': dbid, 'genes': genes, 'groups' : '',
              'param': 'type:thr'}
        response = requests.post(url, opt)
        data = StringIO(response.text)
        df = pd.read_csv(data, sep="\t", header=None)
        df.columns=['ProbeID', 'thr1', 'stat', 'thr0', 'thr2']
        return df
    @staticmethod
    def getHegemonGeneIDs(dbid, genelist=None, urlbase=urlbase):
        genes =''
        if genelist is not None:
              genes = ' '.join(genelist)
        url = urlbase
        opt = {'go': 'dataDownload', 'id': dbid, 'genes': genes, 'groups' : '',
              'param': 'type:geneids'}
        response = requests.post(url, opt)
        data = StringIO(response.text)
        df = pd.read_csv(data, sep="\t")
        return df
    @staticmethod
    def getHegemonPlots(dbid, gA, gB, urlbase=urlbase):
      url = urlbase + "?go=getplotsjson&id=" + dbid + \
              "&A=" + str(gA) + "&B=" + str(gB)
      response = requests.get(url)
      obj = json.loads(response.text)
      return  obj
    @staticmethod
    def getHegemonData(dbid, gA, gB, urlbase=urlbase):
      url = urlbase + "?go=getdatajson&id=" + dbid + \
              "&A=" + str(gA) + "&B=" + str(gB)
      response = requests.get(url)
      obj = json.loads(response.text)
      return  obj
    @staticmethod
    def censor(time, status, ct):
      t = [time[i] if i < 2 or time[i] == '' or float(time[i]) < ct \
              else ct for i in range(len(time))]
      s = [status[i] if i < 2 or time[i] == '' or float(time[i]) < ct \
              else 0 for i in range(len(time))]
      return t,s
    @staticmethod
    def survival(time, status, pGroups=None, ax=None, atrisk=1):
      from lifelines import KaplanMeierFitter
      if pGroups is None:
        order = [i for i in range(2, len(time)) 
                    if time[i] != "" and status[i] != ""]
        t = [float(time[i]) for i in order]
        s = [int(status[i]) for i in order]
        kmf = KaplanMeierFitter()
        kmf.fit(t, s)
        ax = kmf.plot(ax = ax, color='red')
        return ax
      else:
        t1 = []
        s1 = []
        g1 = []
        kmfs = []
        for k in range(len(pGroups)):
          df = pd.DataFrame()
          order = [i for i in pGroups[k][2]
                   if time[i] != "" and status[i] != ""]
          if len(order) <= 0:
              continue
          t = [float(time[i]) for i in order]
          s = [int(status[i]) for i in order]
          g = [k for i in order]
          t1 += t
          s1 += s
          g1 += g
          kmf = KaplanMeierFitter()
          kmf.fit(t, s, label = pGroups[k][0])
          if ax is None:
            ax = kmf.plot(color=pGroups[k][1], ci_show=False, show_censors=True,
                    figsize=(4,4))
          else:
            ax = kmf.plot(ax = ax, color=pGroups[k][1], ci_show=False,
                    show_censors=True)
          kmfs += [kmf]
        if atrisk == 1:
            from lifelines.plotting import add_at_risk_counts
            add_at_risk_counts(*kmfs, ax=ax)
        if len(t1) > 0:
          from lifelines.statistics import multivariate_logrank_test
          from matplotlib.legend import Legend
          res = multivariate_logrank_test(t1, g1, s1)
          leg = Legend(ax, [], [], title = "p = %.2g" % res.p_value,
                       loc='lower left', frameon=False)
          ax.add_artist(leg);
        return ax
    @staticmethod
    def logrank(time, status, pGroups):
        groups = [ "" for i in time]
        for k in range(len(pGroups)):
          df = pd.DataFrame()
          order = [i for i in pGroups[k][2]
                   if time[i] != "" and status[i] != ""]
          if len(order) <= 0:
              continue
          for i in order:
            groups[i] = k
        order = [i for i in range(len(groups)) if groups[i] != ""]
        if len(order) > 0:
          t = [float(time[i]) for i in order]
          s = [int(status[i]) for i in order]
          g = [int(groups[i]) for i in order]
          from lifelines.statistics import multivariate_logrank_test
          from matplotlib.legend import Legend
          res = multivariate_logrank_test(t, g, s)
          return res.p_value
        return 1.0
    @staticmethod
    def getBestThr(time, status, value, order, vrange=None, ct=None, tn=3):
        if ct is not None:
            time, status = HegemonUtil.censor(time, status, ct)
        order = [i for i in order if value[i] is not None and
                value[i] != '' and value[i] != 'NA' and value[i] != "null"]
        vals = sorted([float(value[i]) for i in order])
        p = []
        thr = []
        for fthr in vals:
            if vrange is not None and len(vrange) > 1 and fthr < vrange[0]:
                continue
            if vrange is not None and len(vrange) > 1 and fthr > vrange[1]:
                continue
            g1 = [i for i in order if float(value[i]) < fthr]
            g2 = [i for i in order if float(value[i]) >= fthr]
            pG = [ ["Low", "red", g1], ["High", "green", g2]]
            pv = HegemonUtil.logrank(time, status, pG)
            if not np.isnan(pv):
                p.append(pv)
                thr.append(fthr)
        iord = np.argsort(p)[0:tn]
        return [[p[i], thr[i]] for i in iord]

    @staticmethod
    def getThrData(arr, start = None, length = None):
      if start is None:
        start = 0
      if length is None:
        length = len(arr)
      s_thr = getStepMinerThr(arr, start, start+length-1)
      if s_thr is None:
          return None
      thr = s_thr["threshold"]
      stat = s_thr["statistic"]
      return [thr, stat, thr-0.5, thr+0.5]

    @staticmethod
    def getThrCode (thr_step, value, code = None):
      if (code is None):
          return value 
      thr = value;
      if (code == "thr1"):
        thr = thr_step[0];
      elif (code == "thr0"):
        thr = thr_step[2];
      elif (code == "thr2"):
        thr = thr_step[3];
      else:
        thr = code;
      return float(thr);


def sumf(arr):
  s = 0
  for i in arr:
    s += i
  return s

def meanf(arr):
  if len(arr) == 0:
    return 0
  else:
    return sum(arr)/len(arr)

def variancef(arr):
  sq = 0.0
  m = meanf(arr)
  sumsq = 0
  for item in arr:
    sumsq += item**2
  return (sumsq/len(arr) - m*m)

def stdevf(arr):
  return variancef(arr) ** 0.5

def msef(arr):
  result = 0
  mean = meanf(arr)
  for item in arr:
    result += (item - mean) ** 2
  return result

def fitstep(arr):
  if len(arr) <= 0:
    return None
  #start = 0    # start and end are indices in arr  
  #end = count - 1
  sseArray = [0 for i in range(len(arr))] 
  sum = sumf(arr)
  mean = meanf(arr)
  sstot = msef(arr)
  count = len(arr)
  count1 = 0
  count2 = len(arr)
  sum1 = 0.0
  sum2 = sum
  sum1sq = 0.0
  sum2sq = sstot
  m1 = 0.0
  m2 = mean
  sse = sum1sq + sum2sq
  
  # loops through the array where index is an integer
  for index in range(count):
    entry = arr[index]
    # checks if element in array exists
    if entry is None:
      sseArray[index] = sse 

    count1 += 1
    count2 -= 1
    
    # checking if the division reaches the beginning so if the end counter reaches the beginning counter
    if count2 == 0:
      sseArray[index] = sstot
      continue;

    tmp = (mean - (entry + sum1)/count1)
    sum1sq = sum1sq + (entry - mean)**2 - tmp**2 * count1 + (count1 - 1) * (mean - m1)**2
    tmp = (mean - (sum2 - entry)/count2)
    sum2sq = sum2sq - (entry - mean)**2 - tmp**2 * count2 + (count2 + 1) * (mean - m2)**2
    sum1 += entry
    sum2 -= entry
    m1 = sum1/count1
    m2 = sum2/count2
    sse = sum1sq + sum2sq
    sseArray[index] = sse
  
  # find the minimum sumsq and its index
  bestSse = min(sseArray)
  bestIndex = sseArray.index(bestSse)

  # find mean of the first part and second part
  m1 = meanf(arr[:bestIndex+1])
  m2 = meanf(arr[bestIndex+1:])
  
  #threshold
  thr = (m1 + m2) /2

  # list reversed or not
  label = 0
  if m1 < m2:
    label = 1
  else:
    label = 2

  statistic = 0
  if bestSse > 0 :
    if count > 4:
      statistic = (sstot - bestSse)/3/(bestSse/(count - 4))
    else:
      statistic = (sstot - bestSse)/2/bestSse

  return {"cutoff": bestIndex+1, "bestSse": bestSse, "sstot": sstot, "statistic" : statistic, "threshold": thr, "label":label, "m1": m1, "m2": m2}

def getX(filename, x, debug):
  if not os.path.isfile(filename):
    print("Can't open file {0} <br>".format(filename));
    exit()
  fp = open(filename, "r")
  header = fp.readline().strip()
  fp.seek(int(x), 0)
  in_x = fp.readline()
  if (debug == 1):
    print("Line 1:<br/>",in_x,":<br>");
  fp.close()
  x_arr = in_x.split("\t")
  h_arr = header.split("\t")
  x_arr[-1] = x_arr[-1].strip()
  return (x_arr, h_arr);

def getHash(filename, index=None):
  if index is None:
      index = 0
  if not os.path.isfile(filename):
    print("Can't open file {0} <br>".format(filename));
    exit()
  res = {}
  fp = open(filename, "r")
  for line in fp:
      line = line.strip()
      ll = line.split("\t");
      res[ll[index]] = ll
  return res

def getStepMinerThr(data, start=None, end=None):
  if (start is None):
    start = 0;
  if end is None:
    end = len(data) - 1;
  array = []
  for i in range(start, len(data)):
    if (data[i] is None or data[i] is ""):
      continue;
    array.append(float(data[i]))
  array.sort()
  return fitstep(array);


hu = HegemonUtil

def reactome(idlist):
    import requests
    reactomeURI = 'http://www.reactome.org/AnalysisService/identifiers/projection?pageSize=100&page=1';
    response = requests.post(reactomeURI, data = idlist, \
                headers = { "Content-Type": "text/plain",  "dataType" : "json" })
    obj = json.loads(response.text)
    df = pd.DataFrame()
    df['name'] = [p["name"] for  p in obj["pathways"]]
    df['count'] = [p["entities"]["found"] for  p in obj["pathways"]]
    df['pValue'] = [p["entities"]["pValue"] for  p in obj["pathways"]]
    df['fdr'] = [p["entities"]["fdr"] for  p in obj["pathways"]]
    return df

def getPDF(cfile):
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(cfile)
    return pdf

def closePDF(pdf):
    import datetime
    d = pdf.infodict()
    d['Title'] = 'Plots'
    d['Author'] = 'Debashis Sahoo'
    d['Subject'] = "BoNE"
    d['Keywords'] = 'disease training validation ROC'
    d['CreationDate'] = datetime.datetime(2021, 10, 18)
    d['ModDate'] = datetime.datetime.today()
    pdf.close()

def asciiNorm(ah):
    if sys.version_info[0] >= 3:
        keys = list(ah.keys())
        for k in keys:
            ah[bytes(k, encoding='latin-1').decode('utf-8')] = ah[k]
    return ah

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

def getSimpleName(gene_groups):
    res = []
    for s in gene_groups:
        s1 = set()
        df = hu.getHegemonGeneIDs("G16", s)
        for name in df['Gene Name'].unique():
            if name != "" and name != "---":
                s1.add(name)
        res.append(s1)
    return res
def getGeneGroups(order = None, weight = None):
    data_item = []
    dir1 = basedir
    with open(dir1 + 'path-1.json') as data_file:
        data_item += json.load(data_file)
    with open(dir1 + 'path-0.json') as data_file:
        l1 = json.load(data_file)
        data_item[5] = l1[5]
        data_item[6] = l1[6]
    with open(dir1 + 'path-2.json') as data_file:
        data_item += json.load(data_file)
    with open(dir1 + 'path-3.json') as data_file:
        data_item += json.load(data_file)
    with open(dir1 + 'path-4.json') as data_file:
        data_item += json.load(data_file)
    cfile = dir1 + "mac-net-cls.txt"
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
                    gene_groups[i].add(k)
    print([len(s) for s in gene_groups])
    if order is None:
        order = [8, 9, 10]
    gene_groups = [gene_groups[i] for i in order]
    print([len(s) for s in gene_groups])
    gene_groups = getSimpleName(gene_groups)
    print([len(s) for s in gene_groups])
    if weight is None:
        weight = [-1, -2, -3]
    print(weight)
    genes = []
    return genes, weight, gene_groups

def getCls13a14a3():
    #order = [13, 14, 3]
    #wt1 = [-1, 1, 2]
    #nx = [0, 1, 4, 5, 6, 8, 9, 10, 16, 17, 19, 20, 21, 25, 28]
    #genes, wt1, l1 = getGeneGroups([nx[i] for i in order], wt1)
    wt1, l1 = [-1, 1, 2], [readGenes(basedir +  f"node-{i}.txt") for i in [13, 14, 3]]
    return wt1, l1

def getCls13():
    #order = [13]
    #wt1 = [-1]
    #nx = [0, 1, 4, 5, 6, 8, 9, 10, 16, 17, 19, 20, 21, 25, 28]
    #genes, wt1, l1 = getGeneGroups([nx[i] for i in order], wt1)
    wt1, l1 = [-1], [readGenes(basedir + f"node-{i}.txt") for i in [13]]
    return wt1, l1

def getCls14a3():
    #order = [14, 3]
    #wt1 = [1, 2]
    #nx = [0, 1, 4, 5, 6, 8, 9, 10, 16, 17, 19, 20, 21, 25, 28]
    #genes, wt1, l1 = getGeneGroups([nx[i] for i in order], wt1)
    wt1, l1 = [1, 2], [readGenes(basedir + f"node-{i}.txt") for i in [14, 3]]
    return wt1, l1

def saveNodes():
    nx = [0, 1, 4, 5, 6, 8, 9, 10, 16, 17, 19, 20, 21, 25, 28]
    order = range(1, len(nx))
    wt1 = [1]
    genes, wt1, l1 = getGeneGroups([nx[i] for i in order], wt1)
    for i in range(len(l1)):
        ofh = open(basedir + "node-" + str(i+1) + ".txt", "w")
        for g in l1[i]:
            ofh.write(g + "\n")
        ofh.close()

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
    from matplotlib.transforms import blended_transform_factory
    
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

def mergeRanks(group, start, exp, weight):
    X = np.array([[e[k-start] for e in exp] for k in group])
    arr = np.dot(X, np.array(weight))
    return arr

def getOrder(group, start, exp, weight):
    arr = mergeRanks(group, start, exp, weight)
    return [group[i] for i in np.argsort(arr)]

def getSName(name):
    l1 = re.split(": ", name)
    l2 = re.split(" /// ", l1[0])
    return l2[0]

def getRanksDf2(gene_groups, df_g, df_e, df_t):
    expr = []
    row_labels = []
    row_ids = []
    row_numhi = []
    ranks = []
    g_ind = 0
    counts = []
    noisemargins = []
    for k in range(len(df_e)):
        count = 0
        order = range(2, df_e[k].shape[1])
        avgrank = [0 for i in order]
        noisemargin = 0
        for j in df_g[k]['idx']:
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
            nv1 = 0.5/3
            if sd > 0:
                nv1 = nv1 / sd
            noisemargin += nv1 * nv1
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
        noisemargins.append(noisemargin)
        g_ind += 1
        counts += [count]
    print(counts)
    return ranks, noisemargins, row_labels, row_ids, row_numhi, expr

def getGroupsHs(gene_groups):
    cfile = basedir + "database/ensembl-GRCm38.p6-100-mm-hs.txt"
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

def getGroupsMm(gene_groups):
    cfile = basedir + "database/ensembl-GRCh38.p13-100-hs-mm.txt"
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

def getGroupsMmv1(gene_groups):
    cfile = basedir + "database/Homo_sapiens.GRCh38.95.chr_patch_hapl_scaff.len.txt"
    fp = open(cfile, "r")
    hsdict = {}
    for line in fp:
        line = line.strip();
        ll = re.split("\t", line);
        hsdict[ll[0]] = ll[1]
    fp.close();

    cfile = basedir + "database/mart-export-hs-mm.txt"
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

def plotHeatMapReactome():
    pdf = getPDF("heatmap-1.pdf")
    wt1, l1 = getCls13a14a3()
    ana = MacAnalysis()
    ana.getGEOMacAnn()
    ana.orderDataDf(l1, wt1)
    genes = readList(basedir + "selected-genes.txt")
    ofile = "heatmap-2.pdf"
    params = {'dx': 50, 'dy': 150, 'spaceAnn': len(ana.order)/len(ana.atypes),
              'tAnn': 0.5, 'widthAnn':0.5, 'lfs':8,
              'sy': 8, 'thr': 1, 'w': 3, 'h': 24, 'rowlabels':1,
              'genes': genes, 'atypes': ana.atypes,
              'tl': 6, 'tw': 0.25, 'ts': 10, 'tsi': -6000}
    ana.printHeatmap(ofile, genes, params)
    pdf.savefig(transparent=True,bbox_inches = 'tight')

    def plotReactome(df1):
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.25,4))
        ax = sns.barplot(x='lfdr', y = 'name', color=adj_light('orange', 1.5, 1) , data=df1)
        ax.set_xlabel('-log10(FDR)')
        ax.set_ylabel('')

    df = pd.read_csv(basedir + "reactome-13.txt", sep="\t")
    df['lfdr'] = -np.log10(df['fdr'])
    df1 = df.loc[ [2, 8, 9, 10, 14, 18], :]
    plotReactome(df1)
    pdf.savefig(transparent=True,bbox_inches = 'tight')

    df = pd.read_csv(basedir + "reactome-14.txt", sep="\t")
    df['lfdr'] = -np.log10(df['fdr'])
    df1 = df.loc[ [1, 6, 8, 9, 14, 18], :]
    plotReactome(df1)
    pdf.savefig(transparent=True,bbox_inches = 'tight')

    df = pd.read_csv(basedir + "reactome-3.txt", sep="\t")
    df['lfdr'] = -np.log10(df['fdr'])
    df1 = df.loc[ [0, 6, 8, 12, 14, 15], :]
    plotReactome(df1)
    pdf.savefig(transparent=True,bbox_inches = 'tight')

    l1 = readGenes(basedir + "node-14.txt")
    l2 = readGenes(basedir + "node-3.txt")
    df = reactome(" ".join(l1 + l2))
    df['lfdr'] = -np.log10(df['fdr'])
    df1 = df.loc[ [1, 6, 8, 9, 16, 18], :]
    plotReactome(df1)
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    closePDF(pdf)

def plotM1M2list(l1, wt1, ax, desc):
    ana = MacAnalysis()
    ana.getGEOMacAnn()
    ana.orderDataDf(l1, wt1)
    params = {'spaceAnn': len(ana.order)/len(ana.atypes), 
             'tAnn': 1, 'widthAnn':1, 'acolor': acolor,
             'w': 5, 'h': 0.8, 'ax': ax}
    ax = ana.printTitleBar(params)
    res = ana.getMacMetrics(ana.cval[0])
    print(desc, res)
    ax.text(len(ana.cval[0]), 4, ",".join(res))
    ax.text(-1, 2, desc, horizontalalignment='right',
            verticalalignment='center')

def plotM1M2(order, wt1, ax, desc):
    nx = [0, 1, 4, 5, 6, 8, 9, 10, 16, 17, 19, 20, 21, 25, 28]
    genes, wt1, l1 = getGeneGroups([nx[i] for i in order], wt1)
    plotM1M2list(l1, wt1, ax, desc)

def computeROCAUC(order, wt1):
    nx = [0, 1, 4, 5, 6, 8, 9, 10, 16, 17, 19, 20, 21, 25, 28]
    genes, wt1, l1 = getGeneGroups([nx[i] for i in order], wt1)
    ana = MacAnalysis()
    ana.getGEOMacAnn()
    ana.orderDataDf(l1, wt1)
    res = ana.getMacMetrics(ana.cval[0])
    return res

def learningAlgorithm():
    cfile = basedir + 'graph.txt'
    fp = open(cfile, "r")
    nodes = {}
    edges = {}
    for line in fp:
        line = line.strip();
        ll = line.split(" ");
        ll[0] = int(ll[0])
        ll[2] = int(ll[2])
        nodes[ll[0]] = 1
        nodes[ll[2]] = 1
        if ll[0] not in edges:
            edges[ll[0]] = {}
        if ll[1] not in edges[ll[0]]:
            edges[ll[0]][ll[1]] = {}
        edges[ll[0]][ll[1]][ll[2]] = 1
    fp.close();
    res = []
    for k in edges:
        if "hilo" in edges[k]:
            for l in edges[k]['hilo']:
                res += [[k, l]]
    res2 = []
    for k in res:
        if k[1] in edges and "lolo" in edges[k[1]]:
            for l in edges[k[1]]['lolo']:
                res2 += [ [[-1, 1, 2], [k[0], k[1], l]] ]
        if k[0] in edges and "lolo" in edges[k[0]]:
            for l in edges[k[0]]['lolo']:
                res2 += [ [[-2, -1, 1], [l, k[0], k[1]]] ]
    for i in range(len(res2)):
        k = res2[i]
        res2[i] += [computeROCAUC(k[1], k[0])]
        print(res2[i])
    res3 = max(res2, key=lambda k: float(k[2][1]) * float(k[2][2]))
    return res3, res2


def trainingAlgorithm():
    fig = plt.figure(figsize=(5,7), dpi=100)
    n1 = 7
    axlist = []
    for i in range(n1):
        ax = plt.subplot2grid((n1, 1), (i, 0))
        axlist.extend([ax])
    plotM1M2([1, 2, 3], [-1, 1, 2], axlist[0], "Path 1-2-3")
    plotM1M2([1, 14, 3], [-1, 1, 2], axlist[1], "Path 1-14-3")
    plotM1M2([12, 5, 6], [-1, 1, 2], axlist[2], "Path 12-5-6")
    plotM1M2([11, 10, 8], [-2, -1, 1], axlist[3], "Path 11-10-8")
    plotM1M2([13, 14, 3], [-1, 1, 2], axlist[4], "Path 13-14-3")
    plotM1M2([14, 3], [1, 2], axlist[5], "Path 14-3")
    plotM1M2([13], [-1], axlist[6], "Cluster 13")
    fig.savefig("training.pdf", dpi=100)

def getCoates():
    #PMID: 18199539 DOI: 10.1158/0008-5472.CAN-07-3050
    cfile = basedir + "database/c2.all.v6.2.symbols.gmt"
    fp = open(cfile, "r")
    l1 = l2 = None
    for line in fp:
        line = line.strip();
        ll = line.split("\t");
        if ll[0] == "COATES_MACROPHAGE_M1_VS_M2_UP":
           l1 = ll[2:]
        if ll[0] == "COATES_MACROPHAGE_M1_VS_M2_DN":
           l2 = ll[2:]
    fp.close();
    res = [set(l1), set(l2)]
    return res

def getMartinez():
    #PMID: 17082649 DOI: 10.4049/jimmunol.177.10.7303
    cfile = basedir + "database/c7.all.v6.2.symbols.gmt"
    fp = open(cfile, "r")
    l1 = l2 = None
    for line in fp:
        line = line.strip();
        ll = line.split("\t");
        if ll[0] == "GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP":
           l2 = ll[2:]
        if ll[0] == "GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN":
           l1 = ll[2:]
    fp.close();
    res = [set(l1), set(l2)]
    return res

def getMartinezM2hsmm():
    cfile = basedir + "database/Martinez-Table3.txt"
    res = getEntries(cfile, 0)
    return [res[1:]]

def getMartinezM0hsmm():
    cfile = basedir + "database/Martinez-TableS1.txt"
    res = getEntries(cfile, 0)
    return [res[1:]]

def getBell2016():
# Citation: Bell LCK, Pollara G, Pascoe M, Tomlinson GS, Lehloenya RJ, Roe J, et
# al. (2016) In Vivo Molecular Dissection of the Effects of HIV-1 in Active
# Tuberculosis. PLoS Pathog 12(3): e1005469.
# https://doi.org/10.1371/journal.ppat.1005469
#PMID: 26986567 PMCID: PMC4795555 DOI: 10.1371/journal.ppat.1005469
    cfile = basedir + "database/journal.ppat.1005469.s018.txt"
    fp = open(cfile, "r")
    l1 = l2 = l3 = None
    for line in fp:
        line = line.strip();
        ll = line.split("\t");
        if ll[0] == "LPS":
           l1 = ll[2].replace('"', "").split(", ")
        if ll[0] == "IFNg":
           l2 = ll[2].replace('"', "").split(", ")
        if ll[0] == "IL-4 and IL-13":
           l3 = ll[2].replace('"', "").split(", ")
    fp.close();
    res = [set(l1), set(l2), set(l3)]
    return res

def getBecker():
    #Citation:
    #Becker M, De Bastiani MA, Parisi MM, Guma FT, Markoski MM, Castro
    #MA, Kaplan MH, Barbé-Tuana FM, Klamt F. Integrated Transcriptomics 
    #Establish Macrophage Polarization Signatures and have Potential 
    #Applications for Clinical Health and Disease. Sci Rep. 2015 Aug 25;
    #5:13351. doi: 10.1038/srep13351.
    #PMID: 26302899 PMCID: PMC4548187 DOI: 10.1038/srep13351
    cfile = basedir + "database/srep13351-s1-1.txt"
    l1 = getEntries(cfile, 0)
    cfile = basedir + "database/srep13351-s1-2.txt"
    l2 = getEntries(cfile, 0)
    return [l1[1:], l2[2:]]

def getSViP():
    l1 = [readList(basedir + "database/svip-signature.txt")] # 20 gene signature
    wt1 = [1]
    return wt1, l1

def getViP():
    l1 = [readList(basedir + "database/vip-signature.txt")] # 166 gene signature
    wt1 = [1]
    return wt1, l1

def comparativeAnalysis():
    fig = plt.figure(figsize=(5,3), dpi=100)
    n1 = 5
    axlist = []
    for i in range(n1):
        ax = plt.subplot2grid((n1, 1), (i, 0))
        axlist.extend([ax])

    order = [13, 14, 3]
    wt1 = [-1, 1, 2]
    nx = [0, 1, 4, 5, 6, 8, 9, 10, 16, 17, 19, 20, 21, 25, 28]
    genes, wt1, l1 = getGeneGroups([nx[i] for i in order], wt1)
    plotM1M2list(l1, wt1, axlist[0], "BoNE")

    l1 = getBecker()
    wt1 = [-1, 1]
    plotM1M2list(l1, wt1, axlist[1], "PMID: 26302899")

    l1 = getBell2016()
    wt1 = [-1, -1, 1]
    plotM1M2list(l1, wt1, axlist[2], "PMID: 26986567")

    l1 = getCoates()
    wt1 = [-1, 1]
    plotM1M2list(l1, wt1, axlist[3], "PMID: 18199539")

    l1 = getMartinez()
    wt1 = [-1, 1]
    plotM1M2list(l1, wt1, axlist[4], "PMID: 17082649")

    fig.savefig("comp-1.pdf", dpi=100)

def T2(order, wt1, ofile):
    pathdesc = "-".join([str(i) for i in order])
    acolor = ["#00CC00", "#D8A03D","#EC008C",
	      'cyan', "#B741DC", "#808285",
	      'blue', 'black', 'green', 'red']
    fig = plt.figure(figsize=(5,7), dpi=100)
    n1 = 6
    axlist = []
    for i in range(n1):
        ax = plt.subplot2grid((n1, 1), (i, 0))
        axlist.extend([ax])
    nx = [0, 1, 4, 5, 6, 8, 9, 10, 16, 17, 19, 20, 21, 25, 28]
    genes, wt1, l1 = getGeneGroups([nx[i] for i in order], wt1)
    ana = MacAnalysis()
    ana.getGEOMacAnn()
    desc = "Pooled GSMs"
    ana.orderDataDf(l1, wt1)
    params = {'spaceAnn': len(ana.order)/len(ana.atypes), 
             'tAnn': 1, 'widthAnn':1, 'acolor': acolor,
             'w': 5, 'h': 0.8, 'ax': axlist[0]}
    ax = ana.printTitleBar(params)
    res = ana.getMacMetrics(ana.cval[0])
    print(pathdesc, desc, res)
    ax.text(len(ana.cval[0]), 4, ",".join(res))
    ax.text(-1, 2, desc, horizontalalignment='right',
            verticalalignment='center')

    ana.getBeyer2012()
    desc = "GSE35449"
    ana.orderDataDf(l1, wt1)
    params = {'spaceAnn': len(ana.order)/len(ana.atypes), 
             'tAnn': 1, 'widthAnn':1, 'acolor': acolor,
             'w': 5, 'h': 0.8, 'ax': axlist[1]}
    ax = ana.printTitleBar(params)
    res = ana.getMacMetrics(ana.cval[0])
    print(pathdesc, desc, res)
    ax.text(len(ana.cval[0]), 4, ",".join(res))
    ax.text(-1, 2, desc, horizontalalignment='right',
            verticalalignment='center')

    ana.getXue2014(4)
    desc = "GSE46903"
    ana.orderDataDf(l1, wt1)
    params = {'spaceAnn': len(ana.order)/len(ana.atypes), 
             'tAnn': 1, 'widthAnn':1, 'acolor': acolor,
             'w': 5, 'h': 0.8, 'ax': axlist[2]}
    ax = ana.printTitleBar(params)
    res = ana.getMacMetrics(ana.cval[0])
    print(pathdesc, desc, res)
    ax.text(len(ana.cval[0]), 4, ",".join(res))
    ax.text(-1, 2, desc, horizontalalignment='right',
            verticalalignment='center')

    ana.getOhradanovaRepic2018(2)
    desc = "GSE61298"
    ana.orderDataDf(l1, wt1)
    params = {'spaceAnn': len(ana.order)/len(ana.atypes), 
             'tAnn': 1, 'widthAnn':1, 'acolor': acolor,
             'w': 5, 'h': 0.8, 'ax': axlist[3]}
    ax = ana.printTitleBar(params)
    res = ana.getMacMetrics(ana.cval[0])
    print(pathdesc, desc, res)
    ax.text(len(ana.cval[0]), 4, ",".join(res))
    ax.text(-1, 2, desc, horizontalalignment='right',
            verticalalignment='center')

    ana.getZhang2015(2)
    desc = "GSE55536"
    ana.orderDataDf(l1, wt1)
    params = {'spaceAnn': len(ana.order)/len(ana.atypes), 
             'tAnn': 1, 'widthAnn':1, 'acolor': acolor,
             'w': 5, 'h': 0.8, 'ax': axlist[4]}
    ax = ana.printTitleBar(params)
    res = ana.getMacMetrics(ana.cval[0])
    print(pathdesc, desc, res)
    ax.text(len(ana.cval[0]), 4, ",".join(res))
    ax.text(-1, 2, desc, horizontalalignment='right',
            verticalalignment='center')

    ana.getZhang2015(3)
    desc = "GSE55536"
    ana.orderDataDf(l1, wt1)
    params = {'spaceAnn': len(ana.order)/len(ana.atypes), 
             'tAnn': 1, 'widthAnn':1, 'acolor': acolor,
             'w': 5, 'h': 0.8, 'ax': axlist[5]}
    ax = ana.printTitleBar(params)
    res = ana.getMacMetrics(ana.cval[0])
    print(pathdesc, desc, res)
    ax.text(len(ana.cval[0]), 4, ",".join(res))
    ax.text(-1, 2, desc, horizontalalignment='right',
            verticalalignment='center')

    fig.savefig(ofile, dpi=100)

def validationAlgorithm():
    T2([13, 14, 3], [-1, 1, 2], "validation-1.pdf")
    T2([13], [-1], "validation-2.pdf")
    T2([14, 3], [1, 2], "validation-3.pdf")
    return

def MacData(adata, desc, Org, tn=1):
    ana, wt1, l1, pathdesc = adata
    if Org == "Mm":
        if tn == 5 or tn == 6 or tn == 7:
            l1 = getGroupsMmv1(l1)
        else:
            l1 = getGroupsMm(l1)
    ana.orderDataDf(l1, wt1)
    params = {'spaceAnn': len(ana.order)/len(ana.atypes), 
             'tAnn': 1, 'widthAnn':1, 'acolor': acolor,
             'w': 5, 'h': 0.8}
    ax = ana.printTitleBar(params)
    if tn == 1 or tn == 6:
        res = ana.getROCAUCspecific(1, 0)
        print(pathdesc, desc, res)
        ax.text(len(ana.cval[0]), 4, res)
    if tn == 2 or tn == 7:
        res = ana.getMacMetrics(ana.cval[0])
        print(pathdesc, desc, " ".join(res[1:]))
        ax.text(len(ana.cval[0]), 4, " ".join(res[1:]))
    if tn == 3 or tn == 5:
        res = ana.getROCAUCspecific(0, 2)
        print(pathdesc, desc, res)
        ax.text(len(ana.cval[0]), 4, res)
    if tn == 4:
        res = ana.getROCAUCspecific(1, 2)
        print(pathdesc, desc, res)
        ax.text(len(ana.cval[0]), 4, res)
    ax.text(-1, 2, desc, horizontalalignment='right',
            verticalalignment='center')
    return

def figure2a():
    wt1, l1 = [-1, 1, 2], [readGenes(basedir + f"node-{i}.txt") for i in [13, 14, 3]]
    pathdesc = "13-14-3"
    ta = 1
    ana = MacAnalysis()
    adata = [ana, wt1, l1, pathdesc]
    res = []
    ana.getRock2005()
    res += [MacData(adata, desc="Microglia (Brain)", Org='Hs', tn=ta)]
    ana.getPolak2014(2)
    res += [MacData(adata, desc="Langerhan’s (Skin)", Org='Hs', tn=ta)]
    ana.getQu2016(2)
    res += [MacData(adata, desc="Colon, Sigmoid", Org='Hs', tn=ta)]
    ana.getImmGenULI(1) # Broncho-alveolar lavage Macrophage in 2ug LPS
    res += [MacData(adata, desc="Lung Lavage", Org='Mm', tn=ta)]
    ana.getImmGenULI(2) # Lung Macrophage in 2ug LPS 3 days
    res += [MacData(adata, desc="Lung MHCII+", Org='Mm', tn=ta)]
    ana.getImmGenULI(3) # Lung Macrophage in 2ug LPS 3 days
    res += [MacData(adata, desc="Lung MHCII-", Org='Mm', tn=ta)]
    ana.getImmGenULI(4) # Female Peritoneal macrophages-10kIFN treatment
    res += [MacData(adata, desc="Female Peritoneal", Org='Mm', tn=6)]
    ana.getImmGenULI(5) # Male/Female Peritoneal macrophages-10kIFN treatment
    res += [MacData(adata, desc="Male Peritoneal", Org='Mm', tn=ta)]
    return

def figure2b():
    wt1, l1 = [-1, 1, 2], [readGenes(basedir + f"node-{i}.txt") for i in [13, 14, 3]]
    pathdesc = "13-14-3"
    ta = 2
    ana = MacAnalysis()
    adata = [ana, wt1, l1, pathdesc]
    res = []
    ana.getSvensson2011()
    res += [MacData(adata, desc="Blood,Decidua,Preg", Org='Hs', tn=ta)]
    ana.getChandriani2014(2)
    res += [MacData(adata, desc="Blood,Monocytes", Org='Hs', tn=ta)]
    ana.getMartinez2015()
    res += [MacData(adata, desc="PBMCs Macs", Org='Hs', tn=ta)]
    ana.getGuler2015(2) # mouse bone marrow derive macrophage
    res += [MacData(adata, desc="BMDMs", Org='Mm', tn=ta)]
    ana.getDas2018()
    res += [MacData(adata, desc="BMDMs", Org='Mm', tn=7)]
    ana.getPiccolo2017()
    res += [MacData(adata, desc="BMDMs", Org='Mm', tn=ta)]
    ana.getOstuni2013()
    res += [MacData(adata, desc="BMDMs", Org='Mm', tn=ta)]
    ana.getRochaResende()
    res += [MacData(adata, desc="Peritoneal mac", Org='Mm', tn=7)]
    ana.getHill2018(2)
    res += [MacData(adata, desc="BMDMs", Org='Mm', tn=ta)]
    ana.getDaniel2018()
    res += [MacData(adata, desc="BMDMs", Org='Mm', tn=7)]
    ana.getFreemerman2019(2)
    res += [MacData(adata, desc="BMDMs", Org='Mm', tn=ta)]
    ana.getLi2015()
    res += [MacData(adata, desc="BMDMs", Org='Mm', tn=ta)]
    ana.getHan2017(2)
    res += [MacData(adata, desc="Kupffer cells (Liver)", Org='Mm', tn=ta)]
    return

def figure2c():
    wt1, l1 = [-1, 1, 2], [readGenes(basedir + f"node-{i}.txt") for i in [13, 14, 3]]
    wt2, l2 = [-1], [readGenes(basedir + f"node-{i}.txt") for i in [13]]
    wt3, l3 = [1, 2], [readGenes(basedir + f"node-{i}.txt") for i in [14, 3]]
    pathdesc = "13-14-3"
    ta = 1
    ana = MacAnalysis()
    adata = [ana, wt1, l1, "13-14-3"]
    adata2 = [ana, wt2, l2, "13"]
    adata3 = [ana, wt3, l3, "14-3"]
    res = []
    ana.getJeffrey2006(2)
    desc = "Mac, NK, DC, Eos, Bas"
    res += [MacData(adata2, desc=desc, Org='Hs', tn=ta)]
    ana.getSippel2018(2)
    desc = "Eosinophil"
    res += [MacData(adata2, desc=desc, Org='Hs', tn=ta)]
    ana.getPuan2017(2)
    desc = "Basophil"
    ana.orderDataDf(l1, wt1)
    for i in range(3):
        ana.f_ranks[ana.pair[i][1] - ana.start] -= ana.f_ranks[ana.pair[i][0] - ana.start]
        ana.f_ranks[ana.pair[i][0] - ana.start] = 0
    arr = [ana.f_ranks[i - ana.start] for i in ana.order]
    ana.i1 = [ana.order[i] for i in np.argsort(arr)]
    ana.cval = np.array([[ana.aval[i] for i in ana.i1]])
    params = {'spaceAnn': len(ana.order)/len(ana.atypes),
            'tAnn': 1, 'widthAnn':1, 'w': 5, 'h':0.75,
            'acolor': acolor}
    ax = ana.printTitleBar(params)
    res1 = ana.getROCAUCspecific(1, 0)
    print(adata[3], desc, res1)
    ax.text(len(ana.cval[0]), 4, res1)
    ax.text(-1, 2, desc, horizontalalignment='right',
                verticalalignment='center')

    ana.getKlocperk2020(3)
    desc = "Neutrophil"
    res += [MacData(adata, desc=desc, Org='Hs', tn=ta)]
    ana.getLenardo2020()
    desc = "T Cell"
    res += [MacData(adata3, desc=desc, Org='Hs', tn=3)]
    ana.getNair2015()
    desc = "NKT"
    res += [MacData(adata, desc=desc, Org='Hs', tn=ta)]
    ana.getAbbas2005()
    desc = "B,T,Pl,NK,Mo,Neu"
    res += [MacData(adata, desc=desc, Org='Hs', tn=ta)]
    ana.getMetcalf2015()
    desc = "PBMC"
    res += [MacData(adata2, desc=desc, Org='Hs', tn=ta)]
    ana.getBanchereau2014I(2)
    desc = "DC"
    res += [MacData(adata, desc=desc, Org='Hs', tn=4)]
    ana.getBanchereau2014I(3)
    desc = "DC"
    res += [MacData(adata, desc=desc, Org='Hs', tn=ta)]
    ana.getIampietro2017(2)
    desc = "T Cell"
    res += [MacData(adata, desc=desc, Org='Hs', tn=ta)]
    #ana.getJohnson2020()
    #desc = "DC"
    #res += [MacData(adata2, desc=desc, Org='Mm', tn=ta)]

    ana.getImmGenULI(7) # Female Splenic NKT Cell
    res += [MacData(adata, desc="Splenic NKT", Org='Mm', tn=ta)]
    ana.getImmGenULI(8) # Splenic CD4 T cell-10kIFN treatment
    res += [MacData(adata, desc="Splenic CD4 T", Org='Mm', tn=ta)]
    ana.getImmGenULI(9) # Female Splenic B cells-1kIFN treatment
    res += [MacData(adata, desc="Splenic B Cell", Org='Mm', tn=ta)]
    ana.getImmGenULI(10) # Female Splenic Neutrophils-10kIFN treatment
    res += [MacData(adata, desc="Splenic Neutrophils", Org='Mm', tn=ta)]

    ta = 2
    ana.getHutchins2015TPM(2) # Bone marrow neutrophil
    res += [MacData(adata, desc="BM Neutrophils", Org='Mm', tn=ta)]
    ana.getHutchins2015TPM(3) # spleen-purified dendritic cells
    res += [MacData(adata, desc="Spleen DC", Org='Mm', tn=ta)]
    ana.getHutchins2015TPM(4) # Peritoneal exudate cells
    res += [MacData(adata, desc="Peritoneal Macrophages", Org='Mm', tn=7)]
    ana.getHutchins2015TPM(6) # Bone marrow-derived mast cell
    res += [MacData(adata, desc="BM mast cell", Org='Mm', tn=ta)]

    return

def getMacSignatures(index=0):
    wt1, l1 = [1], [['TYROBP', 'FCER1G']]
    if (index == 0 or index == 'C13'):
        wt1, l1 = getCls13()
    if (index == 1 or index == 'C14-3'):
        wt1, l1 = getCls14a3()
    if (index == 2 or index == 'C13-14-3'):
        wt1, l1 = getCls13a14a3()
    if (index == 5 or index == 'SAM A'):
        cfile = basedir + 'database/SAM_DAM_LAM signatures.xlsx'
        df = pd.read_excel(cfile, sheet_name="SAMac")
        wt1, l1 = [1], [list(df['Signature A'].dropna())]
    if (index == 6 or index == 'SAM B'):
        cfile = basedir + 'database/SAM_DAM_LAM signatures.xlsx'
        df = pd.read_excel(cfile, sheet_name="SAMac")
        wt1, l1 = [1], [list(df['Signature B'].dropna())]
    if (index == 7 or index == 'SAM i'):
        cfile = basedir + 'database/SAM_DAM_LAM signatures.xlsx'
        df = pd.read_excel(cfile, sheet_name="SAMac_inflam", skiprows=4)
        wt1, l1 = [1], [list(df['Inflammatory macs'].dropna())]
    if (index == 8 or index == 'SAM ni'):
        cfile = basedir + 'database/SAM_DAM_LAM signatures.xlsx'
        df = pd.read_excel(cfile, sheet_name="SAMac_NonInf", skiprows=3)
        wt1, l1 = [1], [list(df['Non-Inflammatory Macs'].dropna())]
    if (index == 9 or index == 'DAM'):
        cfile = basedir + 'database/SAM_DAM_LAM signatures.xlsx'
        df = pd.read_excel(cfile, sheet_name="DAM signature")
        up = list(df[df['up/down'] == 1]['Unnamed: 0'])
        down = list(df[df['up/down'] == -1]['Unnamed: 0'])
        wt1, l1 = [1, -1], [up, down]
        l1 = getGroupsHs(l1)
    if (index == 10 or index == 'LAM'):
        cfile = basedir + 'database/SAM_DAM_LAM signatures.xlsx'
        df = pd.read_excel(cfile, sheet_name="LAM signature")
        up = list(df[df['fcLAMvsMacHuman'] >= 3]['gene'])
        down = list(df[df['fcLAMvsMacHuman'] <= -2]['gene'])
        wt1, l1 = [1, 0], [up, down]
    if (index == 11 or index == 'ViP'):
        wt1, l1 = getViP()
        wt1 = [-1]
    if (index == 12 or index == 'sViP'):
        wt1, l1 = getSViP()
        wt1 = [-1]
    if (index == 13 or index == 'Becker' or index == "PMID: 26302899"):
        l1 = getBecker()
        wt1 = [-1, 1]
    if (index == 14 or index == 'Bell2016' or index == "PMID: 26986567"):
        l1 = getBell2016()
        wt1 = [-1, -1, 1]
    if (index == 15 or index == 'Coates' or index == "PMID: 18199539"):
        l1 = getCoates()
        wt1 = [-1, 1]
    if (index == 16 or index == 'Martinez' or index == "PMID: 17082649"):
        l1 = getMartinez()
        wt1 = [-1, 1]
    if (index == 17 or index == 'LM22-M1'):
        df = pd.read_csv(basedir + "database/LM22-signature.csv")
        wt1, l1 = [-1], [list(df['Macrophages_M1'].dropna())]
    if (index == 18 or index == 'LM22-M2'):
        df = pd.read_csv(basedir + "database/LM22-signature.csv")
        wt1, l1 = [1], [list(df['Macrophages_M2'].dropna())]
    if (index == 18 or index == 'LM22-M0'):
        df = pd.read_csv(basedir + "database/LM22-signature.csv")
        wt1, l1 = [1], [list(df['Macrophages_M0'].dropna())]
    if (index == 19 or index == 'Murray-M1'):        
        df = pd.read_csv(basedir + "database/Murray-PMID27813830.txt", sep="\t")
        wt1, l1 = [-1], [list(df['M1'].dropna())]
        l1 = getGroupsHs(l1)
    if (index == 20 or index == 'Murray-M2'):
        df = pd.read_csv(basedir + "database/Murray-PMID27813830.txt", sep="\t")
        wt1, l1 = [1], [list(df['M2'].dropna())]
        l1 = getGroupsHs(l1)
    if (index == 21 or index == 'DAM 1'):
        wt1, l1 = [1, -1],[['Tyrobp', 'Ctsb', 'Cstd', 'Apoe', 'B2m', 'Fth1', 'Lyz2'],
                           ['Cx3cr1',  'P2ry12', 'Tmem119']]
    if (index == 22 or index == 'DAM 2'):
        wt1, l1 = [1],[['Trem2', 'Axl', 'Cst7', 'Cst1', 'Lpl', 'Cd9', 'Csf1',
            'Ccl6', 'Itgax', 'Clec7', 'Lilrb4', 'Timp4']]
    return wt1, l1

def figure2e():
    from scipy.stats import hypergeom
    list1 = ['SAM A', 'SAM B', 'SAM i', 'SAM ni', 'DAM', 'LAM', 'ViP']
    res = []
    for j in range(1, 15):
        k = readGenes(basedir + f"node-{j}.txt")
        v = []
        for i in list1:
            wt1, l1 = getMacSignatures(i)
            v += [1 - hypergeom.cdf(len(set(k).intersection(l1[0])), 13037, len(k), len(l1[0]))]
        res += [v]

    df = pd.DataFrame(res)
    df.index = [i+1 for i in df.index]
    df.columns = list1
    #df = df[list1[5:]]
    #sns.heatmap(df, cmap="YlGnBu")
    #sorted([i for k in res for i in k])
    df[df == 0] = 1e-13
    df = -np.log(df)/np.log(10)
    df[df > 4] = 4
    df[ df < -np.log(0.05)/np.log(10) ] = 0
    sns.heatmap(df, cmap="Blues")
    return

def MacModel(ana, desc='', Org='Hs', id1=None):
    def getL(l1):
        return '(' + ",".join([str(len(k)) for k in l1]) +')'

    list1 = ['C13', 'C14-3', 'C13-14-3']
    res = []

    for k in list1:
        wt1, l1 = getMacSignatures(k)
        ann = [re.sub(" .*", "", ana.source),Org,k, getL(l1)]
        res += [ana.getStats(l1, wt1, ann, id1)]
    
    cols = ['GSEID', 'ROC-AUC', 'pvalue', '#Cont', '#Expt',
            'Series', 'Species', 'Signature', '#Genes']
    df = pd.DataFrame(res, columns=cols)
    df['Condition'] = desc
    return df

def MacDiseaseModel():
    ana = MacAnalysis()
    res = []
    ana.getPeters(3)
    res += [MacModel(ana, desc="UC " + ana.source, Org='Hs')]
    ana.getPeters(4)
    res += [MacModel(ana, desc="CD " + ana.source, Org='Hs')]
    ana.getQu2016(2)
    res += [MacModel(ana, desc="Normal Surface " + ana.source, Org='Hs')]
    ana.getQu2016(3)
    res += [MacModel(ana, desc="Adenoma " + ana.source, Org='Hs')]
    ana.getQu2016(4)
    res += [MacModel(ana, desc="CRC " + ana.source, Org='Hs')]
    ana.getGkouskou2016(2)
    res += [MacModel(ana, desc="Distal/Proximal " + ana.source, Org='Mm')]
    ana.getGkouskou2016ProxDis()
    res += [MacModel(ana, desc="Distal/Proximal " + ana.source, Org='Mm')]
    ana.getWoetzel2014(2)
    res += [MacModel(ana, desc="OA " + ana.source, Org='Hs')]
    ana.getWoetzel2014(3)
    res += [MacModel(ana, desc="RA " + ana.source, Org='Hs')]
    ana.getLefebvre2017(3)
    res += [MacModel(ana, desc="NASH " + ana.source, Org='Hs')]
    ana.getSuppli2019(3)
    res += [MacModel(ana, desc="NASH " + ana.source, Org='Hs')]
    ana.getSuppli2019(2)
    res += [MacModel(ana, desc="NAFLD " + ana.source, Org='Hs')]
    ana.getSuppli2019(4)
    res += [MacModel(ana, desc="Obese " + ana.source, Org='Hs')]
    ana.getWoodruff2005(3)
    res += [MacModel(ana, desc="Smoker " + ana.source, Org='Hs')]
    ana.getWS2009(4)
    res += [MacModel(ana, desc="Smoker " + ana.source, Org='Hs')]
    ana.getWS2009(2)
    res += [MacModel(ana, desc="Asthma " + ana.source, Org='Hs')]
    ana.getWS2009(3)
    res += [MacModel(ana, desc="COPD " + ana.source, Org='Hs')]
    ana.getLissner(4)
    res += [MacModel(ana, desc="Newborn " + ana.source, Org='Hs')]
    ana.getLissner(5)
    res += [MacModel(ana, desc="Old " + ana.source, Org='Hs')]
    ana.getBadawi2019(2)
    res += [MacModel(ana, desc="ICM " + ana.source, Org='Hs')]
    ana.getBondar2017(2)
    res += [MacModel(ana, desc="ICM " + ana.source, Org='Hs')]
    ana.getPatel2019(6)
    res += [MacModel(ana, desc="AD " + ana.source, Org='Hs')]
    ana.getBerchtold2014RMA(6)
    res += [MacModel(ana, desc="AD " + ana.source, Org='Hs')]
    ana.getGelman2012()
    res += [MacModel(ana, desc="HAND " + ana.source, Org='Hs')]
    ana.getChenPlotkin2008(2)
    res += [MacModel(ana, desc="FTD " + ana.source, Org='Hs')]
    ana.getOlmosSerrano2016()
    res += [MacModel(ana, desc="DS " + ana.source, Org='Hs')]
    ana.getBartolettiStella2019()
    res += [MacModel(ana, desc="CJD " + ana.source, Org='Hs')]
    ana.getTsalik2015(2)
    res += [MacModel(ana, desc="Sepsis " + ana.source, Org='Hs')]
    ana.getBarcella2018(2)
    res += [MacModel(ana, desc="Sepsis " + ana.source, Org='Hs')]
    ana.getStenvers2019(6)
    res += [MacModel(ana, desc="T2D " + ana.source, Org='Hs')]
    ana.getWu2007IS(3)
    res += [MacModel(ana, desc="Diabetes " + ana.source, Org='Hs')]
    ana.getDuPlessis2015(2)
    res += [MacModel(ana, desc="MetS " + ana.source, Org='Hs')]
    ana.getWatson2017()
    res += [MacModel(ana, desc="Sleep " + ana.source, Org='Hs', id1="DBP")]
    ana.getUyhelji2018()
    res += [MacModel(ana, desc="Sleep " + ana.source, Org='Hs', id1="ARNTL")]
    ana.getMaret2007(3)
    res += [MacModel(ana, desc="Sleep " + ana.source, Org='Hs')]
    ana.getResuehr2019()
    res += [MacModel(ana, desc="Sleep " + ana.source, Org='Hs')]
    ana.getDAmore2018()
    res += [MacModel(ana, desc="Sleep " + ana.source, Org='Hs', id1="DBP")]
    ana.getChristou2019(3)
    res += [MacModel(ana, desc="Night " + ana.source, Org='Hs')]
    ana.getWagstaffe2020(2)
    res += [MacModel(ana, desc="Ebola " + ana.source, Org='Hs')]
    ana.getPrice2020(2)
    res += [MacModel(ana, desc="Ebola " + ana.source, Org='Hs')]
    ana.getReynard2019(3)
    res += [MacModel(ana, desc="Ebola " + ana.source, Org='Hs')]
    ana.getCameron2007(2)
    res += [MacModel(ana, desc="SARS " + ana.source, Org='Hs')]
    ana.getMitchell2013(2)
    res += [MacModel(ana, desc="SARS " + ana.source, Org='Hs')]
    ana.getGuan2018()
    res += [MacModel(ana, desc="H7N9 " + ana.source, Org='Hs')]
    ana.getZhai2015(2)
    res += [MacModel(ana, desc="H1N1 " + ana.source, Org='Hs')]
    ana.getJosset2014(2)
    res += [MacModel(ana, desc="InfA " + ana.source, Org='Hs')]
    ana.getJosset2014(3)
    res += [MacModel(ana, desc="Cov1 " + ana.source, Org='Hs')]
    ana.getJones2019(3)
    res += [MacModel(ana, desc="Flu " + ana.source, Org='Hs')]
    ana.getDunning2018(2)
    res += [MacModel(ana, desc="Flu " + ana.source, Org='Hs')]

    df = pd.concat(res, sort=True)
    return df

def MacComparison(ana, desc='', Org='Hs'):
    def getL(l1):
        return '(' + ",".join([str(len(k)) for k in l1]) +')'

    list1 = ['C13', 'C14-3', 'C13-14-3',
            'Becker', 'Bell2016', 'Coates', 'Martinez']
    res = []

    for k in list1:
        wt1, l1 = getMacSignatures(k)
        ann = [re.sub(" .*", "", ana.source),Org,k, getL(l1)]
        res += [ana.getStats(l1, wt1, ann)]
    
    cols = ['GSEID', 'ROC-AUC', 'pvalue', '#Cont', '#Expt',
            'Series', 'Species', 'Signature', '#Genes']
    df = pd.DataFrame(res, columns=cols)
    df['Condition'] = desc
    return df

def MacComparisonRes():
    ana = MacAnalysis()
    res = []
    ana.getGEOMacAnn()
    res += [MacComparison(ana, desc="" + ana.source, Org='Hs')]
    ana.getBeyer2012()
    res += [MacComparison(ana, desc="" + ana.source, Org='Hs')]
    ana.getXue2014(4)
    res += [MacComparison(ana, desc="" + ana.source, Org='Hs')]
    ana.getOhradanovaRepic2018(2)
    res += [MacComparison(ana, desc="" + ana.source, Org='Hs')]
    df = pd.concat(res, sort=True)
    return df

def DP1(dfs, pal=None): # VERTICAL

    if len(dfs) <= 0:
        return None
    df1 = dfs[0]
    df1['Name'] = list(df1['Signature'])
    df1['Xl'] = list(df1['#Genes'])
    labels = [k['GSEID'][0] + ' ' + k['Condition'][0] for k in dfs]
    n1 = df1.shape[0]
    rocauc = list(df1['ROC-AUC'])
    p = list(df1['pvalue'])
    y = [1] * n1
    for i in range(1, len(dfs)):
        rocauc += list(dfs[i]['ROC-AUC'])
        p += list(dfs[i]['pvalue'])
        y += [i+1] * n1
    df = pd.DataFrame()
    df['ROC-AUC'] = rocauc
    df['pvalue'] = p
    df['ROC-AUC'] = df['ROC-AUC'].apply(
               lambda x: max([float(k) for k in str(x).split(",")]))
    df['pvalue'] = df['pvalue'].apply(
               lambda x: min([float(k) for k in str(x).split(",")]))
    df['Y'] = y
    df['R'] = [(i - 0.5) if i !=200  else 200  for i in df['ROC-AUC']]      #df['ROC-AUC'] - 0.5
    df['Ra'] = [abs(i)+ 0.5  if i<1  else 0  for i in df['R'] ]     # abs(df['R']) + 0.5
    df['Ra1'] = [abs(i)+ 0.5  if i<1  else 0.51 if i==0 else 0.48  for i in df['R'] ]

    df['AUC'] = ['Up' if i > 0 else 'Down' for i in df['R']]
    #['Up' if i > 0 else 'Dif i!=0 else 'NP' for i in df['R']]
    df['code'] = [getCode(i) for i in df['pvalue']]
    df['ROC-AUC'] = df['Ra1']
    sns.set()
    sns.set_style("white")
    sns.set_style({'xtick.color':'.5', 'ytick.color':'.5', 'axes.labelcolor': '.5'})
    sns.set_context("notebook")
    if pal is None:
        sns.set_palette([adj_light(c, 0.7, 1) for c in ['red', 'blue']])
    else:
        sns.set_palette(pal)
    x = [i + 1 for i in range(n1)] * len(labels)
    y = df['Y']
    fig, ax = plt.subplots(figsize=(12, len(dfs)*0.7+1), dpi=100)
    ax = sns.scatterplot(x=x, y=y, size="ROC-AUC", hue='AUC',
                         sizes = (0, 100), size_norm = (0.49, 1),
                         hue_order = ['Up', 'Down'], ax=ax, data=df);
    roc = list(df['Ra'])
    code = list(df['code'])
    for line in range(n1):
        ax.text(line + 1, len(labels) + .5, df1['Xl'][line],
                horizontalalignment='center', size='small', color='0.8',
                verticalalignment='bottom', rotation=90)
        for i in range(len(labels)):
            ax.text(line + 1, i + 0.9, "%.2f" % roc[line + n1 * i],
                    horizontalalignment='right', size='small', color='0.8',
                    verticalalignment='top',  rotation=90)
            ax.text(line + 1.5, i + 0.9, code[line + n1 * i],
                     horizontalalignment='right', size='small', color='0.8',
                     verticalalignment='top',  rotation=90)

    x1 = [i + 1 for i in range(n1)]
    ax.set_yticks(range(1, len(labels) + 1))
    ax.set_yticklabels(labels)
    ax.set_xlim([0, len(x1)+1])
    ax.set_ylim([0, len(labels) + 2])
    ax.set_xticks(x1)
    ax.set_xticklabels(df1['Name'], rotation=90)
    ax.set_ylabel("")
    ax.grid(False)
    handles, labels = ax.get_legend_handles_labels()
#     labels[4] = '0.5'
    ax.legend(handles, labels, bbox_to_anchor=(1.3, 1))

    return df,ax,fig

def MacComparisonDiverse():
    ana = MacAnalysis()
    res = []
    ana.getSvensson2011()
    res += [MacComparison(ana, desc="" + ana.source, Org='Hs')]
    ana.getMartinez2015()
    res += [MacComparison(ana, desc="" + ana.source, Org='Hs')]
    ana.getGuler2015(2) # mouse bone marrow derive macrophage
    res += [MacComparison(ana, desc="" + ana.source, Org='Mm')]
    ana.getHutchins2015TPM(2) # Bone marrow neutrophil
    res += [MacComparison(ana, desc="" + ana.source, Org='Mm')]
    ana.getHutchins2015TPM(3) # spleen-purified dendritic cells
    res += [MacComparison(ana, desc="" + ana.source, Org='Mm')]
    ana.getHutchins2015TPM(4) # Peritoneal exudate cells
    res += [MacComparison(ana, desc="" + ana.source, Org='Mm')]
    ana.getHutchins2015TPM(6) # Bone marrow-derived mast cell
    res += [MacComparison(ana, desc="" + ana.source, Org='Mm')]
    ana.getDas2018()
    res += [MacComparison(ana, desc="" + ana.source, Org='Mm')]
    ana.getPiccolo2017()
    res += [MacComparison(ana, desc="" + ana.source, Org='Mm')]
    ana.getOstuni2013()
    res += [MacComparison(ana, desc="" + ana.source, Org='Mm')]
    ana.getRochaResende()
    res += [MacComparison(ana, desc="" + ana.source, Org='Mm')]
    ana.getHill2018(2)
    res += [MacComparison(ana, desc="" + ana.source, Org='Mm')]
    ana.getFreemerman2019(2)
    res += [MacComparison(ana, desc="" + ana.source, Org='Mm')]
    ana.getLi2015()
    res += [MacComparison(ana, desc="" + ana.source, Org='Mm')]

    df = pd.concat(res, sort=True)
    return df

def MacComparisonTransplantOutcome():
    ana = MacAnalysis()
    res = []
    ana.getMorgun2006I(2)
    res += [MacComparison(ana, desc="Heart " + ana.source, Org='Hs')]
    ana.getMorgun2006II(2)
    res += [MacComparison(ana, desc="Heart " + ana.source, Org='Hs')]
    ana.getBohne2012II(3)
    res += [MacComparison(ana, desc="Liver " + ana.source, Org='Hs')]
    ana.getVanLoon2019(3)
    res += [MacComparison(ana, desc="Kidny " + ana.source, Org='Hs')]
    ana.getPalmer2019(2)
    res += [MacComparison(ana, desc="IBD blood " + ana.source, Org='Hs')]
    ana.getGharib2019Alv(2)
    res += [MacComparison(ana, desc="ARDS " + ana.source, Org='Hs')]
    ana.getSimpson2019(2)
    res += [MacComparison(ana, desc="Liver " + ana.source, Org='Hs')]

    df = pd.concat(res, sort=True)
    return df

def figure3p():
    #wt1, l1 = [-1, 1, 2], [readGenes(basedir + f"node-{i}.txt") for i in [13, 14, 3]]
    wt1, l1 = [-1], [readGenes(basedir + f"node-{i}.txt") for i in [13]]
    #wt1, l1 = [1, 2], [readGenes(basedir + f"node-{i}.txt") for i in [14, 3]]
    pathdesc = "13"
    ta = 1
    ana = MacAnalysis()
    adata = [ana, wt1, l1, pathdesc]
    res = []
    ana.getTang2020(4)
    res += [MacData(adata, desc="GSE116105", Org='Mm', tn=ta)]
    ana.getLink2018(5)
    res += [MacData(adata, desc="GSE109965", Org='Mm', tn=ta)]
    ana.getHowes2016(5)
    res += [MacData(adata, desc="GSE79809", Org='Mm', tn=ta)]
    return

def printTestSurvival(l1, wt1, fthr, ax = None):
    ana = MacAnalysis()
    ana.getJSTOM()
    ana.orderDataDf(l1, wt1)
    sax = ana.printSurvival(fthr, None, None, None)
    sax.set_title(ana.getTitle())
    return sax

def printTestSurvival2(l1, wt1, fthr, ax = None):
    ana = MacAnalysis()
    ana.getSurvival("LIV15.2")
    ana.orderDataDf(l1, wt1)
    sax = ana.printSurvival(fthr, None, None, None)
    sax.set_title(ana.getTitle())
    return sax

def printTestSurvival3(l1, wt1, fthr, ax = None):
    ana = MacAnalysis()
    ana.getSurvival("MAC104")
    ana.orderDataDf(l1, wt1)
    sax = ana.printSurvival(fthr, None, None, None)
    sax.set_title(ana.getTitle())
    return sax

def printTestSurvival4(l1, wt1, i1, i2, i3, ax = None):
    ana = MacAnalysis()
    ana.getSurvival("MACV74")
    ana.orderDataDf(l1, wt1)
    time = ana.getSurvName('time')
    status = ana.getSurvName('status')
    time1 = [time[i] for i in ana.order]
    status1 = [status[i] for i in ana.order]
    thr = hu.getThrData(ana.f_ranks)
    nm = (np.max(ana.f_ranks) - np.min(ana.f_ranks))/16
    print(thr, nm)
    fthr = -800
    fthr = thr[0]
    print(fthr)
    group1 = [1 if ana.f_ranks[i - ana.start] > fthr else 0 for i in ana.order]

    ana = MacAnalysis()
    ana.getSurvival("MACV75")
    ana.orderDataDf(l1, wt1)
    time = ana.getSurvName('time')
    status = ana.getSurvName('status')
    time2 = [time[i] for i in ana.order]
    status2 = [status[i] for i in ana.order]
    v1 = sorted([ana.f_ranks[i - ana.start] for i in ana.order])
    thr = hu.getThrData(v1[:-1])
    nm = ana.noisemargin
    #nm = (np.max(ana.f_ranks) - np.min(ana.f_ranks))/16
    print(thr, nm)
    fthr = -342
    fthr = thr[0] - int(nm)
    print(fthr)
    group2 = [1 if ana.f_ranks[i - ana.start] > fthr else 0 for i in ana.order]

    df = pd.DataFrame()
    df['time'] = time1 + time2
    df['status'] = status1 + status2
    df['group'] = group1 + group2
    g1 = [i for i in df.index if df['group'][i] == 0]
    g2 = [i for i in df.index if df['group'][i] == 1]
    pG = [ ["Low", "red", g1], ["High", "green", g2]]
    time, status = hu.censor(df['time'], df['status'], 3)
    ax = hu.survival(time, status, pG, ax)
    return ax

def printTestSurvival5(l1, wt1, fthr, ax = None):
    ana = MacAnalysis()
    ana.getSurvival("MACV82")
    ana.orderDataDf(l1, wt1)
    sax = ana.printSurvival(fthr, None, None, None)
    sax.set_title(ana.getTitle())
    return sax

def printTestSurvival6(l1, wt1, fthr, ax = None):
    ana = MacAnalysis()
    ana.getSurvival("MACV80")
    ana.orderDataDf(l1, wt1)
    thr = np.percentile(ana.f_ranks, [2, 100])
    v1 = [ana.f_ranks[i - ana.start] for i in ana.order
            if ana.f_ranks[i - ana.start] >= thr[0] and
            ana.f_ranks[i - ana.start] < thr[1]]
    fthr = hu.getThrData(v1)[0]
    fthr = fthr + ana.noisemargin
    sax = ana.printSurvival(fthr, None, 1500, ax)
    sax.set_title(ana.getTitle())
    return sax

def printTestSurvival7(l1, wt1, fthr, ax = None):
    ana = MacAnalysis()
    ana.getSurvival("MACV120.2")
    sex = ana.getSurvName("c Sex")
    st1 = [i for i in ana.aRange() if sex[i] == 'male']
    st2 = [i for i in ana.aRange() if sex[i] == 'female']
    ana.order = st1
    ana.orderDataDf(l1, wt1)
    fthr = hu.getThrData(sorted(ana.f_ranks)[1:-2])[0]
    sax = ana.printSurvival(fthr, None, 3.5, ax)
    sax.set_title(ana.getTitle())
    return sax

def figure4aOld():
    wt1, l1 = getCls13a14a3()
    ana = MacAnalysis()
    ana.getSurvival("CRC35.3")
    ana.orderDataDf(l1, wt1)
    nm = ana.noisemargin
    v1 = sorted([ana.f_ranks[i - ana.start] for i in ana.order])
    thr = hu.getThrData(v1)
    v2 = sorted([ana.f_ranks[i - ana.start] for i in ana.order
          if ana.f_ranks[i - ana.start] >= thr[0]])
    thr2 = hu.getThrData(v2[6:])
    fthr = thr2[0] - int(nm)
    sax = ana.printSurvival(fthr, None, None, None)
    sax.set_title(ana.getTitle())
    return sax

def figure4a():
    wt1, l1 = [-1, 1, 2], [readGenes(basedir + f"node-{i}.txt") for i in [13, 14, 3]]
    sax = printTestSurvival(l1, wt1, "thr+thr1")
    return sax

def figure4bOld():
    wt1, l1 = getCls13a14a3()
    ana = MacAnalysis()
    ana.getSurvival("LIV15.2")
    ana.orderDataDf(l1, wt1)
    nm = ana.noisemargin
    v1 = sorted([ana.f_ranks[i - ana.start] for i in ana.order])
    thr = hu.getThrData(v1[4:])
    fthr = thr[0] - int(nm)
    sax = ana.printSurvival(fthr, None, None, None)
    sax.set_title(ana.getTitle())
    return sax

def figure4b():
    wt1, l1 = [-1, 1, 2], [readGenes(basedir + f"node-{i}.txt") for i in [13, 14, 3]]
    sax = printTestSurvival2(l1, wt1, "thr1")
    return sax

def figure4cOld():
    wt1, l1 = getCls13()
    ana = MacAnalysis()
    ana.getSurvival("MAC104")
    ana.orderDataDf(l1, wt1)
    nm = ana.noisemargin
    v1 = sorted([ana.f_ranks[i - ana.start] for i in ana.order])
    thr = hu.getThrData(v1)
    v2 = sorted([ana.f_ranks[i - ana.start] for i in ana.order
          if ana.f_ranks[i - ana.start] >= thr[0]])
    thr2 = hu.getThrData(v2[15:])
    fthr = thr2[0] + nm
    sax = ana.printSurvival(fthr, None, None, None)
    sax.set_title(ana.getTitle())
    return sax

def figure4c():
    wt1, l1 = [-1], [readGenes(basedir + f"node-{i}.txt") for i in [13]]
    sax = printTestSurvival3(l1, wt1, "thr+thr1")
    return sax

def figure4d():
    wt1, l1 = [-1, 1, 2], [readGenes(basedir + f"node-{i}.txt") for i in [13, 14, 3]]
    sax = printTestSurvival4(l1, wt1, 0, 1.5, 0)
    sax.set_title("GSE28221")
    return sax

def figure4eOld():
    wt1, l1 = getCls13()
    ana = MacAnalysis()
    ana.getSurvival("MACV82")
    ana.orderDataDf(l1, wt1)
    nm = ana.noisemargin
    v1 = sorted([ana.f_ranks[i - ana.start] for i in ana.order])
    thr = hu.getThrData(v1[5:])
    fthr = thr[0] + int(nm)
    sax = ana.printSurvival(fthr, None, 1500, None)
    sax.set_title(ana.getTitle())
    return sax

def figure4fOld():
    wt1, l1 = getCls13()
    ana = MacAnalysis()
    ana.getSurvival("MACV80")
    ana.orderDataDf(l1, wt1)
    nm = ana.noisemargin
    v1 = sorted([ana.f_ranks[i - ana.start] for i in ana.order])
    thr = hu.getThrData(v1[5:])
    fthr = thr[0] + nm
    sax = ana.printSurvival(fthr, None, 1500, None)
    sax.set_title(ana.getTitle())
    return sax

def figure4e():
    wt1, l1 = [-1], [readGenes(basedir + f"node-{i}.txt") for i in [13]]
    sax = printTestSurvival5(l1, wt1, "thr1")
    return sax

def figure4f():
    wt1, l1 = [-1], [readGenes(basedir + f"node-{i}.txt") for i in [13]]
    sax = printTestSurvival6(l1, wt1, None)
    return sax

def figure4g():
    wt1, l1 = [-1, 1, 2], [readGenes(basedir + f"node-{i}.txt") for i in [13, 14, 3]]
    sax = printTestSurvival7(l1, wt1, None)

    ana = MacAnalysis()
    ana.getSurvival("MACV120.2")
    sex = ana.getSurvName("c Sex")
    st1 = [i for i in ana.aRange() if sex[i] == 'male']
    st2 = [i for i in ana.aRange() if sex[i] == 'female']
    ana.order = st1 + st2
    ana.orderDataDf(l1, wt1)
    time = ana.getSurvName('time')
    status = ana.getSurvName('status')
    ana.order = st1
    time2 = [time[i] for i in ana.order]
    status2 = [status[i] for i in ana.order]
    tlist2 = [ana.f_ranks[i - ana.start] for i in ana.order]
    res1 = hu.getBestThr(time2, status2, tlist2, range(len(time2)), None,
            3.5, len(time2))
    ana.order = st2
    time2 = [time[i] for i in ana.order]
    status2 = [status[i] for i in ana.order]
    tlist2 = [ana.f_ranks[i - ana.start] for i in ana.order]
    res2 = hu.getBestThr(time2, status2, tlist2, range(len(time2)), None,
            3.5, len(time2))
    df1 = pd.DataFrame(res1)
    df2 = pd.DataFrame(res2)
    df1[2] = -np.log(df1[0])/np.log(10)
    df2[2] = -np.log(df2[0])/np.log(10)
    ax = df1.plot.scatter(1, 2, color='blue')
    ax = df2.plot.scatter(1, 2, color='pink', ax=ax)
    ax.legend(["Male", "Female"])
    ax.axhline(y=1.3, color='red')

    return sax, ax

def figure4Old():
    pdf = getPDF("surv-old.pdf")
    figure4aOld()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    figure4bOld()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    figure4cOld()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    figure4d()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    figure4eOld()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    figure4fOld()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    sax, ax = figure4g()
    pdf.savefig(sax.get_figure(), transparent=True,bbox_inches = 'tight')
    pdf.savefig(ax.get_figure(), transparent=True,bbox_inches = 'tight')
    closePDF(pdf)

def figure4():
    pdf = getPDF("surv-1.pdf")
    figure4a()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    figure4b()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    figure4c()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    figure4d()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    figure4e()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    figure4f()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    sax, ax = figure4g()
    pdf.savefig(sax.get_figure(), transparent=True,bbox_inches = 'tight')
    pdf.savefig(ax.get_figure(), transparent=True,bbox_inches = 'tight')
    closePDF(pdf)

def figure6():
    wt1, l1 = [-1, 1, 2], [readGenes(basedir + f"node-{i}.txt") for i in [13, 14, 3]]
    wt2, l2 = [-1], [readGenes(basedir + f"node-{i}.txt") for i in [13]]
    wt3, l3 = [1, 2], [readGenes(basedir + f"node-{i}.txt") for i in [14, 3]]
    pathdesc = "13-14-3"
    ta = 1
    ana = MacAnalysis()
    adata = [ana, wt1, l1, "13-14-3"]
    adata2 = [ana, wt2, l2, "13"]
    adata3 = [ana, wt3, l3, "14-3"]
    res = []
    ana.getLee2012()
    res += [MacData(adata, desc="AIM2 OvExp", Org='Hs', tn=ta)]
    ana.getOakes2017Hs(2)
    res += [MacData(adata, desc="Oas2 OvExp", Org='Hs', tn=ta)]
    ana.getWang2019Mac(2)
    res += [MacData(adata, desc="Il15/IL15R", Org='Hs', tn=ta)]
    ta = 3
    ana.getMan2015(2)
    res += [MacData(adata, desc="Irf1 -/-", Org='Mm', tn=ta)]
    ana.getFensterl2012(2)
    res += [MacData(adata, desc="Ifit2 -/-", Org='Mm', tn=ta)]
    ana.getIrey2019Hs(2)
    res += [MacData(adata, desc="STAT3 inh", Org='Hs', tn=ta)]
    ana.getIrey2019Mm()
    res += [MacData(adata, desc="Stat3 -/-", Org='Mm', tn=ta)]
    ana.getOakes2017Mm(2)
    res += [MacData(adata3, desc="Oas2 (I405N)", Org='Mm', tn=ta)]
    ana.getGoldmann2015(3)
    res += [MacData(adata, desc="Usp18 -/-", Org='Mm', tn=5)]

    ta = 1
    ana.getPG2019lpsRep(2)
    res += [MacData(adata, desc="Ccdc88a -/-", Org='Mm', tn=ta)]
    ana.getLinke2017()
    res += [MacData(adata, desc="Tsc2 -/-", Org='Mm', tn=ta)]
    ana.getLi2019()
    res += [MacData(adata2, desc="Rnf5 -/-", Org='Mm', tn=ta)]
    ana.getUckelmann2020(2)
    res += [MacData(adata, desc="MLL inhibitor", Org='Hs', tn=ta)]
    ana.getUckelmann2020(3)
    res += [MacData(adata, desc="MLL inhibitor", Org='Hs', tn=ta)]
    ana.getUckelmann2020Mm(2)
    res += [MacData(adata, desc="MLL inhibitor", Org='Mm', tn=6)]
    ana.getElDahr2019(2)
    res += [MacData(adata, desc="Ezh1 -/-", Org='Mm', tn=6)]
    ana.getEncode2020(2)
    res += [MacData(adata2, desc="NFX1; PCBP2; EEF2; HNRNPA1 (shRNA)", Org='Hs', tn=ta)]
    return


def getIDhash(adata):
    idhash = {}
    for i in range(len(adata.raw.var['gene_symbols'])):
        k = adata.raw.var['gene_symbols'][i]
        v = adata.raw.var['gene_ids'][i]
        if k not in idhash:
            idhash[k] = []
        idhash[k] += [v]
    return idhash

def getRanks3(gene_groups, adata):
    idhash = getIDhash(adata)
    expr = []
    row_labels = []
    row_ids = []
    row_numhi = []
    ranks = []
    g_ind = 0
    counts = []
    for s in gene_groups:
        count = 0
        avgrank = [0] * adata.raw.n_obs
        for gn in s:
            if gn not in idhash:
                continue
            e = adata.raw.obs_vector(gn)
            v = e
            if (np.max(v) - np.min(v)) <= 0:
                continue
            t = hu.getThrData(v)
            te = []
            for i in range(len(e)):
                if e[i] == "":
                    v1 = - t[3] / 3;
                else:
                    v1 = (float(e[i]) - t[3]) / 3;
                if np.std(v) > 0:
                    v1 = v1 / np.std(v)
                avgrank[i] += v1
                te.append(v1)
            expr.append(te)
            row_labels.append(gn)
            row_ids.append(idhash[gn][0])
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
    print (counts)
    return ranks, row_labels, row_ids, row_numhi, expr
def computeSMART(adata, org='Hs'):
    wt1, l1 = getCls13()
    if org is 'Mm':
        l1 = getGroupsMm(l1)
    ranks, row_labels, row_ids, row_numhi, expr = getRanks3(l1, adata)
    f_ranks = mergeRanks(range(adata.raw.n_obs), 0, ranks, wt1)
    adata.obs['c13'] = f_ranks

    wt1, l1 = getCls14a3()
    if org is 'Mm':
        l1 = getGroupsMm(l1)
    ranks, row_labels, row_ids, row_numhi, expr = getRanks3(l1, adata)
    f_ranks = mergeRanks(range(adata.raw.n_obs), 0, ranks, wt1)
    adata.obs['c14_3'] = f_ranks
    adata.obs['c14'] = ranks[0]
    adata.obs['c3'] = ranks[1]
    return adata
def convertString(data):
    for k in data.obs.columns:
        data.obs[k] = [k.decode('utf-8') if type(k) == bytes else k for k in data.obs[k]]
    for k in data.var.columns:
        data.var[k] = [k.decode('utf-8') if type(k) == bytes else k for k in data.var[k]]
    data.var_names = list(data.var['gene_symbols'])
    data.var_names_make_unique()
    data.obs_names = [k.decode('utf-8') if type(k) == bytes else k for k in data.obs_names]
    return data
def scatterPlot(data, gA, gB, col="red"):
    import scanpy as sc
    plotdf = sc.get.obs_df(data, keys=[gA, gB])
    return plotdf.plot.scatter(gA, gB, c=col)
def computePCAandUMAP(adata):
    import scanpy as sc
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pl.highest_expr_genes(adata, n_top=20)
    sc.pp.highly_variable_genes(adata, min_mean=0.01, max_mean=5, min_disp=0.5)
    sc.pl.highly_variable_genes(adata)
    adata = adata[:, adata.var.highly_variable]
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata, log=True)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    return adata

def figure2f():
    import scanpy as sc
    from collections import Counter
    dir1 = "/booleanfs2/sahoo/Data/Macrophage/CoV/GSE168710/"
    cfile = dir1 + "GSE168710_QC_normalized.h5ad"
    data = sc.read(cfile)
    data.var['gene_symbols'] = data.var['name']
    data.var['gene_ids'] = data.var['name']
    convertString(data)
    data = data[data.obs['stim'].isin(['Control', 'IFNg+','IL4+']),:]
    scatterPlot(data, "TYROBP", "FCER1G")
    data.raw = data
    adata = data.raw.to_adata()
    adata.raw = data
    adata = computePCAandUMAP(adata)
    computeSMART(adata)
    plotdf = sc.get.obs_df(adata, keys=["c13", "c14_3", 'stim'])
    #plotdf['c desc'] = [re.sub(".*_", "", k) for k in plotdf['c desc']]
    ahash = {'Control':acolor[0], 'IFNg+':acolor[1],'IL4+':acolor[2]}
    plotdf['color'] = [ahash[k] for k in plotdf['stim']]
    ax = plotdf.plot.scatter("c13", "c14_3", c="color", alpha=0.5,
            rasterized=True)
    thr1 = hu.getThrData(plotdf['c13'])
    thr2 = hu.getThrData(plotdf['c14_3'])
    ax.axhline(y=thr2[0], color='cyan')
    ax.axvline(x=thr1[0], color='cyan')
    c1 = (plotdf['c13'] <= thr1[0])
    c2 = (plotdf['c14_3'] <= thr2[0])
    print(Counter(plotdf[c1 & c2]['stim']), Counter(plotdf['stim']))
    from statsmodels.stats.proportion import proportions_ztest
    print(proportions_ztest(2149, 2966, 62/5396), 2149/2966, 62/5936, 1/2829)
    print(thr1[0], thr2[0])

def figures9f1():
    wt1, l1 = getCls13a14a3()
    ana = MacAnalysis()
    ana.getBos()
    ana.orderDataDf(l1, wt1)
    ESR1 = "205225_at"
    id1 = ESR1
    expr = ana.getExprData(id1)
    low = ana.getArraysAll(id1, "thr0", "lo")
    high = ana.getArraysAll(id1, "thr0", "hi")
    v1 = sorted([ana.f_ranks[i - ana.start] for i in ana.order])
    thr = hu.getThrData(v1)
    nm = ana.noisemargin
    v2 = sorted([ana.f_ranks[i - ana.start] for i in ana.order
          if ana.f_ranks[i - ana.start] >= thr[0] - nm])
    thr2 = hu.getThrData(v2[:-4])
    fthr = thr2[0]
    #print(v2)
    print(thr)
    print(nm, fthr)
    g1 = [i for i in high if
            ana.f_ranks[i - ana.start] < fthr]
    g2 = [i for i in high if
            ana.f_ranks[i - ana.start] >= fthr]
    pG = [ ["Low", "red", g1], ["High", "green", g2]]
    time = ana.getSurvName('time')
    status = ana.getSurvName('status')
    sax = hu.survival(time, status, pG)
    sax.set_title(ana.getTitle())
    return sax

def figures9f2():
    wt1, l1 = getCls13()
    ana = MacAnalysis()
    ana.getSurvival("PANC2")
    ana.orderDataDf(l1, wt1)
    v1 = sorted([ana.f_ranks[i - ana.start] for i in ana.order])
    thr = hu.getThrData(v1[:-3])
    fthr = thr[0]
    sax = ana.printSurvival(fthr, None, None)
    sax.set_title(ana.getTitle())
    return sax

def figures9f3():
    wt1, l1 = getCls13a14a3()
    ana = MacAnalysis()
    ana.getSurvival("P9")
    ana.orderDataDf(l1, wt1)
    sax = ana.printSurvival(None, None, 60, None)
    sax.set_title(ana.getTitle())
    return sax

def figures9f4():
    wt1, l1 = getCls13a14a3()
    ana = MacAnalysis()
    ana.getSurvival("G2")
    ana.orderDataDf(l1, wt1)
    nm = ana.noisemargin
    v1 = sorted([ana.f_ranks[i - ana.start] for i in ana.order])
    thr = hu.getThrData(v1[13:])
    fthr = thr[0] + nm
    sax = ana.printSurvival(fthr, None, None)
    sax.set_title(ana.getTitle())
    return sax

def figures9f5():
    wt1, l1 = getCls13a14a3()
    ana = MacAnalysis()
    ana.getSurvival("B2.1")
    ana.orderDataDf(l1, wt1)
    nm = ana.noisemargin
    v1 = sorted([ana.f_ranks[i - ana.start] for i in ana.order])
    thr = hu.getThrData(v1)
    v2 = sorted([ana.f_ranks[i - ana.start] for i in ana.order
              if ana.f_ranks[i - ana.start] < thr[0] + nm])
    thr2 = hu.getThrData(v2[1:-1])
    fthr = int(thr2[0] + nm)
    sax = ana.printSurvival(fthr, None, None, None)
    sax.set_title(ana.getTitle())
    return sax

def figures9f():
    pdf = getPDF("surv-sup.pdf")
    figures9f1()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    figures9f2()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    figures9f3()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    figures9f4()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    figures9f5()
    pdf.savefig(transparent=True,bbox_inches = 'tight')
    closePDF(pdf)

class MacAnalysis:
    def __init__(self, urlbase=urlbase):
        self.state = []
        self.params = {}
        self.start = 2
        self.end = 2
        self.urlbase = urlbase
        return
    
    def aRange(self):
        return range(self.start, self.end + 1)

    def getTitle(self):
        title = self.name + " (" + self.source + "; n = " + str(self.num) + ")"
        return title
    
    def printInfo(self):
        print(self.name + " (n = " + str(self.num) + ")")
        url = "http://hegemon.ucsd.edu/Tools/explore.php?key=polyps&id="
        if self.dbid.startswith("MAC"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=mac&id="
        if self.dbid.startswith("MACV"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=macv&id="
        if self.dbid == "G16":
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=gbm&id="
        if self.dbid.startswith("GL"):
            url = "http://hegemon.ucsd.edu/Tools/explore.php?key=global&id="
        print(self.source + " " + url + self.dbid)
        print(len(self.order), [len(i) for i in self.state], \
              self.source, url + self.dbid, self.dbid)
        return
    
    def prepareDataDf(self, dbid, urlbase=urlbase):
        self.dbid = dbid
        self.dataset = hu.getHegemonDataset(self.dbid, urlbase)
        self.num = self.dataset[2]
        self.name = self.dataset[1]
        self.source = self.dataset[3]
        obj = hu.getHegemonPatientData(self.dbid, 'time', urlbase)
        self.headers = obj[0]
        self.hhash = {}
        self.start = 2;
        self.end = len(self.headers) - 1
        for i in range(len(self.headers)):
            self.hhash[self.headers[i]] = i
        return

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
        return
    
    def getSurvName(self, name):
        return hu.getHegemonPatientData(self.dbid, name, self.urlbase)[1]

    def getExprData(self, name):
        return hu.getHegemonData(self.dbid, name, "", self.urlbase)[1]
    
    def getThrData(self, name):
        obj = hu.getHegemonThrFrame(self.dbid, [name])
        thr = [obj['thr1'][0], obj['stat'][0], obj['thr0'][0], obj['thr2'][0]]
        return thr
    
    def getArraysThr (self, id1, thr = None, type1 = None):
        res = []
        expr = self.getExprData(id1);
        thr_step = self.getThrData(id1);
        thr = hu.getThrCode(thr_step, thr_step[0], thr);
        for i in self.aRange():
          if (thr is None):
             res.append(i)
          elif (expr[i] == ""):
             continue
          elif (type1 == "hi" and float(expr[i]) >= thr):
             res.append(i)
          elif (type1 == "lo" and float(expr[i]) < thr):
             res.append(i)
          elif (type1 is not None and type1 != "lo" and type1 != "hi" \
                  and float(expr[i]) >= thr and float(expr[i]) <= float(type1)): 
             res.append(i)
        return res

    def getArraysAll (self, *data):
        res = self.aRange()
        for i in range(0, len(data), 3):
          r = self.getArraysThr(data[i], data[i+1], data[i+2])
          res = list(set(res) & set(r))
        return res;

    def orderDataDf(self, gene_groups, weight):
        data_g = []
        data_e = []
        data_t = []
        for k in gene_groups:
            df_g = hu.getHegemonGeneIDs(self.dbid, k, self.urlbase)
            df_e = hu.getHegemonDataFrame(self.dbid, k, None, self.urlbase)
            df_t = hu.getHegemonThrFrame(self.dbid, k, self.urlbase)
            df_e.fillna(0,inplace=True)
            rhash = {}
            for i in range(df_t.shape[0]):
                rhash[df_t.iloc[i,0]] = i
            order = [rhash[df_e.iloc[i,0]] for i in range(df_e.shape[0])]
            df_t = df_t.reindex(order)
            df_t.reset_index(inplace=True)
            rhash = {}
            for i in df_e.index:
                rhash[df_e.iloc[i,0]] = i
            df_g['idx'] = [rhash[df_g.iloc[i,0]] if df_g.iloc[i,0] in rhash
                    else None for i in df_g.index]
            df_g = df_g.dropna()
            df_g['idx'] = df_g['idx'].astype(np.int64)
            data_g.append(df_g)
            data_e.append(df_e)
            data_t.append(df_t)
        self.col_labels = self.headers[self.start:]
        if len(gene_groups) > 0:
            self.col_labels = data_e[0].columns[self.start:]
        self.chash = {}
        for i in range(len(self.col_labels)):
            self.chash[self.col_labels[i]] = i
        compositres = getRanksDf2(gene_groups, data_g, data_e, data_t)
        ranks, noisemargins, row_labels, row_ids, row_numhi, expr = compositres
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
        f_nm = 0
        for i in range(len(gene_groups)):
            f_nm += abs(weight[i]) * noisemargins[i]
        self.noisemargin = 0.5/3
        if f_nm > 0:
            self.noisemargin = np.sqrt(f_nm)
        self.ranks = ranks
        self.row_labels = row_labels
        self.row_ids = row_ids
        self.row_numhi = row_numhi
        self.expr = expr
        self.i1 = i1
        self.index = index
        return

    def getScores(ana, ahash = None):
        lval = [[] for i in ana.atypes]
        cval = ana.cval[0]
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

    def getMacMetrics(ana, actual, ahash = None):
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
        return res

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

    def getTPvalSpecific(ana, m=0, n=1):
        res = []
        lval, score = ana.getScores()
        if m < 0 or n < 0 or m >= len(lval) or n >= len(lval):
            return 0, 1
        t, p = ttest_ind(lval[m],lval[n], equal_var=False)
        return t, p

    def getTPvals(ana):
        res = []
        lval, score = ana.getScores()
        for k in range(1, len(ana.atypes)):
            t, p = ttest_ind(lval[0],lval[k], equal_var=False)
            res += ["%.3g" % p]
        return ",".join(res)

    def getStats(self, l1, wt1, annotation=[], id1=None):
        src = re.sub(" .*", "", self.source)
        species = annotation[1]
        if species == 'Hs' or species == 'Rm' :
            self.orderDataDf(l1, wt1)
        else:
            l1 = getGroupsMm(l1)
            self.orderDataDf(l1, wt1)
        if id1 is not None:
            self.normGene(id1)
        roc = self.getROCAUC()
        p = self.getTPvals()
        lval, score = self.getScores()
        n1 = np.sum([len(lval[i])for i in range(1, len(lval))])
        return [src, roc, p, len(lval[0]), n1] + annotation
        
    def normGene(ana, id1):
        expr = ana.getExprData(id1)
        df = pd.DataFrame()
        df['x'] = [float(expr[i]) for i in ana.order]
        df['y'] = [ana.f_ranks[i - ana.start] for i in ana.order]
        df['a'] = [ana.aval[i] for i in ana.order]
        for i in range(len(ana.atypes)):
            cond1 = (df['a'] == i)
            if (sum(cond1) > 0):
                s1 = np.max(df[cond1]['y']) - np.min(df[cond1]['y'])
                s2 = np.max(df[cond1]['x']) - np.min(df[cond1]['x'])
                df.loc[cond1, 'y'] += (np.mean(df[cond1]['x']) - df.loc[cond1, 'x']) * (s1+1) / (s2+1)
                df.loc[cond1, 'x'] -= (np.mean(df[cond1]['y']) - df.loc[cond1, 'y']) * (s2+1) / (s1+1)
        from sklearn.linear_model import LinearRegression
        #linreg = LinearRegression(normalize=True)
        linreg = LinearRegression()
        linreg.fit(np.array(df['x']).reshape(-1, 1),df['y'])
        y_pred = linreg.predict(np.array(df['x']).reshape(-1, 1))
        df['y1'] = (df['y'] - y_pred)
        ana.f_ranks = df['y1']
        ana.i1 = [ana.order[i] for i in np.argsort(ana.f_ranks)]
        ana.f_ranks = [0 for i in ana.aRange()]
        for i in range(len(ana.order)):
            ana.f_ranks[ana.order[i] - ana.start] = df['y1'][i]
        index = np.array([i - ana.start for i in ana.i1])
        ana.cval = np.array([[ana.aval[i] for i in ana.i1]])
        ana.data = np.array([ana.expr[i] for i in ana.ind_r])[:,index]
        return

    def printTitleBar(self, params):
        self.params = {'atypes': self.atypes,'cval': self.cval}
        self.params.update(params)
        ax = plotTitleBar(self.params['cval'], \
                self.params['atypes'], self.params)
        return ax

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

    def printHeatmap(self, ofile, genes, params):
        i1 = self.i1
        f_ranks = self.f_ranks
        self.params = {'genes': genes,'atypes': self.atypes,'cval': self.cval}
        self.params.update(params)
        ax, divider = plotHeatmap(ofile, self.data, self.col_labels, 
                self.row_labels, self.params)
        return

    def printSurvival(self, fthr = None, pG = None, ct = None, ax = None):
        f_ranks = self.f_ranks
        order = self.order
        thr = hu.getThrData(f_ranks)
        #thr = hu.getThrData([f_ranks[i - self.start] for i in order])
        nm1 = (np.max(f_ranks) - np.min(f_ranks))/16
        nm = self.noisemargin
        if fthr is None:
            fthr = thr[0]
        if fthr == "thr0":
            fthr = thr[0] - nm
        if fthr == "thr1":
            fthr = thr[0]
        if fthr == "thr2":
            fthr = thr[0] + nm
        if fthr == "thr.5":
            fthr = thr[0] + 0.5 * nm
        if fthr == "thr-.5":
            fthr = thr[0] - 0.5 * nm
        if fthr == "thr1.5":
            fthr = thr[0] + 1.5 * nm
        if fthr == "thr2.5":
            fthr = thr[0] + 2.5 * nm
        if fthr == "thr3":
            fthr = thr[0] + 3 * nm
        if fthr == "thr+thr1":
            v1 = [f_ranks[i - self.start] for i in order
                    if f_ranks[i - self.start] >= thr[0]]
            thr2 = hu.getThrData(v1)
            print(thr2)
            fthr = thr2[0]
        if fthr == "thr-thr1":
            v1 = [f_ranks[i - self.start] for i in order
                    if f_ranks[i - self.start] < thr[0]]
            thr2 = hu.getThrData(v1)
            print(thr2)
            fthr = thr2[0]
        if fthr == "avg":
            fthr = np.mean([f_ranks[i - self.start] for i in order])
        if fthr == "med":
            fthr = np.median([f_ranks[i - self.start] for i in order])
        print(thr)
        print(fthr, nm, nm1, thr[0] + nm1)
        g1 = [i for i in order if f_ranks[i - self.start] < fthr]
        g2 = [i for i in order if f_ranks[i - self.start] >= fthr]
        if pG is None:
            pG = [ ["Low", "red", g1], ["High", "green", g2]]
        time = self.getSurvName('time')
        status = self.getSurvName('status')
        if ct is not None:
            time, status = hu.censor(time, status, ct)
        sax = hu.survival(time, status, pG, ax)
        return sax

    def getDataset(self, dbid):
        self.prepareDataDf(dbid)
        atype = self.getSurvName("time")
        atypes = ['All']
        atype = [atypes[0] for i in atype]
        ahash = {}
        self.initData(atype, atypes, ahash)
        return

    def getGEOMacAnn(self):
        self.prepareDataDf("G16")
        atype = self.getSurvName("c Type")
        atypes = ['M0', 'M1', "M2"]
        ahash = {}
        self.initData(atype, atypes, ahash)
        return

    def getBeyer2012(self):
        self.prepareDataDf("MAC1")
        atype = self.getSurvName("c cell type")
        atypes = ['M0', 'M1', 'M2']
        ahash = {'M1 macrophages':1, 'M0 macrophages':0, 'M2 macrophages':2}
        self.initData(atype, atypes, ahash)
        return

    def getXue2014(self, mt=1):
        self.prepareDataDf("MAC2")
        atype = self.getSurvName("c Type")
        atypes = ['M0', 'M1', 'M2']
        ahash = {'M_GMCSF_IL4_72h':2, 'M_GMCSF_IFNg_72h':1, 'M0_GMCSF_72h':0}
        if mt == 2:
            ahash = {'M_GMCSF_IL4_24h':2,'M_GMCSF_IFNg_24h':1,'M0_GMCSF_24h':0}
        if mt == 3:
            ahash = {'M_GMCSF_IL4_12h':2,'M_GMCSF_IFNg_12h':1,'M0_GMCSF_12h':0}
        if mt == 4:
            ahash = {\
                    'M0_GMCSF_0h':0, 'M0_GMCSF_12h':0, 'M0_GMCSF_24h':0,
                    'M0_GMCSF_48h':0, 'M0_GMCSF_6h':0, 'M0_GMCSF_72h':0,
                    'M0_MCSF_0h':0, 'M1/2_GMCSF_24h':0, 'M_GMCSF_IFNg_30min':0,
                    'M_GMCSF_IFNg_1h':0, 'M_GMCSF_IFNg_2h':0,
                    'M_GMCSF_IFNg_4h':1, 'M_GMCSF_IFNg_6h':1,
                    'M_GMCSF_IFNg_12h':1, 'M_GMCSF_IFNg_24h':1,
                    'M_GMCSF_IFNg_72h':1,
                    'M_GMCSF_IL4_30min':2, 'M_GMCSF_IL4_1h':2,
                    'M_GMCSF_IL4_2h':2, 'M_GMCSF_IL4_4h':2, 'M_GMCSF_IL4_6h':2,
                    'M_GMCSF_IL4_12h':2, 'M_GMCSF_IL4_24h':2,
                    'M_GMCSF_IL4_72h':2, 'M_MCSF_IL4_72h':2}
        self.initData(atype, atypes, ahash)
        return

    def getOhradanovaRepic2018(self, tn=1):
        self.prepareDataDf("MAC14")
        atype = self.getSurvName("c treatment")
        atypes = ['M0', 'M1', 'M2', 'IL10']
        ahash = {'mock-activated (medium only; control) for 2d':0,
                'activated with 100 ng/ml LPS + 25ng/ml IFN\xce\xb3 for 2d':1,
                'activated with 20 ng/ml IL-4 for 2d':2,
                'activated with 20 ng/ml IL-10 for 2d':3}
        if (tn == 2):
            atypes = ['M0', 'M1', 'M2']
            ahash = {'mock-activated (medium only; control) for 2d':0,
                    'activated with 100 ng/ml LPS + 25ng/ml IFN\xce\xb3 for 2d':1,
                    'activated with 20 ng/ml IL-4 for 2d':2}
        if (tn == 3):
            atypes = ['M0', 'IL10']
            ahash = {'mock-activated (medium only; control) for 2d':0,
                    'activated with 20 ng/ml IL-10 for 2d':1}
        ahash = asciiNorm(ahash)
        self.initData(atype, atypes, ahash)
        return

    def getZhang2015(self, tn=1):
        self.prepareDataDf("MAC3")
        atype = self.getSurvName("c Type")
        atypes = ['M', 'iM', 'M1', 'iM1', 'M2', 'iM2', 'i']
        ahash = {'IPSDM M2':5, 'IPSDM MAC':1, 'IPSDM M1':3,
                'iPS':6, 'HMDM M1':2, 'HMDM MAC':0, 'HMDM M2':4}
        if (tn == 2):
            atypes = ['M0', 'M1', 'M2']
            ahash = {'IPSDM M2':2, 'IPSDM MAC':0, 'IPSDM M1':1}
        if (tn == 3):
            atypes = ['M0', 'M1', 'M2']
            ahash = {'HMDM M1':1, 'HMDM MAC':0, 'HMDM M2':2}
        self.initData(atype, atypes, ahash)
        return

    def getRock2005(self):
        self.prepareDataDf("MAC73")
        atype = self.getSurvName('c Group')
        atypes = ['Control', 'IFNG']
        ahash = {}
        self.initData(atype, atypes, ahash)
        return

    def getPolak2014(self, tn=1):
        self.prepareDataDf("MAC69")
        atype = self.getSurvName('c treatment')
        ahash = {'TNF-alpha':0, 'TSLP':1}
        rval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName('c time')
        atypes = ['0', '2', '8', '24']
        ahash = {'8h':2, '2h':1, '24h':3, '0h':0}
        if (tn == 2):
            atypes = ['C', 'TNF-a']
            ahash = {'8h':1, '24h':1, '0h':0}
            atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)
        return

    def getQu2016(self, tn=1):
        self.prepareDataDf("PLP50")
        atype = self.getSurvName("c tissue type")
        atypes = ['NS', 'NC', 'A', 'C', 'M']
        ahash = {'Metastasis':4,
                'Carcinoma':3,
                'Normal crypt epithelium':1,
                'Adenoma':2,
                'Normal surface epithelium':0}
        if (tn == 2):
            atypes = ['NC', 'NS']
            ahash = {'Normal crypt epithelium':0,
                    'Normal surface epithelium':1}
        if (tn == 3):
            atypes = ['NC', 'A']
            ahash = {'Normal crypt epithelium':0, 'Adenoma':1}
        if (tn == 4):
            atypes = ['NC', 'C']
            ahash = {'Normal crypt epithelium':0, 'Carcinoma':1}
        self.initData(atype, atypes, ahash)
        return

    def getImmGenULI(self, tn=1):
        self.prepareDataDf("MAC81")
        atype = self.getSurvName('c Title')
        atypes = ['M0', 'M1']
        ahash = {}
        if (tn == 1):
            ahash = {'MF.11cpSigFp.BAL.1':0,
                    'MF.11cpSigFp.BAL.2':0,
                    'MF.SSChipSigFn.LPS.d3.BAL.1':1,
                    'MF.SSChipSigFn.LPS.d3.BAL.2':1,
                    'MF.11cpSigFp.LPS.d6.BAL.1':1,
                    'MF.11cpSigFp.LPS.d6.BAL.2':1,
                    'MF.SSChipSigFn.LPS.d6.BAL.1':1,
                    'MF.SSChipSigFn.LPS.d6.BAL.2':1}
            ahash = {'MF.11cpSigFp.BAL.1':0,
                    'MF.11cpSigFp.BAL.2':0,
                    'MF.SSChipSigFn.LPS.d3.BAL.1':1,
                    'MF.SSChipSigFn.LPS.d3.BAL.2':1}
        if (tn == 2):
            atypes = ['M0', 'M1']
            ahash = {'MF.64p6Cn206nIIp.LPS.d3.Lu.1':1,
                    'MF.64p6Cn206nIIp.LPS.d3.Lu.2':1,
                    'MF.64p6Cn206nIIp.LPS.d6.Lu.1':1,
                    'MF.64p6Cn206nIIp.LPS.d6.Lu.2':1,
                    'MF.64p6Cn206nIIp.Lu.1':0,
                    'MF.64p6Cn206nIIp.Lu.2':0}
        if (tn == 3):
            atypes = ['M0', 'M1']
            ahash = {'MF.64p6Cn206pIIn.LPS.d3.Lu.1':1,
                    'MF.64p6Cn206pIIn.LPS.d3.Lu.2':1,
                    'MF.64p6Cn206pIIn.LPS.d6.Lu.1':1,
                    'MF.64p6Cn206pIIn.LPS.d6.Lu.2':1,
                    'MF.64p6Cn206pIIn.Lu.1':0,
                    'MF.64p6Cn206pIIn.Lu.2':0}
        if (tn == 4):
            atypes = ['M0', 'M1']
            ahash = {'MF.F.10kIFN.PC#1':1,
                    'MF.F.10kIFN.PC#2':1,
                    'MF.F.10kIFN.PC#3':1,
                    'MF.F.PC#1':0,
                    'MF.F.PC#1_':0,
                    'MF.F.PC#2':0,
                    'MF.F.PC#2_':0,
                    'MF.F.PC#3':0,
                    'MF.F.PC.1':0,
                    'MF.F.PC.2':0,
                    'MF.F.PC.3':0,
                    'MF.Fem.PC#1_RNA-seq':0,
                    'MF.Fem.PC#2_RNA-seq':0}
        if (tn == 5):
            atypes = ['M0', 'M1']
            ahash = {'MF.M.10kIFN.PC#1':1,
                    'MF.M.10kIFN.PC#2':1,
                    'MF.M.10kIFN.PC#3':1,
                    'MF.M.PC#1':0,
                    'MF.M.PC#2':0,
                    'MF.M.PC#3':0}
            for k in ['MF.PC#1.1', 'MF.PC#1.10', 'MF.PC#1.11', 'MF.PC#1.12',
                    'MF.PC#1.13', 'MF.PC#1.14', 'MF.PC#1.15', 'MF.PC#1.16',
                    'MF.PC#1.2', 'MF.PC#1.3', 'MF.PC#1.4', 'MF.PC#1.5', 'MF.PC#1.6',
                    'MF.PC#1.7', 'MF.PC#1.8', 'MF.PC#1.9', 'MF.PC#1_RNA-seq',
                    'MF.PC#2.1', 'MF.PC#2.2', 'MF.PC#2.3', 'MF.PC#2.4', 'MF.PC#2.5',
                    'MF.PC#2.6', 'MF.PC#2.7', 'MF.PC#2.8', 'MF.PC#2_RNA-seq',
                    'MF.PC#3', 'MF.PC#3_RNA-seq', 'MF.PC#4', 'MF.PC#4_RNA-seq',
                    'MF.PC.01', 'MF.PC.02', 'MF.PC.03', 'MF.PC.04', 'MF.PC.05',
                    'MF.PC.06', 'MF.PC.07', 'MF.PC.08', 'MF.PC.09', 'MF.PC.10',
                    'MF.PC.11', 'MF.PC.12', 'MF.PC.13', 'MF.PC.14', 'MF.PC.15',
                    'MF.PC.17', 'MF.PC.18', 'MF.PC.19', 'MF.PC.20', 'MF.PC.21',
                    'MF.PC.23', 'MF.PC.24', 'MF.PC.25', 'MF.PC.26', 'MF.PC.37'
                    'MF.PC.38', 'MF.PC.39', 'MF.PC.40']:
                ahash[k] = 0
        if (tn == 6):
            atypes = ['M0', 'M1']
            ahash = {'Mo.6Chi11bp.APAP.36h.Lv.1':1,
                    'Mo.6Chi11bp.APAP.36h.Lv.2':1,
                    'Mo.6Chi11bp.APAP.36h.Lv.3':1,
                    'Mo.6Chi11bp.APAP.36h.Lv.4':1,
                    'Mo.6Chi11bp.PBS.Lv.1':0,
                    'Mo.6Chi11bp.PBS.Lv.2':0,
                    'Mo.6Chi11bp.PBS.Lv.3':0,
                    'Mo.6Chi11bp.PBS.Lv.4':0}
        if (tn == 7):
            atypes = ['M0', 'M1']
            ahash = {
                    'NKT.F.Sp#1':0,
                    'NKT.F.Sp#2':0,
                    'NKT.F.Sp#3':0,
                    'NKT.M.Sp#1':0,
                    'NKT.M.Sp#2':0,
                    'NKT.M.Sp#3':0,
                    'NKT.Sp#3_RNA-seq':0,
                    'NKT.Sp.LPS.18hr#1_RNA-seq':1,
                    'NKT.Sp.LPS.18hr#2_RNA-seq':1,
                    'NKT.Sp.LPS.3hr#1_RNA-seq':1,
                    'NKT.Sp.LPS.3hr#2_RNA-seq':1}
        if (tn == 8):
            atypes = ['M0', 'M1']
            ahash = {
                    'T4.F.10kIFN.Sp#1':1,
                    'T4.F.10kIFN.Sp#2':1,
                    'T4.F.10kIFN.Sp#3':1,
                    'T4.F.Sp#1':0,
                    'T4.F.Sp#2':0,
                    'T4.F.Sp#3':0,
                    'T4.M.10kIFN.Sp#1':1,
                    'T4.M.10kIFN.Sp#2':1,
                    'T4.M.10kIFN.Sp#3':1,
                    'T4.M.Sp#1':0,
                    'T4.M.Sp#2':0,
                    'T4.M.Sp#3':0}
        if (tn == 9):
            atypes = ['M0', 'M1']
            ahash = {
                    'B.17m.F.Sp#1':0,
                    'B.20m.Sp#1':0,
                    'B.2m.F.Sp#1':0,
                    'B.2m.Sp#1':0,
                    'B.6m.F.Sp#1':0,
                    'B.6m.F.Sp#2':0,
                    'B.6m.Sp#1':0,
                    'B.6m.Sp#2':0,
                    'B.F.10kIFN.Sp#1':1,
                    'B.F.10kIFN.Sp#2':1,
                    'B.F.10kIFN.Sp#3':1,
                    'B.F.1kIFN.Sp#1':1,
                    'B.F.1kIFN.Sp#2':1,
                    'B.F.1kIFN.Sp#3':1,
                    'B.F.Sp#1':0,
                    'B.F.Sp#1_':0,
                    'B.F.Sp#2':0,
                    'B.F.Sp#2_':0,
                    'B.F.Sp#3':0,
                    'B.Fem.Sp#1_RNA-seq':0,
                    'B.Fo.Sp#1_RNA-seq':0,
                    'B.Fo.Sp#2_RNA-seq':0,
                    'B.Fo.Sp#3_RNA-seq':0,
                    'B.Fo.Sp#4_RNA-seq':0}
        if (tn == 10):
            atypes = ['M0', 'M1']
            ahash = {
                    'GN.17m.F.Sp#1':0,
                    'GN.20m.Sp#1':0,
                    'GN.F.10kIFN.Sp#1':1,
                    'GN.F.10kIFN.Sp#2':1,
                    'GN.F.10kIFN.Sp#3':1,
                    'GN.F.Sp#1':0,
                    'GN.F.Sp#2':0,
                    'GN.F.Sp#3':0,
                    'GN.M.10kIFN.Sp#1':1,
                    'GN.M.10kIFN.Sp#2':1,
                    'GN.M.10kIFN.Sp#3':1,
                    'GN.M.Sp#1':0,
                    'GN.M.Sp#2':0,
                    'GN.M.Sp#3':0,
                    'GN.Sp#3_RNA-seq':0,
                    'GN.Sp#4_RNA-seq':0}
        if (tn == 11):
            atypes = ['M0', 'M1', 'M2']
            ahash = {
                    'MF.KC.Clec4FpTim4p64p.APAP.12h.Lv.1':2,
                    'MF.KC.Clec4FpTim4p64p.APAP.12h.Lv.2':2,
                    'MF.KC.Clec4FpTim4p64p.APAP.12h.Lv.4':2,
                    'MF.KC.Clec4FpTim4p64p.APAP.36h.Lv.1':2,
                    'MF.KC.Clec4FpTim4p64p.APAP.36h.Lv.2':2,
                    'MF.KC.Clec4FpTim4p64p.APAP.36h.Lv.3':2,
                    'MF.KC.Clec4FpTim4p64p.APAP.36h.Lv.4':2,
                    'MF.KC.Clec4FpTim4p64p.Lv.2':0,
                    'MF.KC.Clec4FpTim4p64p.Lv.3':0,
                    'MF.KC.Clec4FpTim4p64p.Lv.4':0,
                    'MF.KC.Clec4FpTim4p64p.PBS.Lv.1':0,
                    'MF.KC.Clec4FpTim4p64p.PBS.Lv.2':0,
                    'MF.KC.Clec4FpTim4p64p.PBS.Lv.3':0,
                    'MF.KC.Clec4FpTim4p64p.PBS.Lv.4':0}
        self.initData(atype, atypes, ahash)
        return

    def getSvensson2011(self, tn=1):
        self.prepareDataDf("MAC87")
        atype = self.getSurvName('c factors')
        atypes = ['M0', 'M1', 'M2']
        ahash = {'M-CSF + E/P/IL-10/4/13':2, 'M-CSF + IL-10':2,
                'GM-CSF':0, 'M-CSF':0, 'M-CSF + IL4/13':2,
                'GM-CSF + M-CSF':0, 'GM-CSF/LPS/IFN':1, 'GM-CSF + IL4/13':2,
                'GM-CSF + IL-10':2, 'Decidual macrophages':0, 'Blood monocytes':0,
                'GM-CSF + M-CSF +IL-10':2, 'M-CSF + M-CSF':0, 'GM-CSF/LPS/IFN D6':1,
                'GM-CSF + E/P/IL-10/4/13':2}
        if (tn == 2):
            atypes = ['GMCSF', 'MCSF']
            ahash = {'GM-CSF':0, 'GM-CSF + M-CSF':0, 'M-CSF':1, 'M-CSF + M-CSF':1}
        self.initData(atype, atypes, ahash)
        return

    def getChandriani2014(self, tn=1):
        self.prepareDataDf("MAC88")
        source = self.getSurvName('c src1')
        atype = self.getSurvName('c treatment')
        atypes = ['M0', 'M1', 'M2']
        ahash = {'unstimulated':0, 'IL13':2, 'TGFb':1, 'IL10':2, 'IL4':2,
                'Dex':1}
        if (tn == 2):
            atype = source
            ahash = {'Monocytes, unstimulated, 24h':0, 'Monocytes, IL13, 24h':2,
                    'Monocytes, unstimulated, 6h':0, 'Monocytes, IL4, 24h':2,
                    'Monocytes, IL4, 6h':2, 'Monocytes, IL13, 6h':2}
        if (tn == 3):
            atype = source
            ahash = {'Macrophages, unstimulated, 24h':0, 'Macrophages, IL13, 24h':2,
                    'Macrophages, IL10, 24h':2, 'Macrophages, TGFb, 24h':1,
                    'Macrophages, IL4, 24h':2, 'Macrophages, Dex, 24h':1}
        if (tn == 4):
            atype = source
            ahash = {'Normal lung fibroblasts, TGFb, 24h':1,
                    'Normal lung fibroblasts, IL13, 24h':2,
                    'Normal lung fibroblasts, IL4, 24h':2,
                    'Normal lung fibroblasts, unstimulated, 24h':0}
        if (tn == 5):
            atype = source
            atypes = ['Mono', 'Mac']
            ahash = {'Monocytes, unstimulated, 24h':0,
                    'Monocytes, unstimulated, 6h':0,
                    'Macrophages, unstimulated, 24h':1}
        self.initData(atype, atypes, ahash)
        return

    def getMartinez2015(self, tn=1):
        self.prepareDataDf("MAC89")
        atype = self.getSurvName('c src1')
        atypes = ['M0', 'M1', 'M2']
        ahash = { 'Monocyte-derived macrophages polarized with IL-4 for 5 days':2,
                'Monocyte-derived macrophages polarized with IL-10 for 5 days':2,
                'Monocyte-derived macrophages polarized with IFNgTNFa for 5 days':1,
                'Monocyte-derived macrophages':0}
        self.initData(atype, atypes, ahash)
        return

    def getGuler2015(self, tn=1):
        self.prepareDataDf("MAC84.3")
        atype = self.getSurvName('c Title')
        atype = ["_".join(str(k).split("_")[0:-2]) 
                if len(str(k).split("_")) > 2 else None for k in atype]
        atypes = ['M0', 'M1', 'M2']
        ahash = {'IFNg':1,
                'IL4IL13':2,
                'M.tb_IL4IL13':0,
                'M.tb':1,
                'Ust':0,
                'M.tb_IFNg':1,
                'M.tb_IL41L13':0}
        if (tn == 2):
            atype = self.getSurvName('c Title')
            atype = ["_".join(str(k).split("_")[0:-1]) 
                    if len(str(k).split("_")) > 2 else None for k in atype]
            ahash = {
                    'Ust_28h':0,
                    'IL4IL13_28h':2,
                    'M.tb_IFNg_28h':1,
                    'IFNg_28h':1,
                    'M.tb_28h':1}
        self.initData(atype, atypes, ahash)
        return

    def getHutchins2015TPM(self, tn = 1):
        self.prepareDataDf("MAC41.2")
        atype = self.getSurvName("c source_name")
        ahash = {'Bone marrow neutrophil':2,
                'spleen-purified dendritic cells':3,
                'Peritoneal exudate cells (adherent cells)':4,
                'Bone marrow-derived Eosinophils':5,
                'Bone marrow-derived mast cell':6,
                'CD4+ na\xc3\xafve T cells':7}
        ahash = asciiNorm(ahash)
        rval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c treatment")
        atypes = ['M0', 'M1', 'M2']
        ahash = {'LPS':1, 'IL10':2,'IL10\\, LPS':1, 'None':0}
        aval = [ahash[i] if i in ahash else None for i in atype]
        if (tn >= 2):
            #atype = [atype[i] if rval[i] == tn or aval[i] == 0
            #        else None for i in range(len(atype))]
            atype = [atype[i] if rval[i] == tn
                    else None for i in range(len(atype))]
        self.rval = rval
        self.initData(atype, atypes, ahash)

    def getDas2018(self, tn=1):
        self.prepareDataDf("MAC90")
        atype = self.getSurvName('c src1')
        atypes = ['M0', 'M1', 'M2']
        ahash = {'cultured for 4 hrs':0,
                'treated with IFN-\xce\xb3 (100 U/ml) and LPS (100 ng/ml) for 12 hrs':1,
                'treated with IFN-\xce\xb3 (100 U/ml) and LPS (100 ng/ml) for 4 hrs':1,
                'treated with IFN-\xce\xb3 (100 U/ml) and LPS (100 ng/ml) for 24 hrs':1,
                'treated with IFN-\xce\xb3 (100 U/ml) and LPS (100 ng/ml) for 1 hr':1,
                'treated with LPS (100 ng/ml) for 4 hrs':1,
                'treated with IL-13 (10\xe2\x80\x89ng/ml) for 12 hrs':2,
                'treated with IL-4 (10\xe2\x80\x89ng/ml) for 12 hrs':2}
        ahash = asciiNorm(ahash)
        self.initData(atype, atypes, ahash)
        return

    def getDaniel2018(self, tn=1):
        self.prepareDataDf("MAC91")
        atype = self.getSurvName('c Title')
        atype = [re.sub(".rep.*", "", str(k)) for k in atype]
        atype = [re.sub("mm_BMDM_", "", str(k)) for k in atype]
        atypes = ['M0', 'M1', 'M2']
        ahash = {'Wt_2nd_stim_ctrl_RNA':0,
                'Wt_1st_stim_3hIL4_RNA':2,
                'PpargKO_2nd_stim_3hIL4_RNA':2,
                'ctrl_24hVeh_RNA':0,
                'PpargKO_1st_stim_3hIL4_RNA':2,
                'PpargKO_1st_stim_ctrl_RNA':0,
                'Wt_1st_stim_ctrl_RNA':0,
                'Wt_2nd_stim_3hIL4_RNA':2,
                'PpargKO_2nd_stim_ctrl_RNA':0}
        self.initData(atype, atypes, ahash)
        return

    def getPiccolo2017(self, tn=1):
        self.prepareDataDf("MAC92")
        atype = self.getSurvName('c Name')
        atype = [re.sub("_R.*", "", str(k)) for k in atype]
        atype = [re.sub("^_", "", str(k)) for k in atype]
        atypes = ['M0', 'M1', 'M2']
        ahash = {'IL4_4h':2,
                'IFNy_2h':1,
                'UT':0,
                'IFNy_IL4_2h':1,
                'IL4_2h':2,
                'IFNy_4h':1,
                'IFNy_IL4_4h':1}
        if (tn == 2):
            atype = self.getSurvName('c Name')
            atype = [re.sub("_R.*", "", str(k)) for k in atype]
            ahash = {'_shMYC_UT':0,
                    '_scramble_IL-4_4h':2,
                    '_scramble_IL-4_2h':2,
                    '_scramble_UT':0,
                    '_shMYC_IL-4_4h':2,
                    '_shMYC_IL-4_2h':2}
            if (tn == 3):
                atype = self.getSurvName('c Name')
            atype = [re.sub("_R.*", "", str(k)) for k in atype]
            ahash = {'shCEBP-beta_UT':0,
                    'shJunB_UT':0,
                    'shJunB_IFNy_4h':1,
                    'scramble_UT':0,
                    'scramble_IFNy_4h':1,
                    'shCEBP-beta_IFNy_4h':1}
        self.initData(atype, atypes, ahash)
        return

    def getOstuni2013(self, tn=1):
        self.prepareDataDf("MAC93")
        atype = self.getSurvName('c treatment')
        atypes = ['M0', 'M1', 'M2']
        ahash = {'TGFb1 (1 ng/ml) for 4hrs':1,
                'IL4 (10 ng/ml) for 4hrs':2,
                'No treatment':0,
                'IFNg (100 ng/ml) for 4hrs':1,
                'TNFa (10 ng/ml) for 4hrs':1}
        self.initData(atype, atypes, ahash)
        return

    def getRochaResende(self, tn=1):
        self.prepareDataDf("MAC94.2")
        atype = self.getSurvName('c treatment')
        atypes = ['M0', 'M1', 'M2']
        ahash = {'none':0, 'LPS':1, 'IL-4':2}
        self.initData(atype, atypes, ahash)
        return

    def getHill2018(self, tn=1):
        self.prepareDataDf("MAC95")
        atype = self.getSurvName('c condition')
        atypes = ['M0', 'M1', 'M2']
        ahash = {'Cd9+ macrophage transfer':1,
                'Ly6c+ macrophage transfer':1,
                'PBS transfer':0,
                'High fat diet':2,
                'IL4':2,
                'LPS':1,
                'Veh':0}
        if (tn == 2):
            ahash = {'PBS transfer':0,
                    'IL4':2,
                    'LPS':1,
                    'Veh':0}
        self.initData(atype, atypes, ahash)
        return

    def getFreemerman2019(self, tn=1):
        self.prepareDataDf("MAC96")
        atype = self.getSurvName('c Title')
        atype = [re.sub("\..*", "", str(k)) for k in atype]
        atypes = ['M0', 'M1', 'M2']
        ahash = {'GLUT1 WT M1':1,
                'GLUT1 KO M2':2,
                'GLUT1 WT M0':0,
                'GLUT1 WT M2':2,
                'GLUT1 KO M1':1,
                'GLUT1 KO M0':0}
        if (tn == 2):
            ahash = {'GLUT1 WT M1':1,
                    'GLUT1 WT M0':0,
                    'GLUT1 WT M2':2}
            if (tn == 3):
                ahash = {'GLUT1 KO M2':2,
                        'GLUT1 KO M1':1,
                        'GLUT1 KO M0':0}
        self.initData(atype, atypes, ahash)
        return

    def getEl2010(self, tn=1):
        self.prepareDataDf("MAC97")
        atype = self.getSurvName('c Title')
        atype = [re.sub(" [A-D]$", "", str(k)) for k in atype]
        atypes = ['M0', 'M1', 'M2']
        ahash = {'Irf4+/- IL4 4h':2,
                'Irf4-/- IL4 18h':2,
                'Irf4-/- Mock':0,
                'Irf4+/- IL4 18h':2,
                'Irf4+/- Mock':0,
                'Irf4-/- IL4 4h':2}
        if (tn == 2):
            ahash = {'Irf4+/- IL4 4h':2,
                    'Irf4+/- IL4 18h':2,
                    'Irf4+/- Mock':0}
            if (tn == 3):
                ahash = {'Irf4-/- IL4 18h':2,
                        'Irf4-/- Mock':0,
                        'Irf4-/- IL4 4h':2}
        self.initData(atype, atypes, ahash)
        return

    def getLi2015(self, tn=1):
        self.prepareDataDf("MAC98")
        atype = self.getSurvName('c src1')
        atypes = ['M0', 'M1', 'M2']
        ahash = {'Mouse macrophage at M2 treated with nutlin-3a':2,
                'Mouse macrophage at M2':2,
                'Mouse macrophage at M2 treated with 10058F4':2,
                'Mouse macrophage at M1':1,
                'Mouse macrophage at M0':0}
        self.initData(atype, atypes, ahash)
        return

    def getHan2017(self, tn=1):
        self.prepareDataDf("MAC25")
        atype = self.getSurvName("c Title")
        atype = [str(i).split(" ")[2] if len(str(i).split(" ")) > 3 else i \
                for i in atype]
        atypes = ['SR1078', 'M1', 'M0', 'Veh', 'M2', 'SR3335']
        ahash = {}
        if tn == 2:
            atypes = ['M0', 'M1', 'M2']
        self.initData(atype, atypes, ahash)
        return

    def getJeffrey2006(self, tn=1):
        self.prepareDataDf("MACV52")
        src = self.getSurvName('c src1')
        ahash = {'Cord blood':0, 'Peripheral blood':1}
        rval = [ahash[i] if i in ahash else None for i in src]
        title = self.getSurvName('c Title')
        atype = [1 if str(k).find("unstimulated") >= 0 or \
                str(k).find("control") >= 0 or \
                str(k).find("Immature") >= 0 else 0 for k in title]
        atypes = ['C', 'S']
        ahash = {0:1, 1:0}
        if (tn == 2):
            atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)
        return

    def getSippel2018(self, tn=1):
        self.prepareDataDf("MACV53")
        atype = self.getSurvName('c treatment')
        atypes = ['V', 'IL5', 'PGD2']
        ahash = {'IL5':1, 'Vehicle':0, 'dkPGD2':2}
        if (tn == 2):
            atypes = ['V', 'T']
            ahash = {'IL5':1, 'Vehicle':0, 'dkPGD2':1}
        self.initData(atype, atypes, ahash)
        return

    def getPuan2017(self, tn=1):
        self.prepareDataDf("MACV54")
        state = self.getSurvName('c donor_state')
        ahash = {'reactive':0, 'anergic':1}
        rval = [ahash[i] if i in ahash else None for i in state]
        atype = self.getSurvName('c stimulation')
        atypes = ['C', 'S']
        ahash = {'unstimulated':0, 'Fc-epsilon receptor-crosslinking':1}
        self.pair = [ [2, 11], [3, 6], [7, 5], [4, 9], [10, 8]]
        if (tn == 2):
            atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)
        return

    def getKlocperk2020(self, tn=1):
        self.prepareDataDf("MACV55")
        src = self.getSurvName('c src1')
        src = [re.sub("mo.*of ", "", str(k)) for k in src]
        treat = [re.sub(".* [1-5] *", "", str(k)) for k in src]
        treat = [re.sub("c.* with ", "", str(k)) for k in treat]
        ahash = {'':0,
                'autologous healthy NETs':1,
                'healthy NETs':2,
                'T1D NETs':3,
                'autologous T1D NETs':4}
        rval = [ahash[i] if i in ahash else None for i in treat]
        disease = [k.split(" ")[0] for k in src]
        atype = disease
        atypes = ['healthy', 'T1D']
        ahash = {}
        if (tn == 2):
            atypes = ['C', 'S']
            atype = rval
            ahash = {0:0, 1:1, 2:1, 3:1, 4:1}
            atype = [atype[i] if disease[i] == 'healthy' \
                    else None for i in range(len(atype))]
        if (tn == 3):
            atypes = ['C', 'S']
            atype = rval
            ahash = {0:0, 1:1, 2:1, 3:1, 4:1}
            atype = [atype[i] if disease[i] == 'T1D' \
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)
        return

    def getLenardo2020(self, tn=1):
        self.prepareDataDf("MACV58")
        atype = self.getSurvName('c time')
        atypes = ['D2', ' ', 'D5']
        ahash = {'Day2':0, 'Day5':2}
        self.initData(atype, atypes, ahash)
        return

    def getNair2015(self, tn=1):
        self.prepareDataDf("MACV59")
        atype = self.getSurvName('c protocol')
        atypes = ['C', 'S']
        ahash = {'unstimulated (control)':0, 'stimulated':1}
        self.initData(atype, atypes, ahash)
        return

    def getAbbas2005(self, tn=1):
        self.prepareDataDf("MACV61")
        atype = self.getSurvName("c src1")
        ahash = {'NK cells from PBMC':0,
                'Plasma cells from bone marrow':1,
                'Monocytes from PBMC':2,
                'CD4+ T cells from PBMC':3,
                'B cells from PBMC':4,
                'Neutrophils from PBMC':5,
                'CD14+ cells from PBMC':6,
                'CD4+ CD45RO+ CD45RA- T cells from PBMC':7,
                'CD8+ T cells from PBMC':8,
                'Plasma cells from PBMC':9}
        rval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName('c treatment agent')
        ahash = {'NA':0,
                'macrophage differentiation medium':1,
                'IL2':2,
                'LPS':3,
                'IL15':4,
                'aCD3/aCD28':5}
        tval = [ahash[i] if i in ahash else None for i in atype]
        atype = tval
        atypes = ['C', 'S']
        ahash = {0:0, 1:1, 2:1, 3:1, 4:1, 5:1}
        self.initData(atype, atypes, ahash)
        return

    def getMetcalf2015(self, tn=1):
        self.prepareDataDf("MACV62")
        atype = self.getSurvName('c treatment')
        atypes = ['C', 'S']
        ahash = {'Rig I':1, 'PolyIC':1, 'NoTx':0, 'LyoVec_only':1, 'LPS':1, 'CLO_97':1}
        self.initData(atype, atypes, ahash)
        return

    def getBanchereau2014I(self, tn=1):
        self.prepareDataDf("MACV63")
        atype = self.getSurvName('c cell population')
        ahash = {'IL4 DC':1, 'IFNa DC':0}
        rval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName('c stimulation')
        atypes = ['C', 'S']
        ahash = {'MDP':1, 'LPS':1, 'Poly I:C-LMW-Lyovec':1,
                'TNFa':1, 'CpG 2216':0, 'Poly I:C':1, 'R837':1,
                'CL097':1, 'IFNa':1, 'IL10':1, 'CpG 2006':0, 'Flagellin':1,
                'PAM3':1, 'IL15':1, 'IL1b':1}
        aval = [ahash[i] if i in ahash else None for i in atype]
        if (tn == 2):
            atypes = ['M0', 'M1', 'M2']
            atype = rval
            ahash = {0:1, 1:2}
            atype = [atype[i] if aval[i] == 0 \
                    else None for i in range(len(atype))]
        if (tn == 3):
            ahash['A'] = 1
            atype = [atype[i] if rval[i] == 1 \
                    else 'A' for i in range(len(atype))]
        if (tn == 4):
            atypes = ['C', 'S']
            ahash = {'CpG 2216':0, 'CpG 2006':0, 'IL1b':1}
        self.initData(atype, atypes, ahash)
        return

    def getBanchereau2014II(self, tn=1):
        self.prepareDataDf("MACV64")
        atype1 = self.getSurvName('c culture conditions')
        atype2 = self.getSurvName('c culture condition')
        atype3 = self.getSurvName('c vaccine abbr.')
        atype = [ " ".join([str(k) for k in [atype1[i], atype2[i], atype3[i]]]) 
                         for i in range(len(atype1))]
        atypes = ['C', 'S']
        ahash = {'  RAB':1, 'LPS  ':1, ' HKSE ':1,
                ' Media ':0, '  Medium':0, 'Medium1  ':0,
                '  FZ':1, '  TDAP':1, ' H1N1 Brisbane ':1, ' HKSA ':1, 'Medium2  ':1,
                'HEPB  ':1, 'HPV  ':1, 'HIB  ':1, '  POL':1, '  VAR':1, 'VAR  ':1,
                'PVX  ':1, 'RAB  ':1, 'MGL  ':1, '  HIB':1, '  HER':1, '  HPV':1,
                'HER  ':1, '  PVX':1, '  JPE':1, '  HEPB':1, 'HEPA  ':1, 'JPE  ':1,
                '  HEPA':1, 'FZ  ':1, 'POL  ':1, '  MGL':1, 'TDAP  ':1}
        aval = [ahash[i] if i in ahash else None for i in atype]
        if (tn == 2):
            atype = self.getSurvName('c cell population')
            atypes = ['M0', 'M1', 'M2']
            ahash = {'IL4 DC':2, 'IFNa DC':1, 'Monocytes':0}
            atype = [atype[i] if aval[i] == 0 \
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)
        return

    def getIampietro2017(self, tn=1):
        self.prepareDataDf("MACV66")
        atype = self.getSurvName('c infection')
        atypes = ['Mock', 'EBOV', 'LPS']
        ahash = {}
        if (tn == 2):
            atypes = ['C', 'S']
            ahash = {'Mock':0, 'EBOV':1, 'LPS':1}
        self.initData(atype, atypes, ahash)
        return

    def getJohnson2020(self, tn=1):
        self.prepareDataDf("MACV60")
        atype = self.getSurvName("c treatment_code1")
        atypes = ['C', 'S']
        ahash = {'NS':0}
        for k in atype:
            if k != 'NS':
                ahash[k] = 1
        if (tn == 2):
            atypes = ['C', 'HIV']
            ahash = {'NS': 0, 'HIV2_WT':1, 'HIV2_P86HA':1}
        self.initData(atype, atypes, ahash)
        return

    def getPeters(self, tn=1):
        self.prepareDataDf("PLP7")
        atype = self.getSurvName("c clinical condition")
        atypes = ['Normal', 'UC', 'CD']
        ahash = {"control": 0, "Ulcerative Colitis":1, "Crohn's disease":2}
        aval = [ahash[i] if i in ahash else None for i in atype]
        if (tn == 2):
            atype = self.getSurvName("c gender")
            atypes = ['F', 'M']
            ahash = {'female':0, 'male':1}
            atype = [atype[i] if aval[i] == 0 else None
                    for i in range(len(atype))]
        if (tn == 3):
            atypes = ['Normal', 'UC']
            ahash = {"control": 0, "Ulcerative Colitis":1}
        if (tn == 4):
            atypes = ['Normal', 'CD']
            ahash = {"control": 0, "Crohn's disease":1}
        self.initData(atype, atypes, ahash)
        return

    def getGkouskou2016(self, tn=1):
        self.prepareDataDf("PLP84")
        tissue = self.getSurvName("c tissue")
        ahash = {'proximal colon':3, 'distal colon':4}
        rval = [ahash[i] if i in ahash else None for i in tissue]
        atype = self.getSurvName("c src1")
        atypes = ['normal', 'AD2', 'AD4', 'proximal', 'distal']
        ahash = {'UNTREATED':0, 'AOM, 4 DSS CYCLES':2, 'AOM, 2 DSS CYCLES':1}
        if (tn == 2):
            atype = [str(atype[i])+ " " + str(tissue[i]) for i in
                    range(len(atype))]
            atypes = ['proximal', 'distal']
            ahash = {'UNTREATED proximal colon':0,
                    'UNTREATED distal colon':1}
        if (tn == 3):
            atypes = ['normal', 'colitis']
            ahash = {'UNTREATED':0, 'AOM, 4 DSS CYCLES':1, 'AOM, 2 DSS CYCLES':1}
        self.initData(atype, atypes, ahash)

    def getGkouskou2016ProxDis(self):
        self.prepareDataDf("PLP85")
        atype = self.getSurvName("c tissue")
        atypes = ['proximal', 'distal']
        ahash = {'PROXIMAL COLON':0, 'DISTAL COLON':1}
        self.initData(atype, atypes, ahash)

    def getWoetzel2014(self, tn=1):
        self.prepareDataDf("MAC31")
        atype = self.getSurvName("c disease state")
        atype2 = self.getSurvName("c clinical status")
        atype = [atype[i] + atype2[i] for i in range(len(atype))]
        atypes = ['HC', 'RA', 'OA']
        ahash = {'healthy control':0,
                'rheumatoid arthritis':1,
                'synovial tissue isolated from osteoarthritic joint':2,
                'osteoarthritis':2,
                'normal control':0}
        if (tn == 2):
            atypes = ['HC', 'OA']
            ahash = {'healthy control':0,
                    'synovial tissue isolated from osteoarthritic joint':1,
                    'osteoarthritis':1, 'normal control':0}
        if (tn == 3):
            atypes = ['HC', 'RA']
            ahash = {'healthy control':0, 'normal control':0,
                    'rheumatoid arthritis':1}
        self.initData(atype, atypes, ahash)

    def getLefebvre2017(self, tn=1):
        self.prepareDataDf("LIV3")
        atype = self.getSurvName("c Title")
        atype = [str(k).replace("liver biopsy ", "") for k in atype]
        self.patient = [str(k).split(" ")[0] for k in atype]
        atype = [str(k).split(" ")[1] if len(str(k).split(" ")) > 1 else "" for
                k in atype]
        self.paired = [re.sub("\((.*)\)", "\\1", str(k)) for k in atype]
        self.time = self.getSurvName("c time")
        self.treatment = self.getSurvName("c type of intervention")
        atype = self.getSurvName("c src1")
        atypes = ['b', 'f', 'Nb', 'Nf', 'ub', 'uf']
        ahash = {'no NASH liver baseline':0,
                'no NASH liver follow-up':1,
                'NASH liver follow-up':3,
                'NASH liver baseline':2,
                'undefined liver baseline':4,
                'undefined liver follow-up':5}
        self.aval = [ahash[i] if i in ahash else None for i in atype]
        phash = {}
        for i in self.aRange():
            phash[self.patient[i]] = i
        self.rtype = [None, None] + [str(self.aval[phash[self.paired[i]]])+" "+\
                str(self.aval[i])+" " + str(self.treatment[i]) \
                if self.paired[i] in phash and self.paired[i] != "" \
                else str(self.aval[i])+" " + str(self.treatment[i]) \
                for i in self.aRange()]
        fhash = {}
        for i in self.aRange():
            if self.paired[i] in phash and self.paired[i] != "":
                fhash[self.paired[i]] = i
        self.ftype = [None, None] + [str(self.aval[i])+" "+str(self.treatment[i])\
                + " " + str(self.aval[fhash[self.patient[i]]])\
                if self.patient[i] in fhash\
                else str(self.aval[i])+" " + str(self.treatment[i]) \
                for i in self.aRange()]
        if (tn == 2):
            atypes = ['b', 'Nb', 'ub']
            ahash = {'no NASH liver baseline':0,
                    'NASH liver baseline':1,
                    'undefined liver baseline':2}
        if (tn == 3):
            atypes = ['b', 'Nb']
            ahash = {'no NASH liver baseline':0,
                    'NASH liver baseline':1}
        if (tn == 4):
            atypes = ['f', 'Nf']
            ahash = {'no NASH liver follow-up':0,
                    'NASH liver follow-up':1}
        if (tn == 5):
            atypes = ['R', 'NR']
            atype = self.ftype
            ahash = {'2 Diet 1': 0, '2 Diet 3': 1}
        self.initData(atype, atypes, ahash)

    def getSuppli2019(self, tn=1):
        self.prepareDataDf("LIV14")
        atype = self.getSurvName("c disease")
        atypes = ['H', 'O', 'FL', 'SH']
        ahash = {'NAFLD':2, 'NASH':3, 'obese':1, 'healthy':0}
        if (tn == 2):
            atypes = ['H', 'FL']
            ahash = {'NAFLD':1, 'healthy':0}
        if (tn == 3):
            atypes = ['H', 'SH']
            ahash = {'NASH':1, 'healthy':0}
        if (tn == 4):
            atypes = ['H', 'O']
            ahash = {'obese':1, 'healthy':0}
        self.initData(atype, atypes, ahash)

    def getWoodruff2005(self, tn=1):
        self.prepareDataDf("MAC12")
        atype = self.getSurvName("c status")
        atypes = ['NS', 'S', 'A']
        ahash = {'Asthmatic':2, 'Smoker':1, 'Nonsmoker':0}
        if (tn == 2):
            atypes = ['H', 'A']
            ahash = {'Asthmatic':1, 'Smoker':0, 'Nonsmoker':0}
        if (tn == 3):
            atypes = ['NS', 'S']
            ahash = {'Smoker':1, 'Nonsmoker':0}
        self.initData(atype, atypes, ahash)

    def getWS2009(self, tn=1):
        self.prepareDataDf("MAC12.2")
        atype = self.getSurvName("c Type")
        atypes = ['NS', 'S', 'A', 'C']
        ahash = {'Nonsmoker':0, 'Asthmatic':2, 'Smoker':1, 'smoker':1,
                'non-smoker':0, 'COPD':3}
        if (tn == 2):
            atypes = ['NS', 'A']
            ahash = {'Nonsmoker':0, 'Asthmatic':1, 'non-smoker':0}
        if (tn == 3):
            atypes = ['NS', 'COPD']
            ahash = {'Nonsmoker':0, 'non-smoker':0, 'COPD':1}
        if (tn == 4):
            atypes = ['NS', 'S']
            ahash = {'Nonsmoker':0, 'non-smoker':0, 'Smoker':1, 'smoker':1}
        self.initData(atype, atypes, ahash)

    def getLissner(self, tn=1):
        self.prepareDataDf("MAC56")
        atype = self.getSurvName('c Agent')
        ahash = {'Lm':0, 'LPS':1}
        bval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName('c Timepoint')
        ahash = {'6hr':6, '2hr':2, '1hr':1, '0hr':0}
        rval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName('c Type')
        atypes = ['N', 'A', 'O']
        ahash = {'Neonate':0, 'Adult':1, 'OlderAdult':2}
        if (tn == 2):
            atype = [atype[i] if bval[i] == 0 and rval[i] == 0 else None
                    for i in range(len(atype))]
        if (tn == 3):
            atype = [atype[i] if bval[i] == 1 and rval[i] <= 1 else None
                    for i in range(len(atype))]
        if (tn == 4):
            atype = [atype[i] if bval[i] == 0 and rval[i] == 0 else None
                    for i in range(len(atype))]
            atypes = ['A', 'N']
            ahash = {'Neonate':1, 'Adult':0}
        if (tn == 5):
            atype = [atype[i] if bval[i] == 0 and rval[i] == 0 else None
                    for i in range(len(atype))]
            atypes = ['A', 'O']
            ahash = {'OlderAdult':1, 'Adult':0}
        self.initData(atype, atypes, ahash)

    def getBadawi2019(self, tn=1):
        self.prepareDataDf("MHP5")
        atype = self.getSurvName('c time post surgery')
        atypes = ['0', '45m', '24h']
        ahash = {'45 minutes':1, '24 hours':2, '0 minutes':0}
        if (tn == 2):
            atypes = ['0', '24h']
            ahash = {'45 minutes':0, '24 hours':1, '0 minutes':0}
        self.initData(atype, atypes, ahash)

    def getBondar2017(self, tn=1):
        self.prepareDataDf("MAC76")
        atype = self.getSurvName('c gender')
        ahash = {'Female':0, 'Male':1}
        tval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName('c underlying disease')
        atypes = ['NICM', 'ICM', 'PPCM', 'NCIM', 'ChemoCM']
        if (tn == 2):
            atypes = ['NICM', 'ICM']
        if (tn == 3):
            atypes = ['NICM', 'ICM']
            atype = [atype[i] if tval[i] == 1 else None
                for i in range(len(atype))]
        ahash = {}
        self.initData(atype, atypes, ahash)

    def getPatel2019(self, tn = 1):
        self.prepareDataDf("AD8")
        atype = self.getSurvName('c tissue');
        ahash = {'Temporal_Cortex':3,
                'Cerebellum':4,
                'Frontal_Cortex':5,
                'Entorhinal_Cortex':6}
        rval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName('c disease state');
        atypes = ['N', 'A', 'AD']
        ahash = {'AsymAD':1, 'AD':2, 'control':0}
        if (tn >= 2):
            atypes = ['N', 'AD']
            ahash = {'AD':1, 'control':0}
        if (tn > 2):
            atype = [atype[i] if rval[i] == tn else None for i in range(len(atype))]
        self.rval = rval
        self.initData(atype, atypes, ahash)

    def getBerchtold2014RMA(self, tn=1):
        self.prepareDataDf("AD11")
        atype = self.getSurvName('c AD specific');
        ahash = {'0':0, '1':1}
        rval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c src1")
        atype = [str(i) for i in atype]
        res = []
        for k in atype:
            l1 = k.split(",")
            if (len(l1) != 4):
                res.extend([k])
            else:
                res.extend([l1[1].strip() + "_" + l1[2].strip()])
        atype = res
        atypes = ['N', 'AD']
        ahash = {'entorhinal cortex_male':0,
                'entorhinal cortex_male_AD':1,
                'entorhinal cortex_female':0,
                'entorhinal cortex_female_AD':1,
                'superior frontal gyrus_male':0,
                'superior frontal gyrus_male_AD':1,
                'superior frontal gyrus_female':0,
                'superior frontal gyrus_female_AD':1,
                'postcentral gyrus_male':0,
                'post-central gyrus_male_AD':1,
                'postcentral gyrus_female':0,
                'post-central gyrus_female_AD':1,
                'hippocampus_male':0,
                'hippocampus_male_AD':1,
                'hippocampus_female':0,
                'hippocampus_female_AD':1}
        if (tn == 4):
            ahash = {'entorhinal cortex_male':0,
                    'entorhinal cortex_male_AD':1,
                    'superior frontal gyrus_male':0,
                    'superior frontal gyrus_male_AD':1,
                    'postcentral gyrus_male':0,
                    'post-central gyrus_male_AD':1,
                    'hippocampus_male':0,
                    'hippocampus_male_AD':1}
        if (tn == 5):
            ahash = {'entorhinal cortex_female':0,
                    'entorhinal cortex_female_AD':1,
                    'superior frontal gyrus_female':0,
                    'superior frontal gyrus_female_AD':1,
                    'postcentral gyrus_female':0,
                    'post-central gyrus_female_AD':1,
                    'hippocampus_female':0,
                    'hippocampus_female_AD':1}
        if (tn == 6):
            ahash = {'entorhinal cortex_male':0,
                    'entorhinal cortex_male_AD':1,
                    'entorhinal cortex_female':0,
                    'entorhinal cortex_female_AD':1}
        if (tn == 7):
            ahash = {'superior frontal gyrus_male':0,
                    'superior frontal gyrus_male_AD':1,
                    'superior frontal gyrus_female':0,
                    'superior frontal gyrus_female_AD':1}
        if (tn == 8):
            ahash = {'postcentral gyrus_male':0,
                    'post-central gyrus_male_AD':1,
                    'postcentral gyrus_female':0,
                    'post-central gyrus_female_AD':1}
        if (tn == 9):
            ahash = {'hippocampus_male':0,
                    'hippocampus_male_AD':1,
                    'hippocampus_female':0,
                    'hippocampus_female_AD':1}
        if (tn == 2):
            atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
        if (tn == 3):
            atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getGelman2012(self, tn = 1):
        self.prepareDataDf("AD51")
        atype = self.getSurvName('c phenotype');
        atypes = ['N', 'HAND']
        ahash = {'HIV Infected':0,
                'HIV Infected with neurocognitive impairment (HAD: HIV-associated dementia)':1,
                'HIV Infected with HAD and HIV encephalitis (HIVE)':1,
                'normal (control)':0,
                'HIV Infected with HAD':1,
                'HIV Infected with HAD and encephalitis':1}
        self.initData(atype, atypes, ahash)
        
    def getChenPlotkin2008(self, tn = 1):
        self.prepareDataDf("AD52")
        atype = self.getSurvName("c src1")
        atype = [str(k).split(" ")[1] if len(str(k).split(" ")) > 1 else None
                for k in atype]
        atype = [str(k).split("-")[0] for k in atype]
        atypes = ['N', 'P', 'S']
        ahash = {'Normal':0,
                'Progranulin':1,
                'Sporadic':2}
        if tn == 2:
            atypes = ['N', 'FTD']
            ahash = {'Normal':0,
                    'Progranulin':1}
        if tn == 3:
            atype = self.getSurvName('c disease and tissue')
            atypes = ['N', 'FTD']
            ahash = {'Normal hippocampus':0,
                'Progranulin hippocampus':1,
                'Sporadic hippocampus':1}
        self.initData(atype, atypes, ahash)
        
    def getOlmosSerrano2016(self, tn = 1):
        self.prepareDataDf("AD53")
        atype = self.getSurvName('c disease status');
        atypes = ['N', 'DS']
        ahash = {'CTL':0,
                'DS':1}
        if tn == 2:
            atype = self.getSurvName('c disease and tissue');
            atypes = ['N', 'DS']
            ahash = {'CTL ITC':0,
                'CTL STC':0,
                'CTL HIP':0,
                'DS ITC':1, 'DS STC':1, 'DS HIP':1}
        self.initData(atype, atypes, ahash)
        
    def getBartolettiStella2019 (self, tn = 1):
        self.prepareDataDf("AD54")
        atype = self.getSurvName('c condition');
        atypes = ['N', 'CJD']
        ahash = {'Control':0, 'sCJD affected':1}
        self.initData(atype, atypes, ahash)
        
    def getTsalik2015(self, tn=1):
        self.prepareDataDf("MACV261")
        atype = self.getSurvName("c sirs outcomes")
        atypes = ['SHK', 'SS', 'SIRS', 'US', 'SD']
        ahash = {'Septic shock':0,
                'severe sepsis':1,
                'SIRS':2,
                'Uncomplicated sepsis':3,
                'sepsis death':4}
        if (tn == 2):
            atype = self.getSurvName("c sirs vs sepsis")
            atypes = ['Sepsis', 'SIRS']
            ahash = {}
        self.initData(atype, atypes, ahash)

    def getBarcella2018(self, tn=1):
        self.prepareDataDf("MAC115.1")
        ctype = self.getSurvName("c clinical classification")
        ttype = self.getSurvName("c timepoint")
        atype = [ str(ctype[i]) + " " + str(ttype[i]) for i in
                                range(len(ctype))]
        atypes = ['R T1', 'R T2', 'NR T1', 'NR T2']
        ahash = {}
        if (tn == 2):
            atypes = ['R', 'NR']
            ahash = {'R T2':0, 'NR T1':1}
        self.initData(atype, atypes, ahash)

    def getStenvers2019(self, tn=1):
        self.prepareDataDf("MAC113")
        atype = self.getSurvName("c subject status")
        atypes = ['H', 'T2D']
        ahash = {'Type 2 diabetes':1, 'Healthy':0}
        if (tn >= 2):
            stype = self.getSurvName("c subject status")
            ttype = self.getSurvName("c timepoint")
            atype = [ str(stype[i]) + " " + str(ttype[i]) for i in
                    range(len(stype))]
            atypes = ['H1', 'H2', 'H3', 'H4', 'D1', 'D2', 'D3', 'D4']
            ahash = {'Type 2 diabetes D2_ZT_15:30':4,
                    'Type 2 diabetes D3_ZT_0:15':5,
                    'Type 2 diabetes D3_ZT_5:45':6,
                    'Type 2 diabetes D3_ZT_11:15':7,
                    'Healthy D2_ZT_15:30':0,
                    'Healthy D3_ZT_0:15':1,
                    'Healthy D3_ZT_5:45':2,
                    'Healthy D3_ZT_11:15':3}
        if (tn == 3):
            atypes = ['H1', 'H2', 'H3', 'H4']
            ahash = {'Healthy D2_ZT_15:30':0,
                    'Healthy D3_ZT_0:15':1,
                    'Healthy D3_ZT_5:45':2,
                    'Healthy D3_ZT_11:15':3}
        if (tn == 4):
            atypes = ['D1', 'D2', 'D3', 'D4']
            ahash = {'Type 2 diabetes D2_ZT_15:30':0,
                    'Type 2 diabetes D3_ZT_0:15':1,
                    'Type 2 diabetes D3_ZT_5:45':2,
                    'Type 2 diabetes D3_ZT_11:15':3}
        if (tn == 5):
            atypes = ['H1', 'H2', 'D1', 'D2']
            ahash = {'Type 2 diabetes D2_ZT_15:30':2,
                    'Type 2 diabetes D3_ZT_0:15':3,
                    'Healthy D2_ZT_15:30':0,
                    'Healthy D3_ZT_0:15':1}
        if (tn == 6):
            atypes = ['H', 'D']
            ahash = {'Type 2 diabetes D2_ZT_15:30':1,
                    'Type 2 diabetes D3_ZT_0:15':1,
                    'Healthy D2_ZT_15:30':0,
                    'Healthy D3_ZT_0:15':0}
        self.initData(atype, atypes, ahash)

    def getWu2007IS(self, tn=1):
        self.prepareDataDf("MAC112")
        atype = self.getSurvName("c status")
        atypes = ['IS', 'IR', 'D']
        ahash = {'insulin sensitive':0, 'diabetic':2, 'insulin resistant':1}
        if (tn == 2):
            atypes = ['IS', 'IR']
            ahash = {'insulin sensitive':0, 'insulin resistant':1}
        if (tn == 3):
            atypes = ['IS', 'D']
            ahash = {'insulin sensitive':0, 'diabetic':1}
        self.initData(atype, atypes, ahash)

    def getDuPlessis2015(self, tn=1):
        self.prepareDataDf("MAC111")
        atype = self.getSurvName("c src1")
        atypes = ['S1', 'V1', 'S2', 'V2', 'S3', 'V3', 'S4', 'V4']
        ahash = {'Subc Fat, Histology class 1':0,
                'Visceral Fat, Histology class 1':1,
                'Subc Fat, Histology class 2':2,
                'Visceral Fat, Histology class 2':3,
                'Subc Fat, Histology class 3':4,
                'Visceral Fat, Histology class 3':5,
                'Subc Fat, Histology class 4':6,
                'Visceral Fat, Histology class 4':7}
        if (tn == 2):
            atypes = ['N', 'F']
            ahash = {'Subc Fat, Histology class 1':0,
                    'Visceral Fat, Histology class 1':0,
                    'Subc Fat, Histology class 2':0,
                    'Visceral Fat, Histology class 2':0,
                    'Subc Fat, Histology class 3':0,
                    'Visceral Fat, Histology class 3':0,
                    'Subc Fat, Histology class 4':1,
                    'Visceral Fat, Histology class 4':1}
        self.initData(atype, atypes, ahash)

    def getWatson2017(self, tn=1):
        self.prepareDataDf("MAC120")
        atype = self.getSurvName("c sleep duration")
        atypes = ['long', 'short']
        if (tn == 2):
            atypes = ['short', 'long']
        ahash = {}
        self.initData(atype, atypes, ahash)

    def getUyhelji2018(self, tn=1):
        self.prepareDataDf("MAC119")
        atype = self.getSurvName("c subject group")
        atypes = ['Sleep Deprived', 'Control']
        ahash = {}
        self.initData(atype, atypes, ahash)

    def getMaret2007(self, tn=1):
        self.prepareDataDf("MAC121")
        atype = self.getSurvName("c src1")
        strain = [re.split(", *", str(i))[0] for i in atype]
        expt = [re.split(", *", str(i))[1] if len(re.split(", *", str(i))) > 1
                else None for i in atype]
        time = [re.split(", *", str(i))[2] if len(re.split(", *", str(i))) > 2
                else None for i in atype]
        rep = [re.split(", *", str(i))[3] if len(re.split(", *", str(i))) > 3
                else None for i in atype]
        tissue = [re.split(", *", str(i))[4] if len(re.split(", *", str(i))) >
                4 else None for i in atype]
        atype = expt
        atypes = ['C', 'SD']
        ahash = {'sleep deprived':1,
                'control':0,
                '6 hrs sleep deprivation':1,
                '6hrs sleep deprivation':1}
        if (tn == 2):
            atype = [ ",".join([str(k[i]) for k in [strain, expt, time,
                tissue]]) for i in range(len(atype))]
            atypes = ['C', 'SD']
            ahash = {'C57BL/6J,sleep deprived,time of sacrifice ZT 0,None':1,
                    'C57BL/6J,control,time of sacrifice ZT 0,None':0}
        if (tn == 3):
            atype = [ ",".join([str(k[i]) for k in [strain, expt, time,
                tissue]]) for i in range(len(atype))]
            atypes = ['C', 'SD']
            ahash = {'AKR/J,control,time of sacrifice ZT 0,None':0,
                    'AKR/J,sleep deprived,time of sacrifice ZT 0,None':1}
        self.initData(atype, atypes, ahash)

    def getResuehr2019(self, tn=1):
        self.prepareDataDf("MAC123")
        atype = self.getSurvName("c work shift")
        atypes = ['D', 'N']
        ahash = {'Day-Shift':0, 'Night-Shift':1}
        if (tn == 2):
            shift = self.getSurvName("c work shift")
            time = self.getSurvName("c time of sample")
            atype = [time[i] if shift[i] == 'Day-Shift' else None
                    for i in range(len(time))]
            atypes = ['9', '12', '15', '18', '21', '24', '27', '30']
            ahash = {}
        if (tn == 3):
            shift = self.getSurvName("c work shift")
            time = self.getSurvName("c time of sample")
            atype = [time[i] if shift[i] == 'Night-Shift' else None
                    for i in range(len(time))]
            atypes = ['9', '12', '15', '18', '21', '24', '27', '30']
            ahash = {}
        self.initData(atype, atypes, ahash)

    def getDAmore2018(self, tn=1):
        self.prepareDataDf("MAC114")
        atype = self.getSurvName("c disease state")
        atypes = ['C', 'MetS']
        ahash = {'control':0, 'Metabolic Syndrome':1}
        self.initData(atype, atypes, ahash)

    def getChristou2019(self, tn=1):
        self.prepareDataDf("GL29")
        atype = self.getSurvName('c circadian phase')
        atypes = ['D1', 'D2', 'D3', 'N1', 'N2', 'N3']
        v1 = [-4, 0, 6, 12, 18]
        ahash = {}
        for i in self.aRange():
            t = atype[i]
            if t is None:
                continue
            ahash[t] = np.searchsorted(v1, float(t), 'left')
        if (tn == 2):
            atypes = ['N', 'D']
            atype = [ahash[i] if i in ahash else None for i in atype]
            ahash = { 1: 0, 2: 0, 3:1 }
        if (tn == 3):
            atypes = ['D', 'N']
            atype = [ahash[i] if i in ahash else None for i in atype]
            ahash = { 1: 1, 2: 1, 3:0 }
        self.initData(atype, atypes, ahash)

    def getWagstaffe2020(self, tn=1):
        self.prepareDataDf("MACV121")
        atype = self.getSurvName("c Title")
        atype = [re.sub("^._", "", str(k)) for k in atype]
        atypes = ['U-', 'U+', 'E-', 'E+']
        ahash = {'Med_CD14-':0, 'EBOV_CD14+':3, 'Med_CD14+':1, 'EBOV_CD14-':2}
        if (tn == 2):
            atypes = ['U', 'E']
            ahash = {'Med_CD14-':0, 'EBOV_CD14+':1, 'Med_CD14+':0, 'EBOV_CD14-':1}
        if (tn == 3):
            atypes = ['U', 'E']
            ahash = {'Med_CD14-':0, 'EBOV_CD14-':1}
        self.initData(atype, atypes, ahash)

    def getPrice2020(self, tn=1):
        self.prepareDataDf("MACV111")
        mtype = self.getSurvName("c tissue")
        atype = self.getSurvName("c infection condtion")
        atype = [ " ".join([str(atype[i]), str(mtype[i])]) for i in
                range(len(atype))]
        atypes = ['S', 'SI', 'L', 'LI'];
        ahash = {'Mock Spleen':0, 'Mock Liver':2,
                'MA-EBOV infected Spleen':1, 'MA-EBOV infected Liver':3}
        if (tn == 2):
            atypes = ['C', 'I'];
            ahash = {'Mock Liver':0, 'MA-EBOV infected Liver':1}
        self.initData(atype, atypes, ahash)

    def getReynard2019(self, tn=1):
        self.prepareDataDf("MACV112")
        atype = self.getSurvName("c group")
        atypes = ['HC', 'SR', 'VR', 'F'];
        ahash = {'Fatalities':3, 'Healthy controls':0,
                'Survivors in recovery phase':1, 'Viremic survivors':2}
        if (tn == 2):
            atypes = ['HC', 'R', 'I'];
            ahash = {'Fatalities':2, 'Healthy controls':0,
                    'Survivors in recovery phase':1, 'Viremic survivors':2}
        if (tn == 3):
            atypes = ['C', 'I'];
            ahash = {'Fatalities':1, 'Healthy controls':0,
                    'Survivors in recovery phase':0, 'Viremic survivors':1}
        self.initData(atype, atypes, ahash)

    def getCameron2007(self, tn=1):
        self.prepareDataDf("MACV109")
        atype = self.getSurvName("c Status");
        atypes = ['HC', 'C', 'Pre', 'Post']
        ahash = {'pre-pO2 nadir':2, 'post-pO2 nadir':3, 'healthy control':0,
                'convalescent':1}
        if (tn == 2):
            atypes = ['HC', 'I']
            ahash = {'pre-pO2 nadir':1, 'post-pO2 nadir':1, 'convalescent':1,
                    'healthy control':0}
        self.initData(atype, atypes, ahash)

    def getMitchell2013(self, tn=1):
        self.prepareDataDf("MACV104")
        time = self.getSurvName("c timepoint")
        atype = [re.sub("h.*", "", str(k)) for k in time]
        ahash = {'0':0, '18':3, '12':2, '6':1}
        tval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c infection code")
        atypes = ['C', 'I']
        ahash = {'BatSRBD':1, 'icSARS':1, 'Mock':0, 'dORF6':1, 'H1N1':1, 'mock':0}
        aval = [ahash[i] if i in ahash else None for i in atype]
        if (tn == 2):
            atype = self.getSurvName("c infection code")
            ahash = {'BatSRBD':1, 'icSARS':1, 'Mock':0, 'mock':0}
            aval = [ahash[i] if i in ahash else None for i in atype]
            atype = ['C' if aval[i] == 1 and tval[i] == 0 else atype[i]
                    for i in range(len(atype))]
        if (tn == 3):
            atype = self.getSurvName("c infection code")
            ahash = {'H1N1':1, 'Mock':0, 'mock':0}
            aval = [ahash[i] if i in ahash else None for i in atype]
            atype = ['C' if aval[i] == 1 and tval[i] == 0 else atype[i]
                    for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getGuan2018(self, tn=1):
        self.prepareDataDf("MACV117")
        atype = self.getSurvName("c Title")
        atype = [re.sub(", .*", "", str(k)) for k in atype]
        atype = [ re.split(" ", str(k))[1] if len(re.split(" ", str(k))) > 1
                                else None for k in atype]
        ahash = {'2':2, 'control':0, '3':3, '8':8, '7':7, '4':4,
                '9':9, '6':6, '10':10}
        pval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c number of day post-infection")
        atype = [re.sub("NA", "0", str(k)) for k in atype]
        phash = {}
        for i in range(2, len(atype)):
            if pval[i] not in phash:
                phash[pval[i]] = []
            phash[pval[i]].append([i, int(atype[i])])
        before = []
        after = []
        for i in phash.keys():
            if i == 0:
                before += [k[0] for k in phash[i]]
                after += [k[0] for k in phash[i]]
            else:
                ll = sorted(phash[i], key = lambda k: k[1])
                before.append(ll[0][0])
                before.append(ll[1][0])
                after.append(ll[-1][0])
        beforehash = set(before)
        afterhash = set(after)
        ahash = {'2':0, '3':1, '8':1, '7':0, '4':1,
                '9':0, '6':0, '10':1}
        gender = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c subject category")
        atypes = ['C', 'P']
        ahash = {'Patient':1, 'Control':0}
        if (tn == 2):
            atype = [atype[i] if gender[i] == 0 or pval[i] == 0
                    else None for i in range(len(atype))]
        if (tn == 3):
            #Categories: 1) Control; 2) Mild - non-invasive ventilation (#4,#5);
            #3) Moderate- MV, but discharged within 2 mo (#1, 2, 7, 8);
            #4) Severe-  MV+prolonged hospitalization (#6 and 9);
            #5) Death after MV + ECMO = (#3 and 10)
            atype = pval
            atypes = ['C', 'Mi', 'Mo', 'S', 'D']
            ahash = {0:0, 4:1, 5:1, 1:2, 2:2, 7:2, 8:2, 6:3, 9:3, 3:4, 10:4}
            #atype = [atype[i] if gender[i] == 1 or pval[i] == 0
            #        else None for i in range(len(atype))]
        if (tn == 4):
            atype = pval
            atypes = ['C', 'Mi', 'MV', 'D']
            ahash = {0:0, 4:1, 5:1, 1:2, 2:2, 7:2, 8:2, 6:2, 9:2, 3:3, 10:3}
            atype = [atype[i] if i in beforehash
                    else None for i in range(len(atype))]
        if (tn == 5):
            atype = pval
            atypes = ['C', 'Mi', 'Mo', 'S', 'D']
            ahash = {0:0, 4:1, 5:1, 1:2, 2:2, 7:2, 8:2, 6:3, 9:3, 3:4, 10:4}
            atype = [atype[i] if i in afterhash
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getZhai2015(self, tn=1):
        self.prepareDataDf("MACV103")
        atype = self.getSurvName("c time point")
        atypes = ['C', 'I']
        ahash = {'Baseline':0, 'Spring':0, 'Day0':1, 'Day6':1, 'Day21':1,
                'Day2':1, 'Day4':1}
        if (tn == 2):
            ahash = {'Baseline':0, 'Day0':1}
        if (tn == 3):
            atypes = ['CV', 'AV']
            ahash = {'Day0':1, 'Day21':0}
        self.initData(atype, atypes, ahash)

    def getJosset2014(self, tn=1):
        self.prepareDataDf("MACV110")
        atype = self.getSurvName("c virus");
        atypes = ['C', 'InfA', 'CoV'];
        ahash = {'MA15':2, 'MOCK':0, 'PR8':1}
        if (tn == 2):
            atypes = ['C', 'InfA']
            ahash = {'MOCK':0, 'PR8':1}
        if (tn == 3):
            atypes = ['C', 'CoV'];
            ahash = {'MA15':1, 'MOCK':0}
        self.initData(atype, atypes, ahash)

    def getJones2019(self, tn=1):
        self.prepareDataDf("MACV107")
        atype = self.getSurvName("c visit")
        ahash = {'AV':0, 'CV':1}
        gval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c src1")
        ahash = {'NMS':0, 'PBMC':1}
        rval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c virus positive at av (1=yes, 0=no, 9=not measured)")
        atypes = ['0', '1']
        ahash = {}
        if (tn == 2):
            atype = self.getSurvName("c human coronavirus at av (1=yes, 0=no, 9=not measured)")
            atype = [atype[i] if rval[i] == 0 and gval[i] == 0
                    else None for i in range(len(atype))]
        if (tn >= 3):
            atype = self.getSurvName("c visit")
            ahash = {'AV':1, 'CV':0}
            atypes = ['CV', 'AV']
        if (tn == 4):
            atype = [atype[i] if rval[i] == 0
                    else None for i in range(len(atype))]
        if (tn == 5):
            atype = [atype[i] if rval[i] == 1
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getDunning2018(self, tn=1):
        self.prepareDataDf("MACV113")
        atype = self.getSurvName("c t1severity")
        atypes = ['HC', '1', '2', '3'];
        ahash = {'HC':0, '1':1, '2':2, '3':3}
        aval = [ahash[i] if i in ahash else None for i in atype]
        if (tn == 2):
            atypes = ['HC', 'I']
            ahash = {'1':1, '2':1, '3':1}
        if (tn == 3):
            atype = self.getSurvName("c ethnicity")
            atypes = ['W', 'B', 'A', 'O']
            ahash = {'White':0, 'Other':3, 'Black':1, 'Asian':2}
            atype = [atype[i] if aval[i] == 0
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getMorgun2006I(self, tn=1):
        self.prepareDataDf("MACV89")
        atype = self.getSurvName("c Final Clinical diagnosis ")
        atypes = ['N', 'Ch', 'Pre-Ch', 'R', 'Pre-R', 'Tox']
        ahash = {'toxoplasma myocarditis':5}
        if (tn == 2):
            atypes = ['S', 'R']
            ahash = {'R':1, 'N':0}
        self.initData(atype, atypes, ahash)

    def getMorgun2006II(self, tn=1):
        self.prepareDataDf("MACV90")
        atype = self.getSurvName("c Final Clinical diagnosis")
        atypes = ['N', 'R', 'Pre-R']
        ahash = {}
        if (tn == 2):
            atypes = ['S', 'R']
            ahash = {'R':1, 'N':0}
        self.initData(atype, atypes, ahash)

    def getVanLoon2019(self, tn=1):
        self.prepareDataDf("MACV88")
        atype = self.getSurvName("c tissue")
        ahash = {'BLOOD':0, 'kidney allograft biopsy':1}
        rval = [ahash[i] if i in ahash else None for i in atype]
        tcmr = self.getSurvName("c tcmr (no:0_borderline:1_TCMR:2)")
        abmr = self.getSurvName("c abmr (no:0_Yes:1)")
        atype = [str(tcmr[i]) + " " + str(abmr[i]) for i in range(len(tcmr))]
        atypes = ['S', 'R']
        ahash = {'2 0':1, '0 1':1, '0 0':0, '1 0':1, '2 1':1, '1 1':1}
        if (tn == 2):
            atype = [atype[i] if rval[i] == 0 \
                    else None for i in range(len(atype))]
        if (tn == 3):
            atype = [atype[i] if rval[i] == 1 \
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getPalmer2019(self, tn=1, tb=0):
        self.prepareDataDf("MACV50")
        atype = self.getSurvName('c src1')
        ahash = {'colon':0, 'blood':1}
        tval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName('c disease type')
        atypes = ['C', 'IBD', 'UC', 'CD']
        ahash = {'Control':0,
                'Control - infect':0,
                'Control - celiac':0,
                'Control (E.coli) --> CD':0,
                'Control--> CD':0,
                'IBD w/ oral':1,
                'IBDU':1,
                'UC - pancolitis':2,
                'UC - L colitis':2,
                'UC - proctitis':2,
                "Crohn's Disease":3}
        if (tn == 2):
            atypes = ['C', 'IBD']
            ahash = {'Control':0,
                    'Control - infect':0,
                    'Control - celiac':0,
                    'Control (E.coli) --> CD':0,
                    'Control--> CD':0,
                    'UC - pancolitis':1,
                    'UC - L colitis':1,
                    'UC - proctitis':1,
                    "Crohn's Disease":1}
            atype = [atype[i] if tval[i] == tb
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getGharib2019Alv(self, tn=1, ri=0):
        self.prepareDataDf("MAC16")
        atype = self.getSurvName("c time point")
        atypes = ['D1', 'D4', 'D8']
        ahash = {'Day 1':0, 'Day 4':1, 'Day 8':2}
        rval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("n ventilator-free days (vfd)")
        atypes = ['VFD-Extubated/Alive', 'VFD-Intubated/Dead']
        ahash = {'0':1, '7':0, '18':0, '19':0, '21':0,
                '22':0, '23':0, '24':0, '25':0}
        atype = [atype[i] if rval[i] == ri else None for i in range(len(atype))]
        if (tn == 2):
            atypes = ['VFD-Extubated/Alive', 'VFD-Intubated/Dead']
            ahash = {'0':1, '7':1, '18':1, '19':0, '21':0,
                    '22':0, '23':0, '24':0, '25':0}
        self.initData(atype, atypes, ahash)

    def getSimpson2019(self, tn=1):
        self.prepareDataDf("MACV45")
        atype = self.getSurvName('c patient outcome')
        atypes = ['C', 'I', 'F']
        ahash = {'control':0,
                'acute liver failure':2,
                'acute liver injury':1}
        if (tn == 2):
            atype = self.getSurvName('c survival')
            atypes = ['S', 'D']
            ahash = {'spontaneously survived':0, 'dead or transplanted':1}
        self.initData(atype, atypes, ahash)

    def getBohne2012II(self, tn=1):
        self.prepareDataDf("MACV171")
        atype = self.getSurvName("c liver sample group")
        atypes = ['Cont', 'Non-TOL', 'HEPC', 'Non-TOL REJ', 'TOL', 'Cont-Tx', 'REJ']
        ahash = {}
        if (tn == 2):
            atypes = ['Cont-Tx', 'REJ']
        if (tn == 3):
            atypes = ['TOL', 'REJ']
        self.initData(atype, atypes, ahash)

    def getArijs2009uc(self, tn=1):
        self.prepareDataDf("PLP142")
        atype = self.getSurvName("c WK8RSPHM")
        atypes = ["R", "NR"]
        ahash = {"Yes": 0, "No": 1}
        self.initData(atype, atypes, ahash)

    def getArijs2018(self, tn=1):
        self.prepareDataDf("PLP10")
        atype = self.getSurvName("c Response")
        atype = [str(i).strip() for i in atype]
        atypes = ['Control', 'UC R', 'UC NR', 'Active UC', 'UC other']
        ahash = {}
        if (tn == 2):
            atypes = ['UC R', 'UC NR']
        if (tn == 3):
            pid = self.getSurvName("c study individual number")
            res = self.getSurvName("c Response")
            res = [str(i).strip() for i in res]
            phash = {}
            for i in range(len(atype)):
                if res[i] == 'UC R':
                    phash[pid[i]] = 'R'
                if res[i] == 'UC NR':
                    phash[pid[i]] = 'NR'
            time = self.getSurvName("c week (w)")
            atype = [phash[pid[i]] if pid[i] in phash and time[i] == 'W0'
                    else None for i in range(len(atype))]
            atypes = ['R', 'NR']
            ahash = {}
        self.initData(atype, atypes, ahash)

    def getArijs2009(self, tn=1):
        self.prepareDataDf("PLP27")
        atype = self.getSurvName("c tissue")
        ahash = {'Colon':0, 'Ileum':1}
        tissue = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c before or after first infliximab treatment")
        atype = [re.sub(" .*", "", str(k)) for k in atype]
        ahash = {'Before':1, 'After':2, 'Not':0}
        treatment = [ahash[i] if i in ahash else None for i in atype]
        response = self.getSurvName("c response to infliximab")
        atype = self.getSurvName("c disease")
        atypes = ['Control', 'UC', 'CD']
        ahash = {}
        if (tn == 2):
            atypes = ["R", "NR"]
            ahash = {"Yes": 0, "No": 1}
            atype = response
            atype = [atype[i] if tissue[i] == 0 and treatment[i] == 1
                    else None for i in range(len(atype))]
        if (tn == 3):
            atypes = ["R", "NR"]
            ahash = {"Yes": 0, "No": 1}
            atype = response
            atype = [atype[i] if tissue[i] == 0 and treatment[i] == 2
                    else None for i in range(len(atype))]
        self.initData(atype, atypes, ahash)

    def getVerstockt2019(self):
        self.prepareDataDf("PLP73")
        atype = self.getSurvName("c clinical history")
        atypes = ['R', 'NR']
        ahash = {'responder':0, 'non-responder':1}
        self.initData(atype, atypes, ahash)

    def getTang2020(self, tn=1):
        self.prepareDataDf("MACV95")
        atype = self.getSurvName("c macrophage phenotype")
        ahash = {'Tissue resident, F4/80hi CD206-':1,
                'Monocyte-derived, F4/80int CD206+':0}
        mval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName("c Title")
        strain = [str(k).split("_")[0] if len(str(k).split("_")) > 0 else None
                for k in atype]
        ahash = {'B6':1, 'Balbc':0}
        rval = [ahash[i] if i in ahash else None for i in strain]
        atype = [str(k).split("_")[1] if len(str(k).split("_")) > 1 else None
                for k in atype]
        atypes = ['IL4c', 'ThioIL4c']
        ahash = {}
        if (tn == 2):
            atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
        if (tn == 3):
            atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
        if (tn == 4):
            atype = strain
            atype = [atype[i] if mval[i] == 0 else None for i in range(len(atype))]
            atypes = ['Balbc', 'B6']
            ahash = {}
        if (tn == 5):
            atype = strain
            atype = [atype[i] if mval[i] == 1 else None for i in range(len(atype))]
            atypes = ['Balbc', 'B6']
            ahash = {}
        self.initData(atype, atypes, ahash)

    def getLink2018(self, tn=1):
        self.prepareDataDf("MACV98")
        atype = self.getSurvName("c Title")
        atype = [re.sub("_rep.*", "", str(k)) for k in atype]
        group = [re.sub("BMDM_RNA_polyA_", "", str(k)) for k in atype]
        strain = self.getSurvName("c strain")
        ahash = {'C57BL/6J':1, 'SPRET/EiJ':2, 'BALB/cJ':0, 'NOD/ShiLtJ':3, 'PWK/PhJ':4}
        rval = [ahash[i] if i in ahash else None for i in strain]
        atype = self.getSurvName("c ligands in culture")
        atypes = ['M0', 'M1']
        ahash = {'no treatment':0, 'KLA 6h':1, 'KLA 1h':1}
        if (tn == 2):
            atype = group
            ahash = {'BALB_notx':0, 'BALB_KLA_1h':1}
        if (tn == 3):
            atype = group
            ahash = {'C57_notx_6h':0, 'C57_KLA_6h':1}
        if (tn == 4):
            atype = group
            atypes = ['Bl6', 'Bl6t', 'Bc', 'Bct']
            ahash = {'C57_notx_6h':0, 'C57_KLA_6h':1,
                    'BALB_notx':2, 'BALB_KLA_1h':3}
        if (tn == 5):
            atype = group
            atypes = ['Bc', 'Bl6']
            ahash = {'C57_notx_6h':1, 'C57_KLA_6h':1,
                    'BALB_notx':0, 'BALB_KLA_1h':0}
        self.initData(atype, atypes, ahash)

    def getHowes2016(self, tn=1):
        self.prepareDataDf("MACV101")
        strain = self.getSurvName("c strain")
        ahash = {'C57BL/6':1, 'BALB/c':0}
        rval = [ahash[i] if i in ahash else None for i in strain]
        gtype = self.getSurvName("c genotype/variation")
        ahash = {'WT':0, 'IFNabRKO':1}
        gval = [ahash[i] if i in ahash else None for i in gtype]
        atype = self.getSurvName("c treatment")
        atypes = ['M0', 'M1']
        ahash = {'HkBps':1, 'media':0}
        aval = [ahash[i] if i in ahash else None for i in atype]
        if (tn == 2):
            atype = [atype[i] if rval[i] == 0 else None for i in range(len(atype))]
        if (tn == 3):
            atype = [atype[i] if rval[i] == 1 else None for i in range(len(atype))]
        if (tn == 4):
            atype = strain
            atypes = ['Bc', 'B6']
            ahash = {'C57BL/6':1, 'BALB/c':0}
        if (tn == 5):
            atype = strain
            atype = [atype[i] if aval[i] == 0 and gval[i] == 0 \
                    else None for i in range(len(atype))]
            atypes = ['Bc', 'B6']
            ahash = {'C57BL/6':1, 'BALB/c':0}
        self.initData(atype, atypes, ahash)

    def getSurvival(self, dbid = "CRC35.3"):
        self.prepareDataDf(dbid)
        atype = self.getSurvName("status")
        atypes = ['Censor', 'Relapse']
        ahash = {"0": 0, "1":1}
        self.initData(atype, atypes, ahash)

    def getJSTOM(self):
        self.getSurvival("CRC35.3")

    def getBos(self):
        self.getSurvival("BC20")

    def getLee2012(self, tn=1):
        self.prepareDataDf("MACV16")
        atype = self.getSurvName('c aim2 expression')
        atypes = ['Wt', 'Mut']
        ahash = {'persistent expression':1, 'absent expression':0}
        self.initData(atype, atypes, ahash)

    def getOakes2017Hs(self, tn=1):
        self.prepareDataDf("MACV19")
        atype = self.getSurvName('c src1')
        atypes = ['W-D', 'W+D', 'M-D', 'M+D']
        ahash = {'T47D cell line, WT mouse Oas2, -DOX':0,
                'T47D cell line, MUT mouse Oas2, +DOX':3,
                 'T47D cell line, WT mouse Oas2, +DOX':1,
                 'T47D cell line, MUT mouse Oas2, -DOX':2}
        if (tn == 2):
            atypes = ['W-D', 'W+D']
            ahash = {'T47D cell line, WT mouse Oas2, -DOX':0,
                     'T47D cell line, WT mouse Oas2, +DOX':1}
        self.initData(atype, atypes, ahash)

    def getWang2019Mac(self, tn=1):
        self.prepareDataDf("MACV26")
        atype = self.getSurvName('c treatment')
        atypes = ['None', 'IL-15', 'IL-4', 'Media']
        ahash = {}
        if (tn == 2):
            atypes = ['None', 'IL-15']
        self.initData(atype, atypes, ahash)

    def getMan2015(self, tn=1):
        self.prepareDataDf("MACV4")
        atype = self.getSurvName('c genotype/variation')
        atypes = ['Wt', 'Irf1', 'Aim2', 'Ifnar1']
        ahash = {'wild-type':0, 'Irf1-/-':1, 'Aim2-/-':2, 'Ifnar1-/-':3}
        if (tn == 2):
            atypes = ['Wt', '', 'Irf1']
            ahash = {'wild-type':0, 'Irf1-/-':2}
        self.initData(atype, atypes, ahash)

    def getFensterl2012(self, tn=1):
        self.prepareDataDf("MACV14")
        atype = self.getSurvName('c Title')
        atypes = ['WV6', 'IV6', 'WV2', 'IV2']
        ahash = {'wtbrain-VSV-6d-rep1':0,
                'Ifit2KObrain-VSV-6d-rep1':1,
                'wtbrain-VSV-2d-rep1':2,
                'Ifit2KObrain-VSV-2d-rep1':3}
        if (tn == 2):
            atypes = ['W', '', 'Ifit2']
            ahash = {'wtbrain-VSV-2d-rep1':0,
                    'Ifit2KObrain-VSV-2d-rep1':2}
        self.initData(atype, atypes, ahash)

    def getIrey2019Hs(self, tn=1):
        self.prepareDataDf("MACV18")
        atype = self.getSurvName('c Title')
        atype = [ re.split(",", str(i))[0] for i in atype]
        atypes = ['V', '7', '231']
        ahash = {'MCF7-conditioned-medium':1,
                'vehicle-conditioned-medium':0,
                'MDA-MB-231-conditioned-medium':2}
        media = [ahash[i] if i in ahash else None for i in atype]
        if (tn == 2):
            atype = self.getSurvName('c Title')
            atype = [ re.split(", ", str(atype[i]))[1] \
                    if len(re.split(", ", str(atype[i]))) > 1 and \
                    media[i] != 2 \
                    else "" for i in range(len(atype))]
            atypes = ['V', 'I', 'I']
            ahash = {'vehicle':0, 'ruxolitinib':2}
        self.initData(atype, atypes, ahash)

    def getIrey2019Mm(self, tn=1):
        self.prepareDataDf("MACV18.2")
        atype = self.getSurvName('c phenotype')
        atypes = ['Wt', '', 'STAT3']
        ahash = {'STAT3 wild type':0, 'STAT3 knockout':2}
        self.initData(atype, atypes, ahash)

    def getOakes2017Mm(self, tn=1):
        self.prepareDataDf("MACV19.2")
        atype = self.getSurvName('c src1')
        atypes = ['W2', 'W18', 'M2', 'M18']
        ahash = {'Mammary gland, MUT Oas2, 2dpp':2,
                'Mammary gland, WT Oas2, 2dpp':0,
                'Mammary gland, WT Oas2, 18dpc':1,
                'Mammary gland, MUT Oas2, 18dpc':3}
        if (tn == 2):
            atypes = ['W2', '', 'M2']
            ahash = {'Mammary gland, MUT Oas2, 2dpp':2,
                    'Mammary gland, WT Oas2, 2dpp':0}
        self.initData(atype, atypes, ahash)

    def getGoldmann2015(self, tn=1):
        self.prepareDataDf("MACV24")
        atype = self.getSurvName('c genotype/variation')
        atypes = ['Wt', 'Ifnar', 'Usp18', 'C61A', 'KD', 'NC', 'DKO']
        ahash = {'WT':0,
                'IFNARko':1,
                'Usp18_C61A':3,
                'Usp18 ko':2,
                'si RNA Usp18':4,
                'si RNA control':5,
                'USP18ko:IFNARko':6}
        if (tn == 2):
            atypes = ['Wt', '', 'DKO']
            ahash = {'WT':0,
                    'USP18ko:IFNARko':2}
        if (tn == 3):
            atypes = ['Ifnar', '', 'DKO']
            ahash = {'IFNARko':0,
                    'USP18ko:IFNARko':2}
        self.initData(atype, atypes, ahash)

    def getPG2019lpsRep(self, tn=1):
        self.prepareDataDf("MAC125")
        ttype = self.getSurvName("c type")
        mtype = self.getSurvName("c times")
        atype = [ str(ttype[i]) + " " + str(mtype[i]) for i in
                range(len(ttype))]
        atypes = ['KO 0', 'WT 0', 'WT 6hr', 'KO 6hr']
        ahash = {}
        if (tn == 2):
            atypes = ['WT 6hr', 'KO 6hr']
        self.initData(atype, atypes, ahash)

    def getLinke2017(self, tn=1):
        self.prepareDataDf("MACV20")
        atype = self.getSurvName('c genotype/variation')
        atypes = ['Wt', 'Mut']
        ahash = {'Tsc2fl/fl LysM+/+':0, 'Tsc2fl/fl LysM+/cre':1}
        self.initData(atype, atypes, ahash)

    def getLi2019(self, tn=1):
        self.prepareDataDf("MACV22")
        atype = self.getSurvName('c mouse genotype')
        atypes = ['Wt', 'Rnf5']
        ahash = {'WT':0, 'Rnf5 KO':1}
        self.initData(atype, atypes, ahash)

    def getUckelmann2020(self, tn=1):
        self.prepareDataDf("MACV27")
        rtype = self.getSurvName('c cell line')
        ahash = {'IMS-M2':0, 'OCI-AML-3':1}
        line = [ahash[i] if i in ahash else None for i in rtype]
        ttype = self.getSurvName('c timepoint')
        ahash = {'day 3':3, 'day 5':5, 'day 7':7}
        time = [ahash[i] if i in ahash else None for i in ttype]
        atype = self.getSurvName('c treatment')
        ahash = {'':0, '330nM VTP50469':1}
        treat = [ahash[i] if i in ahash else None for i in atype]
        atype = ["-".join([str(line[i]), str(treat[i]), str(time[i])])
                                    for i in range(len(atype))]
        atypes = ['N', 'T']
        ahash = {'0-0-3':0, '0-0-5':0, '0-0-7':0,
                '1-0-3':0, '1-0-5':0, '1-0-7':0,
                '0-1-3':1, '0-1-5':1, '0-1-7':1,
                '1-1-3':1, '1-1-5':1, '1-1-7':1}
        if (tn == 2):
            ahash = {'1-0-3':0, '1-0-5':0, '1-0-7':0, '1-1-5':1}
        if (tn == 3):
            ahash = {'0-0-3':0, '0-0-5':0, '0-0-7':0,
                    '0-1-3':1, '0-1-5':1, '0-1-7':1}
        if (tn == 4):
            ahash = {'1-0-3':0, '1-0-5':0, '1-0-7':0,
                    '1-1-3':1, '1-1-5':1, '1-1-7':1}
        self.initData(atype, atypes, ahash)

    def getUckelmann2020Mm(self, tn=1):
        self.prepareDataDf("MACV28")
        rtype = self.getSurvName('c cell type')
        rtype = [str(i).replace('mouse ', '').replace(' cells', '') \
                for i in rtype]
        ahash = {'':0, 'LSK':1, 'GMP':2, 'LSK-derived GMP':3,
                'GMP-derived GMP':4, 'LSK-derived LSK':5,
                'GMP-derived LSK':6, 'long-term GMP-derived GMP':7,
                'LSK-derived GMP-like':8, 'GMP-derived GMP-like':9}
        line = [ahash[i] if i in ahash else None for i in rtype]
        ttype = self.getSurvName('c timepoint')
        ahash = {'4 weeks':28, '':0, '9 months post transplant':270,
                'day 5':5, 'day 3':3}
        time = [ahash[i] if i in ahash else None for i in ttype]
        atype = self.getSurvName('c treatment')
        ahash = {'pIpC induction':1, '':0, '1% VTP50469 chow':2, 
                '30nM VTP50469':3}
        treat = [ahash[i] if i in ahash else None for i in atype]
        ctype = ["-".join([str(line[i]), str(treat[i]), str(time[i])])
                                    for i in range(len(atype))]
        atype = treat
        atypes = [0, 1, 2, 3]
        ahash = {}
        if (tn == 2):
            atype = [treat[i] if (line[i] == 0) else None
                    for i in range(len(atype))]
            atypes = [0, 3]
        if (tn == 3):
            atype = [treat[i] if (line[i] == 2) else None
                    for i in range(len(atype))]
            atypes = [0, 1, 2, 3]
        self.initData(atype, atypes, ahash)

    def getElDahr2019 (self, tn=1):
        self.prepareDataDf("MACV31")
        atype = self.getSurvName('c src1')
        atype = [re.sub("\s*[0-9]$", "", str(i)) for i in atype]
        atypes = ['W', '2', '21het', '21']
        ahash = {'WT':0,
                'EZH2-KO & EZH1-WT':1,
                'EZH2-KO & EZH1-hetKO':2,
                'EZH2-KO & EZH1-homoKO':3}
        if (tn == 2):
            atypes = ['W', 'M']
            ahash = {'WT':0,
                    'EZH2-KO & EZH1-WT':0,
                    'EZH2-KO & EZH1-hetKO':1,
                    'EZH2-KO & EZH1-homoKO':1}
        self.initData(atype, atypes, ahash)

    def getEncode2020 (self, tn=1):
        self.prepareDataDf("MACV32")
        atype = self.getSurvName('c Type')
        ahash = {'K562':1, 'HepG2':2, '':0}
        rval = [ahash[i] if i in ahash else None for i in atype]
        atype = self.getSurvName('c Target')
        atypes = ['W', 'M']
        ahash = {'':0, 'EEF2':1, 'NFX1':1, 'HMBOX1':1, 
                'HNRNPA1':1, 'NFYB':1, 'PCBP2':1}
        for k in atype:
            if k != '':
                ahash[k] = 1
        if (tn == 2):
            atype = [atype[i] if rval[i] == 1 else None
                    for i in range(len(atype))]
            ahash = {'':0, 'EEF2':1, 'NFX1':1, 'HMBOX1':1, 
                    'HNRNPA1':1, 'NFYB':1, 'PCBP2':1}
        self.initData(atype, atypes, ahash)

