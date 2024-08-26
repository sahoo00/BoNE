import bone
pd = bone.pd
re = bone.re
hu = bone.hu
np = bone.np

def getMSigDB(gs):
    url = "https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=" + gs + "&fileType=txt"
    df = pd.read_csv(url, sep="\t")
    df.columns.values[0] = 'ID'
    l1 = [list(df.ID[1:])]
    wt1 = [1]
    return wt1, l1
bone.getMSigDB = getMSigDB

def getYang2021(self, tn=1):
    self.prepareData("MACV412")
    atype = self.h.getSurvName("c treatment")
    atypes = ['C', 'LPS']
    ahash = {'saline':0, 'challenge':1}
    if tn == 2:
        bhash = {'GSM4291101', 'GSM4291098', 'GSM4291105', 'GSM4291100'}
        atype = [atype[i] if self.h.headers[i] not in bhash
                else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getYang2021 = getYang2021
def getSharma2022Mac(self, tn=1, ta=None):
    self.prepareData("MACV365")
    pdom = self.h.getSurvName('c Placental_Domain')
    infl = self.h.getSurvName('c Inflam_vasc_none')
    cell = self.h.getSurvName('c cell population')
    labor = self.h.getSurvName('c Preterm Labor')
    pterm = self.h.getSurvName('c Preterm_0no1yes')
    vpterm = self.h.getSurvName('c Very Preterm0no1yes')
    bpd = self.h.getSurvName('c BPD Classify')
    atype = cell
    atypes = ['cMNC', 'iMNC', 'ncMNC']
    ahash = {}
    if (tn == 2):
        atype = infl
        atypes = ['I', 'V']
        if ta == None:
            atype = [atype[i] if labor[i] == 'Yes'
                     else None for i in range(len(atype))]
        else:
            atype = [atype[i] if labor[i] == 'Yes' and cell[i] == ta
                     else None for i in range(len(atype))]
    if (tn == 3):
        atype = infl
        atypes = ['I', 'V']
        if ta == None:
            atype = [atype[i] if labor[i] == 'No'
                     else None for i in range(len(atype))]
        else:
            atype = [atype[i] if labor[i] == 'No' and cell[i] == ta
                     else None for i in range(len(atype))]
    if (tn == 4):
        atype = labor
        atypes = ['No', 'Yes']
    if (tn == 5):
        atype = bpd
        atypes = ['no', 'yes']
        atype = [atype[i] if pterm[i] == '1' and cell[i] == 'iMNC' and
                 labor[i] == 'No'
                 else None for i in range(len(atype))]
    if (tn == 6):
        atype = bpd
        atypes = ['no', 'yes']
        atype = [atype[i] if pterm[i] == '1' and cell[i] == 'cMNC' and
                 labor[i] == 'No'
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getSharma2022Mac = getSharma2022Mac
def getLance2007(self, tn=1):
    self.prepareData("LP3")
    treatment = self.h.getSurvName('c Treatment')
    place = self.h.getSurvName('c Collection            Place ')
    bpd = self.h.getSurvName('c BPD                  Yes/No')
    time = self.h.getSurvName('c Time           (Letter)')
    sid = self.h.getSurvName('c Original               Sample Name')
    sid = [re.sub("[A-Z]\s*$", "", str(k)) for k in sid]
    sid = [re.sub("\s", "", str(k)) for k in sid]
    shash = {'M01':'I', 'M07':'I', 'M10': 'V', 'M11':'I',
             'M13': 'I', 'M14': 'I', 'M17': 'I', 'M18': 'I',
             'GM02': 'I', 'GM03': 'I', 'GM08': 'V', 'GM14': 'V', 'GM16': 'V',
             'GM17': 'N', 'GM19': 'I', 'GM21': 'I', 'GM22': 'I', 'GM26': 'I'}
    atype = treatment
    atypes = ['CTRL', 'LPS']
    ahash = {}
    if (tn == 2):
        atype = place
        atypes = ['VAND', 'GM', 'GMR', 'UCSD']
    if (tn == 3):
        atype = [atype[i] if place[i] == 'GMR'
                 else None for i in range(len(atype))]
    if (tn == 4):
        atype = bpd
        atypes = ['No', 'Yes']
        atype = [atype[i] if time[i] == 'A'
                 else None for i in range(len(atype))]
        #atype = [atype[i] if time[i] != 'A'
        #         else None for i in range(len(atype))]
    if (tn == 5):
        atype = [treatment[i] if sid[i] in shash
                 else None for i in range(len(sid))]
        atypes = ['CTRL', 'LPS']
    if (tn == 6):
        atype = [shash[sid[i]] if sid[i] in shash
                 else None for i in range(len(sid))]
        atypes = ['I', 'V']
    if (tn == 7):
        atype = [shash[sid[i]] if sid[i] in shash
                 else None for i in range(len(sid))]
        atypes = ['I', 'V']
        atype = [atype[i] if time[i] != 'A' and time[i] != 'B'
                 else None for i in range(len(atype))]
    if (tn == 8):
        atype = [shash[sid[i]] if sid[i] in shash
                 else None for i in range(len(sid))]
        atypes = ['I', 'V']
        atype = [atype[i] if time[i] == 'A' or time[i] == 'B'
                 else None for i in range(len(atype))]
    if (tn == 9):
        atype = [shash[sid[i]] if sid[i] in shash
                 else None for i in range(len(sid))]
        atypes = ['I', 'V', 'N']
        atype = [atype[i] if time[i] != 'A'
                 else None for i in range(len(atype))]
    if (tn == 10):
        atype = ['AB' if time[i] == 'A' or time[i] == 'B'
                 else 'C' for i in range(len(sid))]
        atypes = ['AB', 'C']
    if (tn == 11):
        atype = bpd
        atypes = ['No', 'Yes']
        atype = [atype[i] if time[i] == 'A' or time[i] == 'B'
                 else None for i in range(len(atype))]
    if (tn == 12):
        atype = bpd
        atypes = ['No', 'Yes']
        atype = [atype[i] if time[i] != 'A' and time[i] != 'B'
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLance2007 = getLance2007

def getBell2023(self, tn=1, ta = 'CordBlood'):
    self.prepareData("MACV367")
    atype = self.h.getSurvName("c variable-disease")
    atypes = ['nonBPD', 'BPD']
    ahash = {}
    if (tn == 2):
        atype = self.h.getSurvName("c variable-time")
        atypes = ['CB', 'D14', 'D28']
        ahash = {'Day14':1, 'CordBlood':0, 'Day28':2}
    if (tn == 3):
        timepoint = self.h.getSurvName("c variable-time")
        atype = [atype[i] if timepoint[i] == ta
                 else None for i in range(len(atype))]
    if (tn == 4):
        o2 = self.h.getSurvName("c variable-o2 treatment")
        atype = [str(atype[i])+'-'+str(o2[i])for i in range(len(atype))]
        timepoint = self.h.getSurvName("c variable-time")
        atype = [atype[i] if timepoint[i] == ta
                 else None for i in range(len(atype))]
        atypes = ['nonBPD-no', 'BPD-yes', 'nonBPD-yes']
        ahash = {}
    if (tn == 5):
        o2 = self.h.getSurvName("c variable-o2 treatment")
        atype = [str(atype[i])+'-'+str(o2[i])for i in range(len(atype))]
        timepoint = self.h.getSurvName("c variable-time")
        atype = [atype[i] if timepoint[i] == ta
                 else None for i in range(len(atype))]
        atypes = ['nonBPD-no', 'BPD-yes']
        ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getBell2023 = getBell2023
def getWindhorst2023(self, tn=1):
    self.prepareData("MACV368")
    atype = self.h.getSurvName("c disease state")
    atypes = ['No', 'Mild', 'M/S']
    ahash = {'No BPD':0, 'Mild BPD':1, 'Moderate/severe BPD':2}
    if (tn == 2):
        weight = [0,0] + ana.h.getSurvName("c birth weight")[2:]
        atypes = ['No', 'M/S']
        ahash = {'No BPD':0, 'Moderate/severe BPD':1}
        atype = [atype[i] if ((atype[i] == 'No BPD') and (float(weight[i]) < 1300)) or
                 ((atype[i] == 'Moderate/severe BPD') and (float(weight[i]) < 1000))
                 else None for i in range(len(atype))]
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWindhorst2023 = getWindhorst2023
def getWang2022(self, tn=1):
    self.prepareData("MACV369")
    atype = self.h.getSurvName("c disease status")
    atypes = ['nonBPD', 'BPD']
    ahash = {}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getWang2022 = getWang2022
def getRyan2019(self, tn=1):
    self.prepareData("MACV370")
    atype = self.h.getSurvName("c src1")
    atypes = ['NC', 'BC', 'ND', 'BD']
    ahash = {'No_BPD_DHA_peripheral blood':2,
             'BPD_DHA_peripheral blood':3,
             'No_BPD_Control_peripheral blood':0,
             'BPD_Control_peripheral blood':1}
    if (tn == 2):
        atypes = ['noBPD', 'BPD']
        ahash = {'No_BPD_Control_peripheral blood':0,
                 'BPD_Control_peripheral blood':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getRyan2019 = getRyan2019
def getCheng2018(self, tn=1):
    self.prepareData("MACV372")
    atype = self.h.getSurvName("c diagnosis")
    atypes = ['C', 'BPD']
    ahash = {'normal':0, 'Bronchopulmonary dysplasia (BPD)':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getCheng2018 = getCheng2018

def getEldredge2024(self, tn=1, ta=0):
    self.prepareData("MACV506")
    atype = self.h.getSurvName("c cell type")
    atypes = ['E', 'M']
    ahash = {'cd14+cd16- monocytes':1, 'cd14+cd16+ monocytes':1,
              'epithelial':0}
    if tn == 2:
        cell = self.h.getSurvName("c cell type")
        ahash = {'cd14+cd16- monocytes':1, 'cd14+cd16+ monocytes':2,
                  'epithelial':0}
        ctype = [ahash[i] if i in ahash else None for i in cell]
        atype = self.h.getSurvName("c severity")
        atypes = ['M', 'S']
        ahash = {'severe I':1, 'severe II':1, 'moderate':0, 'mild':0}
        atype = [atype[i] if ctype[i] == ta
                 else None for i in range(len(atype))]
    if tn == 3:
        atypes = ['E', 'M']
        ahash = {'cd14+cd16- monocytes':1,
                  'epithelial':0}
    if tn == 4:
        atypes = ['E', 'M']
        ahash = {'cd14+cd16+ monocytes':1,
                  'epithelial':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getEldredge2024 = getEldredge2024

def getZec2023mm(self, tn=1, ta=0):
    self.prepareData("MACV515")
    atype = self.h.getSurvName("c src1")
    atypes = ['M0', 'M1']
    ahash = {'wildtype':1, 'wildtype_mock':0}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getZec2023mm = getZec2023mm

def getLi2024mm(self, tn=1, ta=0):
    self.prepareData("MACV516")
    atype = self.h.getSurvName("c treatment")
    atypes = ['M0', 'M1']
    ahash = {'no':0, 'IL-4':1, 'IL-4/ARV825':1,
            'IL-4/DMSO':1, 'IL-4/JQ1':1, 'IL-4/ZL0420':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2024mm = getLi2024mm

def getLi2024mmII(self, tn=1, ta=0):
    self.prepareData("MACV516.2")
    atype = self.h.getSurvName("c treatment")
    atypes = ['M0', 'M2']
    ahash = {'none':0, 'IL-4/ARV825':1, 'IL-4/ZL0420':1,
            'IL-4':1, 'IL-4/DMSO':1, 'IL-4/JQ1':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2024mmII = getLi2024mmII

def getLi2024mmIII(self, tn=1, ta=0):
    self.prepareData("MACV516.3")
    atype = self.h.getSurvName("c Title")
    atype = [hu.re.sub("[-_].*", "", str(k)) for k in atype]
    atypes = ['M0', 'M1', 'M2']
    ahash = {}
    if tn == 2:
        atype = self.h.getSurvName("c Title")
        atypes = ['M0', 'M1', 'M2']
        ahash = {'M1-1':1, 'M1-2':1, 'M1-3':1,
                'M2-1':2, 'M2-2':2, 'M2-3':2,
                'M0-a':0, 'M0-b':0, 'M0-c':0}
    if tn == 3:
        atype = self.h.getSurvName("c Title")
        atypes = ['M0', 'M1', 'M2']
        ahash = {'M1-1':1, 'M1-2':1, 'M1-3':1,
                'M2-1':2, 'M2-2':2, 'M2-3':2}
    if tn == 4:
        atype = self.h.getSurvName("c Title")
        atypes = ['M1', 'M2']
        ahash = {'M1-1':0, 'M1-2':0, 'M1-3':0,
                'M2-1':1, 'M2-2':1, 'M2-3':1}
    self.initData(atype, atypes, ahash)
    return
bone.IBDAnalysis.getLi2024mmIII = getLi2024mmIII

