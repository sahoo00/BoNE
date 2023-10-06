# IMPORT STATEMENTS
import cv2
import re
import numpy as np
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import matplotlib.patches as patches
import matplotlib.colors as colors
from matplotlib.transforms import *
import PIL
import math
import pandas as pd
import seaborn as sns
import json
from sklearn.metrics import *
from scipy.stats import fisher_exact, ttest_ind
from collections import Counter
from pprint import pprint
import os
import pickle
import sys
sys.path.append("/booleanfs2/sahoo/Hegemon/")
sys.path.append("BoNE")
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
        
import bone
reload(bone)

class PolypAnalysis(bone.IBDAnalysis):

    def __init__(self):
        bone.IBDAnalysis.__init__(self)
        
    def getqPCR_WT_APC_Inf(self, tn=1):
        self.prepareData("PG14", cfile="/Users/dtv004/public_html/Hegemon/explore.conf")
        atype = self.h.getSurvName('c Description')
        atypes = ['FN_WT_Un', 'FN_WT_Inf']
        ahash = {'FN_WT_Un':0, 'FN_WT_Inf':1}
        
        if tn == 2:
            atypes = ['K12_APC_Un', 'K12_APC_Inf']
            ahash = {'K12_APC_Un':0, 'K12_APC_Inf':1}
            
        if tn == 3:
            atypes = ['LF82_APC_Un', 'LF82_APC_Inf']
            ahash = {'LF82_APC_Un':0, 'LF82_APC_Inf':1}
            
        if tn == 4:
            atypes = ['NC101_APC_Un', 'NC101_APC_Inf']
            ahash = {'NC101_APC_Un':0, 'NC101_APC_Inf':1}
            
        if tn == 5:
            atypes = ['BF_APC_Un', 'BF_APC_Inf']
            ahash = {'BF_APC_Un':0, 'BF_APC_Inf':1}
            
        if tn == 6:
            atypes = ['FN_APC_Un', 'FN_APC_Inf']
            ahash = {'FN_APC_Un':0, 'FN_APC_Inf':1}
            
        if tn == 7:
            atypes = ['WT_UN', 'WT_FN', 'APC_UN', 'APC_FN']
            ahash = {'FN_WT_Un':0, 'FN_WT_Inf':1, 'FN_APC_Un':2, 'FN_APC_Inf':3}
        
        self.initData(atype, atypes, ahash)
        return
        
        
    def getqPCR_Hs_tissue_NI_P(self, tn=1):
        self.prepareData("PG15", cfile="/Users/dtv004/public_html/Hegemon/explore.conf")
        atype = self.h.getSurvName('c Description')
        atypes = ['H', 'NI', 'P']
        ahash = {'Healthy':0, 
                 'FAP NI':1, 'Lynch NI':1, 'PJS NI':1, 'JPS NI':1,
                 'FAP P':2, 'JPS P':2}
        
        if tn == 2:
            atypes = ['NI', 'P']
            ahash = {'FAP NI':0, 'Lynch NI':0, 'PJS NI':0, 'JPS NI':0,
                     'FAP P':1, 'JPS P':1}
            
        if tn == 3:
            atypes = ['PJS NI', 'JPS NI', 'Lynch NI', 'FAP NI']
            ahash = {'FAP NI':3, 'Lynch NI':2, 'PJS NI':0, 'JPS NI':1}
            
        if tn == 4:
            atypes = ['JPS NI', 'JPS P']
            ahash = {'JPS NI':0, 'JPS P':1}
            
        if tn == 5:
            atypes = ['FAP NI', 'FAP P']
            ahash = {'FAP NI':0, 'FAP P':1}
            
        self.initData(atype, atypes, ahash)
        return
        
    def getqPCR_Hs_org_NI_P(self, tn=1):
        self.prepareData("PG16", cfile="/Users/dtv004/public_html/Hegemon/explore.conf")
        atype = self.h.getSurvName('c Description')
            
        if tn == 1: 
            atypes = ['JPS H', 'JPS NI', 'JPS P']
            ahash = {'JPS_H':0, 'JPS_NI':1, 'JPS_P':2}
        
        if tn == 2:
            atypes = ['PJS H', 'PJS NI', 'PJS P']
            ahash = {'PJS_H':0, 'PJS_NI':1, 'PJS_P':2}
            
        if tn == 3:
            atypes = ['Lynch H', 'Lynch NI']
            ahash = {'Lynch_H':0, 'Lynch_NI':1}
            
        if tn == 4:
            atypes = ['FAP1 H', 'FAP1 NI', 'FAP1 P']
            ahash = {'FAP1_H':0, 'FAP1_NI':1, 'FAP1_P':2}
        
        if tn == 5:
            atypes = ['FAP4 H', 'FAP4 NI', 'FAP4 P']
            ahash = {'FAP4_H':0, 'FAP4_NI':1, 'FAP4_P':2}
            
        if tn == 6:
            atypes = ['FAP7 H', 'FAP7 NI', 'FAP7 P']
            ahash = {'FAP7_H':0, 'FAP7_NI':1, 'FAP7_P':2}
            
        if tn == 7:
            atypes = ['FAP5 H', 'FAP5 NI', 'FAP5 P']
            ahash = {'FAP5_H':0, 'FAP5_NI':1, 'FAP5_P':2}
            
        if tn == 8:
            atypes = ['FAP H', 'FAP NI', 'FAP P']
            ahash = {'FAP5_H':0, 'FAP1_H':0, 'FAP7_H':0, 'FAP4_H':0,
                     'FAP5_NI':1,'FAP1_NI':1, 'FAP7_NI':1, 'FAP4_NI':1,
                     'FAP5_P':2, 'FAP1_P':2, 'FAP7_P':2, 'FAP4_P':2}
            
        if tn == 9:
            atype = self.h.getSurvName('c Type')
            atypes = ['H', 'NI', 'P']
            ahash = {'H':0, 'NI':1, 'P':2}
            
        if tn == 9.5:
            atype = self.h.getSurvName('c Type')
            atypes = ['NI', 'P']
            ahash = {'NI':0, 'P':1}
            
        if tn == 10:
            atypes = ['JPS NI', 'PJS NI', 'Lynch NI', 'FAP NI']
            ahash = {'JPS_NI':0, 
                     'PJS_NI':1, 
                     'Lynch_NI':2,
                     'FAP5_NI':3,'FAP1_NI':3, 'FAP7_NI':3, 'FAP4_NI':3}
        
        if tn == 11:
            atypes = ['JPS P', 'PJS P', 'FAP P']
            ahash = {'JPS_P':0, 
                     'PJS_P':1, 
                     'FAP5_P':2,'FAP1_P':2, 'FAP7_P':2, 'FAP4_P':2}
        
        if tn == 12:
            atypes = ['H', 'JPS NI', 'PJS NI', 'Lynch NI', 'FAP NI']
            ahash = {'JPS_H':0, 'PJS_H':0, 'Lynch_H':0, 'FAP5_H':0, 'FAP1_H':0, 'FAP7_H':0, 'FAP4_H':0,
                     'JPS_NI':1, 
                     'PJS_NI':2, 
                     'Lynch_NI':3,
                     'FAP5_NI':4,'FAP1_NI':4, 'FAP7_NI':4, 'FAP4_NI':4}
            
        if tn == 13:
            atypes = ['H', 'JPS P', 'PJS P', 'FAP P']
            ahash = {'JPS_H':0, 'PJS_H':0, 'Lynch_H':0, 'FAP5_H':0, 'FAP1_H':0, 'FAP7_H':0, 'FAP4_H':0,
                     'JPS_P':1, 
                     'PJS_P':2, 
                     'FAP5_P':3,'FAP1_P':3, 'FAP7_P':3, 'FAP4_P':3}
        
        if tn == 14:
            atypes = ['H', 'JPS NI', 'PJS NI', 'Lynch NI', 'FAP NI', 'JPS P', 'PJS P', 'FAP P']
            ahash = {'JPS_H':0, 'PJS_H':0, 'Lynch_H':0, 'FAP5_H':0, 'FAP1_H':0, 'FAP7_H':0, 'FAP4_H':0,
                     'JPS_NI':1, 
                     'PJS_NI':2, 
                     'Lynch_NI':3,
                     'FAP5_NI':4,'FAP1_NI':4, 'FAP7_NI':4, 'FAP4_NI':4,
                     'JPS_P':5, 
                     'PJS_P':6, 
                     'FAP5_P':7,'FAP1_P':7, 'FAP7_P':7, 'FAP4_P':7}
            
        if tn == 15:
            atypes = ['JPS NI', 'JPS P']
            ahash = {'JPS_NI':0, 'JPS_P':1}
        
        if tn == 16:
            atypes = ['PJS NI', 'PJS P']
            ahash = {'PJS_NI':0, 'PJS_P':1}
            
        if tn == 18:
            atypes = ['FAP1 NI', 'FAP1 P']
            ahash = {'FAP1_NI':0, 'FAP1_P':1}
        
        if tn == 19:
            atypes = ['FAP4 NI', 'FAP4 P']
            ahash = {'FAP4_NI':0, 'FAP4_P':1}
            
        if tn == 20:
            atypes = ['FAP7 NI', 'FAP7 P']
            ahash = {'FAP7_NI':0, 'FAP7_P':1}
            
        if tn == 21:
            atypes = ['FAP5 NI', 'FAP5 P']
            ahash = {'FAP5_NI':0, 'FAP5_P':1}
            
        if tn == 22:
            atypes = ['FAP NI', 'FAP P']
            ahash = {'FAP5_NI':0,'FAP1_NI':0, 'FAP7_NI':0, 'FAP4_NI':0,
                     'FAP5_P':1, 'FAP1_P':1, 'FAP7_P':1, 'FAP4_P':1}
    
        self.initData(atype, atypes, ahash)
        return

    def getRNASeq_Hs_tissue_NI_P(self, tn=1):
        self.prepareData("PG0.5", cfile="/Users/dtv004/public_html/Hegemon/explore.conf")
        atype = self.h.getSurvName('c Sample ID')
        
        atypes = ['H', 'Adult NI', 'Adult P']
        ahash = {'H9_S73':0, 'org_H14_Colon_p13_S1':0,
                 'NI_ADULT_S72':1, 'C112216_S47':2}
            
        if tn == 3:
            atypes = ['NI AD', 'P AD']
            ahash = {'NI_ADULT_S72':0, 'C112216_S47':1}
            
        if tn == 4:
            atypes = [ 'H', 'NI PED', 'P PED']
            ahash = {'H9_S73':0, 'org_H14_Colon_p13_S1':0,
                     'FAP1_NI_S40':1, 'FAP_7_NI_S45':1, 'Lynch_2_TCNI_S43':1, 'Lynch1_no_PF_S51':1,
                     'FAP1_polyp_S41':2, 'FAP4_polyp_S42':2, 'FAP_7_polyp_S46':2, 'PJS1_TC_polyp_S44':2}
        
        if tn == 5:
            atypes = ['FAP NI', 'FAP P']
            ahash = {'FAP1_NI_S40':0, 'FAP_7_NI_S45':0,
                 'FAP1_polyp_S41':1, 'FAP4_polyp_S42':1, 'FAP_7_polyp_S46':1}
            
            
        if tn == 6:
            atypes = ['NI PED', 'P PED']
            ahash = {'FAP1_NI_S40':0, 'Lynch_2_TCNI_S43':0, 'Lynch1_no_PF_S51':0,
                     'FAP1_polyp_S41':1, 'FAP4_polyp_S42':1, 'FAP_7_polyp_S46':1, 'PJS1_TC_polyp_S44':1}
            
        if tn == 7:
            atypes = ['FAP NI', 'FAP P']
            ahash = {'FAP1_NI_S40':0,
                 'FAP1_polyp_S41':1, 'FAP4_polyp_S42':1, 'FAP_7_polyp_S46':1}
            
        self.initData(atype, atypes, ahash)   
        return
        
    def getPGFusoMM2021(self, tn=1):
        self.prepareData("PG0.75", cfile="/Users/dtv004/public_html/Hegemon/explore.conf")
        atype = self.h.getSurvName('c Sample ID')
        
        atypes = ['APC UN', 'APC FN']
        ahash = {'APC_UN1_S67':0, 'APC_UN2_S68':0, 'APC_Uninfecetd_1_S74':0, 'APC_Uninfecetd_2_S75':0,
                 'APC_FN1_S69':1, 'APC_FN2_S70':1, 'APC_FN3_S71':1, 'APC_FN_1_S78':1, 'APC__FN_2_S79':1}
        
        if tn == 0:
            atypes = ['APC UN', 'APC FN']
            ahash = {'APC_UN1_S67':0, 'APC_UN2_S68':0,
                     'APC_FN1_S69':1, 'APC_FN2_S70':1, 'APC_FN3_S71':1}

        if tn == 2:
            atypes = ['WT UN', 'WT FN']
            ahash = {'WT3_UN_S55':0, 'WT3_un1_S59':0, 'WT3_un2_S60':0, 
                     'WT3_FN_S57':1, 'WT3_FN1_S61':1, 'WT3_FN2_S62':1}
            
        if tn == 3:
            atypes = ['WT UN', 'WT FN', 'APC UN', 'APC FN']
            ahash = {'WT3_UN_S55':0, 'WT3_un1_S59':0, 'WT3_un2_S60':0, 
                     'WT3_FN_S57':1, 'WT3_FN1_S61':1, 'WT3_FN2_S62':1,
                    'APC_UN1_S67':2, 'APC_UN2_S68':2, 'APC_Uninfecetd_1_S74':2, 'APC_Uninfecetd_2_S75':2,
                     'APC_FN1_S69':3, 'APC_FN2_S70':3, 'APC_FN3_S71':3, 'APC_FN_1_S78':3, 'APC__FN_2_S79':3}
            
        if tn == 4:
            atypes = ['WT UN', 'WT FN', 'APC UN', 'APC FN']
            ahash = {'WT3_UN_S55':0, 'WT3_un1_S59':0, 'WT3_un2_S60':0, 
                     'WT3_FN_S57':1, 'WT3_FN1_S61':1, 'WT3_FN2_S62':1,
                    'APC_UN1_S67':2, 'APC_UN2_S68':2, 
                     'APC_FN1_S69':3, 'APC_FN2_S70':3, 'APC_FN3_S71':3}
            
        if tn == 5:
            atypes = ['UN', 'FN']
            ahash = {'WT3_UN_S55':0, 'WT3_un1_S59':0, 'WT3_un2_S60':0, 
                     'WT3_FN_S57':1, 'WT3_FN1_S61':1, 'WT3_FN2_S62':1,
                    'APC_UN1_S67':0, 'APC_UN2_S68':0, 
                     'APC_FN1_S69':1, 'APC_FN2_S70':1, 'APC_FN3_S71':1}
            
        if tn == 6:
            atypes = ['UN', 'FN']
            ahash = {'WT3_UN_S55':0, 'WT3_un1_S59':0, 'WT3_un2_S60':0, 
                     'WT3_FN_S57':1, 'WT3_FN1_S61':1, 'WT3_FN2_S62':1,
                    'APC_UN1_S67':0, 'APC_UN2_S68':0, 'APC_Uninfecetd_1_S74':0, 'APC_Uninfecetd_2_S75':0,
                 'APC_FN1_S69':1, 'APC_FN2_S70':1, 'APC_FN3_S71':1, 'APC_FN_1_S78':1, 'APC__FN_2_S79':1}
        
        if tn == 7:
            atypes = ['APC UN', 'WT_UN']
            ahash = {'APC_UN1_S67':0, 'APC_UN2_S68':0,
                     'WT3_UN_S55':1, 'WT3_un1_S59':1, 'WT3_un2_S60':1}
            
        if tn == 8:
            atypes = ['APC UN', 'WT_UN']
            ahash = {'APC_UN1_S67':0, 'APC_UN2_S68':0, 'APC_Uninfecetd_1_S74':0, 'APC_Uninfecetd_2_S75':0,
                     'WT3_UN_S55':1, 'WT3_un1_S59':1, 'WT3_un2_S60':1}
        
        self.initData(atype, atypes, ahash)
        return
        
    def getRNASeq_Hs_tissue_NI_P_data2(self, tn=1):
        self.prepareData("PG3", cfile="/Users/dtv004/public_html/Hegemon/explore.conf")
        atype = self.h.getSurvName('c Sample ID')
        
        atypes = ['H', 'NI']
        ahash = {'VA-DC_S45':0, 'VA-TC_S44':0,
                 'FAP-1_sig_S48':1, 'FAP-4_Cecum_NI_S49':1, 'FAP-7_DC_NI_S51':1, 'FAP-7_TC_NI_S52':1, 
                 'JPS1_REC_NI_S55':1, 'Lyn1-sig_colon_S46':1, 'Lynch2-sig_colon_S47':1,
                 'PJS2_SIG_NI_S54':1, 'PJS2__TC_NI_S53':1}
        
        if tn == 2:
            atypes = ['H', 'P']
            ahash = {'VA-DC_S45':0, 'VA-TC_S44':0,
                     'FAP-6_TC_polyp_S50':1, 'JPS_1_TC_POLYP_S56':1}
            
        if tn == 3:
            atypes = ['FAP NI', 'FAP P']
            ahash = {'FAP-1_sig_S48':0, 'FAP-4_Cecum_NI_S49':0, 'FAP-7_DC_NI_S51':0, 'FAP-7_TC_NI_S52':0,
                     'FAP-6_TC_polyp_S50':1}
            
        if tn == 4:
            atypes = ['JPS NI', 'JPS P']
            ahash = {'JPS1_REC_NI_S55':0,
                     'JPS_1_TC_POLYP_S56':1}
            
        if tn == 5:
            atypes = ['PJS NI', 'JPS NI', 'Lynch NI','FAP NI']
            ahash = {'FAP-1_sig_S48':3, 'FAP-4_Cecum_NI_S49':3, 'FAP-7_DC_NI_S51':3, 'FAP-7_TC_NI_S52':3, 
                     'JPS1_REC_NI_S55':1, 'Lyn1-sig_colon_S46':2, 'Lynch2-sig_colon_S47':2,
                     'PJS2_SIG_NI_S54':0, 'PJS2__TC_NI_S53':0}
            
        if tn == 6:
            atypes = [ 'JPS P','FAP P']
            ahash = {'FAP-6_TC_polyp_S50':1, 'JPS_1_TC_POLYP_S56':0}
        
        self.initData(atype, atypes, ahash)
        return
        
    def getJia2017(self, tn=1):
        self.prepareData("CRC112")
        atype = self.h.getSurvName('c src1')
        
        atypes = ['NI', 'FN']
        ahash = {'Caco-2 cells non-infected 24h':0, 'Caco-2 cells F.nucleatum 24h':1}
        
        self.initData(atype, atypes, ahash)
        return
        
    def getJia2017RMA(self, tn=1):
        self.prepareData("CRC112.2")
        atype = self.h.getSurvName('c treatment')
        
        atypes = ['NI', 'FN']
        ahash = {'non-infected':0, 'F. nucleatum':1}
        
        self.initData(atype, atypes, ahash)
        return
        
    def getPooled(self):
        self.prepareData("PLP33")
        atype = self.h.getSurvName("c Histology")
        atypes = ['Normal', 'Adenoma']
        ahash = {}
        self.initData(atype, atypes, ahash)
        return
        
    def getCAP(self):
        self.prepareData("PLP62.3")
        atype = self.h.getSurvName("c Subtype")
        btype = self.h.getSurvName("c Type")
        atypes = ['CFP', 'CAP']
        ahash = {}
        self.initData(atype, atypes, ahash)
        return
        
    def getPekow(self):
        self.prepareData("PLP59")
        atype = self.h.getSurvName("c src1")
        atypes = ['C', 'qUC', 'nUC']
        ahash = {'normal control': 0, 
                'quiescent ulcerative colitis': 1,
                'ulcerative colitis with neoplasia': 2}
        self.initData(atype, atypes, ahash)
        return
    
    def getPooled(self):
        self.prepareData("PLP33")
        atype = self.h.getSurvName("c Histology")
        atypes = ['Normal', 'Adenoma']
        ahash = {}
        self.initData(atype, atypes, ahash)
        return
    
    def getKanth2016(self, tn=1):
        self.prepareData("PLP9")
        atype = self.h.getSurvName("c Tissue")
        atypes = ['C', 'U', 'H', 'A', 'SSA', 'C']
        ahash = {'Control left colon':0, 'Control right colon':0,'Uninvolved left colon':1,'Uninvolved right colon':1,
                 'Hyperplastic polyp':3,'Adenomatous polyp':4,'Sessile serrated adenoma/polyp':5,'Colon adenocarcinoma':6}
        if (tn == 2):
            atypes = ['N', 'A', 'C']
            ahash = {'Control left colon':0,
                    'Control right colon':0,
                    'Uninvolved left colon':0,
                    'Uninvolved right colon':0,
                    'Hyperplastic polyp':1,
                    'Adenomatous polyp':1,
                    'Sessile serrated adenoma/polyp':1,
                    'Colon adenocarcinoma':2}
        if (tn == 3):
            atypes = ['N', 'A']
            ahash = {'Control left colon':0,
                    'Control right colon':0,
                    'Uninvolved left colon':0,
                    'Uninvolved right colon':0,
                    'Hyperplastic polyp':1,
                    'Adenomatous polyp':1,
                    'Sessile serrated adenoma/polyp':1}
        self.initData(atype, atypes, ahash)
        return
    
    def getFuso1(self):
        self.prepareData("CRC112.2")
        atype = self.h.getSurvName("c treatment")
        atypes = ['Uninfected', 'FN-infected']
        ahash = {"non-infected": 0, "F. nucleatum":1}
        self.initData(atype, atypes, ahash)
        return
    
    def getFuso2(self):
        self.prepareData("CRC114")
        atype = self.h.getSurvName("c Type")
        atypes = ['Normal', 'Tumor']
        ahash = {}
        self.initData(atype, atypes, ahash)
        return
    
    def getQu(self):
        self.prepareData("PLP50")
        atype = self.h.getSurvName("c tissue type")
        atypes = ['Normal', 'Adenoma']
        ahash = { 'Normal crypt epithelium': 0, 'Normal surface epithelium': 0}
        self.initData(atype, atypes, ahash)
        return

    def getReumers(self):
        self.prepareData("PLP65")
        atype = self.h.getSurvName("c tissue")
        atypes = ['Normal', 'Adenoma']
        ahash = { 'Adjacent Mucosa': 0, 'SSA': 1, 'colorectal adenoma': 1}
        self.initData(atype, atypes, ahash)
        return

    def getDruliner(self, tn=1):
        self.prepareData("PLP62.3")
        atype = self.h.getSurvName("c Type")
        atypes = ['Normal', 'Tubular', 'Villous']
        ahash = {'NORMAL':0, 'TUBULAR':1, "VILLOUS":2}
        if (tn == 2):
            atypes = ['Normal', 'Adenoma']
            ahash = {'NORMAL':0, 'TUBULAR':1, "VILLOUS":1}
        self.initData(atype, atypes, ahash)
        return

    def getPeters(self):
        self.prepareData("PLP67")
        atype = self.h.getSurvName("c clinical condition")
        atypes = ['Normal', 'UC', 'CD']
        ahash = {"control": 0, "Ulcerative Colitis":1, "Crohn's disease":2}
        self.initData(atype, atypes, ahash)
        return

    def getPaoni(self):
        self.prepareData("PLP51")
        atype = self.h.getSurvName("c Title")
        atypes = ['Normal', 'Adenoma', 'Cancer']
        atypes = ['Normal', 'Adenoma']
        ahash = {}
        aval = [ahash[i] if i in ahash else None for i in atype]
        normal = [ i for i in self.h.aRange() if atype[i].find("-WT") > 0]
        for i in normal:
            aval[i] = 0
        ad = [ i for i in self.h.aRange() if atype[i].find("-adenoma") > 0]
        for i in ad:
            aval[i] = 1
        ca = [ i for i in self.h.aRange() if atype[i].find("-carc") > 0]
        for i in ca:
            aval[i] = 2
        self.aval = aval
        self.atype = atype
        self.atypes = atypes
        self.order = normal + ad
        self.normal = normal
        self.ad = ad
        self.ca = ca
        self.printInfo()
        return


    def getBelmont(self):
        self.prepareData("PLP41")
        atype = self.h.getSurvName("c genotype")
        atypes = ['Normal', 'Adenoma', 'Cancer']
        atypes = ['Normal', 'Adenoma']
        ahash = {'wild type':0, 'APC mutant': 1, 'APC KRAS mutant': 2}
        ahash = {'wild type':0, 'APC mutant': 1}
        self.initData(atype, atypes, ahash)
        return

    def getLeclerc(self):
        self.prepareData("PLP45")
        atype = self.h.getSurvName("c src1")
        atypes = ['Normal', 'Adenoma']
        ahash = {'Normal intestine':0, 'Tumor': 1}
        self.initData(atype, atypes, ahash)
        return
    
    def getGalamb2008(self, tn=1):
        self.prepareData("PLP83")
        atype = self.h.getSurvName("c desc")
        atype = [ " ".join(str(i).split(" ")[-2:]) for i in atype]
        atypes = ['N', 'IBD', 'A', 'C']
        ahash = {'bowel disease':1,
                'colorectal cancer':3,
                'colon adenoma':2,
                'healthy control':0}
        if (tn == 2):
            atypes = ['N', 'A']
            ahash = {'colon adenoma':1,
                    'healthy control':0}
        self.initData(atype, atypes, ahash)
        return

    def getSabatesBellver2007(self, tn=1):
        self.prepareData("CRC19")
        atype = self.h.getSurvName("c Type")
        atypes = ['N', 'A']
        ahash = {'normal mucosa':0, 'colorectal adenoma':1}
        self.initData(atype, atypes, ahash)
        return

    def getLaPointe2010(self, tn=1):
        self.prepareData("CRC41")
        atype = self.h.getSurvName("c phenotype")
        atypes = ['N', 'A', 'C']
        ahash = {'Cancer':2, 'Normal':0, 'Adenoma':1}
        if (tn == 2):
            atypes = ['N', 'A']
            ahash = {'Normal':0, 'Adenoma':1}
        self.initData(atype, atypes, ahash)
        return

    def getSheffer2012(self, tn=1):
        self.prepareData("CRC49")
        atype = self.h.getSurvName("c tissue")
        atypes = ['NC', 'NV', 'NL', 'A', 'C', 'VM', 'LM']
        ahash = {'Primary Tumor':4,
                'Liver Metastasis':5,
                'Lung Metastasis':6,
                'Normal Lung':2,
                'Normal Colon':0,
                'Normal Liver':1,
                'Polyp':3,
                'Microadenoma':3,
                'Polyp, high grade':3}
        if (tn == 2):
            atypes = ['N', 'A']
            ahash = {'Normal Colon':0,
                    'Polyp':1,
                    'Microadenoma':1,
                    'Polyp, high grade':1}
        self.initData(atype, atypes, ahash)
        return

    def getSkrzypczak2010(self, tn=1):
        self.prepareData("CRC137")
        atype = self.h.getSurvName("c tissue")
        atypes = ['N', 'A', 'C']
        ahash = {'adenoma':1, 'adenocarcinoma':2, 'normal colon':0,
                'colon tumor':2}
        if (tn == 2):
            atypes = ['N', 'A']
            ahash = {'adenoma':1, 'normal colon':0}
        self.initData(atype, atypes, ahash)
        return

    def getFujii2017(self, tn=1):
        self.prepareData("ORG3")
        atype = self.h.getSurvName("c Type")
        atypes = ['N', 'A', 'C']
        ahash = {'CRC':2, 'Normal':0, 'Adenoma':1}
        if (tn == 2):
            atypes = ['N', 'A']
            ahash = {'Normal':0, 'Adenoma':1}
        self.initData(atype, atypes, ahash)
        return

    def getMatano2015(self, tn=1):
        self.prepareData("ORG4")
        gtype = self.h.getSurvName("c genetic modification")
        atype = self.h.getSurvName("c src1")
        atype = [str(i).split(" ")[1] if len(str(i).split(" ")) > 1 else None
                for i in atype]
        atype = [ str(atype[i]) + " " + str(gtype[i]) for i in
                range(len(atype))]
        atypes = ['KRAS', 'PIK3CA', 'Smad4', 'A', 'AKSTP', 
                'BKST', 'AKST', 'BKSTP', 'oA', 'C', 'M']
        ahash = {
                'epithelial KRAS G12V knock in':0,
                'epithelial PIK3CA E545K overexpression':1,
                'epithelial Smad4 deletion':2,
                'adenoma ':3,
                'organoids AKSTP (CRISPR/Cas9)':4,
                'organoids BKST (Lenti virus vector)':5,
                'organoids AKST (CRISPR/Cas9)':6,
                'organoids BKSTP (Lenti virus vector)':7,
                'organoids A (CRISPR/Cas9)':8,
                'cancer primary CRC':9,
                'colorectal LM CRC':10}
        self.initData(atype, atypes, ahash)
        return

    def getFessler2016(self, tn=1):
        self.prepareData("PLP3")
        atype = self.h.getSurvName("c tissue")
        atypes = ['N', 'TA', 'SSA', 'OC', 'OCTA']
        ahash = {
                'organoid culture':3,
                'organoid culture of tubular adenoma':4,
                'normal colon':0,
                'tubular adenomas':1,
                'sessile serrated adenomas':2}
        if (tn == 2):
            atypes = ['N', 'A']
            ahash = {'organoid culture':0,
                    'organoid culture of tubular adenoma':1,
                    'normal colon':0}

        self.initData(atype, atypes, ahash)
        return
    
    def getThiruvengadam2018(self, tn=1):
        self.prepareData("PLP28")
        atype = self.h.getSurvName("c tissue")
        atypes = ['N', 'A', 'C']
        ahash = {'adenoma':1, 'normal':0, 'adenocarcinoma':2}
        if (tn == 2):
            atypes = ['N', 'A']
            ahash = {'adenoma':1, 'normal':0}
        self.initData(atype, atypes, ahash)
        return

    def getDelker2018(self, tn=1):
        self.prepareData("PLP29")
        atype = self.h.getSurvName("c tissue")
        ttype = self.h.getSurvName("c timepoint")
        rtype = self.h.getSurvName("c treatment")
        atype = [ " ".join([str(atype[i]), str(ttype[i]), str(rtype[i])]) \
            for i in range(len(atype))]
        atypes = ['UED', 'UEP', 'UBN', 'PED', 'PEP']
        ahash = {'Uninvolved Endpoint Drug':0,
            'Uninvolved Endpoint Placebo':1,
            'Uninvolved Baseline None':2,
            'Polyp Endpoint Drug':3,
            'Polyp Endpoint Placebo':4}
        if (tn == 2):
            atype = self.h.getSurvName("c tissue")
        atypes = ['N', 'A']
        ahash = {'Uninvolved':0, 'Polyp':1}
        if (tn == 3):
            atypes = ['N', 'A']
            ahash = {'Uninvolved Baseline None':0,
                'Uninvolved Endpoint Placebo':0,
                'Polyp Endpoint Drug':1,
                'Polyp Endpoint Placebo':1}
        if (tn == 4):
            atypes = ['N', 'A']
            ahash = {'Uninvolved Endpoint Drug':0,
                'Polyp Endpoint Drug':1,
                'Polyp Endpoint Placebo':1}
        if (tn == 5):
            atypes = ['N', 'A']
            ahash = {'Uninvolved Endpoint Placebo':0,
                'Polyp Endpoint Drug':1,
                'Polyp Endpoint Placebo':1}
        self.initData(atype, atypes, ahash)
        return
    
    def getArbibe(self):
        self.prepareData("CRC143")
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
        self.initData(atype, atypes, ahash)
        return
    
    def getMatsuki(self):
        self.prepareData("CRC142")
        atype = self.h.getSurvName("c src1")
        atypes = ['C', 'Lc', 'Bb']
        ahash = {'Caco-2 cells cultured with B. breve':2,
            'Caco-2 cells alone':0,
            'Caco-2 cells cultured with L. casei': 1}
        self.initData(atype, atypes, ahash)
        return
    
    def getEColi(self):
        self.prepareData("CRC141")
        series = self.h.getSurvName("c Series")
        atype = self.h.getSurvName("c treatment")
        atypes = ['C', 'K12', 'O157']
        ahash = {'control, 60min':0, 'K-12, 60min':1,
                'O157:H7, 120min':2, 'control, 90min':0,
                'O157:H7, 60min':2, 'O157:H7, 90min':2,
                'control, 120min':0, 'K12, 120min':1, 'K-12, 90min':1}
        self.initData(atype, atypes, ahash)
        return
    
    def getIftekhar2020(self, tn=1):
        self.prepareData("PG19", cfile="/Users/dtv004/public_html/Hegemon/explore.conf")
        atype = self.h.getSurvName('c condition (ch1)')
        atypes = ['UN', 'INF']
        ahash = {'Uninfected control pool':0, 'Wnt-independent pool':1}
        
        if tn == 2:
            atypes = ['UN', 'INF']
            ahash = {'Uninfected control pool':0, 'Uninfected control clone':0,
                     'Wnt-independent pool':1, 'Wnt-independent clone':1}
        
        if tn == 3:
            atypes = ['UN', 'INF', 'UN C', 'INF C']
            ahash = {'Uninfected control pool':0, 'Uninfected control clone':2,
                     'Wnt-independent pool':1, 'Wnt-independent clone':3}
        
        self.initData(atype, atypes, ahash)    
        return
    
    def printHeatmap_test(self, ofile, genes, params):
        i1 = self.i1
        f_ranks = self.f_ranks
        self.params = {'genes': genes,'atypes': self.atypes,'cval': self.cval}
        self.params.update(params)
        ax = bone.plotHeatmap(ofile, self.data, self.col_labels, 
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
        score=bone.convertScore(predicted)
        print(list(actual))
        print(list(predicted))
        ax1 = bone.printReport(actual, predicted, score, target_names)
        tab = pd.crosstab(df.y > 0, df.x > 0)
        print(tab)
        print(fisher_exact(tab))
        print('Fisher Exact pvalue =', fisher_exact(tab)[1])
        return ax, ax1
        
def readGenes(cfile):
    if not os.path.isfile(cfile):
        print("Can't open file")
        exit()
    fp = open(cfile, "r")
    nodelist = []
    for line in fp:
        line = line.strip();
        ll = re.split("[\[\]()\s]", line);
        nodelist += ll
    fp.close();
    return [i for i in hu.uniq(nodelist) if i != '']

def getGeneGroups():
    reload(hu)
    db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
    h = hu.Hegemon(db.getDataset("PLP33"))
    h.init()
    h.initPlatform()
    h.initSurv()
    with open('Supplementary/path-8.json') as data_file:
        data_item = json.load(data_file)
    cfile = "Supplementary/colon-network-clusters.txt"
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
        print(i, gn, h.getSimpleName(gn), data_item[i][0], len(gene_groups[i]))
    print([len(s) for s in gene_groups])
    order = [8, 7, 6, 1];
    gene_groups = [gene_groups[i] for i in order]
    print([len(s) for s in gene_groups])
    weight = [0, 0, 0, 1]
    weight = [-3, -2, -0.5, 2]
    genes = readGenes("Supplementary/cluster-names.txt")
    return genes, weight, gene_groups

def getGeneGroups2(order = None, weight = None, debug = 1):
    db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
    h = hu.Hegemon(db.getDataset("PLP33"))
    h.init()
    h.initPlatform()
    h.initSurv()
    dir1 = "/booleanfs2/sahoo/Data/Colon/Adenoma/"
    cfile = dir1 + "Supplementary/colon-network-clusters.txt"
    if not os.path.isfile(cfile):
        print("Can't open file {0}".format(cfile));
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
    gene_groups = [[]]
    for i in range(1, 16):
        cfile = dir1 + "Supplementary/node-" + str(i) + ".txt"
        res = readGenes(cfile)
        s1 = set()
        for g in res:
            s1.add(g)
            if g in nodelist:
                for k in nodelist[g]:
                    s1.add(k)
        gene_groups.append(s1);
        if debug == 1:
            print(i, res[0], h.getSimpleName(res[0]), len(s1))
    print([len(s) for s in gene_groups])
    if order is None:
        order = [1, 2, 3, 4]
    gene_groups = [gene_groups[i] for i in order]
    print([len(s) for s in gene_groups])
    #gene_groups = getSimpleName(gene_groups, h)
    #print([len(s) for s in gene_groups])
    if weight is None:
        weight = [-3, -2, -0.5, 2]
    print(weight)
    genes = readGenes(dir1 + "cluster-names.txt")
    return genes, weight, gene_groups

def processGeneGroupsMm(ana, l1, wt1, debug = 0, fthr = None):
    ana.convertMm(l1, [])
    ana.orderData(ana.gene_groups, wt1)
    print("ROC-AUC", ana.getMetrics())
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

def plotqPCRHeatmap(ofile, data, col_labels, row_labels, params):
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
    if 'cval2' in params:
        cval = params['cval2']
    
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
    
    bone.barTop(lax, atypes, color_sch1, params)

    fig.savefig(ofile, dpi=dpi)
    return ax, divider    

def processData(ana, l1, wt1, desc=None, violin=1):
    ana.orderData(l1, wt1)
    if (violin == 1):
        return plotViolinBar(ana, desc)
    return plotDensityBar(ana, desc)

def processComp(ana, cap=0):
	fig = plt.figure(figsize=(8,2), dpi=100)
	res = []

	order = [1, 2, 3, 4, 5]
	wt1 = [-5, -0.3, 0.1, 2.9, -4]
	genes, wt1, l1 = getGeneGroups2(order, wt1, 0)
	c_dict, fpr, tpr, roc_auc = bone.processGeneGroups(ana, l1, wt1)
	res += [[c_dict, fpr, tpr, roc_auc]]
	ax = plt.subplot2grid((4, 5), (0, 2), colspan=3)
	params = {'spaceAnn': len(ana.order)/len(ana.atypes), 'tAnn': 1,
		'widthAnn':1, 'genes': genes, 'acolor': acolor, 'ax': ax}
	ax = ana.printTitleBar(params)
	desc = "Boolean"
	ax.text(-1, 2, desc, horizontalalignment='right',
		verticalalignment='center')

	l2 = [["SLC1A2", "TCF7L2", "COL18A1", "KLF6", "OAT"],
	  ['CTNNB1','AXIN2','TCF4','LEF1','EPHB2','EPHB3','HNF1A',
	   'MYC','CCND1','FAT1','SP5','ZNRF3','LGR5','KIAA1199',
	   'ASCL2','PPARG','RNF43','CD44','BMP7','ADRA2C','FZD7',
	   'IL8','TBX3','NKD1','DKK1','DEFA6','GLUL','REG1B','SOX9',
	   'TDGF1']]
	wt2 = [-1, 1]
	c_dict, fpr, tpr, roc_auc = bone.processGeneGroups(ana, l2, wt2)
	res += [[c_dict, fpr, tpr, roc_auc]]
	ax = plt.subplot2grid((4, 5), (1, 2), colspan=3)
	params = {'spaceAnn': len(ana.order)/len(ana.atypes), 'tAnn': 1,
		'widthAnn':1, 'genes': genes, 'acolor': acolor, 'ax': ax}
	ax = ana.printTitleBar(params)
	desc = "Bayesian"
	ax.text(-1, 2, desc, horizontalalignment='right',
		verticalalignment='center')

	ifile = "Supplementary/sabates-bellver-diff.txt"
	high, low = bone.getFdrStats(ifile, 1e-6, 0)
	print(len(high), len(low))
	l2 = [high, low]
	df = pd.read_csv(ifile, sep='\t',
		    names=["transcript_id", "name", "stat", "pval", "diff"])
	l2 = [set(df.nlargest(857, 'stat')['name']),
		       set(df.nsmallest(808, 'stat')['name'])]
	wt2 = [-1, 1]
	c_dict, fpr, tpr, roc_auc = bone.processGeneGroups(ana, l2, wt2)
	res += [[c_dict, fpr, tpr, roc_auc]]
	ax = plt.subplot2grid((4, 5), (2, 2), colspan=3)
	params = {'spaceAnn': len(ana.order)/len(ana.atypes), 'tAnn': 1,
		'widthAnn':1, 'genes': genes, 'acolor': acolor, 'ax': ax}
	ax = ana.printTitleBar(params)
	desc = "Differential"
	ax.text(-1, 2, desc, horizontalalignment='right',
		verticalalignment='center')

	l2 = [bone.getEntries("Supplementary/lee-2006-down.txt", 0),
	bone.getEntries("Supplementary/lee-2006-up.txt", 0)]
	wt2 = [-1, 1]
	c_dict, fpr, tpr, roc_auc = bone.processGeneGroups(ana, l2, wt2)
	res += [[c_dict, fpr, tpr, roc_auc]]
	ax = plt.subplot2grid((4, 5), (3, 2), colspan=3)
	params = {'spaceAnn': len(ana.order)/len(ana.atypes), 'tAnn': 1,
		'widthAnn':1, 'genes': genes, 'acolor': acolor, 'ax': ax}
	ax = ana.printTitleBar(params)
	desc = "Differential"
	ax.text(-1, 2, desc, horizontalalignment='right',
		verticalalignment='center')

	actual = [1 if ana.aval[i] >= 1 else 0 for i in ana.i1]
	data_list = { 'y': actual }
	for i in range(len(res)):
	    id1 = 'c' + str(i)
	    c = [res[i][0][j] for j in ana.i1]
	    data_list[id1] = c
	df1 = pd.DataFrame(data_list)

	df = bone.printOLS('y ~ c0 + c1 + c2 + c3', df1)
	df = df.drop(['Intercept'])
	if (cap == 1):
	    df = bone.printUOLS('y', ['c0', 'c1', 'c2', 'c3'], df1)
	df["Name"] = ["Boolean", "Bayesian", "Differential", "Differential"]
	df['ROC-AUC-1'] = [k[3] for k in res]
	df = df.reindex(index=df.index[::-1])
	ax = plt.subplot2grid((1, 5), (0, 1))
	ax.errorbar(df["coeff"], range(len(df.index)),
	    yerr=0,
	    xerr=[list(df["coeff"] - df["lower 0.95"]), list(df["upper 0.95"] - df["coeff"])],
	    fmt='o', capsize=3)
	ax.set_yticks(range(len(df.index)))
	ax.set_yticklabels(df["Name"])
	ax.set_xlabel("coefficient")
	ax.axvline(x=0)
	#ax.set_xlim([0, 5.5])
	ax.set_ylim([-0.5, len(df.index) - 0.5])
	for i in range(len(df.index)):
	    ax.text(df["upper 0.95"][i] + 0.02,i, df['codes'][i], verticalalignment=
	'center')
	ax.set_title("Normal vs Adenoma (GSE76987; N = 41, AD = 41)")
	ax = plt.subplot2grid((1, 5), (0, 0))
	ax = df.plot.barh(x="Name", y="ROC-AUC-1", color='orange', ax=ax)
	if (cap == 0):
		ax.set_xlim([0.59, 0.9])
	else:
		ax.set_xlim([0.2, 0.9])
	return fig

def printUOLS(target, list1, df1):
    import statsmodels.formula.api as smf
    res = []
    for v in list1:
        lm1 = smf.ols(formula= target + " ~ " + v, data=df1).fit()
        idx = lm1.params.index
        ci = lm1.conf_int()
        ci_1 = [ ci[0][i] for i in range(len(idx))]
        ci_2 = [ ci[1][i] for i in range(len(idx))]
        c_1 = [ bone.getCode(p) for p in lm1.pvalues]
        res += [[idx[1], lm1.params[1], ci_1[1], ci_2[1], lm1.pvalues[1], c_1[1]]]
    df = pd.DataFrame(res, columns=['Name', 'coeff', 'lower 0.95',
        'upper 0.95', 'pvalues', 'codes'])
    df.index = list1
    return df

def getPDF(cfile):
    import bone
    reload(bone)
    from matplotlib.backends.backend_pdf import PdfPages

    pdf = PdfPages(cfile)
    return pdf

def closePDF(pdf):
    import datetime
    d = pdf.infodict()
    d['Title'] = 'Plots'
    d['Author'] = 'Daniella Vo'
    d['Subject'] = "Microbe Polyp"
    d['Keywords'] = 'disease training validation ROC'
    d['CreationDate'] = datetime.datetime(2021, 1, 18)
    d['ModDate'] = datetime.datetime.today()
    pdf.close()