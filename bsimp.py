##
#---------------------------------------------------------------------
# SERVER only input all files (.bam and .fa) output MeH matrix in .csv
# Oct 19, 2021 ML after imputation test
# github
#---------------------------------------------------------------------

import random
import math
import pysam
import csv
import sys
import pandas as pd
import numpy as np
import datetime
import time as t
from collections import Counter, defaultdict, OrderedDict


#---------------------------------------

# Functions definition

#---------------------------------------

    
def open_log(fname):
    open_log.logfile = open(fname, 'w', 1)
    

def logm(message):
    log_message = "[%s] %s\n" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), message)
    print(log_message),
    open_log.logfile.write(log_message)

def close_log():
    open_log.logfile.close()

# Check whether a window has enough reads for complete/impute
def enough_reads(window,w,complete):
    temp=np.isnan(window).sum(axis=1)==0
    if complete: # For heterogeneity estimation
        return temp.sum()>=3
    else:  # for imputation
        tempw1=np.isnan(window).sum(axis=1)==1
        #return temp.sum()>=2**(w-2) and tempw1.sum()>0
        return temp.sum()>=2 and tempw1.sum()>0
    
def impute(window,w):
    full_ind=np.where(np.isnan(window).sum(axis=1)==0)[0]
    part_ind=np.where(np.isnan(window).sum(axis=1)==1)[0]
    for i in range(len(part_ind)):
        sam = []
        # which column is nan
        pos=np.where(np.isnan(window[part_ind[i],:]))[0]
        if np.unique(window[np.where(np.invert(np.isnan(window[:,pos])))[0],pos]).shape[0]==1:
            window[part_ind[i],pos]=window[np.where(np.invert(np.isnan(window[:,pos])))[0],pos][0]
        else:
            #print("win_part i pos =",window[part_ind[i],pos])
            for j in range(len(full_ind)):
                if (window[part_ind[i],:]==window[full_ind[j],:]).sum()==w-1:
                    sam.append(j)
            if len(sam)>0:
                s1=random.sample(sam, 1)
                s=window[full_ind[s1],pos]
            else:
                s=random.sample(window[np.where(np.invert(np.isnan(window[:,pos])))[0],pos].tolist(), k=1)[0]
            window[part_ind[i],pos]=np.float64(s)
            #print("win_part i =",window[part_ind[i],pos])
            #print("s = ",np.float64(s))
    return window 
   
def outwindow(pat,pos,chrom,w,M,UM,mC=4,strand='f'): 
    
    # get complete reads
    temp=np.isnan(pat).sum(axis=1)==0
    pat=pat[np.where(temp)[0],:]
    
    count=np.zeros((2**w,1))
    m=np.shape(pat)[0]
    
    pat=np.array(pat)
    if w==2:
        pat = Counter([str(i[0])+str(i[1]) for i in pat.astype(int).tolist()])
        count=np.array([float(pat[i]) for i in ['00','10','01','11']])
    if w==3:
        pat = Counter([str(i[0])+str(i[1])+str(i[2]) for i in pat.astype(int).tolist()])
        count=np.array([float(pat[i]) for i in ['000','100','010','110','001','101','011','111']])
    if w==4:
        pat = Counter([str(i[0])+str(i[1])+str(i[2])+str(i[3]) for i in pat.astype(int).tolist()])
        count=np.array([float(pat[i]) for i in ['0000','1000','0100','1100','0010','1010','0110','1110','0001',\
                                '1001','0101','1101','0011','1011','0111','1111']])
    if w==5:
        pat = Counter([str(i[0])+str(i[1])+str(i[2])+str(i[3])+str(i[4]) for i in pat.astype(int).tolist()])
        count=np.array([float(pat[i]) for i in ['00000','10000','01000','11000','00100','10100','01100','11100','00010',\
                                '10010','01010','11010','00110','10110','01110','11110','00001','10001','01001','11001','00101',\
                                '10101','01101','11101','00011','10011','01011','11011','00111','10111','01111','11111']])
    if w==6:
        pat = Counter([str(i[0])+str(i[1])+str(i[2])+str(i[3])+str(i[4])+str(i[5]) for i in pat.astype(int).tolist()])
        count=np.array([float(pat[i]) for i in ['000000','100000','010000','110000','001000','101000','011000','111000','000100',\
                                '100100','010100','110100','001100','101100','011100','111100','000010','100010','010010','110010','001010',\
                                '101010','011010','111010','000110', '100110','010110','110110','001110','101110','011110','111110',\
                                '000001','100001','010001','110001','001001','101001','011001','111001','000101',\
                                '100101','010101','110101','001101','101101','011101','111101','000011','100011','010011','110011','001011',\
                                '101011','011011','111011','000111', '100111','010111','110111','001111','101111','011111','111111']])
    
    count=count.reshape(2**w)
    count=np.concatenate((count[[0]],count))
    if w==3:
        opt=pd.DataFrame({'chrom':chrom,'pos':pos,'p01':count[1],'p02':count[2],'p03':count[3],'p04':count[4],\
                        'p05':count[5],'p06':count[6],'p07':count[7],'p08':count[8],'M':M,'UM':UM,'strand':strand}, index=[0])     
    if w==4:
        opt=pd.DataFrame({'chrom':chrom,'pos':pos,'p01':count[1],'p02':count[2],'p03':count[3],'p04':count[4],\
                        'p05':count[5],'p06':count[6],'p07':count[7],'p08':count[8],'p09':count[9],'p10':count[10],\
                        'p11':count[11],'p12':count[12],'p13':count[13],'p14':count[14],'p15':count[15],\
                        'p16':count[16],'M':M,'UM':UM,'strand':strand}, index=[0])   
    if w==5:
        opt=pd.DataFrame({'chrom':chrom,'pos':pos,'p01':count[1],'p02':count[2],'p03':count[3],'p04':count[4],\
                        'p05':count[5],'p06':count[6],'p07':count[7],'p08':count[8],'p09':count[9],'p10':count[10],\
                        'p11':count[11],'p12':count[12],'p13':count[13],'p14':count[14],'p15':count[15],\
                        'p16':count[16],'p17':count[17],'p18':count[18],'p19':count[19],'p20':count[20],\
                        'p21':count[21],'p22':count[22],'p23':count[23],'p24':count[24],'p25':count[25],\
                        'p26':count[26],'p27':count[27],'p28':count[28],'p29':count[29],'p30':count[30],\
                        'p31':count[31],'p32':count[32],'M':M,'UM':UM,'strand':strand}, index=[0])    
    if w==6:
        opt=pd.DataFrame({'chrom':chrom,'pos':pos,'p01':count[1],'p02':count[2],'p03':count[3],'p04':count[4],\
                        'p05':count[5],'p06':count[6],'p07':count[7],'p08':count[8],'p09':count[9],'p10':count[10],\
                        'p11':count[11],'p12':count[12],'p13':count[13],'p14':count[14],'p15':count[15],\
                        'p16':count[16],'p17':count[17],'p18':count[18],'p19':count[19],'p20':count[20],\
                        'p21':count[21],'p22':count[22],'p23':count[23],'p24':count[24],'p25':count[25],\
                        'p26':count[26],'p27':count[27],'p28':count[28],'p29':count[29],'p30':count[30],\
                        'p31':count[31],'p32':count[32],'p33':count[33],'p34':count[34],'p35':count[35],\
                        'p36':count[36],'p37':count[37],'p38':count[38],'p39':count[39],'p40':count[40],\
                        'p41':count[41],'p42':count[42],'p43':count[43],'p44':count[44],'p45':count[45],\
                        'p46':count[46],'p47':count[47],'p48':count[48],'p49':count[49],'p50':count[50],\
                        'p51':count[51],'p52':count[52],'p53':count[53],'p54':count[54],'p55':count[55],\
                        'p56':count[56],'p57':count[57],'p58':count[58],'p59':count[59],'p60':count[60],\
                        'p61':count[61],'p62':count[62],'p63':count[63],'p64':count[64],'M':M,'UM':UM,\
                       'strand':strand}, index=[0])    
    
    return opt


def impute(window,w):
    full_ind=np.where(np.isnan(window).sum(axis=1)==0)[0]
    part_ind=np.where(np.isnan(window).sum(axis=1)==1)[0]
    for i in range(len(part_ind)):
        sam = []
        # which column is nan
        pos=np.where(np.isnan(window[part_ind[i],:]))[0]
        if np.unique(window[np.where(np.invert(np.isnan(window[:,pos])))[0],pos]).shape[0]==1:
            window[part_ind[i],pos]=window[np.where(np.invert(np.isnan(window[:,pos])))[0],pos][0]
        else:
            #print("win_part i pos =",window[part_ind[i],pos])
            for j in range(len(full_ind)):
                if (window[part_ind[i],:]==window[full_ind[j],:]).sum()==w-1:
                    sam.append(j)
            if len(sam)>0:
                s1=random.sample(sam, 1)
                s=window[full_ind[s1],pos]
            else:
                s=random.sample(window[np.where(np.invert(np.isnan(window[:,pos])))[0],pos].tolist(), k=1)[0]
            window[part_ind[i],pos]=np.float64(s)
            #print("win_part i =",window[part_ind[i],pos])
            #print("s = ",np.float64(s))
    return window 


def CGgenome_scr(bamfile,chrom,w,fa,optional,mC=4,silence=False):
    filename, file_extension = os.path.splitext(bamfile)
    coverage = cov_context = 0
    # load bamfile
    samfile = pysam.AlignmentFile("MeHdata/%s.bam" % (filename), "rb")
    # load reference genome
    fastafile = pysam.FastaFile('MeHdata/%s.fa' % fa)
    
    # initialise data frame for genome screening (load C from bam file)
    aggreR = aggreC = pd.DataFrame(columns=['Qname'])
    
    # if user wants to output compositions of methylation patterns at every eligible window, initialise data frame

    if w==3:
        ResultPW=pd.DataFrame(columns=['chrom','pos','p01','p02','p03','p04',\
                        'p05','p06','p07','p08','M','UM','strand'])     
    if w==4:
        ResultPW=pd.DataFrame(columns=['chrom','pos','p01','p02','p03','p04',\
                        'p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16','M','UM','strand'])   
    if w==5:
        ResultPW=pd.DataFrame(columns=['chrom','pos','p01','p02','p03','p04',\
                        'p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16','p17','p18',\
                        'p19','p20','p21','p22','p23','p24','p25','p26','p27','p28','p29','p30','p31','p32',\
                        'M','UM','strand'])    
    if w==6:
        ResultPW=pd.DataFrame(columns=['chrom','pos','p01','p02','p03','p04',\
                        'p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16','p17','p18',\
                        'p19','p20','p21','p22','p23','p24','p25','p26','p27','p28','p29','p30','p31','p32',\
                        'p33','p34','p35','p36','p37','p38','p39','p40','p41','p42','p43','p44','p45','p46',\
                        'p47','p48','p49','p50','p51','p52','p53','p54','p55','p56','p57','p58','p59','p60',\
                        'p61','p62','p63','p64','M','UM','strand'])  
        
    if w==7:
        ResultPW = pd.DataFrame(columns=\
                        ['chrom','pos','M','UM','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','p33','p34','p35','p36','p37','p38','p39','p40','p41','p42','p43','p44','p45','p46'\
                         ,'p47','p48','p49','p50','p51','p52','p53','p54','p55','p56','p57','p58','p59','p60','p61','p62','p63','p64'\
                         ,'p65','p66','p67','p68','p69','p70','p71','p72','p73','p74','p75','p76','p77','p78','p79','p80','p81','p82','p83','p84','p85','p86'\
                         ,'p87','p88','p89','p90','p91','p92','p93','p94','p95','p96','p97','p98','p99','p100','p101','p102','p103','p104'\
                        ,'p105','p106','p107','p108','p109','p120','p121','p122','p123','p124','p125','p126','p127','p128','strand'])

    

    neverr = never = True
    
    chrom_list = []
    # all samples' bam files
    for i in samfile.get_index_statistics():
        chrom_list.append(i.contig)

    if chrom in chrom_list:
        # screen bamfile by column
        for pileupcolumn in samfile.pileup(chrom):
            coverage += 1
            if not silence:
                if (pileupcolumn.pos % 2000000 == 1):
                    print("CG %s s %s w %s %s pos %s Result %s" % (datetime.datetime.now(),filename,w,chrom,pileupcolumn.pos,ResultPW.shape[0]))

            # Forward strand, check if 'CG' in reference genome 
            if (fastafile.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos+2)=='CG'):        
                cov_context += 1
                temp = pd.DataFrame(columns=['Qname',pileupcolumn.pos+1])
                pileupcolumn.set_min_base_quality(0)
                # append reads in the column
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and not pileupread.alignment.is_reverse:  # C
                        d = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos+1: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                        df2 = pd.DataFrame(data=d)
                        temp=temp.append(df2, ignore_index=True)

                # merge with other columns
                if (not temp.empty):
                    aggreC = pd.merge(aggreC,temp,how='outer',on=['Qname'])
                    aggreC = aggreC.drop_duplicates()

            # Reverse strand, check if 'CG' in reference genome 
            if pileupcolumn.pos>1:
                if (fastafile.fetch(chrom,pileupcolumn.pos-1,pileupcolumn.pos+1)=='CG'):
                    cov_context += 1
                    tempr = pd.DataFrame(columns=['Qname',pileupcolumn.pos+1])
                    pileupcolumn.set_min_base_quality(0)
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.is_reverse:  # C
                            dr = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos+1: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                            dfr2 = pd.DataFrame(data=dr)
                            tempr=tempr.append(dfr2, ignore_index=True)

                    if (not tempr.empty):
                        aggreR = pd.merge(aggreR,tempr,how='outer',on=['Qname'])
                        aggreR = aggreR.drop_duplicates()

            # Impute and estimate, if there are 2w-1 columns
            if never and aggreC.shape[1] == (2*w):
                # C/G to 1, rest to 0, N to NA
                never = False
                aggreC = aggreC.replace(['C'],1)
                aggreC = aggreC.replace(['T'],0)
                aggreC = aggreC.replace(['A','N','G'],np.nan)
                methbin = aggreC 
                meth = methbin.copy()
                # remove read ID
                meth = meth.drop('Qname',axis=1)
                # back up for imputation
                methtemp = meth.copy()
                # imputation by sliding windows of w C by 1 C
                for i in range(0,meth.shape[1]-w+1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    # save methylation statuses before imputation
                    # check if eligible for imputation, impute
                    if enough_reads(window,w,complete=False):
                        window=pd.DataFrame(data=impute(window,w))
                        ind=np.where(window.notnull().sum(axis=1)==w)[0]
                        methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
                # overwrite imputed window
                meth = methtemp.copy()
                # Evaluate methylation level and methylation heterogeneity and append to result
                for i in range(0,w,1): # w windows
                    window = meth.iloc[:,range(i,i+w)].values
                    M=(window==1).sum(axis=0)[0]
                    UM=(window==0).sum(axis=0)[0]
                    if depth>=mC:
                    # check if enough complete patterns for evaluating MeH
                        toappend=outwindow(window,w=w,pos=meth.iloc[:,range(i,i+w)].columns[0],\
                                            chrom=chrom,strand='f',mC=mC)
                        ResultPW=ResultPW.append(toappend)

                # remove 1 column
                aggreC = aggreC.drop(meth.columns[0:1],axis=1)
                # drop rows with no values
                aggreC.dropna(axis = 0, thresh=2, inplace = True)
                # total += w

            # Reverse
            if neverr and aggreR.shape[1] == (2*w):
                neverr = False
                aggreR = aggreR.replace(['G'],1)
                aggreR = aggreR.replace(['A'],0)
                aggreR = aggreR.replace(['C','N','T'],np.nan)
                methbin = aggreR # backup
                #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
                meth = methbin.copy()
                meth = meth.drop('Qname',axis=1)
                methtemp = meth.copy()
                # impute once if valid
                for i in range(0,meth.shape[1]-w+1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    # if eligible for imputation
                    if enough_reads(window,w,complete=False):
                        window=pd.DataFrame(data=impute(window,w))
                        ind=np.where(window.notnull().sum(axis=1)==w)[0]
                        methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
                meth = methtemp.copy()
                # compute coverage and output summary
                # for i in range(0,meth.shape[1]-w+1,1):
                # if i<w:
                for i in range(0,w,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    M=(window==1).sum(axis=0)[0]
                    UM=(window==0).sum(axis=0)[0]
                    if depth>=mC:
                        toappend=outwindow(window,w=w,pos=meth.iloc[:,range(i,i+w)].columns[0],\
                                            chrom=chrom,strand='r',mC=mC)
                        ResultPW=ResultPW.append(toappend)
                aggreR = aggreR.drop(meth.columns[0:1],axis=1)
                aggreR.dropna(axis = 0, thresh=2, inplace = True)

            #------------------
            #  SECONDARY CASE
            #------------------

            if (aggreC.shape[1] == (3*w-1)):
                aggreC = aggreC.replace(['C'],1)
                aggreC = aggreC.replace(['T'],0)
                aggreC = aggreC.replace(['A','N','G'],np.nan)
                methbin = aggreC # backup
                #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
                meth = methbin.copy()
                meth = meth.drop('Qname',axis=1)
                methtemp = meth.copy()
                # impute once if valid
                for i in range(0,meth.shape[1]-w+1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    if enough_reads(window,w,complete=False):
                        window=pd.DataFrame(data=impute(window,w))
                        ind=np.where(window.notnull().sum(axis=1)==w)[0]
                        methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
                meth = methtemp.copy()        
                # compute coverage and output summary
                for i in range(w-1,2*w-1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    M=(window==1).sum(axis=0)[0]
                    UM=(window==0).sum(axis=0)[0]
                    if depth>=mC:
                        toappend=outwindow(window,w=w,pos=meth.iloc[:,range(i,i+w)].columns[0],\
                                            chrom=chrom,strand='f',mC=mC)
                        ResultPW=ResultPW.append(toappend)

                        if ResultPW.shape[0] % 100000 == 1:   
                            ResultPW.to_csv(r"MeHdata/CG_%s_%s.csv"%(filename,chrom),index = False, header=True)
                            if not silence: 
                                print("Checkpoint CG. For file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos))
                aggreC = aggreC.drop(meth.columns[0:w],axis=1)
                aggreC.dropna(axis = 0, thresh=2, inplace = True)
                #print(aggreC)
                #total += w

            # reverse
            if (aggreR.shape[1] == (3*w-1)):
                aggreR = aggreR.replace(['G'],1)
                aggreR = aggreR.replace(['A'],0)
                aggreR = aggreR.replace(['C','N','T'],np.nan)
                methbin = aggreR # backup
                #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
                meth = methbin.copy()
                meth = meth.drop('Qname',axis=1)
                methtemp = meth.copy()
                # impute once if valid
                for i in range(0,meth.shape[1]-w+1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    # if eligible for imputation
                    if enough_reads(window,w,complete=False):
                        window=pd.DataFrame(data=impute(window,w))
                        ind=np.where(window.notnull().sum(axis=1)==w)[0]
                        methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
                meth = methtemp.copy()        
                # compute coverage and output summary
                for i in range(w-1,2*w-1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    M=(window==1).sum(axis=0)[0]
                    UM=(window==0).sum(axis=0)[0]
                    if depth>=mC:
                        toappend=outwindow(window,w=w,pos=meth.iloc[:,range(i,i+w)].columns[0],\
                                            chrom=chrom,strand='r',mC=mC)
                        ResultPW=ResultPW.append(toappend)
                        if ResultPW.shape[0] % 100000 == 1:   
                            ResultPW.to_csv(r"MeHdata/CG_%s_%s.csv"%(filename,chrom),index = False, header=True)

                aggreR = aggreR.drop(meth.columns[0:w],axis=1)
                aggreR.dropna(axis = 0, thresh=2, inplace = True)

        if ResultPW.shape[0]>0:   
            ResultPW.to_csv(r"MeHdata/CG_%s_%s.csv"%(filename,chrom),index = False, header=True)

        return filename, coverage, cov_context, 'CG'        
        print("Done CG for file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos))

    #samfile.close()  
    
def CHHgenome_scr(bamfile,chrom,w,fa,mC=4,silence=False):
    filename, file_extension = os.path.splitext(bamfile)
    coverage = cov_context = 0
    # load bamfile
    samfile = pysam.AlignmentFile("MeHdata/%s.bam" % (filename), "rb")
    # load reference genome
    fastafile = pysam.FastaFile('MeHdata/%s.fa' % fa)
    
    # initialise data frame for genome screening (load C from bam file)
    aggreR = aggreC = pd.DataFrame(columns=['Qname'])
    
    # if user wants to output compositions of methylation patterns at every eligible window, initialise data frame
    if w==3:
        ResultPW = pd.DataFrame(columns=\
                        ['chrom','pos','M','UM','p01','p02','p03','p04','p05','p06','p07','p08','strand'])
    if w==4:
        ResultPW = pd.DataFrame(columns=\
                        ['chrom','pos','M','UM','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11',\
                         'p12','p13','p14','p15','p16','strand'])
    if w==5:
        ResultPW = pd.DataFrame(columns=\
                        ['chrom','pos','M','UM','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','strand'])
    if w==6:
        ResultPW = pd.DataFrame(columns=\
                        ['chrom','pos','M','UM','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','p33','p34','p35','p36','p37','p38','p39','p40','p41','p42','p43','p44','p45','p46'\
                         ,'p47','p48','p49','p50','p51','p52','p53','p54','p55','p56','p57','p58','p59','p60','p61','p62','p63','p64'\
                         ,'strand'])
    if w==7:
        ResultPW = pd.DataFrame(columns=\
                        ['chrom','pos','M','UM','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','p33','p34','p35','p36','p37','p38','p39','p40','p41','p42','p43','p44','p45','p46'\
                         ,'p47','p48','p49','p50','p51','p52','p53','p54','p55','p56','p57','p58','p59','p60','p61','p62','p63','p64'\
                         ,'p65','p66','p67','p68','p69','p70','p71','p72','p73','p74','p75','p76','p77','p78','p79','p80','p81','p82','p83','p84','p85','p86'\
                         ,'p87','p88','p89','p90','p91','p92','p93','p94','p95','p96','p97','p98','p99','p100','p101','p102','p103','p104'\
                        ,'p105','p106','p107','p108','p109','p120','p121','p122','p123','p124','p125','p126','p127','p128','strand'])

    
    neverr = never = True
    
    if samfile.is_valid_reference_name(chrom):
        for pileupcolumn in samfile.pileup(chrom):
            coverage += 1
            if not silence:
                if (pileupcolumn.pos % 2000000 == 1):
                    print("CHH %s s %s w %s %s pos %s Result %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),filename,w,chrom,pileupcolumn.pos,ResultPW.shape[0]))

            # forward
            if fastafile.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos+1)=='C' and fastafile.fetch(chrom,pileupcolumn.pos+1,pileupcolumn.pos+2)!='G' and fastafile.fetch(chrom,pileupcolumn.pos+2,pileupcolumn.pos+3)!='G':        
                cov_context += 1
                temp = pd.DataFrame(columns=['Qname',pileupcolumn.pos+1])
                pileupcolumn.set_min_base_quality(0)
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and not pileupread.alignment.is_reverse:  # C
                        d = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos+1: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                        df2 = pd.DataFrame(data=d)
                        #df2.head()
                        temp=temp.append(df2, ignore_index=True)
                if (not temp.empty):
                    #temp.head()
                    aggreC = pd.merge(aggreC,temp,how='outer',on=['Qname'])
                    aggreC = aggreC.drop_duplicates()

            # reverse
            if pileupcolumn.pos>2:
                if fastafile.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos+1)=='G' and fastafile.fetch(chrom,pileupcolumn.pos-1,pileupcolumn.pos)!='C' and fastafile.fetch(chrom,pileupcolumn.pos-2,pileupcolumn.pos-1)!='C':        
                    cov_context += 1
                    tempr = pd.DataFrame(columns=['Qname',pileupcolumn.pos+1])
                    pileupcolumn.set_min_base_quality(0)
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.is_reverse:  # C
                            d = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos+1: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                            df2 = pd.DataFrame(data=d)
                            #df2.head()
                            tempr=tempr.append(df2, ignore_index=True)
                    if (not tempr.empty):
                        aggreR = pd.merge(aggreR,tempr,how='outer',on=['Qname'])
                        aggreR = aggreR.drop_duplicates()   

            if never and aggreC.shape[1] == (2*w):
                never = False
                aggreC = aggreC.replace(['C'],1)
                aggreC = aggreC.replace(['T'],0)
                aggreC = aggreC.replace(['A','N','G'],np.nan)
                methbin = aggreC # backup
                #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
                meth = methbin.copy()
                meth = meth.drop('Qname',axis=1)
                methtemp = meth.copy()
                # impute once if valid
                for i in range(0,meth.shape[1]-w+1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    # if eligible for imputation
                    if enough_reads(window,w,complete=False):
                        window=pd.DataFrame(data=impute(window,w))
                        ind=np.where(window.notnull().sum(axis=1)==w)[0]
                        methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
                meth = methtemp.copy()
                # compute coverage and output summary
                for i in range(0,w,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    M=(window==1).sum(axis=0)[0]
                    UM=(window==0).sum(axis=0)[0]
                    if depth>=mC:
                        toappend=outwindow(window,w=w,pos=meth.iloc[:,range(i,i+w)].columns[0],\
                                            chrom=chrom,strand='f',mC=mC)
                        ResultPW=ResultPW.append(toappend)

                aggreC = aggreC.drop(meth.columns[0:1],axis=1)
                aggreC.dropna(axis = 0, thresh=2, inplace = True)
                #total += w

            # reverse
            if neverr and aggreR.shape[1] == (2*w):
                neverr = False
                aggreR = aggreR.replace(['G'],1)
                aggreR = aggreR.replace(['A'],0)
                aggreR = aggreR.replace(['C','N','T'],np.nan)
                methbin = aggreR # backup
                #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
                meth = methbin.copy()
                meth = meth.drop('Qname',axis=1)
                methtemp = meth.copy()
                # impute once if valid
                for i in range(0,meth.shape[1]-w+1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    # if eligible for imputation
                    if enough_reads(window,w,complete=False):
                        window=pd.DataFrame(data=impute(window,w))
                        ind=np.where(window.notnull().sum(axis=1)==w)[0]
                        methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
                meth = methtemp.copy()
                # compute coverage and output summary
                for i in range(0,w,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    M=(window==1).sum(axis=0)[0]
                    UM=(window==0).sum(axis=0)[0]
                    if depth>=mC:
                        toappend=outwindow(window,w=w,pos=meth.iloc[:,range(i,i+w)].columns[0],\
                                            chrom=chrom,strand='r',mC=mC)
                        ResultPW=ResultPW.append(toappend)
                aggreR = aggreR.drop(meth.columns[0:1],axis=1)
                aggreR.dropna(axis = 0, thresh=2, inplace = True)

            #------------------
            #  SECONDARY CASE
            #------------------

            if (aggreC.shape[1] == (3*w-1)):
                aggreC = aggreC.replace(['C'],1)
                aggreC = aggreC.replace(['T'],0)
                aggreC = aggreC.replace(['N','G','A'],np.nan)
                methbin = aggreC # backup
                #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
                meth = methbin.copy()
                meth = meth.drop('Qname',axis=1)
                methtemp = meth.copy()
                # impute once if valid
                for i in range(0,meth.shape[1]-w+1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    # if eligible for imputation
                    if enough_reads(window,w,complete=False):
                        window=pd.DataFrame(data=impute(window,w))
                        ind=np.where(window.notnull().sum(axis=1)==w)[0]
                        methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
                meth = methtemp.copy()        
                # compute coverage and output summary
                for i in range(w-1,2*w-1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    M=(window==1).sum(axis=0)[0]
                    UM=(window==0).sum(axis=0)[0]
                    if depth>=mC:
                        toappend=outwindow(window,w=w,pos=meth.iloc[:,range(i,i+w)].columns[0],\
                                            chrom=chrom,strand='f',mC=mC)
                        ResultPW=ResultPW.append(toappend)

                        if ResultPW.shape[0] % 100000 == 1:   
                            ResultPW.to_csv(r"MeHdata/CHH_%s_%s.csv"%(filename,chrom),index = False, header=True)
                            if not silence: 
                                print("Checkpoint CHH. For file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos))

                aggreC = aggreC.drop(meth.columns[0:w],axis=1)
                aggreC.dropna(axis = 0, thresh=2, inplace = True)
                #print(aggreC)
                #total += w

            if (aggreR.shape[1] == (3*w-1)):
                aggreR = aggreR.replace(['G'],1)
                aggreR = aggreR.replace(['A'],0)
                aggreR = aggreR.replace(['N','T','C'],np.nan)
                methbin = aggreR # backup
                #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
                meth = methbin.copy()
                meth = meth.drop('Qname',axis=1)
                methtemp = meth.copy()
                # impute once if valid
                for i in range(0,meth.shape[1]-w+1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    # if eligible for imputation
                    if enough_reads(window,w,complete=False):
                        window=pd.DataFrame(data=impute(window,w))
                        ind=np.where(window.notnull().sum(axis=1)==w)[0]
                        methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
                meth = methtemp.copy()        
                # compute coverage and output summary
                for i in range(w-1,2*w-1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    M=(window==1).sum(axis=0)[0]
                    UM=(window==0).sum(axis=0)[0]
                    if depth>=mC:
                        toappend=outwindow(window,w=w,pos=meth.iloc[:,range(i,i+w)].columns[0],\
                                            chrom=chrom,strand='r',mC=mC)
                        ResultPW=ResultPW.append(toappend)
                        if ResultPW.shape[0] % 100000 == 1:
                            ResultPW.to_csv(r"MeHdata/CHH_%s_%s.csv"%(filename,chrom),index = False, header=True)
                            if not silence: 
                                print("Checkpoint CHH. For file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos))

                aggreR = aggreR.drop(meth.columns[0:w],axis=1)
                aggreR.dropna(axis = 0, thresh=2, inplace = True)
                #print(aggreC)
                #total += w 

        if ResultPW.shape[0]>0:
            ResultPW.to_csv(r"MeHdata/CHH_%s_%s.csv"%(filename,chrom),index = False, header=True)
        return sample, coverage, cov_context, 'CHH'                        
        print("Done CHH for file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos))

def CHGgenome_scr(bamfile,chrom,w,fa,mC=4,silence=False):
    filename, file_extension = os.path.splitext(bamfile)
    coverage = cov_context = 0
    # load bamfile
    samfile = pysam.AlignmentFile("MeHdata/%s.bam" % (filename), "rb")
    # load reference genome
    fastafile = pysam.FastaFile('MeHdata/%s.fa' % fa)
    
    # initialise data frame for genome screening (load C from bam file)
    aggreR = aggreC = pd.DataFrame(columns=['Qname'])
    
    # if user wants to output compositions of methylation patterns at every eligible window, initialise data frame
    if w==3:
        ResultPW = pd.DataFrame(columns=\
                        ['chrom','pos','M','UM','p01','p02','p03','p04','p05','p06','p07','p08','strand'])
    if w==4:
        ResultPW = pd.DataFrame(columns=\
                        ['chrom','pos','M','UM','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11',\
                         'p12','p13','p14','p15','p16','strand'])
    if w==5:
        ResultPW = pd.DataFrame(columns=\
                        ['chrom','pos','M','UM','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','strand'])
    if w==6:
        ResultPW = pd.DataFrame(columns=\
                        ['chrom','pos','M','UM','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','p33','p34','p35','p36','p37','p38','p39','p40','p41','p42','p43','p44','p45','p46'\
                         ,'p47','p48','p49','p50','p51','p52','p53','p54','p55','p56','p57','p58','p59','p60','p61','p62','p63','p64'\
                         ,'strand'])
    if w==7:
        ResultPW = pd.DataFrame(columns=\
                        ['chrom','pos','M','UM','p01','p02','p03','p04','p05','p06','p07','p08','p09','p10','p11','p12','p13','p14','p15','p16'\
                        ,'p17','p18','p19','p20','p21','p22','p23','p24','p25','p26','p27','p28',\
                        'p29','p30','p31','p32','p33','p34','p35','p36','p37','p38','p39','p40','p41','p42','p43','p44','p45','p46'\
                         ,'p47','p48','p49','p50','p51','p52','p53','p54','p55','p56','p57','p58','p59','p60','p61','p62','p63','p64'\
                         ,'p65','p66','p67','p68','p69','p70','p71','p72','p73','p74','p75','p76','p77','p78','p79','p80','p81','p82','p83','p84','p85','p86'\
                         ,'p87','p88','p89','p90','p91','p92','p93','p94','p95','p96','p97','p98','p99','p100','p101','p102','p103','p104'\
                        ,'p105','p106','p107','p108','p109','p120','p121','p122','p123','p124','p125','p126','p127','p128','strand'])

    
    neverr = never = True
    
    start=datetime.datetime.now()
    
    if samfile.is_valid_reference_name(chrom):
        for pileupcolumn in samfile.pileup(chrom):
            coverage += 1
            #chrom = pileupcolumn.reference_name
            if not silence:
                if (pileupcolumn.pos % 2000000 == 1):
                    print("CHG %s s %s w %s %s pos %s Result %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),filename,w,chrom,pileupcolumn.pos,ResultPW.shape[0]))

            if fastafile.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos+1)=='C' and fastafile.fetch(chrom,pileupcolumn.pos+1,pileupcolumn.pos+2)!='G' and fastafile.fetch(chrom,pileupcolumn.pos+2,pileupcolumn.pos+3)=='G':        
                cov_context += 1
                temp = pd.DataFrame(columns=['Qname',pileupcolumn.pos+1])
                pileupcolumn.set_min_base_quality(0)
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and not pileupread.alignment.is_reverse:  # C
                        d = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos+1: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                        df2 = pd.DataFrame(data=d)
                        #df2.head()
                        temp=temp.append(df2, ignore_index=True)

                if (not temp.empty):
                    #temp.head()
                    aggreC = pd.merge(aggreC,temp,how='outer',on=['Qname'])
                    aggreC = aggreC.drop_duplicates()

            # reverse
            if pileupcolumn.pos>2:
                if fastafile.fetch(chrom,pileupcolumn.pos,pileupcolumn.pos+1)=='G' and fastafile.fetch(chrom,pileupcolumn.pos-1,pileupcolumn.pos)!='C' and fastafile.fetch(chrom,pileupcolumn.pos-2,pileupcolumn.pos-1)=='C':        
                    cov_context += 1
                    tempr = pd.DataFrame(columns=['Qname',pileupcolumn.pos+1])
                    pileupcolumn.set_min_base_quality(0)
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.is_reverse:  # G
                            dr = {'Qname': [pileupread.alignment.query_name], pileupcolumn.pos+1: [pileupread.alignment.query_sequence[pileupread.query_position]]}
                            df2r = pd.DataFrame(data=dr)
                            #df2.head()
                            tempr=tempr.append(df2r, ignore_index=True)

                    if (not tempr.empty):
                        #temp.head()
                        aggreR = pd.merge(aggreR,tempr,how='outer',on=['Qname'])
                        aggreR = aggreR.drop_duplicates()        

            if never and aggreC.shape[1] == (2*w):
                never = False
                aggreC = aggreC.replace(['C'],1)
                aggreC = aggreC.replace(['T'],0)
                aggreC = aggreC.replace(['A','G','N'],np.nan)
                methbin = aggreC # backup
                #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
                meth = methbin.copy()
                meth = meth.drop('Qname',axis=1)
                methtemp = meth.copy()
                # impute once if valid
                for i in range(0,meth.shape[1]-w+1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    # if eligible for imputation
                    if enough_reads(window,w,complete=False):
                        window=pd.DataFrame(data=impute(window,w))
                        ind=np.where(window.notnull().sum(axis=1)==w)[0]
                        methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
                meth = methtemp.copy()
                # compute coverage and output summary
                for i in range(0,w,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    M=(window==1).sum(axis=0)[0]
                    UM=(window==0).sum(axis=0)[0]
                    if depth>=mC:
                        toappend=outwindow(window,w=w,pos=meth.iloc[:,range(i,i+w)].columns[0],\
                                            chrom=chrom,strand='f',mC=mC)
                        ResultPW=ResultPW.append(toappend)

                aggreC = aggreC.drop(meth.columns[0:1],axis=1)
                aggreC.dropna(axis = 0, thresh=2, inplace = True)
                #total += w

            # reverse
            if neverr and aggreR.shape[1] == (2*w):
                neverr = False
                aggreR = aggreR.replace(['G'],1)
                aggreR = aggreR.replace(['A'],0)
                aggreR = aggreR.replace(['N','C','T'],np.nan)
                methbin = aggreR # backup
                meth = methbin.copy()
                meth = meth.drop('Qname',axis=1)
                methtemp = meth.copy()
                # impute once if valid
                for i in range(0,meth.shape[1]-w+1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    # if eligible for imputation
                    if enough_reads(window,w,complete=False):
                        window=pd.DataFrame(data=impute(window,w))
                        ind=np.where(window.notnull().sum(axis=1)==w)[0]
                        methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
                meth = methtemp.copy()
                # compute coverage and output summary
                for i in range(0,w,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    M=(window==1).sum(axis=0)[0]
                    UM=(window==0).sum(axis=0)[0]
                    if depth>=mC:
                        toappend=outwindow(window,w=w,pos=meth.iloc[:,range(i,i+w)].columns[0],\
                                            chrom=chrom,strand='r',mC=mC)
                        ResultPW=ResultPW.append(toappend)

                aggreR = aggreR.drop(meth.columns[0:1],axis=1)
                aggreR.dropna(axis = 0, thresh=2, inplace = True)
                #total += w
            #------------------
            #  SECONDARY CASE
            #------------------

            if (aggreC.shape[1] == (3*w-1)):
                aggreC = aggreC.replace(['C'],1)
                aggreC = aggreC.replace(['T'],0)
                aggreC = aggreC.replace(['N','A','G'],np.nan)
                methbin = aggreC # backup
                #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
                meth = methbin.copy()
                meth = meth.drop('Qname',axis=1)
                methtemp = meth.copy()
                # impute once if valid
                for i in range(0,meth.shape[1]-w+1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    # if eligible for imputation
                    if enough_reads(window,w,complete=False):
                        window=pd.DataFrame(data=impute(window,w))
                        ind=np.where(window.notnull().sum(axis=1)==w)[0]
                        methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
                meth = methtemp.copy()        
                # compute coverage and output summary
                for i in range(w-1,2*w-1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    M=(window==1).sum(axis=0)[0]
                    UM=(window==0).sum(axis=0)[0]
                    if depth>=mC:
                        toappend=outwindow(window,w=w,pos=meth.iloc[:,range(i,i+w)].columns[0],\
                                            chrom=chrom,strand='f',mC=mC)
                        ResultPW=ResultPW.append(toappend)

                        if ResultPW.shape[0] % 100000:
                            ResultPW.to_csv(r"MeHdata/CHG_%s_%s.csv"%(filename,chrom),index = False, header=True)
                            if not silence: 
                                print("Checkpoint CHG. For file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos+1))

                aggreC = aggreC.drop(meth.columns[0:w],axis=1)
                aggreC.dropna(axis = 0, thresh=2, inplace = True)
                #print(aggreC)
                #total += w
            # reverse
            if (aggreR.shape[1] == (3*w-1)):
                aggreR = aggreR.replace(['G'],1)
                aggreR = aggreR.replace(['A'],0)
                aggreR = aggreR.replace(['N','T','C'],np.nan)
                methbin = aggreR # backup
                #meth = methbin.iloc[:,methbin.columns!='Qname'] # pd to np
                meth = methbin.copy()
                meth = meth.drop('Qname',axis=1)
                methtemp = meth.copy()
                # impute once if valid
                for i in range(0,meth.shape[1]-w+1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    # if eligible for imputation
                    if enough_reads(window,w,complete=False):
                        window=pd.DataFrame(data=impute(window,w))
                        ind=np.where(window.notnull().sum(axis=1)==w)[0]
                        methtemp.loc[methtemp.iloc[ind,:].index,meth.iloc[:,range(i,i+w)].columns]=window.loc[ind,:].values
                meth = methtemp.copy()        
                # compute coverage and output summary
                for i in range(w-1,2*w-1,1):
                    window = meth.iloc[:,range(i,i+w)].values
                    M=(window==1).sum(axis=0)[0]
                    UM=(window==0).sum(axis=0)[0]
                    if depth>=mC:
                        toappend=outwindow(window,w=w,pos=meth.iloc[:,range(i,i+w)].columns[0],\
                                            chrom=chrom,strand='r',mC=mC)
                        ResultPW=ResultPW.append(toappend)

                        if ResultPW.shape[0] % 100000 == 1:   
                            ResultPW.to_csv(r"MeHdata/CHG_%s_%s.csv"%(filename,chrom),index = False, header=True)
                            if not silence: 
                                print("Checkpoint CHG. For file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos+1))

                aggreR = aggreR.drop(meth.columns[0:w],axis=1)
                aggreR.dropna(axis = 0, thresh=2, inplace = True) 

        if ResultPW.shape[0]>0:   
            ResultPW.to_csv(r"MeHdata/CHG_%s_%s.csv"%(filename,chrom),index = False, header=True)

        return filename, coverage, cov_context, 'CHG'
        print("Done CHG for file %s: %s results obtained up to position chr %s: %s." % (filename,ResultPW.shape[0],chrom,pileupcolumn.pos+1))

        
def split_bam(samplenames,Folder): 
    # get bam size
    spbam_list = []
    bamfile = samplenames + '.bam'
    statinfo_out = os.stat(Folder+bamfile)
    bamsize = statinfo_out.st_size
    samfile = pysam.Samfile(Folder+bamfile, "rb")
    fileout_base = os.path.splitext(bamfile)[0] # filename
    ext = '.bam'
    x = 0
    fileout = Folder+fileout_base+"_" + str(x)+ext # filename_x.bam
    print("fileout",fileout)
    header = samfile.header
        
    outfile = pysam.Samfile(fileout, "wb", header = header)
    sum_Outfile_Size=0
    for reads in samfile.fetch():
        outfile.write(reads)
        statinfo_out = os.stat(fileout)
        outfile_Size = statinfo_out.st_size
        if(outfile_Size >=337374182 and sum_Outfile_Size <= bamsize):
            sum_Outfile_Size = sum_Outfile_Size + outfile_Size
            x = x + 1
            spbam_list.append(fileout_base + "_" + str(x)+ext)
            outfile.close()
            pysam.index(fileout)
            fileout = Folder+fileout_base + "_" + str(x)+ext
            print("fileout",fileout)
            outfile = pysam.Samfile(fileout, "wb",header = header)
            
    outfile.close()
    pysam.index(fileout)



import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-w", "--windowsize",type=int, default=4 ,help='number of CGs')
parser.add_argument("-c", "--cores",type=int, default=4, help='number of cores')
parser.add_argument("--CG", default=False, action='store_true', help='Include genomic context CG')
parser.add_argument("--CHG", default=False, action='store_true', help='Include genomic context CHG')
parser.add_argument("--CHH", default=False, action='store_true', help='Include genomic context CHH')
parser.add_argument("-mC", "--mindepth",type=int, default=4, help='Minimum depth per cytosine')


args = parser.parse_args()

import sys
import time
import os
import pandas as pd
import multiprocessing
from joblib import Parallel, delayed

#num_cores = multiprocessing.cpu_count()
                                                
if __name__ == "__main__":
    
    #open_log('MeHscreening.log')
    
    #start = time.time()
    Folder = 'MeHdata/'

    files = os.listdir(Folder)
    bam_list = []
    # all samples' bam files
    for file in files: 
        filename, file_extension = os.path.splitext(file)
        if file_extension == '.fa':
            fa = filename
        if file_extension == '.bam':
            bam_list.append(filename)
   
    fastafile = pysam.FastaFile('%s%s.fa' % (Folder,fa)) 
    chromosomes=[]
    for chrom in fastafile.references:
        chromosomes.append(chrom)
    fastafile.close()    
    
    topp = pd.DataFrame(columns=['sample','coverage','context_coverage','context'])    
    #CG = []
    #start=t.time()
    if args.CG:
        con='CG'
        CG=Parallel(n_jobs=args.cores)(delayed(CGgenome_scr)(bam,chrom=c,w=args.windowsize,fa=fa,optional=args.opt,mC=args.mindepth) for bam in bam_list for c in chromosomes)
        
        if args.opt:
            for file in bam_list:
                for c in chromosomes:
                    res_dir = Folder + con + '_'+ file + '.csv'
                    toapp_dir = Folder + con + '_'+ file + '_'+ c + '.csv'
                    if os.path.exists(res_dir):
                        Tomod = pd.read_csv(res_dir) 
                        Toappend = pd.read_csv(toapp_dir)
                        Tomod = Tomod.append(Toappend)
                        Tomod.to_csv(res_dir,index = False,header=True)
                        os.remove(toapp_dir)
                    else:
                        Toappend = pd.read_csv(toapp_dir)
                        Toappend.to_csv(res_dir,index = False,header=True)
                        os.remove(toapp_dir)
 
        
        logm("All done. "+str(len(bam_list))+" bam files processed and merged for CG.")
        for i in CG:
            toout=pd.DataFrame({'sample':i[0],'coverage':i[1],'context_coverage':i[2],'context':i[3]},index=[0])
            topp=topp.append(toout)
        #topp.groupby(['context','sample']).agg({'coverage': 'sum', 'context_coverage': 'sum'})
        #print(topp)
          
    if args.CHG:
        con='CHG'
        CG=Parallel(n_jobs=args.cores)(delayed(CHGgenome_scr)(bam,chrom=c,w=args.windowsize,fa=fa,optional=args.opt,mC=args.mindepth) for bam in bam_list for c in chromosomes)
        
        logm("Merging within samples for CHG.")  
        # not into bins of 400bp
        if args.opt:
            for file in bam_list:
                for c in chromosomes:
                    res_dir = Folder + con + '_'+ file + '.csv'
                    toapp_dir = Folder + con + '_'+ file + '_'+ c + '.csv'
                    if os.path.exists(res_dir):
                        Tomod = pd.read_csv(res_dir) 
                        Toappend = pd.read_csv(toapp_dir)
                        Tomod = Tomod.append(Toappend)
                        Tomod.to_csv(res_dir,index = False,header=True)
                        os.remove(toapp_dir)
                    else:
                        Toappend = pd.read_csv(toapp_dir)
                        Toappend.to_csv(res_dir,index = False,header=True)
                        os.remove(toapp_dir)
 
 
        logm("All done. "+str(len(bam_list))+" bam files processed and merged for CHG.")
        
        for i in CG:
            toout=pd.DataFrame({'sample':i[0],'coverage':i[1],'context_coverage':i[2],'context':i[3]},index=[0])
            topp=topp.append(toout)
        #topp.groupby(['context','sample']).agg({'coverage': 'sum', 'context_coverage': 'sum'})
        
    if args.CHH:
        con='CHH'
        CG=Parallel(n_jobs=args.cores)(delayed(CHHgenome_scr)(bam,chrom=c,w=args.windowsize,fa=fa,optional=args.opt,mC=args.mindepth) for bam in bam_list for c in chromosomes)
        
        if args.opt:
            for file in bam_list:
                for c in chromosomes:
                    res_dir = Folder + con + '_'+ file + '.csv'
                    toapp_dir = Folder + con + '_'+ file + '_'+ c + '.csv'
                    if os.path.exists(res_dir):
                        Tomod = pd.read_csv(res_dir) 
                        Toappend = pd.read_csv(toapp_dir)
                        Tomod = Tomod.append(Toappend)
                        Tomod.to_csv(res_dir,index = False,header=True)
                        os.remove(toapp_dir)
                    else:
                        Toappend = pd.read_csv(toapp_dir)
                        Toappend.to_csv(res_dir,index = False,header=True)
                        os.remove(toapp_dir)
 
        
        print("All done.",len(bam_list),"bam files processed and merged for CHH.")    
        logm("All done. "+str(len(bam_list))+" bam files processed and merged for CHH.")
        
        for i in CG:
            toout=pd.DataFrame({'sample':i[0],'coverage':i[1],'context_coverage':i[2],'context':i[3]},index=[0])
            topp=topp.append(toout)

    topp=topp.groupby(['context','sample']).agg({'context_coverage': 'sum', 'coverage': 'sum'}).reset_index()
    
    end = time.time()
    
    for i in range(topp.shape[0]):
        #print('i = ',i)
        print('Sample', topp.iloc[i,1],'has coverage',topp.iloc[i,2],'for context',topp.iloc[i,0],'out of data coverage',topp.iloc[i,3])
        logm('Sample '+str(topp.iloc[i,1])+' has coverage '+str(topp.iloc[i,2])+' for context '+str(topp.iloc[i,0])+' out of data coverage '+str(topp.iloc[i,3])+ '.')


# Imp
# /MH/OTU5/
# python impmlonly.py -w 4 -c 60 --CG --opt
# only for CG
# FINAL FINAL Nov3
