import numpy as np
import math
import os

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


from random import seed
from random import random
seed(10)

PI    = 4 * math.atan2(1, 1)

def randn(m,sigma):
  r1  = random()
  r2  = random()
  while r1==0:
      r1 = random()
  value = (sigma* (( -2 * math.log(r1) )**1/2) * math.sin(2 * PI * r2)) + m
  return value


def lograndn(CV):
  sigma = math.log((CV**2)+1)
  sigma **= 1/2

  return math.exp(randn(0,sigma))


def mean_sigma(array):
    mean = np.mean(array)
    sigma = 0

    for n in array:
        sigma += (n-mean)**2;

    sigma /= len(array)
    sigma **= 1/2

    CV = sigma/mean

    return (mean,sigma,CV)



def adjust_num(array,CV,pop_size):
    for n in range(len(array)):
        val = array[n][1]* lograndn(CV)
        array[n][1] = val

    sum = np.sum([i[1] for i in array])

    for n in range(len(array)):
        val = float(array[n][1])*(float(pop_size)/sum)
        array[n][1] = int(val)
    return array




def adjust_num_d(array_d,CV,pop_size):
    for n in range(len(array_d)):
        val = array_d[n]['count']* lograndn(CV)
        array_d[n]['count'] = val

    sum = np.sum([i['count'] for i in array_d])

    for n in range(len(array_d)):
        val = float(array_d[n]['count'])*(float(pop_size)/sum)
        array_d[n]['count'] = int(val)
    return array_d


def adjust_num_d_positions(big_array_d,CV,pop_size):
    sum = 0
    for pos in range(len(big_array_d)):
        pos_sum = float(np.sum([ n['count'] for n in big_array_d[pos] ]))
        val = pos_sum * lograndn(CV)
        fold_change = val/pos_sum
        for n in big_array_d[pos]:
            n['count'] *= fold_change
        sum += pos_sum

    for pos in range(len(big_array_d)):
        for n in big_array_d[pos]:
            val = float(n['count'])*(float(pop_size)/sum)
            n['count'] = int(val)
    return big_array_d

def adjust_num_d_dms(big_array_d,CV,pop_size):
    sum = 0
    for pos in range(len(big_array_d)):
        for n in range(len(big_array_d[pos])):
            val =big_array_d[pos][n]['count']* lograndn(CV)
            big_array_d[pos][n]['count'] = val

            sum += val

    for pos in range(len(big_array_d)):
        for n in range(len(big_array_d[pos])):
            val = float(big_array_d[pos][n]['count'])*(float(pop_size)/sum)
            big_array_d[pos][n]['count'] = int(val)
    return big_array_d




    return big_array_d




def LL2csv(LL,name):
    with open(name,"w") as F:
        for L in LL:
            F.write("%s\n" %  (",").join( [str(i) for i in L]))
    # This script outputs dict as string, so do not print anything
    print("Generated : %s"%(name))
    F.close()

def hap_coverage(x,y,nm):
    all = [i[1] for i in x]
    all += [i[1] for i in y]
    mx = int(math.log(max(all),10)+2)
    cov = [["Threashold".ljust(mx+1),"x_percentage".ljust(mx+1),"y_percentage".ljust(mx+1)]]
    for th in range(0,mx):
        val_x = above_threshold([i[1] for i in x],10**th)*100
        val_y = above_threshold([i[1] for i in y],10**th)*100
        cov.append([str(10**th).ljust(mx+1),str(round(val_x,2)).ljust(mx+1),str(round(val_y,2)).ljust(mx+1)])
    return cov

def dip_coverage(dip,nm):
    all = [i[1] for i in dip]
    mx = int(math.log(max(all),10)+2)
    cov = [["Threashold".ljust(mx+1),"Percentage".ljust(mx+1)]]
    for th in range(0,mx):
        val = above_threshold([i[1] for i in dip],10**th)*100
        cov.append([str(10**th).ljust(mx+1),str(round(val,2)).ljust(mx+1)])
    return cov

def above_threshold(pop,th):
    above = [ strain for strain in pop if (strain > th)]
    value = float(len(above))/len(pop)
    return value

def print_LL(LL):
    for L in LL:
        l = [str(i) for i in L]
        print(("\t").join(l))



def dms_coverage(dms_array):
    all_l = []
    for pos in dms_array:
        for codon in pos:
            all_l.append(codon['count'])
    mx = int(math.log(max(all_l),10)+2)
    cov = [["Threashold".ljust(mx+1),"Percentage".ljust(mx+1)]]
    for th in range(0,mx):
        val = above_threshold(all_l,10**th)*100
        cov.append([str(10**th).ljust(mx+1),str(round(val,2)).ljust(mx+1)])
    return cov


def mating(hap_x,hap_y):
    dip = []
    for strain_x in hap_x:
        for strain_y in hap_y:
            xy_value = strain_x[1] * strain_y[1]
            dip.append([["%s"%(strain_x[0]),"%s"%(strain_y[0])],xy_value])
    return dip

def positive_interaction(array,CV_pos,CV_autoactivity,pos_rate,pop_size):

    strains = [i[0][0] for i in array]
    strains += [i[0][1] for i in array]

    autoactivity = {}
    for strain in strains:
        autoactivity[strain] = lograndn(CV_autoactivity)

    for i in array:
        signal = lograndn(CV_pos)
        val    = i[1] * signal
        i[1]   = val

    array.sort(key = lambda x:x[1],reverse=True)

    positives = len(array)*pos_rate
    for i in range(len(array)):
        if i<positives:
            array[i][0].append(1)
        else:
            array[i][0].append(0)

    for i in array:
        x = i[0][0]
        y = i[0][1]

        aa = 0
        aa = autoactivity[x]*autoactivity[y]

        value = i[1] * aa
        i[1] = value

    sum = np.sum([i[1] for i in array])

    for n in range(len(array)):
        val = float(array[n][1])*(float(pop_size)/sum)
        array[n][1] = val

    return array


####FIX this. Make ir like te BFG pipeline.
def dms_selection(big_array_d,CV_pos,CV_codon,lethality_rate_codon,lethality_rate_pos,pop_size):
    th_pos   =  len(big_array_d)*lethality_rate_pos
    th_codon =  len(big_array_d)*lethality_rate_codon
    POS = {}
    for pos in range(len(big_array_d)):
        POS[pos] = 1 - np.abs(1-lograndn(CV_pos))
    sum = 0
    for pos in range(len(big_array_d)):
        # Only affecting codon which aren;t WT.
        pos_sum = float(np.sum([ n['count'] for n in big_array_d[pos][1:] ]))
        val = pos_sum / (1 -  np.abs(1-lograndn(CV_codon))) / POS[pos]
        fold_change = val/pos_sum
        for n in big_array_d[pos][1:]:
            n['count'] *= fold_change
        pos_sum = float(np.sum([ n['count'] for n in big_array_d[pos] ]))
        sum += pos_sum


    for pos in range(len(big_array_d)):
        for n in big_array_d[pos]:
            val = float(n['count'])*(float(pop_size)/sum)
            n['count'] = int(val)
    return big_array_d









def abundance(array):
    dip = {}
    sum = 0
    for i in array:
        dip_nm = ("-").join([i[0][0],i[0][1]])
        try:
            dip[dip_nm]  += i[1]+1
        except KeyError:
            dip[dip_nm]  = i[1]+1
        sum += i[1]+1

    for i in dip:
        dip[i] /= sum

    return dip



def marginal(array):
    x_hap = {}
    y_hap = {}
    sum = 0

    for i in array:
        x = i[0][0]
        y = i[0][1]
        try:
            x_hap[x]  += i[1]+1
        except KeyError:
            x_hap[x]  = i[1]+1
        try:
            y_hap[y]  += i[1]+1
        except KeyError:
            y_hap[y]  = i[1]+1
        sum += i[1]+1

    for i in x_hap:
        x_hap[i] /= sum
    for i in y_hap:
        y_hap[i] /= sum

    return [x_hap[i] for i in x_hap],[y_hap[i] for i in y_hap]

def format_array(array):
    LL = [["x","y","value"]]
    for i in array:
        l = [i[0][0],i[0][1],i[1]]
        LL.append(l)
    return LL


def plot_DMS_hm(DMS_data,title):
    my_amino_acid_order = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y',"*"]
    codon_table = {
                'A': ['GCT', 'GCC', 'GCA', 'GCG'],
                'C': ['TGT', 'TGC'],
                'D': ['GAT', 'GAC'],
                'E': ['GAA', 'GAG'],
                'F': ['TTT', 'TTC'],
                'G': ['GGT', 'GGC', 'GGA', 'GGG'],
                'H': ['CAT', 'CAC'],
                'I': ['ATT', 'ATC', 'ATA'],
                'K': ['AAA', 'AAG'],
                'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
                'M': ['ATG'],
                'N': ['AAT', 'AAC'],
                'P': ['CCT', 'CCC', 'CCA', 'CCG'],
                'Q': ['CAA', 'CAG'],
                'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                'T': ['ACT', 'ACC', 'ACA', 'ACG'],
                'V': ['GTT', 'GTC', 'GTA', 'GTG'],
                'W': ['TGG'],
                'Y': ['TAT', 'TAC'],
                "*": ['TAG','TAG','TAA']
                }
    array_data = []
    for pos in DMS_data:
        array = [i['count'] for i in pos ]
        array_data.append(array)
    array_data = np.array(array_data).transpose()

    a2 = []
    for pos in DMS_data:
        for codon in pos:
            a2.append(codon['count'] )

    dms_region_len = len(DMS_data)
    #pprint(array_data)
    fig, ax = plt.subplots( figsize=(15,dms_region_len/2))
    im = ax.imshow(array_data,vmin=0,
                        vmax=max(a2)*1.5)

    ylab = []
    for aa in my_amino_acid_order:
        #print(aa)
        for codon in codon_table[aa]:
            #print(codon)
            ylab.append("%s  (%s)"%(codon,aa))

    dms_region_len = len(DMS_data)

    # We want to show all ticks...
    ax.set_yticks(np.arange(len(ylab)))
    ax.set_xticks(np.arange(dms_region_len))
    # ... and label them with the respective list entries
    ax.set_xticklabels([i for i in range(1,dms_region_len+1)])
    ax.set_yticklabels(ylab)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
             rotation_mode="anchor")

    """
    # Loop over data dimensions and create text annotations.
    for i in range(len(vegetables)):
        for j in range(len(farmers)):
            text = ax.text(j, i, harvest[i, j],
                           ha="center", va="center", color="w")
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2.5%", pad=0.5)

    plt.colorbar(im, cax=cax,ticks=range(0,int(max(a2)*1.5),int(int(max(a2)*1.5)/10)))
    ax.set_title(title)
    fig.tight_layout()
    plt.show()
