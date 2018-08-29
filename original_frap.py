# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 16:52:22 2018

@author: Minja
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 12:58:00 2018

@author: Minja
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import math
from scipy.optimize import curve_fit
f1 = 'C:/Users/Minja/Desktop/python petnica/CLSM Proba/Minja Radna verzija 1/Podaci 1/'
f2 = 'C:/Users/Minja/Desktop/python petnica/CLSM Proba/Minja Radna verzija 2/Podaci 2/'
f3 = 'C:/Users/Minja/Desktop/python petnica/CLSM Proba/Minja Radna verzija 3/Podaci 3/'
f4 = 'C:/Users/Minja/Desktop/python petnica/CLSM Proba/Minja Radna verzija 4/Podaci 4/'

def extract_data(fname):
    x = []
    y = []
    cnt = 0
    hlines = 10
    with open(fname,'r') as inf:
        for line in inf:
            if cnt < hlines:
                cnt += 1
                continue
            else:
                parts = line.rstrip().split(',')
                x.append(float(parts[1]))
                y.append(float(parts[2]))
            
    return np.array(x), np.array(y)

fp = [f1, f2, f3, f4]

for i, fname in enumerate(fp):

    plt.figure(figsize=(20,10))

    x,y = extract_data(fname + 'Background.csv')
    plt.plot(x,y,label='Background')

    x,y = extract_data(fp[i] + 'Baseline celija.csv')
    plt.plot(x,y,label='Baseline celija')

    x,y = extract_data(fp[i] + 'Bleachovana celija.csv')
    plt.plot(x,y,label='Bleachovana celija')

    x,y = extract_data(fp[i] + 'ROI.csv')
    plt.plot(x,y,label='ROI')

    plt.title('Prva verzija')

    plt.legend()
    #plt.savefig('prva_verzija.png')


    #Background subtr 1 verzija
    plt.figure(figsize=(20,10))
    x,y1 = extract_data(fp[i] + 'ROI.csv')
    x,y2 = extract_data(fp[i] + 'Background.csv')
    x,y3 = extract_data(fp[i] + 'Bleachovana celija.csv')
    y0 = np.subtract(y1,y2)
    t = np.subtract(y3,y2)
    plt.plot(x,y0,label='ROI(t)-BG(t)')
    plt.plot(x,t,label='TOT(t)-BG(t)')
    plt.legend()
    plt.title('Background subtraction')
    plt.savefig('Background subtraction.png')

    #Correction 1 verzija
    plt.figure(figsize=(20,10))
    d = np.divide(y0,t)
    plt.plot(x,d,label='Correction')
    plt.legend()
    plt.title('Correction')
    plt.savefig('Correction.png')

    #Normalization 1 verzija
    plt.figure(figsize=(20,10))
    y_t = y3[0]
    y_bg = y2[0]
    y_roi = y1[0]
    p = (y_t - y_bg)/(y_roi - y_bg)
    m = np.multiply(d,p)
    plt.plot(x,m,label='Normalization')
    plt.legend()
    plt.title('Normalization')
    plt.savefig('Normalization.png')

    # min tacka na normalizaciji 1 verzija
    dx = np.diff(m)
    loc1 = np.argmin(dx)
    local_max = m[loc1]
    x_data1 = x[loc1:]
    y_data1 = m[loc1:]

    plt.figure(figsize=(20,10))
    plt.plot(x_data1,y_data1,label='min tacke')
    plt.savefig('min_tacke_.png')

    #pretpostavljena asimptotska vrednost
    b = 1.05 * max(y_data1)


    def pretpostavljeni_parametri(x,y,b):
    
        y_p = np.log(b-y)
        #y = k * x + n
        k_p,n_p, _, _, _ = stats.linregress(x,y_p) 
        k = -k_p
        a = b - np.exp(n_p)
    
        return k,a

    k,a = pretpostavljeni_parametri(x_data1, y_data1,b)

    t = x_data1
  
    def curve(x,a,b,k):
        return a + (b-a)*(1-np.exp(-k*t))
 
    popt, pcov = curve_fit(curve, x_data1, y_data1)

    plt.plot(x_data1, y_data1, 'g-')
    plt.plot(x_data1, curve(x_data1,a,b,k), 'r-')
    plt.plot(x_data1, curve(x_data1, *popt), 'k-',
              label='fit: a=%2.3f, b=%2.3f, k=%2.3f' % tuple(popt))
    plt.savefig('eksp.png')