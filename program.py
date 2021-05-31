#The program dont work if we don have the Data
#If we work the program must in writing in terminal ths command: python3 nameofprogram.py 0 76 . The 0 is dt and 76 is final T

import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
import math
import statistics
import glob
from scipy.optimize import minimize_scalar
import sys, os
import fnmatch
from scipy.optimize import curve_fit
from scipy.stats import chisquare
#from scipy.stats import sem

#SIMPLE MASS
def func(x, a, b, c):
    return a * np.exp(-b * x) + c

#CORRELATOR FIT
def func2(x, a, m0 ,b, m1):
    return a * np.exp(-m0 * x)*( 1 + b * np.exp(-m1 * x))

#EFFECTIVE MASS FIT
def func3(x,E0,a,DE):
    return E0 + np.log((1+a*np.exp(-DE*x))/(1+a*np.exp(-DE*(x+1))))

#PLATEAU FIT
def func4(E0):
    return E0


def jackknife(dsetj):
            jackknife_s=[]
            jackknife_err=[]
            l,n=np.shape(dsetj)
            idx = np.arange(n)
            np.where(idx!=1)
            jackknife_samples=[]
            for j in range(l):
                        jacknife_samples = np.average(np.delete(dsetj,j,0),axis=0)
                        jackknife_samples.append(jacknife_samples)

            print(len(jackknife_samples[0]))
            print(len(jackknife_samples))
            j_avg=np.sum(jackknife_samples,axis=0)/l
            j_err=[np.std(jackknife_samples,axis=0)/np.sqrt((l))]
            return j_avg , j_err

def jackknife2(dsetj,dsetj2):
            jackknife_s=[]
            jackknife_err=[]
            l,n=np.shape(dsetj)
            idx = np.arange(n)
            np.where(idx!=1)
            jackknife_samples=[]
            for j in range(l):
                        jacknife_samples = np.log(np.abs(np.average(np.delete(dsetj,j,0),axis=0)))-np.log(np.abs(np.average(np.delete(dsetj2,j,0),axis=0)))
                        jackknife_samples.append(jacknife_samples)

            print(len(jackknife_samples[0]))
            print(len(jackknife_samples))
            j_avg=np.sum(jackknife_samples,axis=0)/l
            #j_err=[np.std(jackknife_samples,axis=0)]
            x=(jackknife_samples-j_avg)**2
            j_err=np.sqrt(l*np.sum(x,axis=0)/(l-1))
            return j_avg , j_err

def objective(x):
    return exp(-x)
#entries = glob.glob('/home/loukmix/jacob/data/*.h5')


print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:', str(sys.argv))
j=0
jB=0
cort=[]
c4=[]
c5=[]
k=0
kB=0
dt=int(sys.argv[1])
L=int(sys.argv[2])-int(sys.argv[1])-1

T=int(sys.argv[2])

i=0

#Reading data
fi=h5py.File('cB211A25_48.h5','r')
fi2=h5py.File('cB211B25_48.h5','r')
filist=list(fi)
filist2=list(fi2)
for k in filist:
        print(k)
        mlist=list(fi[k])
        am1p=np.zeros((L))
        am1p1=np.zeros((L))
        am1pdt=np.zeros((L))
        am1pa=np.zeros((L-1))
        am1pb=np.zeros((L-1))

        am1m=np.zeros((L))
        am1mdt=np.zeros((L))
        am1m1=np.zeros((L))
        for j in mlist:
            am1p=np.array(fi["path for data" % (k,j)])[dt+2:L+dt+1,0]
            am1pdt=np.array(fi["path for data" % (k,j)])[1:L,0]
            am1p1=np.array(fi["path for data" % (k,j)])[dt+1:L+dt,0]

            #am1m=np.array(fi["path for data" % (k,j)])[T-dt-2:0:-1,0]
            #am1mdt=np.array(fi["path for data" % (k,j)])[T+1:dt+1:-1,0]
            #am1m1=np.array(fi["%path for data" % (k,j)])[T-dt-1:1:-1,0]
            #am1m=np.flip(am1m)
            #am1mdt=np.flip(am1mdt)
            #am1m1=np.flip(am1m1)



            if dt>0:
                am1pa=am1p/am1pdt   + am1pa
            else:
                am1pa=am1p  +  am1pa

            if dt>0:
                am1pb=am1p1/am1pdt  + am1pb
            else:
                am1pb=am1p1 + am1pb
        c4.append(am1pa)
        c5.append(am1pb)

c4B=[]
c5B=[]
for kB in filist2:
        print(kB)
        mlistB=list(fi2[kB])
        am1pB=np.zeros((L))
        am1p1B=np.zeros((L))
        am1pdtB=np.zeros((L))
        am1paB=np.zeros((L-1))
        am1pbB=np.zeros((L-1))

        am1mB=np.zeros((L))
        am1mdtB=np.zeros((L))
        am1m1B=np.zeros((L))
        for jB in mlistB:
            am1pB=np.array(fi2["path for data" % (kB,jB)])[dt+2:L+dt+1,0]
            am1pdtB=np.array(fi2["path for data" % (kB,jB)])[1:L,0]
            am1p1B=np.array(fi2["%path for data" % (kB,jB)])[dt+1:L+dt,0]

            #am1m=np.array(fi["path for data" % (k,j)])[T-dt-2:0:-1,0]
            #am1mdt=np.array(fi["path for data" % (k,j)])[T+1:dt+1:-1,0]
            #am1m1=np.array(fi["%path for data" % (k,j)])[T-dt-1:1:-1,0]
            #am1m=np.flip(am1m)
            #am1mdt=np.flip(am1mdt)
            #am1m1=np.flip(am1m1)



            if dt>0:
                am1paB=am1pB/am1pdtB   + am1paB
            else:
                am1paB=am1pB  +  am1paB

            if dt>0:
                am1pbB=am1p1B/am1pdtB  + am1pbB
            else:
                am1pbB=am1p1B + am1pbB
        c4B.append(am1paB)
        c5B.append(am1pbB)

#Average of correlator without jacknife
cort=np.average(c4,axis=0)
dcort=[np.std(c4,axis=0)]
relerr=[]
rel=np.array(dcort)/np.array(cort)

ll=len(c4)
l=L-1
llB=len(c4B)
lB=L-1

#Time for jacknife
c0=np.reshape(c4,(len(c4),l))
c1=np.reshape(c5,(len(c5),l))

c0B=np.reshape(c4B,(len(c4B),lB))
c1B=np.reshape(c5B,(len(c5B),lB))

javg=[]
jerr=[]
javg2=[]
jerr2=[]
javg3=[]
jerr3=[]

javgB=[]
jerrB=[]
javg2B=[]
jerr2B=[]
javg3B=[]
jerr3B=[]

# i jacknife gia 2 config ine kanoniko avarage/!!

javg3,jerr3= jackknife2(c0,c1)
yan=np.abs(javg3)
dyan=np.abs(jerr3)
yerr=np.reshape(dyan,L-1)
print(dyan[:20])
print(javg3)
l=len(yan)
nx=np.arange(1,L,1)

javg3B,jerr3B= jackknife2(c0B,c1B)
yanB=np.abs(javg3B)
dyanB=np.abs(jerr3B)
yerrB=np.reshape(dyanB,L-1)
print(dyanB[:20])
print(javg3B)
lB=len(yanB)
nxB=np.arange(1,L,1)

#Fitting time
popt, pcov= curve_fit (func3, nx,yan )
perr=np.sqrt(np.diag(pcov))
print(popt,pcov)
print(perr)
yn=[]

yn2=[]
x0    = np.array([0.0, 0.0, 0.0])
#params = Parameters()
for x in range(1,L,1):
    yline=func3(x , *popt)
    yn.append(yline)

yn2B=[]
for x in range(tmin,tmax,1):
    yline2=plateau fit
    yn2.append(yline2)
x2=chisquare(yan,yn)
print (x2)

#If we want to do all fits is one program
#popt, pcov= curve_fit (func3, nxB,yanB )
#perr=np.sqrt(np.diag(pcov))
#print(popt,pcov)
#print(perr)

#for x in range(tmin,tmax,1):
#    yline2='plateau fit 1'
#    yn2B.append(yline2)
#yn3=[]
#for x in range(tmin,tmax,1):
#    yline2='plateau fit 2'
#    yn3.append(yline2)

#We choose what do we plot

#Data only
#plt.errorbar(np.arange(1,L), np.reshape(yan,L-1), yerr=np.reshape(dyan,L-1) , fmt='o', label='dt=0 ENS. cB211A25_48')
#plt.errorbar(np.arange(1,L), np.reshape(yanB,L-1), yerr=np.reshape(dyanB,L-1) , fmt='o', label='dt=0 ENS. cB211B25_48')
#(if we want to shift the y axis)plt.errorbar(np.arange(1,L), np.reshape(yanB,L-1)-0.002, yerr=(np.reshape(dyanB,L-1)-0.008)/2 , fmt='o', label='traditional method ENS. cB211B25_48')

#Fit for func3 and func4
#plt.plot(np.arange(16,49),yn2,label='Plateau fit Enseble cB211B25_48',linewidth='3.0')
#plt.plot(np.arange(16,49),yn2B,label='Plateau fit Enseble cB211A25_48',linewidth='3.0')
#plt.plot(np.arange(16,49),yn3,label='Plateau fit(traditional method) Enseble cB211B25_48',linewidth='3.0')
#plt.plot(np.arange(6,49),yn2)

#Fit for dt method
#plt.errorbar(xmdt,ym2dt,errymdt,fmt='o',label='enseble ...',color='orange')
#plt.errorbar(xmdt,ym2dtA,errymdtA,fmt='o',label='enseble ...',color='blue')

#plt.ylim(0.5,1)
#plt.xlim(0,50)
#plt.xlabel('t')
#plt.xlabel('dt')
plt.ylabel('meff')
plt.legend()
plt.savefig('name_picture.png')
