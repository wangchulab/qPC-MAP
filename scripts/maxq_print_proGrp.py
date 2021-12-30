#!/usr/bin/env python
import sys
from math import fabs, isnan
from scipy import stats
import numpy as np
from numpy import log2, log10
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy import interpolate
import matplotlib
matplotlib.use('Agg')
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt

its = "Intensity "
lfq = "LFQ intensity "
rat = "Ratio H/L normalized "

scale = 1.5
if len(sys.argv)<5:
  ratio_cutoff = 1.5
else:
  ratio_cutoff = float(sys.argv[4])

#load fasta to get kD data
lib_fn = sys.argv[1]
map_mass_kD = {}
for record in SeqIO.parse(lib_fn, "fasta"):
  try:
    analysed_seq = ProteinAnalysis(str(record.seq))
    map_mass_kD[record.id] = int(analysed_seq.molecular_weight()/1000)
    #print(record.id, map_mass_kD[record.id])
  except:
    #print "Warning:", record.id, "has Z"
    pass

#markers
markers = open( sys.argv[3], 'r' ).readlines()
ms = []
ls = []
for marker in markers:
  es = marker.split()
  m = float(es[0])
  l = float(es[1])
  ms.append(m)
  ls.append(l)
ms = np.array(ms)
ls = np.array(ls)

def decay( x, p ):
  return p[0] / ( x + p[1] ) ** scale + p[2]

def residuals( p, y, x ):
  return (y - decay(x, p))

#fit N_frac, kD function
p0 = [ 5000.0, 0.0, 0.0 ]
plsq = leastsq(residuals,p0,args=(ls, ms))
a = plsq[0][0]
b = plsq[0][1]
c = plsq[0][2]
#print a, b, c
#B-spline
#tck = interpolate.splrep(ls, ms)
#print tck

#fraction list
experiment_lst = []
map_frac_mass = {} #start from 1
frac_xtics = []
lines = open(sys.argv[2], 'r').readlines()
for n, l in enumerate(lines):
  es = l.strip().split()
  experiment_lst.append(es[0]) 
  l1 = float( es[1] )
  l2 = float( es[2] )
  m1 = a / (b+l1) ** scale + c
  m2 = a / (b+l2) ** scale + c
  if l1<0: m1 = 1000.0
  map_frac_mass[n+1] = ( m2, m1 )
  frac_xtics.append( str(int(m2)) + '-' + str(int(m1)) )
#print map_frac_mass

#fit gaussian?
def gaussian( x, *params ):
  N = len( params ) / 3
  #print N
  sgau = np.zeros(x.size)
  for i in xrange(N):
    a = params[i*3]
    x0 = params[i*3+1]
    sigma = params[i*3+2]
    #print type(a), type(x0), type(sigma)
    #print type(x[0])
    sgau += np.exp( -(x-x0)**2/(2*sigma**2) )
  return sgau

########
#results
########
lines = open("proteinGroups.txt", 'r').readlines()
#tags
tags = {}
es = lines[0].strip().split('\t')
for i, e in enumerate(es):
  tags[e] = i

pdf = matplotlib.backends.backend_pdf.PdfPages("outputs.pdf")

def process(name, Ls, Hs, Rs):
  # real mass
  Nfrac = len(Rs)
  mass_lst = []
  for p in name.split(';'):
    if p in map_mass_kD.keys():
      mass_lst.append( map_mass_kD[p] )
  if len(mass_lst) == 0: return
  real_mass = np.median(mass_lst)
  real_mark = 0
  for i in xrange(1, Nfrac+1):
    m1, m2 = map_frac_mass[i]
    if m2>=real_mass and m1<=real_mass: #m1<r<m2
      real_mark = i
      break
  if real_mark == 0:
    if real_mass < 80: real_mark=16
    if real_mass > 200: real_mark=1
  print "real_mark:", real_mark, real_mass

  p1_l = 0
  p2_l = 0
  p1_h = 0
  p2_h = 0
  params_l = []
  params_h = []
  x = []
  yl = []
  yh = []
  sumA_l = 0.0
  sumA_h = 0.0
  for n, (l, h, r) in enumerate(zip(Ls, Hs, Rs)):
    if p2_l <= p1_l and p1_l > l:
      #print "peakL:", n, p2_l, p1_l, l
      A_l = p2_l + p1_l + l
      params_l.append( (n, p1_l, 1.0, A_l) )
      sumA_l += A_l
    if p2_h <= p1_h and p1_h > h:
      #print "peakH:", n, p2_h, p1_h, h
      A_h = p2_h + p1_h + h
      params_h.append( (n, p1_h, 1.0, A_h) )
      sumA_h += A_h
    if l > p1_l and n+1 == len(Ls):
      #print "peakL:", n+1, p1_l, l, 0
      A_l = p1_l + l
      params_l.append( (n+1, l, 1.0, A_l ) )
      sumA_l += A_l
    #else:
      #print "db:", n, len(Ls), p1_l, l
    if h > p1_h and n+1 == len(Hs):
      #print "peakH:", n+1, p1_h, h, 0
      A_h = p1_h + h
      params_h.append( (n+1, h, 1.0, A_h ) )
      sumA_h += A_h
    p2_l = p1_l
    p1_l = l
    p2_h = p1_h
    p1_h = h
    x.append(n+1)
    yl.append(l)
    yh.append(h)
  #print params_l
  #print params_h

  real_peak_ndx_l = 0
  real_peak_ndx_h = 0
  dM_l = 1000
  dM_h = 1000
  rate_l = 0.0
  rate_h = 0.0
  p_l = 0.0
  p_h = 0.0
  print "#", name, real_mass
  for n, (l, h, r) in enumerate(zip(Ls, Hs, Rs)):
    print n+1, l, h, r

  print "#L:"
  for p in params_l:
    if fabs( map_frac_mass[p[0]][0] - real_mass ) < dM_l:
      dM_l = fabs( map_frac_mass[p[0]][0] - real_mass )
      real_peak_ndx_l = p[0]
      rate_l = real_mass / map_frac_mass[p[0]][0]
      p_l = p[3]/sumA_l
    if fabs( map_frac_mass[p[0]][1] - real_mass ) < dM_l:
      dM_l = fabs( map_frac_mass[p[0]][1] - real_mass )
      real_peak_ndx_l = p[0]
      rate_l = real_mass / map_frac_mass[p[0]][1]
      p_l = p[3]/sumA_l
    print "#", p[0], p[3]/sumA_l, "M=", map_frac_mass[p[0]]
  print "#H:"
  for p in params_h:
    if fabs( map_frac_mass[p[0]][0] - real_mass ) < dM_h:
      dM_h = fabs( map_frac_mass[p[0]][0] - real_mass )
      real_peak_ndx_h = p[0]
      rate_h = real_mass / map_frac_mass[p[0]][0]
      p_h = p[3]/sumA_h
    if fabs( map_frac_mass[p[0]][1] - real_mass ) < dM_h:
      dM_h = fabs( map_frac_mass[p[0]][1] - real_mass )
      real_peak_ndx_h = p[0]
      rate_h = real_mass / map_frac_mass[p[0]][1]
      p_h = p[3]/sumA_h
    print "#", p[0], p[3]/sumA_h, "M=", map_frac_mass[p[0]]
  if real_peak_ndx_l>0:
    print "realL", name, real_peak_ndx_l, map_frac_mass[real_peak_ndx_l], real_mass, rate_l, p_l
  if real_peak_ndx_h>0:
    print "realH", name, real_peak_ndx_h, map_frac_mass[real_peak_ndx_h], real_mass, rate_h, p_h

  #fit gaussian peaks
  #print params
  #popt, pcov = curve_fit( gaussian, x, y, p0=params )
  #print popt
  
  #old version
  #if np.nanmax(Rs) > ratio_cutoff:
  #  #output pdf
  #  fig = plt.figure(21)
  #  plt.xlim([0,16])
  #  #ymax = np.max( [h for (h,r) in zip(Hs, Rs) if r>ratio_cutoff] )
  #  ymax = np.max( [l for (l,r) in zip(Ls, Rs) if r>ratio_cutoff] )

  #  plt.subplot(211)
  #  plt.title(name)
  #  plt.bar( np.array(range(len(Rs)))-0.2, Ls, color="b", width=0.4, edgecolor='b', label="L" )
  #  plt.bar( np.array(range(len(Rs)))+0.2, Hs, color='r', width=0.4, edgecolor='r', label="H" )
  #  plt.subplot(212)
  #  plt.title("Mass: "+str(real_mass))
  #  plt.ylim([0, ymax*1.2])
  #  plt.bar( np.array(range(len(Rs)))-0.2, Ls, color="b", width=0.4, edgecolor='b', label="L" )
  #  plt.bar( np.array(range(len(Rs)))+0.2, Hs, color='r', width=0.4, edgecolor='r', label="H" )

  #  plt.xticks(np.arange(len(Rs)), frac_xtics, rotation=45)

  #  pdf.savefig(fig)
  #  plt.close('all')

  #new version
  xtics = np.arange(Nfrac)+1
  if np.all(np.isnan(Rs)):
    print("Warning: no intensity or ratio")
    return
  if np.nanmax(Rs) > ratio_cutoff:
    #check gap
    maxInt_l = np.argmax(Ls)
    maxInt_h = np.argmax(Hs)
    if maxInt_l != maxInt_h:
      print("Warning: L/H peaks don't match")
    mono_ndx = np.max([maxInt_l, maxInt_h])
    print("Mono peak:", mono_ndx)
    top_band = -1
    for ii in range(mono_ndx):
      if Rs[ii]>ratio_cutoff:
        top_band = ii
        print("Top band:", ii+1)
        break
    if top_band == -1:
      print("Warning: No top band found!")
      return
    with_gap = False
    for ii in range(top_band+1, mono_ndx):
      if isnan(Rs[ii]):
        print("Gap found at:", ii+1)
        with_gap = True
        break
    if not with_gap:
      print("No gap! skipping ...")
      return

    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    plt.gca().invert_yaxis()
    tnames = [ t for t in name.split(';') if "CON_" not in t ]
    tnames = [ "|".join(t.split("|")[1:]) for t in tnames ]
    for t in tnames:
      print "OUT:", t
    ax1.set_title("\n".join(tnames), loc="right")
    ax1.set_yticks(xtics)
    ax1.set_xticks([-1,0,1,2,3,4])
    ax1.axvline(0,0,Nfrac,linestyle='--',linewidth=1,c='black')
    rcut = log2(ratio_cutoff)
    #ax1.axvline(-rcut,0,Nfrac,linestyle='--',linewidth=0.5,c='black')
    ax1.axvline( rcut,0,Nfrac,linestyle='--',linewidth=0.5,c='black')
    for i in xrange(Nfrac+1):
      ax1.axhline(i+0.5, -5, 5, linestyle='--', linewidth=0.4, c='grey', alpha=0.5)
    newRs = []
    for r in Rs:
      if log2(r)>4.0: newRs.append(4.0)
      elif log2(r)<-4.0: newRs.append(-4.0)
      else:
        newRs.append(log2(r))
    ax1.scatter(newRs, xtics, c='orange')
    ax1.axhline( real_mark, -5, 5, linestyle='-', linewidth=1, c='black')
    ax1.set_xlim([-1,4.5])
    ax1.set_xlabel("$log_2$(H/L)")
    
    ax2.set_title("(MW: "+str(real_mass)+" kDa)", loc="left")
    sumL = np.sum(Ls)
    sumH = np.sum(Hs)
    #perctL = Ls / sumL * 100
    #perctH = Hs / sumH * 100
    #txtL = [ str(t)+"%" for t in perctL ]
    #ax2.plot(Ls,xtics, linewidth=2)
    #ax2.plot(Hs,xtics)
    ax2.barh( xtics-0.2, Ls, 0.4, color = "red" )
    ax2.barh( xtics+0.2, Hs, 0.4, color = "blue")
    gap = (sumL+sumH)/2.0/100
    #zoom in
    y_max_l = np.max( [l for (l,r) in zip(Ls, Rs) if r>ratio_cutoff] )
    y_max_h = np.max( [h for (h,r) in zip(Hs, Rs) if r>ratio_cutoff] )
    y_max_1 = np.max( [y_max_l, y_max_h] )
    #new_Ls = []
    #new_Hs = []
    #for l, h in zip(Ls, Hs):
    #  new_Ls.append(l/y_max_1*y_max_2/1.5)
    #  new_Hs.append(h/y_max_1*y_max_2/1.5)
    #ax2.barh( xtics-0.2, new_Ls, 0.2, color='orange' )
    #ax2.barh( xtics,     new_Hs, 0.2, color='red' )

    for i, l, h in zip( xtics, Ls, Hs ):
      pL = l / sumL * 100
      pH = h / sumH * 100
      if pL > 0.05: ax2.text( y_max_1*1.22, i-0.1, "%4.2f%%" % pL, fontsize=7, color="red" )
      if pH > 0.05: ax2.text( y_max_1*1.22, i+0.3, "%4.2f%%" % pH, fontsize=7, color="blue" )

    ax2.set_xlim([0, y_max_1*1.2])
    ax2.set_xlabel("$LFQ$")
    #plt.grid(True, axis='y', linestyle="--", alpha=0.5)
    for i in xrange(Nfrac+1):
      ax2.axhline(i+0.5, 0, y_max_l*1.2, linestyle='--', linewidth=0.4, c='grey', alpha=0.5)

    f.subplots_adjust(wspace=0)
    pdf.savefig(f)
    plt.close('all')
    print("...")

#scan
for l in lines[1:]:
  es = l.strip().split('\t')
  pros = es[tags["Protein IDs"]]
  #if sys.argv[4] not in pros: continue

  #load
  raw_intens_L = [] 
  lfq_intens_L = []
  raw_intens_H = []
  lfq_intens_H = []
  ratios = []
  sum_lfq = 0.0
  sum_its = 0.0
  for exp in experiment_lst:
    raw_intens_L.append( float(es[tags[its+"L "+exp]]) ) 
    lfq_intens_L.append( float(es[tags[lfq+"L "+exp]]) ) 
    raw_intens_H.append( float(es[tags[its+"H "+exp]]) ) 
    lfq_intens_H.append( float(es[tags[lfq+"H "+exp]]) ) 
    tmp_r = float(es[tags[rat+exp]])
    #reverse LH
    #if tmp_r > 0.01:
    #  tmp_r = 1.0 / tmp_r
    #elif tmp_r > 0.001:
    #  tmp_r = 100.0
    #use LFQ result if NaN
    if isnan(tmp_r) and lfq_intens_H[-1]>1000.0 and lfq_intens_L[-1]>1000.0:
      tmp_r = lfq_intens_H[-1] / lfq_intens_L[-1]
    elif isnan(tmp_r) and raw_intens_H[-1]>1000.0 and raw_intens_L[-1]>1000.0:
      tmp_r = raw_intens_H[-1] / raw_intens_L[-1]
    ratios.append( tmp_r )
    sum_lfq += lfq_intens_L[-1] + lfq_intens_H[-1]
    sum_its += raw_intens_L[-1] + raw_intens_H[-1]

  if sum_lfq > 10.0:
    process(pros, lfq_intens_L, lfq_intens_H, ratios)
  elif sum_its > 10.0:
    process(pros, raw_intens_L, raw_intens_H, ratios)

pdf.close()

