#!/usr/bin/env python

import os, sys
from ROOT import gROOT, gStyle, TFile, TCanvas, TGraph
from array import array

path1 = os.path.basename(sys.argv[1])

dirList=os.listdir(path1)
for fname in dirList:
     print fname
     gROOT.Reset()

#     n = 1100
     n = 35
     En, Y1, Y2, diff, pres = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )

     c1 = TCanvas('c1', 'c1',6,6,800,600)
     gStyle.SetOptStat(0)
     c1.SetLogx()
     c1.SetFillColor(0)
     c1.SetBorderMode(0)
     c1.SetBorderSize(0)
     c1.SetFrameBorderMode(0)
   
     infile = open(fname, 'r')

     lines = infile.readlines()
     title = lines[0]

     for line in lines[1:]:
          e = (float(line.split() [0]))
          y1 = (float(line.split() [1]))
          y2 = (float(line.split() [2]))
          di = (float(line.split() [3]))
          if y2 == 0:
               pr = 0
          else:
               pr = ((y1/y2 - 1)*100)     
          En.append (e)
          Y1.append (y1)
          Y2.append (y2)
          diff.append (di)
          pres.append (pr) 

     gr = TGraph(n,En,pres)
     grTit = fname + ' / 100 pts per order, Pb'
     gr.SetTitle(grTit)
     gr.SetMarkerStyle(22)
     gr.SetMarkerSize(0.8)
     gr.GetYaxis().SetLabelFont(132)
     gr.GetYaxis().SetLabelSize(0.04)
     gr.GetXaxis().SetLabelFont(132)
     gr.GetXaxis().SetLabelSize(0.04)
     gr.GetXaxis().SetTitle('E, MeV')
     gr.GetYaxis().SetTitle('(Y1/Y2 - 1), %')
     gr.GetXaxis().SetTitleOffset(1.2)
     gr.Draw('AP')

     c1.Modified()
     c1.cd()
     plNam = fname + '.png'
     c1.Print(plNam)
