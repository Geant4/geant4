#!/usr/bin/env python

import os, sys
from ROOT import gROOT, TFile, TCanvas, TGraph, Riostream

   gROOT.Reset()

n = 1100

En, Y1, Y2, diff, pres = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )

c1 = TCanvas('c1', 'c1',6,6,800,600)
gStyle.SetOptStat(0)
c1.SetLogx()
c1.SetFillColor(0)
c1.SetBorderMode(0)
c1.SetBorderSize(0)
c1.SetFrameBorderMode(0)
   

in = open('DEDX.hIoni.proton.asc_20bins_spl.out', 'r')
#in = open('DEDX.eIoni.e-.asc_spl.out, 'r')

lines = in.readlines()
title = lines[0]

for line in lines[1:]:
     En.append (line.split() [0])
     Y1.append (line.split() [1])
     Y2.append (line.split() [2])
     diff.append (line.split() [3])
     pres.append ((Y1/Y2 - 1)*100) 

gr = new TGraph(n,En,pres)
gr.SetTitle('proton in Pb, table 100 pts / 20 pts_spl per order')
gr.SetMarkerStyle(22)
gr.SetMarkerSize(0.8)
gr.GetYaxis()->SetLabelFont(132)
gr.GetYaxis()->SetLabelSize(0.04)
gr.GetXaxis()->SetLabelFont(132)
gr.GetXaxis()->SetLabelSize(0.04)
gr.GetXaxis()->SetTitle('E, MeV')
gr.GetYaxis()->SetTitle('dE/dx: (100pts/20pts_spl - 1), %')
gr.GetXaxis()->SetTitleOffset(1.2)
gr.Draw('AP')


c1.Modified()
c1.cd()
