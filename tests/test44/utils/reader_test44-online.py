#!/usr/bin/env python

import os, sys
from ROOT import TROOT, gROOT, gStyle, gPad, TCanvas, TStyle
from ROOT import TFile, TH1F, TH1D, TH2F, TGraph, TLegend
from array import array


def main():

    if tyPart == 'p':
        hist_title = 'p 110 MeV in Water, Geant4  ' + refer
        fname2 = 'H-110MeV-endep-EXP-norm-max.txt'
        zmax = 120.
    elif tyPart == 'he4':
        hist_title = '^{4}He 144.3 MeV/u in Water, Geant4  ' + refer
        fname2 = '4He-144.3MeV-endep-EXP-M03-norm-max.txt'
        zmax = 180.
    elif tyPart == 'c12':
        hist_title = '^{12}C 100 MeV/u in Water, Geant4  ' + refer
        fname2 = '12C100MeVen-dep-EXP-norm-max.txt'
        zmax = 40.

    num_opt = '0', '2', '3'
    x_exp, y_exp = array( 'd' ), array ( 'd' )

    gROOT.Reset()
    gROOT.SetStyle('Plain')
    c1 = TCanvas('c1', 'c1',6,6,800,600)
    gStyle.SetOptStat(0)
    c1.SetFillColor(0)
    c1.SetBorderMode(0)
    c1.SetBorderSize(0)
    c1.SetFrameBorderMode(0)

    openFile(fname2)
    
    lines = infile.readlines()
    title = lines[0]

    for line in lines[1:]:
        xx = float(line.split() [0])
        yy = float(line.split() [1])
        xx *= 10.
        x_exp.append(xx)
        y_exp.append(yy)

    infile.close()

    h0 = gPad.DrawFrame(0.0,0.0,zmax,1.2,hist_title)
    h0.GetXaxis().SetTitle('z (mm)')
    h0.GetYaxis().SetTitle('dose (relative unit)')
    h0.Draw('AXIS SAME')

    nn = len(lines) - 1
    gr = TGraph(nn,x_exp,y_exp)
    gr.SetMarkerStyle(22)
    gr.SetMarkerSize(1.2)
    gr.Draw('P SAME')

    global leg
    leg = TLegend(0.2,0.65,0.45,0.86)
    leg.SetTextFont(52)
    leg.SetTextSize(0.035)
    leg.SetLineColor(1)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillStyle(0)
    leg.SetMargin(0.4)
    leg.SetBorderSize(1)
    leg.AddEntry(gr,'Data','p')
    c1.Update()

    global num, j
    j = 0
    for num in num_opt:
        finName = tyPart + '_opt' + num + '.out'
        fillHist(finName)
        j += 1
        
    leg.Draw()
    fout = 'A_' + tyPart + '_water.gif'
    c1.Print(fout)

def openFile(fileName):
    global infile
    try:
        infile = open(fileName, 'r')
    except IOError:
        print 'Input file <',fileName, '> does not exist! Exit'
        sys.exit(2)

def fillHist(fn):
    xmin = 0.0
    xmax = 300.
    nbin = 3000 
    legend = 'QBBC opt' + num
    
    openFile(fn)
    print 'File with MC <', fn, '> is opened'
    lines = infile.readlines()
    title = lines[0]
    
    hhh = 'h' + num
    hh = TH1D(hhh,'',nbin,xmin,xmax)
    hh.SetLineStyle(1)
    hh.SetLineWidth(2)
    hh.SetLineColor(j+2)

    nl = len(lines)
    for line in lines[1:(nl-1)]:
      xt = (float(line.split() [0]))
      yt = (float(line.split() [1]))
      k = lines.index(line) 
      hh.SetBinContent(k, yt)
      hh.SetBinError(k, yt/100.)

    maxJ = float (lines[nl-1].split() [0])
    maxX = float (lines[nl-1].split() [1])
    maxY = float (lines[nl-1].split() [2])

    print 'Histo filled N = ', nbin
    print ' maxJ = ', maxJ
    print ' maxX = ', maxX
    print ' maxY = ', maxY

    norm = maxY
    hh.Scale(1./norm)
    hh.Draw('HISTO SAME')

    infile.close()

    entry=leg.AddEntry(hh, legend, 'l')
    entry.SetLineColor(j+2)
    entry.SetLineStyle(1)
    entry.SetLineWidth(2)
    entry.SetTextColor(1)
    
###______________________________
if __name__ == "__main__":
    main()

