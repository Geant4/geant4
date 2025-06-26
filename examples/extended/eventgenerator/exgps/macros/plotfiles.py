#!/usr/bin/env python

import ROOT
from array import array

def plot_1_file (file):
    ROOT.gROOT.Reset()
    ROOT.gROOT.SetBatch(True)

    input_file_1=ROOT.TFile(file+'.root','READ')
    h11 = input_file_1.Get('histo/h1.1')
    h12 = input_file_1.Get('histo/h1.2')
    h13 = input_file_1.Get('histo/h1.3')
    h14 = input_file_1.Get('histo/h1.4')
    h21 = input_file_1.Get('histo/h2.1')
    h22 = input_file_1.Get('histo/h2.2')
    h23 = input_file_1.Get('histo/h2.3')
    h24 = input_file_1.Get('histo/h2.4')

    c1 = ROOT.TCanvas('c1', file+'.h1', 200, 10, 1200, 700)
    c1.Divide(4,2)
    histos1 = [h11, h12, h13, h14]
    pad = 1
    for h1 in histos1:
        if h1:
            c1.cd(pad)
            h1.Draw()
            pad = pad + 1

    histos2 = [h21, h22, h23, h24]
    pad = 5
    for h2 in histos2:
       if h2:
           c1.cd(pad)
           h2.Draw()
           pad = pad + 1
    c1.Print(file+'.png')

    input_file_1.Close()

def plot_2_files (file):
   ROOT.gROOT.Reset()
   ROOT.gROOT.SetBatch(True)

   input_file_1=ROOT.TFile(file+"a.root",'READ')
   input_file_2=ROOT.TFile(file+"b.root",'READ')

   input_file_1.cd()
   h_1_1 = input_file_1.Get('histo/h1.1')

   c1 = ROOT.TCanvas('c1', file, 200, 10, 700, 500)
   c1.SetGridx()
   c1.SetGridy()
   c1.SetLogx()
   c1.SetLogy()

   #histogram for energy spectra
   n = 41
   bin = array( 'f' )

   for i in range( n ):
       bin.append(pow(10,(-2+0.1*i)))

   h_1 = ROOT.TH1F('unbiased','Source Spectrum',40,bin)
   h_2 = ROOT.TH1F('biased','Source Spectrum',40,bin)

   input_file_1.cd()
   #get the tuple t1
   t1 = input_file_1.Get('ntuple/101')
   for i in range(t1.GetEntries()):
       t1.GetEntry(i)
       h_1.Fill(t1.Ekin,t1.weight)

   input_file_2.cd()
   # get the tuple t1
   t2 =  input_file_2.Get('ntuple/101')
   for i in range(t2.GetEntries()):
       t2.GetEntry(i)
       h_2.Fill(t2.Ekin,t2.weight)

   h_2.Draw();
   h_1.Draw("same") ;
   c1.Update()
   c1.Print(file+".png")

   input_file_1.Close()
   input_file_2.Close()

