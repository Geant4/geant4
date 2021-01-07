#!/usr/bin/python

from ROOT import *
from array import array

def plot_1_file (file):
    gROOT.Reset()
    input_file_1=TFile(file+'.root','READ')
    h1 = input_file_1.Get("histo/1")
    h2 = input_file_1.Get("histo/2")
    h3 = input_file_1.Get("histo/3")
    h4 = input_file_1.Get("histo/4")
    h5 = input_file_1.Get("histo/5")
    h6 = input_file_1.Get("histo/6")

    c1 = TCanvas('c1', file, 200, 10, 700, 900)
    c1.Divide(2,3)

    c1.cd(1)
    h1.Draw()
    c1.cd(2)
    h2.Draw()
    c1.cd(3)
    h3.Draw()
    c1.cd(4)
    h4.Draw()
    c1.cd(5)
    h5.Draw()
    c1.cd(6)
    h6.Draw()
    c1.Update()
    c1.Print("./"+file+".png")
    
    input_file_1.Close()

#    h_gam.SetLineColor(2)

def plot_2_files (file):
    gROOT.Reset()

    input_file_1=TFile(file+'a.root','READ')
    input_file_2=TFile(file+'b.root','READ')

#input_file_1.cd()
#h_1_1 = input_file_1.Get("h16")

    c1 = TCanvas('c1', file, 200, 10, 700, 500)
    c1.SetGridx()
    c1.SetGridy()
    c1.SetLogx()
    c1.SetLogy()

# histogram for energy spectra
    n = 41
    bin = array( 'f' )
    
    for i in range( n ):
        bin.append(pow(10,(-2+0.1*i)))
#
    h_1 = TH1F('unbiased','Source Spectrum',40,bin)
    h_2 = TH1F('biased','Source Spectrum',40,bin)

#
    input_file_1.cd()
# get the tuple t1
    t1 = input_file_1.Get('ntuple/MyTuple')
    print t1
    for i in range(t1.GetEntries()):
        t1.GetEntry(i)
        h_1.Fill(t1.Ekin,t1.weight)
            
    input_file_2.cd()
# get the tuple t1
    t2 =  input_file_2.Get('ntuple/MyTuple')
    for i in range(t2.GetEntries()):
        t2.GetEntry(i)
        h_2.Fill(t2.Ekin,t2.weight)
                
    h_2.SetLineStyle(kDashed);
    h_2.SetLineColor(kBlue);
    h_2.Draw();
    h_1.Draw("same") ;
    c1.Update()
    c1.Print("./"+file+".png")
    
    input_file_1.Close()
    input_file_2.Close()

