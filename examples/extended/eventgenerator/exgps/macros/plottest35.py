#!/usr/bin/python

from ROOT import *
from array import array

gROOT.Reset()

input_file_1=TFile('test35a.root','READ')
input_file_2=TFile('test35b.root','READ')

#input_file_1.cd()
#h_1_1 = input_file_1.Get("h16")

c1 = TCanvas('c1', 'test35', 200, 10, 700, 500)
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
t1 = gROOT.FindObject('MyTuple')
for i in range(t1.GetEntries()):
    t1.GetEntry(i)
    h_1.Fill(t1.Energy,t1.Weight)

input_file_2.cd()
# get the tuple t1
t1 = gROOT.FindObject("MyTuple")
for i in range(t1.GetEntries()):
    t1.GetEntry(i)
    h_2.Fill(t1.Energy,t1.Weight)

h_2.SetLineStyle(kDashed);
h_2.SetLineColor(kBlue);
h_2.Draw();
h_1.Draw("same") ;
c1.Update()
c1.Print("./test35.png")

input_file_1.Close()
input_file_2.Close()




