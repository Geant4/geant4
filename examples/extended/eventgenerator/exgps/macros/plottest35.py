#!/usr/bin/python

from ROOT import gROOT, TCanvas, TH1F, TFile, TNtuple, TDirectory
from array import array

gROOT.Reset()

input_file_1=TFile('test35a.root','READ')
input_file_2=TFile('test35b.root','READ')

#input_file_1.cd()
#h_1_1 = input_file_1.Get("h16")
#input_file_2.cd()
#h_2_1 = input_file_2.Get("h16")

c1 = TCanvas('c1', 'test35', 200, 10, 700, 500)
c1.SetGridx()
c1.SetGridy()
c1.SetLogx()
c1.SetLogy()
# c1.Divide(1,2)
# c1.cd(1)
#h_1_1.Draw()
# c1.cd(2)
#h_2_1.Draw('same')
#c1.Update()
#c1.Print("./exrdm.in.png")
#c1.SetLogx(0)

# histogram for energy spectra
n = 31
bin = array( 'd' )

for i in range( n ):
    bin.append(pow(10,0.1*i))
#
h_1 = TH1F('unbiased','Unbiased',30,bin)
h_2 = TH1F('biased','Biased',30,bin)
input_file_1.cd()
# get the tuple t1
t1 = gROOT.FindObject('MyTuple')
energy = array('f')
energy.append(0.)
weight = array('f')
weight.append(0.)
t1.SetBranchAddress("Energy",energy);
t1.SetBranchAddress("Weight",weight);
for i in range(t1.GetEntries()):
    t1.GetEntry(i)
    h_1.Fill(energy[0],weight[0])
#t1.Draw('Energy>>h_1','','')
#t1.Project('unbiased','Energy')
h_1.Draw()
input_file_2.cd()
# get the tuple t1
t1 = gROOT.FindObject("MyTuple")
t1.Project("biased","Energy")
h_1.Draw()
h_2.Draw('same')
c1.Update()
c1.Print("./test35.png")

input_file_1.Close()
input_file_2.Close()




