# $Id: readBremLPM.py,v 1.1 2008-08-21 15:16:53 schaelic Exp $
# GEANT4 tag $Name: not supported by cvs2svn $
from ROOT import gROOT, gApplication, TFile, TCanvas

# set nice style
gROOT.SetStyle("Plain")

# create 6 pads
c1 = TCanvas('c1','CrossSection')
c1.Divide(3,2)

#open root file
f=TFile('eBremRel01.root')
t=f.Get('info')
zlink=t.GetLeaf('Z')
t.Scan()

#plot histograms
for i in range(t.GetEntries()):
    t.GetEntry(i)
    Z=int(zlink.GetValue())
    c1.cd(i+1)
    c1.GetPad(i+1).SetLogx()

    gR=f.Get('modelR;'+str(i+1))
    gR.GetHistogram().SetTitle('Z='+str(Z))
    gR.Draw('Alp');
    minR=gR.GetHistogram().GetMinimum()
    maxR=gR.GetHistogram().GetMaximum()
    
    g1=f.Get('model1;'+str(i+1))
    g1.Draw('lp');
    min1=g1.GetHistogram().GetMinimum()
    max1=g1.GetHistogram().GetMaximum()

    gR.SetMaximum(1.1*max(max1,maxR))
    gR.SetMinimum(0.9*min(min1,minR))
    print Z

#done
gApplication.Run() 
