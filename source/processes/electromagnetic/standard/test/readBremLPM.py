# $Id: readBremLPM.py,v 1.2 2008-08-22 08:16:17 schaelic Exp $
# GEANT4 tag $Name: not supported by cvs2svn $
#
# works together with BremLPMTest.cc
#
from ROOT import gROOT, gApplication, TFile, TCanvas

# set nice style
gROOT.SetStyle("Plain")

# create 6 pads
c1 = TCanvas('c1','CrossSection')
c1.Divide(3,2)

c2 = TCanvas('c2','DEDX')
c2.Divide(3,2)

#open root file
f=TFile('eBremRel01.root')
t=f.Get('info')
zlink=t.GetLeaf('Z')
t.Scan()

#plot cross section histograms
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

    g2=f.Get('model2;'+str(i+1))
    g2.Draw('lp');
    min2=g2.GetHistogram().GetMinimum()
    max2=g2.GetHistogram().GetMaximum()

    gR.SetMaximum(1.1*max(max1,maxR))
    gR.SetMinimum(0.9*min(min1,minR))
    print Z

#plot dEdx histograms
for i in range(t.GetEntries()):
    t.GetEntry(i)
    Z=int(zlink.GetValue())
    c2.cd(i+1)
    c2.GetPad(i+1).SetLogx()
    c2.GetPad(i+1).SetLogy()

    gR=f.Get('xDEDXR;'+str(i+1))
    gR.GetHistogram().SetTitle('Z='+str(Z))
    gR.Draw('Alp');
    minR=gR.GetHistogram().GetMinimum()
    maxR=gR.GetHistogram().GetMaximum()
    
    g1=f.Get('xDEDX1;'+str(i+1))
    g1.Draw('lp');
    min1=g1.GetHistogram().GetMinimum()
    max1=g1.GetHistogram().GetMaximum()

    g2=f.Get('xDEDX2;'+str(i+1))
    g2.Draw('lp');
    min2=g2.GetHistogram().GetMinimum()
    max2=g2.GetHistogram().GetMaximum()

#    gR.SetMaximum(1.e-3)
    gR.SetMinimum(1.e-6*max(max1,maxR))
    print Z

#done
gApplication.Run() 
