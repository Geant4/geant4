{
gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/Style.C");

iac[0] = 1;
iac[1] = 1;
iac[2] = 1;

iener = 2;
c1 = new TCanvas("c1"," ",0.5, 5, 800, 800);

hh[0] = gPad->DrawFrame(0.,0.,600,0.08,part[iener] + " in ECAL+HCAL");
hh[0]->GetYaxis()->SetTitle(axtit[1]);
//gPad->SetLogy();
hh[0]->GetXaxis()->SetTitle(axtit[0]);

leg[0] = new TLegend(0.6, 0.6, 0.9, 0.9);
 
iplot = 5;

for(idir=0; idir<ndir; idir++) {
if(iac[idir]>0)gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/AddMC1.C");
}
 
leg[0]->Draw("SAME");
c1->Print("AA"+fil[iener] + ".gif");
}
