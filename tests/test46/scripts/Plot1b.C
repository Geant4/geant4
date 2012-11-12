{
gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/Style.C");

iac[0] = 1;
iac[1] = 1;
iac[2] = 1;

iener = 6;
c1 = new TCanvas("c1"," ",0.5, 5, 800, 800);

hh[0] = gPad->DrawFrame(0.,0.,80,0.06,part[iener] + " in ECAL+HCAL");
hh[0]->GetYaxis()->SetTitle(axtit[1]);
//gPad->SetLogy();
hh[0]->GetXaxis()->SetTitle(axtit[0]);

leg[0] = new TLegend(0.7, 0.5, 0.9, 0.8);
 
iplot =  4;
idir  = 17;
 iener = 7;

for(idir=17; idir<19; idir++) {
  gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/AddMC1.C");
  iener++;
}
 
leg[0]->Draw("SAME");
 cout << "Legend is added" << endl;
c1->Print("Abias_"+fil[7] + ".gif");
}
