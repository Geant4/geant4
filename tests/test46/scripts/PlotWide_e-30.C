{
gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/Style.C");
gStyle->SetMarkerSize(2);
//gStyle->SetTextSize(2);
TCanvas c1("c1"," ",0.5, 5, 1000, 600);

iac[7] = 1;
iac[0] = 1;
//iac[2] = 1;

iplot0 = 6;
iener = 6;
c1.Divide(2,1);

c1.cd(1);
gPad->SetGrid();
hh[0] = gPad->DrawFrame(0.55,0.,0.85,0.012,"ECAL e- 30 GeV");
hh[0]->GetYaxis()->SetTitle("");
hh[0]->GetXaxis()->SetTitle("E_{0}/E_{5x5}");

leg[0] = new TLegend(0.65, 0.50, 0.95, 0.95);
leg[0]->SetHeader("FTFP_BERT_EMV");

c1.cd(2);
gPad->SetGrid();
hh[1] = gPad->DrawFrame(0.88,0.,0.98,0.2,"ECAL e- 30 GeV");
hh[1]->GetYaxis()->SetTitle("");
hh[1]->GetXaxis()->SetTitle("E_{3x3}/E_{5x5}");

leg[1] = new TLegend(0.15, 0.45, 0.5, 0.85);
leg[1]->SetHeader("FTFP_BERT_EMV");
 
for(iplot=0; iplot<2; iplot++) {

  c1.cd(iplot + 1);
  for(idir = 13; idir<17; idir++) {
    gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/AddMC.C");
  }
  cout << "PlotSingle " << iplot << " done " << endl;
}
leg[1]->Draw("SAME");
c1.Print("A_e-30d.gif");

}
