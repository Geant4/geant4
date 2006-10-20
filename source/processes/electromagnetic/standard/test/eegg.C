{
gROOT->Reset();
gStyle->SetOptStat(0);
gStyle->SetNdivisions(210, "x");
gStyle->SetNdivisions(210, "y");
gStyle->SetLabelOffset(0.005, "x");
gStyle->SetLabelOffset(0.005, "y");
gStyle->SetLabelSize(0.07, "x");
gStyle->SetLabelSize(0.07, "y");
gStyle->SetTitleOffset(0.8, "x");
gStyle->SetTitleOffset(0.8, "y");
gStyle->SetTitleSize(0.08, "x");
gStyle->SetTitleSize(0.08, "y");
gStyle->SetPadBottomMargin(0.15);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadRightMargin(0.1);
gStyle->SetTickLength(0.04, "x");
gStyle->SetTickLength(0.02, "y");
gStyle->SetPadBorderMode(0);

TCanvas c1("c1","eegg ",0.5, 5, 600, 600);

c1->Divide(1,2);

TString ang[2] = {"",""};

TLegend* leg[2];
TH1F* hh[4];
leg[0] = new TLegend(0.2, 0.5, 0.35, 0.8);
leg[0]->SetTextSize(0.07);
leg[1] = new TLegend(0.2, 0.5, 0.35, 0.8);
leg[1]->SetTextSize(0.07);
leg[0]->SetHeader("100 keV"); 
leg[1]->SetHeader("100 keV"); 
//  c1.cd(i+1);
//  hh[i] = gPad->DrawFrame(0, 0, 6.5, 700, "");
//  hh[i]->GetXaxis()->SetTitle("p (GeV)");
//  hh[i]->GetYaxis()->SetTitle("#sigma (mb/GeV/srad)");
//  hh[i]->Draw("AXIS");

TString file1 = "eegg100kev.root";
//TString file1 = "eegg10kev.root";
//TString file1 = "eegg1mev.root";

TFile ff(file1);
c1.cd(1);
h2->SetLineColor(4); 
h2->Draw("HIST");
leg[0]->AddEntry(h2, "G4 8.2", "l"); 

h1->SetLineColor(2); 
h1->Draw("HIST SAME");
leg[0]->AddEntry(h1, "G4 8.1", "l"); 

c1.cd(2);
h4->SetLineColor(4); 
h4->Draw("HIST");
leg[1]->AddEntry(h4, "G4 8.2", "l"); 

h3->SetLineColor(2); 
h3->Draw("HIST SAME");
leg[1]->AddEntry(h3, "G4 8.1", "l"); 

for(i=0; i<2; i++) {
  c1.cd(i+1);
  leg[i]->Draw();
}

c1.Update();

c1.Print("eegg100kev.eps");
//c1.Print("eegg1mev.eps");

}
