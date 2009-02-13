{ 
gROOT->Reset();
gROOT->SetStyle("Plain"); 
gStyle->SetOptStat(0);
gStyle->SetNdivisions(210, "x");
gStyle->SetNdivisions(10, "y");
gStyle->SetTextFont(1);
gStyle->SetLabelOffset(0.005, "x");
gStyle->SetLabelOffset(0.005, "y");
gStyle->SetLabelSize(0.05, "x");
gStyle->SetLabelSize(0.05, "y");
gStyle->SetTitleOffset(1.2, "x");
gStyle->SetTitleOffset(1.4, "y");
gStyle->SetTitleSize(0.05, "x");
gStyle->SetTitleSize(0.05, "y");
gStyle->SetPadBottomMargin(0.15);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadRightMargin(0.1);
gStyle->SetTickLength(0.05, "x");
gStyle->SetTickLength(0.03, "y");
gStyle->SetPadBorderMode(0);

gStyle->SetMarkerSize(1.5);
TCanvas c1("c1"," ",0.5, 5, 600, 800);

TH1F* hh[6];

Int_t col[6] = {4, 2, 3, 4, 2, 3};
Int_t mar[6] = {21, 20, 22, 21, 20, 22};

c1->Divide(1,2);

TString ahed[3] = {"e- 1 GeV no sat.", "e- 1 GeV Birks sat.","p 100 MeV Birks sat."};
TString tit[2] = {"e^{-} 1 GeV in BGO", ""};
TString axtit[3] = {"cut (mm)","E_{vis}/E_{0}","#sigma_{E}/E_{vis} (%)"};
TLegend* legc = new TLegend(0.4, 0.25, 0.8, 0.45);
Double_t y1[2] = {0.9, 1.02};
Double_t x[2] = {0.00005, 20};
Double_t y2[2] = {0., 0.6};

Int_t i;
for(i=0; i<2; i++) {
  c1.cd(i+1);
  gPad->SetLogx();
  if(i==0) {
    hh[i] = gPad->DrawFrame(x[0],y1[0],x[1],y1[1],tit[i]);
  }
  if(i==1) {
    hh[i] = gPad->DrawFrame(x[0],y2[0],x[1],y2[1],tit[i]);
  }
  hh[i]->GetXaxis()->SetTitle(axtit[0]);
  hh[i]->GetYaxis()->SetTitle(axtit[i+1]);
  gPad->SetGrid();
  hh[i]->Draw("AXIS SAME");
}

c1.Update();

Int_t nd = 6;
Double_t cut[nd]  = {0.0001, 0.001, 0.01, 0.1, 1, 10};
Double_t err0[nd] = {0., 0., 0., 0., 0., 0.};

Double_t e0[nd]  = {0.9965, 0.9965, 0.9964, 0.9964, 0.9963, 0.996};
Double_t e1[nd] = {0.9625, 0.9759, 0.9868, 0.9871, 0.9869, 0.986};
Double_t e2[nd] = {0.9247, 0.9397, 0.9443, 0.9423, 0.9422, 0.942};
Double_t s0[nd]  = {0.3912, 0.4004, 0.3902, 0.3633, 0.3752, 0.4481};
Double_t s1[nd] = {0.3537, 0.3505, 0.3601, 0.3609, 0.3748, 0.4042};
Double_t s2[nd] = {0.5142, 0.5529, 0.6985, 0.4947, 0.695, 0.4181};

TGraphErrors* gr[6];

c1.cd(1);

for(i=0; i<4; i++) {
  if(i == 0) 
    gr[i] = new TGraphErrors(nd,cut,e0,err0,err0);
  if(i == 1) 
    gr[i] = new TGraphErrors(nd,cut,e1,err0,err0);
  if(i == 2) 
    gr[i] = new TGraphErrors(nd,cut,e2,err0,err0);
  if(i == 3) {
    c1.Update();
    c1.cd(2);
    gr[i] = new TGraphErrors(nd,cut,s0,err0,err0);
  }
  if(i == 4) 
    gr[i] = new TGraphErrors(nd,cut,s1,err0,err0);
  //if(i == 5) 
  //  gr[i] = new TGraphErrors(nd,cut,s2,err0,err0);

  cout << "Graph " << i << "  nd= " << nd << endl;

  gr[i]->SetMarkerColor(col[i]);
  gr[i]->SetMarkerStyle(mar[i]);
  gr[i]->Draw("P SAME");
  if(i <= 2) legc->AddEntry(gr[i],ahed[i],"p");
}

legc->Draw();

c1.Print("bgo.gif");

}
