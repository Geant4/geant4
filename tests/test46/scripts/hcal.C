{ 
gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/StyleCut.C");

gStyle->SetTitleSize(0.04, "x");
gStyle->SetTitleSize(0.04, "y");
gStyle->SetMarkerSize(1.2);

TLegend* legc = new TLegend(0.45, 0.6, 0.9, 0.9);
TString bhed[6] = {"QGSP_BERT_EMV pi-","QGSP_BERT pi-","FTFP_BERT pi-","QGSP_BERT_EMV p","QGSP_BERT p","FTFP_BERT p"};

//gPad->SetLogx();
TH1F* hh = gPad->DrawFrame(0,95,310,120,"HCAL response");
hh->GetXaxis()->SetTitle("E (GeV)");
hh->GetYaxis()->SetTitle("E_{HCAL}/E_{vis}");
gPad->SetGrid();
hh->Draw("AXIS SAME");

Int_t nd = 17;
Double_t ener[nd]  = {1,2,3,4,5,6,7,8,9,10,12,20,30,50,100,200,300};
Double_t err0[nd] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

Double_t e0[nd] = {121,121,121,121,121,121,121,114.1,116.9,118.7,113.8,104.8,103.5,101.1,100.9,102.0,102.7};
Double_t e1[nd] = {121,121,121,121,121,121,113.3,114.3,107.4,114.8,108.1,100.1,97.5,97.2,96.3,97.7,95.9};
Double_t e2[nd] = {121,121,121,121,121,110,116.4,110.3,105.1,104.5,105.0,99.5,98.9,97.5,98.0,98.3,98.6};

Double_t e3[nd] = {121,121,121,121,121,121,121,121,121,121,121,120.3,104.1,101.7,102.5,100.1,101.2};
Double_t e4[nd] = {121,121,121,121,121,121,137.6,124.8,121.8,1,1,106.,100.1,97.3,95.4,96.4,96.1};
Double_t e5[nd] = {121,121,121,121,121,121,121,121,122.7,114.5,106.7,105.8,100.7,98.4,99.0,98.4,97.0};

TGraphErrors* gr[6];

for(int i=0; i<3; i++) {
  if(i == 0) {
    gr[i] = new TGraphErrors(nd,ener,e0,err0,err0);
    gr[i+3] = new TGraphErrors(nd,ener,e3,err0,err0);
  } else if(i == 1) { 
    gr[i] = new TGraphErrors(nd,ener,e1,err0,err0);
    gr[i+3] = new TGraphErrors(nd,ener,e4,err0,err0);
  } else if(i == 2) { 
    gr[i] = new TGraphErrors(nd,ener,e2,err0,err0);
    gr[i+3] = new TGraphErrors(nd,ener,e5,err0,err0);
  }  
  gr[i]->SetMarkerColor(col[i]);
  gr[i]->SetLineColor(col[i]);
  //  gr[i]->SetTextColor(col[i]);
  gr[i]->SetMarkerStyle(mar[0]);
  gr[i]->Draw("PL SAME 9");
  legc->AddEntry(gr[i],bhed[i],"p");
  gr[i+3]->SetMarkerColor(col[i]);
  gr[i+3]->SetLineColor(col[i]);
  //  gr[i+3]->SetTextColor(col[i]);
  gr[i+3]->SetMarkerStyle(mar[1]);
  gr[i+3]->Draw("PL SAME 9");
  legc->AddEntry(gr[i+3],bhed[i+3],"p");
}

legc->Draw();

c1.Print("Khcal.gif");

}
