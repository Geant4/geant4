{ 
gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/StyleCut.C");

TString ahed[5] = {"QGSP_BERT", "QGSP_BERT_EMV","QGSP_BERT_EML","QGSP_BERT 4T field","QGSP_BERT_EML 4T field"};

TLegend* legc = new TLegend(0.3, 0.25, 0.6, 0.45);

Int_t i;
gPad->SetLogx();
TH1F* hh = gPad->DrawFrame(0.005,11,200,15,tit[1]);
hh[i]->GetXaxis()->SetTitle(axtit[0]);
hh[i]->GetYaxis()->SetTitle(axtit[4]);
gPad->SetGrid();
hh->Draw("AXIS SAME");

//c1.Update();

Int_t nd = 9;
Double_t cut[nd]  = {0.01, 0.1, 0.3, 0.7, 1, 3, 7, 10, 100};
Double_t err0[nd] = {0., 0., 0., 0., 0., 0.,0., 0., 0.,};
Double_t err1[nd] = {0.035, 0.035, 0.035, 0.035, 0.035, 0.035,0.035, 0.035, 0.035,};

Double_t e0[nd] = {14.49, 14.58, 14.55, 14.46, 14.53, 14.6, 14.53, 14.53, 14.57};
Double_t e1[nd] = {14.50, 14.51, 14.71, 14.53, 14.50, 14.53, 14.61, 14.38, 14.49};
Double_t e2[nd] = {14.40, 14.64, 14.53, 14.51, 14.57, 14.57, 14.52, 14.46, 14.17};

TGraphErrors* gr[6];

//c1.cd(1);

for(i=0; i<3; i++) {
  if(i == 0) 
    gr[i] = new TGraphErrors(nd,cut,e0,err0,err1);
  if(i == 1) 
    gr[i] = new TGraphErrors(nd,cut,e1,err0,err1);
  if(i == 2) 
    gr[i] = new TGraphErrors(nd,cut,e2,err0,err1);

  cout << "Graph " << i << "  nd= " << nd << endl;

  gr[i]->SetMarkerColor(col[i]);
  gr[i]->SetMarkerStyle(mar[i]);
  gr[i]->Draw("PL SAME");
  legc->AddEntry(gr[i],ahed[i],"p");
}

legc->Draw();

c1.Print("hcal_tot_pi-30gev.gif");

}
