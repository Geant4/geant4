{ 
gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/StyleCut.C");

TString ahed[5] = {"QGSP_BERT", "QGSP_BERT_EMV","QGSP_BERT_EML","QGSP_BERT 4T field","QGSP_BERT_EML 4T field"};

TLegend* legc = new TLegend(0.3, 0.25, 0.6, 0.45);

Int_t i;
gPad->SetLogx();
TH1F* hh = gPad->DrawFrame(0.005,20,200,25,tit[1]);
hh[i]->GetXaxis()->SetTitle(axtit[0]);
hh[i]->GetYaxis()->SetTitle(axtit[5]);
gPad->SetGrid();
hh->Draw("AXIS SAME");

//c1.Update();

Int_t nd = 9;
Double_t cut[nd]  = {0.01, 0.1, 0.3, 0.7, 1, 3, 7, 10, 100};
Double_t err0[nd] = {0., 0., 0., 0., 0., 0.,0., 0., 0.,};
Double_t err1[nd] = {0.06, 0.06, 0.06, 0.06, 0.06, 0.06,0.055, 0.052, 0.052,};

Double_t e0[nd] = {24.11, 24.29, 24.36, 24.25, 24.25, 24.10, 23.90, 23.84, 23.94};
Double_t e1[nd] = {23.59, 23.59, 23.59, 23.54, 23.47, 23.42, 23.28, 23.18, 23.27};
Double_t e2[nd] = {23.52, 23.63, 23.48, 23.50, 23.48, 23.03, 22.33, 21.87, 17.06};

TGraphErrors* gr[6];

c1.cd(1);

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

c1.Print("erec_pi-30gev.gif");

}
