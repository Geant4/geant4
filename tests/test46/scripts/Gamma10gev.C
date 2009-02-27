{ 
gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/StyleCut.C");

TLegend* legc = new TLegend(0.3, 0.25, 0.7, 0.45);
TString ahed[5] = {"QGSP_BERT", "QGSP_BERT_EMV","QGSP_BERT 4T field","QGSP_BERT_EMV 4T field","QGSP_BERT_EML 4T field"};

gPad->SetLogx();
TH1F* hh = gPad->DrawFrame(0.005,0.06,200,0.12,tit[0]);
hh->GetXaxis()->SetTitle(axtit[0]);
hh->GetYaxis()->SetTitle(axtit[1]);
gPad->SetGrid();
hh->Draw("AXIS SAME");

//c1.Update();

Int_t nd = 9;
Double_t cut[nd]  = {0.01, 0.1, 0.3, 0.7, 1, 3, 7, 10, 100};
Double_t err0[nd] = {0., 0., 0., 0., 0., 0.,0., 0., 0.,};

Double_t e0[nd] = {0.1007, 0.1014, 0.1017, 0.1016, 0.1015, 0.1003, 0.09947, 0.09911, 0.09956};
Double_t e1[nd] = {0.09659, 0.09642, 0.09664, 0.09629, 0.09676, 0.09561, 0.09550, 0.09500, 0.0942};
Double_t e2[nd] = {0.1153, 0.1175, 0.1181, 0.1177, 0.1175, 0.1165, 0.1154, 0.1144, 0.1100};
Double_t e3[nd] = {0.1114, 0.1119, 0.1115, 0.1120, 0.1119, 0.1109, 0.1100, 0.11, 0.1086};
Double_t e4[nd] = {0.1115, 0.1114, 0.1113, 0.1120, 0.1109, 0.1085, 0.09972, 0.09425, 0.02132};

TGraphErrors* gr[6];

//c1.cd(1);

for(int i=0; i<5; i++) {
  if(i == 0) 
    gr[i] = new TGraphErrors(nd,cut,e0,err0,err0);
  if(i == 1) 
    gr[i] = new TGraphErrors(nd,cut,e1,err0,err0);
  if(i == 2) 
    gr[i] = new TGraphErrors(nd,cut,e2,err0,err0);
  if(i == 3) 
    gr[i] = new TGraphErrors(nd,cut,e3,err0,err0);
  if(i == 4) 
    gr[i] = new TGraphErrors(nd,cut,e4,err0,err0);
  

  cout << "Graph " << i << "  nd= " << nd << endl;

  gr[i]->SetMarkerColor(col[i]);
  gr[i]->SetMarkerStyle(mar[i]);
  gr[i]->Draw("PL SAME");
  legc->AddEntry(gr[i],ahed[i],"p");
}

legc->Draw();

c1.Print("hcal_g10gev.gif");

}
