{ 
gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/StyleCut.C");

TString ahed[5] = {"QGSP_BERT", "QGSP_BERT_EMV","QGSP_BERT_EML","QGSP_BERT 4T field","QGSP_BERT_EML 4T field"};

TLegend* legc = new TLegend(0.3, 0.25, 0.6, 0.45);

Int_t i;
gPad->SetLogx();
TH1F* hh = gPad->DrawFrame(0.005,0.11,200,0.14,tit[1]);
hh[i]->GetXaxis()->SetTitle(axtit[0]);
hh[i]->GetYaxis()->SetTitle(axtit[3]);
gPad->SetGrid();
hh->Draw("AXIS SAME");

//c1.Update();

Int_t nd = 9;
Double_t cut[nd]  = {0.01, 0.1, 0.3, 0.7, 1, 3, 7, 10, 100};
Double_t err0[nd] = {0., 0., 0., 0., 0., 0.,0., 0., 0.,};

//visible energy in HCAL
Double_t e0[nd] = {0.1353, 0.1381, 0.1383, 0.1369, 0.1374, 0.1363, 0.1347, 0.1339, 0.1355};
Double_t e1[nd] = {0.1314, 0.1317, 0.133,  0.1315, 0.1308, 0.1305, 0.1306, 0.1278, 0.1294};
Double_t e2[nd] = {0.1308, 0.1332, 0.1312, 0.1309, 0.1316, 0.1278, 0.1214, 0.1172, 0.0729};

TGraphErrors* gr[6];

for(i=0; i<3; i++) {
  if(i == 0) 
    gr[i] = new TGraphErrors(nd,cut,e0,err0,err0);
  if(i == 1) 
    gr[i] = new TGraphErrors(nd,cut,e1,err0,err0);
  if(i == 2) 
    gr[i] = new TGraphErrors(nd,cut,e2,err0,err0);

  cout << "Graph " << i << "  nd= " << nd << endl;

  gr[i]->SetMarkerColor(col[i]);
  gr[i]->SetMarkerStyle(mar[i]);
  gr[i]->Draw("PL SAME");
  legc->AddEntry(gr[i],ahed[i],"p");
}

legc->Draw();

c1.Print("hcal_vis_pi-30gev.gif");

}
