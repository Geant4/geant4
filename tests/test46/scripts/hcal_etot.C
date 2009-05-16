{ 
gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/StyleCut.C");

gStyle->SetTitleSize(0.04, "x");
gStyle->SetTitleSize(0.04, "y");
gStyle->SetMarkerSize(1.2);

TLegend* legc = new TLegend(0.25, 0.65, 0.5, 0.88);
TString bhed[6] = {"QGSP_BERT_EMV pi-","QGSP_BERT pi-","FTFP_BERT pi-","QGSP_BERT_EMV p","QGSP_BERT p","FTFP_BERT p"};

gPad->SetLogx();
TH1F* hh = gPad->DrawFrame(0.9,0.65,310,0.9,"ECAL+HCAL energy deposition");
hh->GetXaxis()->SetTitle("E (GeV)");
hh->GetYaxis()->SetTitle("(E_{HCAL}+E_{ECAL}/E_0");
gPad->SetGrid();
hh->Draw("AXIS SAME");

Int_t nd = 17;
Double_t ener[nd]  = {1,2,3,4,5,6,7,8,9,10,12,20,30,50,100,200,300};
Double_t err0[nd] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

Double_t e0[nd] = {.799,.799,.787,.791,.767,.767,.765,.767,.767,.708,.710,.781,.811,.818,.831,.843,.848};
Double_t e1[nd] = {.771,.771,.775,.765,.776,.768,.766,.758,.762,.706,.721,.789,.818,.817,.830,.8425,.849};
Double_t e2[nd] = {.681,.756,.729,.760,.748,.792,.771,.785,.793,.781,.804,.817,.826,.840,.851,.863,.867};

Double_t e3[nd] = {.744,.739,.728,.728,.719,.727,.726,.728,.735,.661,.697,.738,.763,.778,.802,.818,.828};
Double_t e4[nd] = {.772,.742,.729,.738,.727,.732,.730,.734,.722,.667,.691,.732,.775,.781,.799,.817,.828};
Double_t e5[nd] = {.762,.762,.769,.738,.713,.748,.744,.737,.757,.758,.760,.778,.789,.808,.824,.839,.845};

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

c1.Print("Khcal_etot.gif");

}
