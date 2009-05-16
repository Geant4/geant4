{ 
gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/StyleCut.C");

gStyle->SetTitleSize(0.04, "x");
gStyle->SetTitleSize(0.04, "y");
gStyle->SetMarkerSize(1.2);

TLegend* legc = new TLegend(0.25, 0.65, 0.5, 0.88);
TString bhed[6] = {"QGSP_BERT pi-","QGSP_BERT_EMV pi-","FTFP_BERT_EML pi-","QGSP_BERT p","QGSP_BERT_EMV p","FTFP_BERT_EML p"};

gPad->SetLogx();
TH1F* hh = gPad->DrawFrame(0.9,0.5,310,0.9,"Response of combined calorimeter");
hh->GetXaxis()->SetTitle("E_{0} (GeV)");
hh->GetYaxis()->SetTitle("E_{rec}/E_{0}");
gPad->SetGrid();
hh->Draw("AXIS SAME");

Int_t nd = 17;
Double_t ener[nd]  = {1,2,3,4,5,6,7,8,9,10,12,20,30,50,100,200,300};
Double_t err0[nd] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

Double_t e0[nd] = {.560,.559,.636,.644,.661,.677,.691,.699,.712,.713,.718,.795,.829,.839,.853,.862,.865};
Double_t e1[nd] = {.560,.589,.619,.635,.651,.659,.672,.679,.689,.694,.702,.771,.805,.812,.824,.830,.819};
Double_t e2[nd] = {.564,.585,.623,.647,.663,.676,.687,.688,.703,.709,.721,.791,.824,.838,.848,.859,.862};

Double_t e3[nd] = {.641,.607,.608,.645,.654,.659,.661,.670,.680,.684,.695,.761,.791,.811,.834,.851,.858};
Double_t e4[nd] = {.646,.599,.613,.613,.645,.640,.657,.663,.670,.670,.681,.747,.772,.787,.804,.822,.826};
Double_t e5[nd] = {.645,.603,.616,.638,.647,.660,.668,.672,.684,.687,.689,.763,.790,.812,.830,.848,.853};

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
  gr[i]->SetMarkerStyle(mar[i]);
  gr[i]->Draw("PL SAME 9");
  legc->AddEntry(gr[i],bhed[i],"p");
  gr[i+3]->SetMarkerColor(col[i]);
  gr[i+3]->SetLineColor(col[i]);
  // gr[i+3]->SetTextColor(col[i]);
  gr[i+3]->SetMarkerStyle(mar[i+3]);
  gr[i+3]->Draw("PL SAME 9");
  legc->AddEntry(gr[i+3],bhed[i+3],"p");
}

legc->Draw();

c1.Print("Erec.gif");

}
