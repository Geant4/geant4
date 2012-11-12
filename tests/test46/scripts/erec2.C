{ 
gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/StyleCut.C");

gStyle->SetTitleSize(0.04, "x");
gStyle->SetTitleSize(0.04, "y");
gStyle->SetMarkerSize(1.2);

TLegend* legc = new TLegend(0.25, 0.65, 0.5, 0.88);
TString bhed[6] = {"pi- g4 9.4p04","pi- g4 9.6beta QBBC","pi- g4 9.6beta QGSP_FTFP_BERT_EMV","p g4 9.4p04","p g4 9.6beta QBBC","p g4 9.6beta QGSP_FTFP_BERT_EML"};

gPad->SetLogx();
TH1F* hh = gPad->DrawFrame(0.9,0.5,310,0.9,"Response of combined calorimeter");
hh->GetXaxis()->SetTitle("E_{0} (GeV)");
hh->GetYaxis()->SetTitle("E_{rec}/E_{0}");
gPad->SetGrid();
hh->Draw("AXIS SAME");

Int_t nd = 17;
Double_t ener[nd]  = {1,2,3,4,5,6,7,8,9,10,12,20,30,50,100,200,300};
Double_t err0[nd] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
/*
Double_t e0[nd] = {.560,.559,.636,.644,.661,.677,.691,.699,.712,.713,.718,.795,.829,.839,.853,.862,.865};
Double_t e1[nd] = {.560,.589,.619,.635,.651,.659,.672,.679,.689,.694,.702,.771,.805,.812,.824,.830,.819};
Double_t e2[nd] = {.564,.585,.623,.647,.663,.676,.687,.688,.703,.709,.721,.791,.824,.838,.848,.859,.862};

Double_t e3[nd] = {.641,.607,.608,.645,.654,.659,.661,.670,.680,.684,.695,.761,.791,.811,.834,.851,.858};
Double_t e4[nd] = {.646,.599,.613,.613,.645,.640,.657,.663,.670,.670,.681,.747,.772,.787,.804,.822,.826};
Double_t e5[nd] = {.645,.603,.616,.638,.647,.660,.668,.672,.684,.687,.689,.763,.790,.812,.830,.848,.853};
*/
Double_t e0[nd] = {.553,.585,.607,.630,.645,.660,.677,.696,.712,.724,.743,.782,.800,.809,.821,-.821,-.865};
Double_t e1[nd] = {.557,.578,.602,.619,.638,.646,.667,.698,.712,.718,.729,.768,.787,.797,.807,.810,.809};
Double_t e2[nd] = {.543,.587,.608,.623,.632,.646,.667,.693,.701,.711,.721,.771,.788,.801,.806,.813,.814};

Double_t e3[nd] = {.631,.615,.615,.629,.635,.645,.652,.666,.677,.688,.700,.746,.773,.789,.809,-.851,-.858};
Double_t e4[nd] = {.637,.611,.610,.615,.630,.634,.655,.671,.680,.683,.698,.739,.765,.779,.797,.806,.812};
Double_t e5[nd] = {.644,.616,.613,.623,.630,.635,.649,.664,.673,.675,.687,.738,.765,.781,.797,.809,-.853};

TGraphErrors* gr[6];

for(int i=0; i<3; i++) {
  if(i == 0) {
    gr[i] = new TGraphErrors(nd-2,ener,e0,err0,err0);
    gr[i+3] = new TGraphErrors(nd-2,ener,e3,err0,err0);
  } else if(i == 1) { 
    gr[i] = new TGraphErrors(nd,ener,e1,err0,err0);
    gr[i+3] = new TGraphErrors(nd,ener,e4,err0,err0);
  } else if(i == 2) { 
    gr[i] = new TGraphErrors(nd,ener,e2,err0,err0);
    gr[i+3] = new TGraphErrors(nd-1,ener,e5,err0,err0);
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
