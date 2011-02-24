{ 
gROOT->ProcessLine(".x $G4INSTALL/tests/test46/scripts/StyleCut.C");

gStyle->SetTitleSize(0.04, "x");
gStyle->SetTitleSize(0.04, "y");
gStyle->SetMarkerSize(1.2);

TLegend* legc = new TLegend(0.45, 0.7, 0.9, 0.9);

Int_t i=0;

Double_t y1[3] = {0.3, 0.3, 0.8};
Double_t y2[3] = {1, 1, 1.2};

TString ahed[3] = {"MIP fraction of #pi^{-} in ECAL","MIP fraction of p in ECAL","MIP fraction ratio"};
TString bhed[6] = {"QGSP_BERT_EMV pi-","QGSP_BERT_EMV p","pi-","QGSP_BERT_EMLGG pi-","QGSP_BERT_EMLGG p","p"};
TString ched[3] = {"MIP fraction","MIP fraction","QGSP_BERT_EMLGG/QGSP_BERT_EMV"};
TString dhed[3] = {"mip_pi-.gif","mip_p.gif","mip_ratio.gif"};

gPad->SetLogx();
TH1F* hh = gPad->DrawFrame(0.8,y1[i],500,y2[i],ahed[i]);
hh->GetXaxis()->SetTitle("E (GeV)");
hh->GetYaxis()->SetTitle(ched[i]);
gPad->SetGrid();
hh->Draw("AXIS SAME");

Int_t nd = 17;
Double_t ener[nd]  = {1,2,3,4,5,6,7,8,9,10,12,20,30,50,100,200,300};
Double_t err0[nd] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

// 9.2ref02  pi-
//Double_t e0[nd] = {.9744,.7693,.5577,.492,.456,.475,.4857,.487,.472,.3964,.3789,.4,.4065,.4042,.3877,.3622,.3528};
// 9.3beta pi-
Double_t e0[nd] = {.9736,.7713,.5278,.4473,.4251,.4186,.4142,.4157,.4097,.3965,.3795,.4054,.4118,.4023,.3916,.3669,.3449};
// 9.3beta GG  pi-
Double_t e1[nd] = {.9749,.7724,.5274,.4571,.426,.4088,.4263,.4202,.4058,.3779,.3874,.3999,.4054,.4125,.3917,.3607,.3353};
// 9.2ref02  p
//Double_t e2[nd] = {.9902,.8095,.5823,.4744,.4449,.4374,.4311,.429,.4138,.3566,.3472,.3406,.3322,.3331,.3295,.3178,.313};
// 9.3beta p
Double_t e2[nd] = {.9917,.8161,.5838,.4696,.4379,.4103,.4048,.4032,.3792,.3448,.3359,.33,.3341,.3379,.317,.3209,.3223};
// 9.3beta GG p
Double_t e3[nd] = {.9913,.8006,.5628,.4718,.4192,.4166,.3987,.3841,.3848,.3472,.3465,.3332,.3423,.3366,.3335,.3313,.3182};
Double_t e4[nd], e5[nd];

for(int j=0; j<nd; j++) {
  e4[j] = e1[j]/e0[j];
  e5[j] = e3[j]/e2[j];
}

TGraphErrors* gr[6];

gr[0] = new TGraphErrors(nd,ener,e0,err0,err0);
gr[1] = new TGraphErrors(nd,ener,e2,err0,err0);
gr[2] = new TGraphErrors(nd,ener,e4,err0,err0);
gr[3] = new TGraphErrors(nd,ener,e1,err0,err0);
gr[4] = new TGraphErrors(nd,ener,e3,err0,err0);
gr[5] = new TGraphErrors(nd,ener,e5,err0,err0);

gr[i]->SetMarkerColor(col[0]);
gr[i]->SetLineColor(col[0]);
  //  gr[i]->SetTextColor(col[i]);
gr[i]->SetMarkerStyle(mar[0]);
gr[i]->Draw("PL SAME 9");
legc->AddEntry(gr[i],bhed[i],"p");
gr[i+3]->SetMarkerColor(col[1]);
gr[i+3]->SetLineColor(col[1]);
  //  gr[i+3]->SetTextColor(col[i]);
gr[i+3]->SetMarkerStyle(mar[1]);
gr[i+3]->Draw("PL SAME 9");
legc->AddEntry(gr[i+3],bhed[i+3],"p");
legc->Draw("SAME");

c1.Print(dhed[i]);


}
