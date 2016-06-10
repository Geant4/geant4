// *********************************************************************
// To execute this macro under ROOT after your simulation ended, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
//
// Author: Sebastien Incerti, CNRS, France
// Date: 25 Feb. 2015
// The Geant4-DNA collaboration
// *********************************************************************

{

gROOT->Reset();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");
gStyle->SetOptStat(00000);

//***************************************
//***************************************
// MAKE YOUR SELECTION OF B, min and max
// for log histograms
//***************************************
//***************************************

Int_t B=1000; // bins per decade
Int_t min=-2; // minimum x-axis value as 10^min
Int_t max=2; // maximum x-axis value as 10^max

//
FILE *fp = fopen("yz.root","r");
if( fp ) {
  // exists
  cout << "*** Notice: the output file yz.root exists ***"<< endl;
  fclose(fp);
} else {
  cout << "*** Notice: the output file yz.root does not exist ***"<< endl;
  cout << "***         it will be created from merged ROOT files ***"<< endl;
  system ("rm -rf yz.root");
  system ("hadd yz.root yz_*.root");}
//

c1 = new TCanvas ("c1","",60,60,800,600);
c1.Divide(3,2);

/*
// for testing only
FILE * fp = fopen("yz.txt","r");
Float_t radius,y,z;
Int_t ncols = 0;
Int_t nlines = 0;

TNtuple *ntuple = new TNtuple("ntuple","micro","radius:y:z");
while (1) 
   {
      ncols = fscanf(fp,"%f %f %f",&radius,&y,&z);
      if (ncols < 0) break;
      ntuple->Fill(radius,y,z);
      nlines++;
   }
fclose(fp);
*/

TFile f("yz.root"); 

TNtuple* ntuple;
ntuple = (TNtuple*)f->Get("yz"); 
   
// --->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->

//plot f(y)
c1.cd(1);

ntuple->Draw("y>>hfy","","");

hfy->Scale (1./(hfy->GetEntries()*hfy->GetBinWidth(1))); // DIVIDE BY BIN WIDTH !!!
hfy->SetTitle("f(y) (um/keV)");
hfy->GetXaxis()->SetTitle("y (keV/um)");
hfy->SetFillColor(1);
hfy->Draw("");

//check normalization

Double_t norm=0;

for (Int_t j=0;j<hfy->GetNbinsX(); j++) 
 norm=norm+hfy->GetBinContent(j)*hfy->GetBinWidth(1); // MULTIPLY BY BIN WIDTH !!!

cout << endl;
cout << "**** I - Results from lin-lin histograms ****" << endl;
cout << endl;
cout << "---> sum of f(y)dy =" << norm << endl;

//plot y*f(y)
c1.cd(2);

TH1F *hyfy = (TH1F*)hfy->Clone("hyfy");

for (Int_t i=0;i<hyfy->GetNbinsX(); i++) 
 hyfy->SetBinContent(i,hyfy->GetBinCenter(i)*hyfy->GetBinContent(i));

hyfy->SetLineColor(2);
hyfy->SetFillColor(2);
hyfy->SetTitle("y*f(y)");
hyfy->GetXaxis()->SetTitle("y (keV/um)");
hyfy->Draw("");

// --->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->

//calculate yF

Double_t yF=0;

for (Int_t j=0;j<hyfy->GetNbinsX(); j++) 
 yF=yF+hyfy->GetBinContent(j)*hyfy->GetBinWidth(1); // MULTIPLY BY BIN WIDTH !!!

cout << "---> yF=" << yF << " keV/um" << endl;

//plot y*f(y)/yF = d(y) (cf. Burigo et al., NIMB 320 (2014) 89-99)
c1.cd(4);

TH1F *hdy = (TH1F*)hyfy->Clone("hdy");

hdy->Scale (1./yF);

hdy->SetLineColor(4);
hdy->SetFillColor(4);
hdy->SetTitle("d(y) (um/keV)");
hdy->GetXaxis()->SetTitle("y (keV/um)");
hdy->Draw("");

//check normalization

Double_t normD=0;

for (Int_t j=0;j<hdy->GetNbinsX(); j++) 
 normD=normD+hdy->GetBinContent(j)*hdy->GetBinWidth(1); // MULTIPLY BY BIN WIDTH !!!

cout << "---> sum of d(y)dy =" << normD << endl;

//plot y*d(y) 
c1.cd(5);

TH1F *hydy = (TH1F*)hdy->Clone("hydy");

for (Int_t k=0;k<hydy->GetNbinsX(); k++) 
 hydy->SetBinContent(k,hydy->GetBinCenter(k)*hdy->GetBinContent(k));

hydy->SetLineColor(3);
hydy->SetFillColor(3);
hydy->SetTitle("y*d(y)");
hydy->GetXaxis()->SetTitle("y (keV/um)");
hydy->Draw("");

//calculate yD

Double_t yD=0;

for (Int_t l=0;l<hydy->GetNbinsX(); l++) 
 yD=yD+hydy->GetBinContent(l)*hydy->GetBinWidth(1); // MULTIPLY BY BIN WIDTH !!!

cout << "---> yD=" << yD << " keV/um" << endl;

// --->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->
// --->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->
// --->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->

//log plot of y*f(y)

//goto end;

c1.cd(3);

cout << endl;
cout << "**** II - Results from log-lin histograms ****" << endl;
cout << endl;
cout << " You have selected "<< B << " bins per decade, from 10^" << min << " (a.u.) to 10^"<< max << " (a.u.) "<< endl;
cout << endl;


Int_t bins=B*(max-min);

TH1F *hlogyfy = new TH1F("1","1",bins,min,max); 
TAxis *axis = hlogyfy->GetXaxis(); 

TH1F *hlogy2fy = new TH1F("2","2",bins,min,max); 
TAxis *axis2 = hlogy2fy->GetXaxis(); 

TH1F *hlogydy = new TH1F("3","3",bins,min,max); 
TAxis *axis3 = hlogydy->GetXaxis(); 

TH1F *hlogy3fy = new TH1F("4","4",bins,min,max); 
TAxis *axis4 = hlogy3fy->GetXaxis(); 


Axis_t from = axis->GetXmin();
Axis_t to = axis->GetXmax();
Axis_t width = (to - from) / bins;
//cout << "*** width=" << width << endl;
//cout << "*** bins=" << bins << endl;
Axis_t *new_bins = new Axis_t[bins + 1];
for (int i = 0; i <= bins; i++) 
{
  new_bins[i] = TMath::Power(10, from + i * width);
} 

axis->Set(bins, new_bins); 
axis2->Set(bins, new_bins); 
axis3->Set(bins, new_bins); 
axis4->Set(bins, new_bins); 

delete new_bins; 

//
/*
//for testing only
FILE * fp2 = fopen("yz.txt","r");
while (1) 
   {
      ncols = fscanf(fp2,"%f %f %f",&radius,&y,&z);
      if (ncols < 0) break;
      hlogyfy->Fill(y);
      hlogy2fy->Fill(y);
      hlogydy->Fill(y);
      hlogy3fy->Fill(y);
      hlogzfz->Fill(z);
      hlogz2fz->Fill(z);
      hlogzdz->Fill(z);
      hlogz3fz->Fill(z);
      nlines++;
   }
fclose(fp2);
*/
//

Double_t radius,y,z;

ntuple->SetBranchAddress("radius",&radius);
ntuple->SetBranchAddress("y",&y);
ntuple->SetBranchAddress("z",&z);
Int_t nentries = (Int_t)ntuple->GetEntries();
for (Int_t i=0; i<nentries; i++) 
{
 ntuple->GetEntry(i);
      hlogyfy->Fill(y);
      hlogy2fy->Fill(y);
      hlogydy->Fill(y);
      hlogy3fy->Fill(y);

}

//

hlogyfy->Scale(1./( hlogyfy->GetEntries() ));
hlogy2fy->Scale(1./( hlogy2fy->GetEntries() ));
hlogydy->Scale(1./( hlogydy->GetEntries() ));
hlogy3fy->Scale(1./( hlogy3fy->GetEntries() ));


//plot y*f(y)

gPad->SetLogx();

for (Int_t i=0;i<hlogyfy->GetNbinsX(); i++) 
{
 hlogyfy->SetBinContent(i,
//(log(10)/B)*
hlogyfy->GetBinCenter(i)
*hlogyfy->GetBinContent(i)
/(TMath::Power(10, from + (i+1) * width)-TMath::Power(10, from + i * width))
);
}

hlogyfy->SetLineColor(2);
hlogyfy->SetFillColor(2);
hlogyfy->SetTitle("y*f(y)");
hlogyfy->GetXaxis()->SetTitle("y (keV/um)");
hlogyfy->Draw("");

//check normalization

Double_t normLogfy=0;

for (Int_t j=0;j<hlogyfy->GetNbinsX(); j++) 
 normLogfy=normLogfy+(log(10)/B)*hlogyfy->GetBinContent(j); 

cout << "---> sum of Log f(y)dy =" << normLogfy << endl;

//calculate yF
//requires plot of log y2fy and integration of this plot

for (Int_t i=0;i<hlogy2fy->GetNbinsX(); i++) 
{
 hlogy2fy->SetBinContent(i,
//(log(10)/B)*
hlogy2fy->GetBinCenter(i)
*hlogy2fy->GetBinCenter(i)
*hlogy2fy->GetBinContent(i)
/(TMath::Power(10, from + (i+1) * width)-TMath::Power(10, from + i * width))
);
}

Double_t logyF=0;

for (Int_t j=0;j<hlogy2fy->GetNbinsX(); j++) 
 logyF=logyF+(log(10)/B)*hlogy2fy->GetBinContent(j); 

cout << "---> yF =" << logyF << " keV/um" << endl;

// --->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->

//plot y*d(y)

c1.cd(6);

for (Int_t i=0;i<hlogydy->GetNbinsX(); i++) 
{
 hlogydy->SetBinContent(i,
//(log(10)/B)*
hlogydy->GetBinCenter(i)
*hlogydy->GetBinCenter(i)
*hlogydy->GetBinContent(i)
/(TMath::Power(10, from + (i+1) * width)-TMath::Power(10, from + i * width))
/logyF
);
}

hlogydy->SetLineColor(3);
hlogydy->SetFillColor(3);
hlogydy->SetTitle("y*d(y)");
hlogydy->GetXaxis()->SetTitle("y (keV/um)");
gPad->SetLogx();
hlogydy->Draw("");

//check normalization

Double_t normLogdy=0;

for (Int_t j=0;j<hlogydy->GetNbinsX(); j++) 
 normLogdy=normLogdy+(log(10)/B)*hlogydy->GetBinContent(j);

cout << "---> sum of Log d(y)dy =" << normLogdy << endl;

//calculate yD
//requires plot of log y3fy and integration of this plot

for (Int_t i=0;i<hlogy3fy->GetNbinsX(); i++) 
{
 hlogy3fy->SetBinContent(i,
//(log(10)/B)*
hlogy3fy->GetBinCenter(i)
*hlogy3fy->GetBinCenter(i)
*hlogy3fy->GetBinCenter(i)
*hlogy3fy->GetBinContent(i)
/(TMath::Power(10, from + (i+1) * width)-TMath::Power(10, from + i * width))
/logyF
);
}

Double_t logyD=0;

for (Int_t j=0;j<hlogy3fy->GetNbinsX(); j++) 
 logyD=logyD+(log(10)/B)*hlogy3fy->GetBinContent(j); 

cout << "---> yD =" << logyD << " keV/um" << endl;

// --->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->--->

//
end:
cout << endl;
}

