// -------------------------------------------------------------------
// $Id: plot.C 70323 2013-05-29 07:57:44Z gcosmo $
// -------------------------------------------------------------------
//
// *********************************************************************
// To execute this macro under ROOT after your simulation ended, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
//   3 - OR type directly 'root plot.C'
// *********************************************************************

{
gROOT->Reset();

gStyle->SetPalette(1);
gROOT->SetStyle("Plain");
gStyle->SetOptStat(000000);
	
c1 = new TCanvas ("c1","",60,60,800,800);
c1->Divide(1,1);

FILE * fp = fopen("s.txt","r");

Float_t radius,E,s,ss;
Int_t ncols = 0;
Int_t nlines = 0;

TNtuple *ntuple = new TNtuple("ntuple","s","radius:E:s:ss");
while (1) 
{
  ncols = fscanf(fp,"%f %f %f %f",&radius,&E,&s,&ss);
  if (ncols < 0) break;
  ntuple->Fill(radius,E,s,ss);
  nlines++;
}
fclose(fp);
 
c1->cd(1);
gPad->SetLogx();
gPad->SetLogy();

TH2F * h2 = new TH2F ("h2","",2,99.999,1e4,2,1e2,1e4);
h2->Draw();
ntuple->SetMarkerStyle(20);
ntuple->SetMarkerSize(1.);
ntuple->Draw("s:E","","LPsame");

h2->GetXaxis()->SetLabelSize(0.025);
h2->GetYaxis()->SetLabelSize(0.025);
h2->GetXaxis()->SetTitleSize(0.035);
h2->GetYaxis()->SetTitleSize(0.035);
h2->GetXaxis()->SetTitleOffset(1.4);
h2->GetYaxis()->SetTitleOffset(1.4);
h2->GetXaxis()->SetTitle("E (eV)");
h2->GetYaxis()->SetTitle("S (Gy/Bq.s)");

}
