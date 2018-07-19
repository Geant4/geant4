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

FILE * fp = fopen("wvalue.txt","r");

Float_t E,nbioni,snbioni,w,sw;
Int_t ncols = 0;
Int_t nlines = 0;

TNtuple *ntuple = new TNtuple("ntuple","w","E:nbioni:snbioni:w:sw");
while (1) 
{
  ncols = fscanf(fp,"%f %f %f %f %f",&E,&nbioni,&snbioni,&w,&sw);
  if (ncols < 0) break;
  ntuple->Fill(E,nbioni,snbioni,w,sw);
  nlines++;
}
fclose(fp);
 
c1->cd(1);
gPad->SetLogx();
gPad->SetLogy();

TH2F * h2 = new TH2F ("h2","",2,9.99,1e3,2,9.99,1e3);
h2->Draw();
ntuple->SetMarkerStyle(20);
ntuple->SetMarkerSize(1.);
ntuple->Draw("w:E","","LPsame");

h2->GetXaxis()->SetLabelSize(0.025);
h2->GetYaxis()->SetLabelSize(0.025);
h2->GetXaxis()->SetTitleSize(0.035);
h2->GetYaxis()->SetTitleSize(0.035);
h2->GetXaxis()->SetTitleOffset(1.4);
h2->GetYaxis()->SetTitleOffset(1.4);
h2->GetXaxis()->SetTitle("E (eV)");
h2->GetYaxis()->SetTitle("W (eV)");

}
