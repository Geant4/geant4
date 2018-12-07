// -------------------------------------------------------------------
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

Float_t radius,E,thick,sc,ssc,sn,ssn;
Int_t ncols = 0;
Int_t nlines = 0;

TNtuple *ntuple = new TNtuple("ntuple","s","radius:thick:E:sc:ssc:sn:ssn");
while (1) 
{
  ncols = fscanf(fp,"%f %f %f %f %f %f %f",&radius,&thick,&E,&sc,&ssc,&sn,&ssn);
  if (ncols < 0) break;
  ntuple->Fill(radius,thick,E,sc,ssc,sn,ssn);
  nlines++;
}
fclose(fp);
 
c1->cd(1);
gPad->SetLogx();
gPad->SetLogy();

TH2F * h2 = new TH2F ("h2","",2,99.999,1e4,2,1e0,1e3);
h2->Draw();

ntuple->SetMarkerStyle(24);
ntuple->SetMarkerSize(1.);
ntuple->Draw("sc:E","","LPsame");

ntuple->SetMarkerStyle(20);
ntuple->SetMarkerSize(1.);
ntuple->Draw("sn:E","","LPsame");

h2->GetXaxis()->SetLabelSize(0.025);
h2->GetYaxis()->SetLabelSize(0.025);
h2->GetXaxis()->SetTitleSize(0.035);
h2->GetYaxis()->SetTitleSize(0.035);
h2->GetXaxis()->SetTitleOffset(1.4);
h2->GetYaxis()->SetTitleOffset(1.4);
h2->GetXaxis()->SetTitle("E (eV)");
h2->GetYaxis()->SetTitle("S (Gy/Bq.s)");

TText *pt1 = new TText(200,200,"Cytoplasm");
pt1->Draw("SAME");

TText *pt2 = new TText(400,15,"Nucleus");
pt2->Draw("SAME");

}
