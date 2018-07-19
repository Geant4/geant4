//
// *********************************************************************
// To execute this macro under ROOT after your simulation ended, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
//   3 - OR directly type 'root plot.C'
// *********************************************************************

{

gROOT->Reset();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");
gStyle->SetOptStat(00000);

c1 = new TCanvas ("c1","Range",60,60,700,700);
c1->Divide(1,1);
c1->cd(1);
gPad->SetLogx();
gPad->SetLogy();

TH2F * h2  = new TH2F("h2","",2,9.99,1e+4,2,1e-1,1e+4);
h2->Draw();
h2->GetXaxis()->SetLabelSize(0.025);
h2->GetYaxis()->SetLabelSize(0.025);
h2->GetXaxis()->SetTitleSize(0.035);
h2->GetYaxis()->SetTitleSize(0.035);
h2->GetXaxis()->SetTitleOffset(1.4);
h2->GetYaxis()->SetTitleOffset(1.4);
h2->GetXaxis()->SetTitle("E (eV)");
h2->GetYaxis()->SetTitle("Distance (nm)");

FILE * fp = fopen("range.txt","r");

Float_t e,track,strack,proj,sproj,pene,spene;
Int_t ncols = 0;
Int_t nlines = 0;

TNtuple *ntuple = new TNtuple("ntuple","range","e:track:strack:proj:sproj:pene:spene");

while (1) 
{
      ncols = fscanf(fp,"%f %f %f %f %f %f %f",&e,&track,&strack,&proj,&sproj,&pene,&spene);
      if (ncols < 0) break;
      ntuple->Fill(e,track,strack,proj,sproj,pene,spene);
      nlines++;
}
   
fclose(fp);

ntuple->SetLineWidth(3);
ntuple->SetLineColor(2);
ntuple->Draw("track:e","","Lsame");
ntuple->SetLineColor(3);
ntuple->Draw("pene:e","","Lsame");
ntuple->SetLineColor(4);
ntuple->Draw("proj:e","","Lsame");
}
