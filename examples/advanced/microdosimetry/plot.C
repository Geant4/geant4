// -------------------------------------------------------------------
// $Id: plot.C,v 1.4 2010-09-13 08:42:59 sincerti Exp $
// -------------------------------------------------------------------
//
// *********************************************************************
// To execute this macro under ROOT, 
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// This macro needs the track.txt file
// *********************************************************************
{
gROOT->Reset();
gStyle->SetPalette(1);
gROOT->SetStyle("Plain");
Double_t scale;
	
c1 = new TCanvas ("c1","",20,20,1000,500);
c1.Divide(2,1);

FILE * fp = fopen("track.txt","r");
Float_t process,part,x,y,z;
Int_t ncols=0;
Int_t nlines = 0;

TNtuple *ntuple = new TNtuple("result","ntuple","process:part:x:y:z");

while (1) 
{
   ncols = fscanf(fp,"%f %f %f %f %f",&part,&process,&x,&y,&z);
   if (ncols < 0) break;
   ntuple->Fill(process,part,x,y,z);
   nlines++;
}
fclose(fp);
      
c1.cd(1);
  gStyle->SetOptStat(000000);
  ntuple->Draw("process","");
  ntuple->SetFillColor(2);
  ntuple->Draw("process","process==12||process==15||process==17||process==22||process==25||process==29","same");
  ntuple->SetFillColor(3);
  ntuple->Draw("process","process==11","same");
  ntuple->SetFillColor(4);
  ntuple->Draw("process","process==13||process==18||process==20||process==23||process==26||process==30||process==32||process==33","same");
  ntuple->SetFillColor(5);
  ntuple->Draw("process","process==19||process==24||process==27","same");
  ntuple->SetFillColor(6);
  ntuple->Draw("process","process==21||process==28||process==31","same");
  
  gPad->SetLogy();
  htemp->GetXaxis()->SetLabelSize(0.025);
  htemp->GetYaxis()->SetLabelSize(0.025);
  htemp->GetZaxis()->SetLabelSize(0.025);
  htemp->GetXaxis()->SetTitleSize(0.035);
  htemp->GetYaxis()->SetTitleSize(0.035);
  htemp->GetXaxis()->SetTitleOffset(1.4);
  htemp->GetYaxis()->SetTitleOffset(1.4);
  htemp->GetXaxis()->SetTitle("Process");
  htemp->GetYaxis()->SetTitle("");
  htemp->SetTitle("Processes");

c1.cd(2);
  ntuple->SetMarkerColor(2);
  ntuple->Draw("x:y:z/1000","part==1");

  ntuple->SetMarkerColor(4);
  ntuple->SetMarkerSize(4);
  ntuple->Draw("x:y:z/1000","part==5","same");

  htemp->GetXaxis()->SetLabelSize(0.025);
  htemp->GetYaxis()->SetLabelSize(0.025);
  htemp->GetZaxis()->SetLabelSize(0.025);
  htemp->GetXaxis()->SetTitleSize(0.035);
  htemp->GetYaxis()->SetTitleSize(0.035);
  htemp->GetZaxis()->SetTitleSize(0.035);
  htemp->GetXaxis()->SetTitleOffset(1.6);
  htemp->GetYaxis()->SetTitleOffset(1.6);
  htemp->GetZaxis()->SetTitleOffset(1.6);
  htemp->GetXaxis()->SetTitle("z (micrometer)");
  htemp->GetYaxis()->SetTitle("x (nanometer)");
  htemp->GetZaxis()->SetTitle("y (nanometer)");
  htemp->SetTitle("Track Structure in liquid water");
  
}
