// Draw Longitudinal Shower Profile

#include "TH1D.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"


TCanvas *cmpL() {

  //TString dataDir("/net/llrdata1.in2p3.fr/data/DATA/data.harpo/BH5D/");
  TString dataDir("./");

  TString dataFile[2];
  dataFile[0]  =  "gflash00.root";
  dataFile[1]  =  "gflash01.root";


  TFile *f[2];
  TProfile *p[2];
  //  UInt_t col[] = { kRed, kBlue, kGreen, kViolet };
  UInt_t col[] = { kRed, kGreen, kBlue, kViolet };
  //  UInt_t mark[] = { 20, 32, 21, 22};
  UInt_t mark[] = { 20, 21, 32, 22};
//
   THStack *hs = new THStack("hs","");

   for (UInt_t i = 0; i < 2; i++) {
     TString tmp = dataDir + dataFile[i];
     f[i] = TFile::Open(tmp);
     if (!f[i]) return 0;
     p[i] = (TProfile *) (f[i])->Get("p0");
     //ShiftUp(p[i],0.0);
   // histogram
     p[i]->SetFillColor(col[i]);
     p[i]->SetLineStyle(i+2);
     p[i]->SetLineColor(col[i]);
     p[i]->SetLineWidth(2);
     p[i]->SetMarkerStyle(mark[i]);
     p[i]->SetMarkerColor(col[i]);
     p[i]->SetMarkerSize(1.5);
     hs->Add(p[i]);
   }
   TCanvas *cst = new TCanvas("cst","stacked hists",10,10,800,700);
   //gPad->SetGrid();
   hs->Draw("p,nostack");
   hs->SetTitle("Longitudinal Profile");
   hs->GetXaxis()->SetTitle("Depth (RadLen)");
   hs->GetYaxis()->SetTitle("E/E_{tot} (%/RadLen)");

   cst->RedrawAxis();
   cst->Update();

   TLegend* legend = new TLegend(0.79,0.84,0.94,0.94);
   legend->AddEntry(p[0],"full ","p");
   legend->AddEntry(p[1],"gflash","p");
   legend->Draw();
   cst->Update();

   //img->FromPad(c, 10, 10, 300, 200);

   // TImage *img = TImage::Create();

   // img->FromPad(cst);

   // img->WriteImage("hist0-640MeV-p.png");

   // for (UInt_t i = 0; i < 4; i++) {
   // f[i]->Close();
   // }
   return cst;
}
