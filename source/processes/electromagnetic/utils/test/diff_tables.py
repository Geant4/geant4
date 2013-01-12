{
   gROOT->Reset();
#include "Riostream.h"
#include "TSystem.h"

   Char_t buffer[256];

   const Int_t n = 1100;

   Double_t En[n], Y1[n], Y2[n], diff[n], pres[n];

   TCanvas *c1 = new TCanvas("c1", "c1",6,6,800,600);
   gStyle->SetOptStat(0);
   c1->SetLogx();
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(0);
   c1->SetFrameBorderMode(0);
   

   ifstream in;

   in.open("DEDX.hIoni.proton.asc_7bins_spl.out");
//   in.open("DEDX.eIoni.e-.asc_spl.out");
   if (!in.good()) break;

// Ignore first blank line

   in.getline(buffer,256);

   for (Int_t i=0; i<n; i++) {
     in >> En[i] >> Y1[i] >> Y2[i] >> diff[i];
     if (!in.good()) break;
     pres[i] = (Y1[i]/Y2[i] - 1)*100; 
  }


/*   Int_t n = 0;
   while(!(in.eof())) {
     n++;
     in >> En[n] >> Y1[n] >> Y2[n] >> diff[n];
     pres[n] = (Y1[n]/Y2[n] - 1)*100;
  }
*/
   gr = new TGraph(n,En,pres);
   gr->SetTitle("proton in Pb, table 100 pts / 7 pts_spl per order");
   gr->SetMarkerStyle(22);
   gr->SetMarkerSize(0.8);
   gr->GetYaxis()->SetLabelFont(132);
   gr->GetYaxis()->SetLabelSize(0.04);
   gr->GetXaxis()->SetLabelFont(132);
   gr->GetXaxis()->SetLabelSize(0.04);
   gr->GetXaxis()->SetTitle("E, MeV");
   gr->GetYaxis()->SetTitle("dE/dx: (100pts/7pts_spl - 1), %");
   gr->GetXaxis()->SetTitleOffset(1.2);
   gr->Draw("AP");

   // draw the legend
/*   TLegend *leg = new TLegend(0.2,0.65,0.61,0.86);
   leg->SetTextFont(52);
   leg->SetTextSize(0.035);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillStyle(0);
   leg->SetMargin(0.4);
   leg->SetBorderSize(1);
   TLegendEntry *entry=leg->AddEntry("NULL","p (110 MeV) in H_{2}O","h");

   entry->SetTextAlign(22);
   entry=leg->AddEntry("NULL","Geant 4.9.1ref01","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetTextAlign(22);
   entry->SetTextColor(1);
   entry=leg->AddEntry(gr,"exp. data","p");  
   leg->Draw();
*/
   c1->Modified();
   c1->cd();
}
