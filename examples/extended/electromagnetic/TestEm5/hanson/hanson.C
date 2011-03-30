
{
   gROOT->Reset();

   // Draw histograms fill by Geant4 TestEm11 simulation
   TFile f1("./9um.opt0.root");
   TH1D* h1 = (TH1D*) f1.Get("12");
   //h1->SetTitle("Angular distribution of 15.7 MeV e- after 19um Au foil");
   h1->SetTitle("Angular distribution of 15.7 MeV e- after 9um Au foil");   
   h1->GetXaxis()->SetTitle("theta (deg)");
   h1->GetYaxis()->SetTitle("dN/dOmega");
   h1->SetStats(kFALSE);  // Eliminate statistics box
   h1->SetLineColor(4);   // blue
   h1->Draw("HIST");
   
   TFile f2("./9um.gs.root");
   TH1D* h2 = (TH1D*) f2.Get("12");
   h2->SetStats(kFALSE);  // Eliminate statistics box
   h2->SetLineColor(3);   // green
   h2->Draw("SAME HIST");
/* data
* angle distribution of  15.7 MeV electrons
* transmitted through thin gold foils.
* A.O.Hanson et al. Phys.Rev.84 (1951) page 634.
*/

   ifstream in;
   //in.open("./19um.ascii");
   in.open("./9um.ascii");
   
   TMarker *pt;
   Double_t x, y;
   // First indicate number of data
   int nbdata = 0;
   in >> nbdata;
   for ( int i = 0 ; i < nbdata ; i++ ) {
      in >> x >> y ;
      if (!in.good()) break;
      pt = new TMarker(x,y,22); // 22 for triangle TMatker
      pt->SetMarkerColor(kRed);
      pt->Draw();
   }
   in.close();

   // Print the histograms legend
   TLegend* legend = new TLegend(0.6,0.55,0.8,0.68);
   legend->AddEntry(h1,"9.4-opt0 ","l");
   legend->AddEntry(h2,"9.4-gs ","l");   
   legend->AddEntry(pt,"Hanson data","P");
   legend->Draw();
}
