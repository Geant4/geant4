
{
   gROOT->Reset();

   // Draw histogram fill by Geant4 TestBruce simulation
   TFile f("./nofoil.13MeV.root");
   TH1D* h1d = (TH1D*) f.Get("4");
   h1d->SetTitle("Fluence distribution of 13 MeV e- in nofoil");
   h1d->GetXaxis()->SetTitle("r (mm)");
   h1d->GetYaxis()->SetTitle("Fluence");
   h1d->SetStats(kFALSE);  // Eliminate statistics box
   h1d->Draw("HIST");

/* data
* Bruce et al.
*/

   ifstream in;
   in.open("../data/nofoil.13MeV.ascii");

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
   legend->AddEntry(h1d,"G4-9.2-cand-01 ","l");
   legend->AddEntry(pt,"NRCC data","P");
   legend->Draw();
}
