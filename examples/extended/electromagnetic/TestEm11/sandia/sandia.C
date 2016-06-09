
{
   gROOT->Reset();

   // Draw histogram fill by Geant4 TestEm11 simulation
   TFile f("./Al.314keV.root");
   TH1D* h1d = (TH1D*) f.Get("8");
   h1d->SetTitle("Depth dose distribution of 314 keV e- in Al");
   h1d->GetXaxis()->SetTitle("Edep (Mev.cm2/g) along x/r0                    x/r0");
   h1d->GetYaxis()->SetTitle("MeV*cm2/g");
   h1d->SetStats(kFALSE);  // Eliminate statistics box
   h1d->Draw("HIST");

/* data
* G.J.Lockwood et al.
*     Sandia report SAND79-0414.UC-34a, February 1987
* O.Kadri et al. NIM B 258 (2007) 381
*/

   ifstream in;
   in.open("./data/Al_314keV.ascii");

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
   legend->AddEntry(h1d,"G4-9.2-ref-04 ","l");
   legend->AddEntry(pt,"Sandia data","P");
   legend->Draw();
}
