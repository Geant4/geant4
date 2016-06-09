{
   gROOT->Reset();

   // Draw histogram fill by Geant4 TestEm5 simulation
   TFile f("./tramu.root");
   TH1D* h1d = (TH1D*) f.Get("10");
   h1d->SetTitle("1 TeV muon in 3 m iron : kinetic energy at exit (GeV)");
   h1d->GetXaxis()->SetTitle("Ekine (GeV)");
   h1d->GetYaxis()->SetTitle("nb/GeV");
   h1d->SetStats(kFALSE);  // Eliminate statistics box
   h1d->Draw("HIST");
   
/*
* muon 1 TeV/c in 3 m Iron
* Particle Data Group. Physics Letters B 592 (2004) page 251
* distribution of the muon kinetic energy
* (from 950 GeV to 1000 GeV by bin of 0.5 GeV --> 100 bins) 
*/

   ifstream in;
   in.open("mars14.ascii");

   // Create a new histogram with mars14.acsii values
   int nb_bins = 100;   
   float x_min = 950.;
   float x_max = 1000.;
   TH1F* h1f = new TH1F("h1f","",nb_bins,x_min,x_max);

   Float_t x, y;
   while (1) {
      in >> x >> y ;
      if (!in.good()) break;
      h1f->Fill(x,y);
   }
   in.close();

   // Draw histogram fill by mars14.acsii values
   h1f->SetLineColor(2);
   h1f->Draw("SAME");

   // Print the histograms legend
   TLegend* legend = new TLegend(0.2,0.55,0.45,0.70);
   legend->AddEntry(h1d,"G4-9.1-ref-08++","l");
   legend->AddEntry(h1f,"Mars14 simul ","L");
   legend->Draw();


}
