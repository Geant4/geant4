{
   gROOT->Reset();

   // Draw histogram fill by Geant4 TestEm5 simulation
   TFile f("./tramu.root");
   TH1D* h1d = (TH1D*) f.Get("10");
   h1d->SetTitle("1 TeV muon in 3 m iron : kinetic energy at exit (GeV)");
   h1d->GetXaxis()->SetTitle("Ekine(GeV)");
   h1d->GetYaxis()->SetTitle("1/GeV");
   h1d->SetStats(kFALSE);  // Eliminate statistics box
   h1d->Draw();

   // data2.ascii came from Mars14
   ifstream in;
   in.open("mars14.ascii");

   // Create a new histogramm which mars14.acsii values
   int x_min = 950;
   int x_max = 1000;
   int nb_bins = 100;
   TH1F *h1f = new TH1F("h1f","",nb_bins,x_min,x_max);

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
   TLegend *legend = new TLegend(0.2,0.55,0.36,0.68);
   legend->AddEntry(h1d,"G4-9.1","l");
   legend->AddEntry(h1f,"Mars14","L");
   legend->Draw();


}
