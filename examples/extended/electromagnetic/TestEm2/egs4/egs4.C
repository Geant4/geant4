{
   gROOT->Reset();

   // Draw histogram fill by Geant4 TestEm2 simulation
   TFile f("./04jan08.root");
   TH1D* h1d = (TH1D*) f.Get("4");
   h1d->SetTitle("30 GeV e- on 20X0 Fe : energy deposit, longitudinal profil");
   h1d->GetXaxis()->SetTitle("depth(X0)");
   h1d->GetYaxis()->SetTitle("(100/E0) (dE/dt)");
   h1d->SetStats(kFALSE);  // Eliminate statistics box
   h1d->Draw();

   // pdg.ascii came from egs4 simulation
   ifstream in;
   in.open("pdg.ascii");

   // Create a new histogramm which egs4.acsii values
   int x_min = 0;
   int x_max = 20;
   int nb_bins = 40;
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
   legend->AddEntry(h1d,"Geant4-9.1","l");
   legend->AddEntry(h1f,"EGS4","L");
   legend->Draw();


}
