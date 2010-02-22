{
   gROOT->Reset();

   // Draw histogram filled by Geant4 TestEm2 simulation
   TFile f("./93ref0.root");
   TH1D* h1d = (TH1D*) f.Get("4");
   h1d->SetTitle("30 GeV e- on 20 X0 Fe : energy dep, longit profil");
   h1d->GetXaxis()->SetTitle("depth (X0)");
   h1d->GetYaxis()->SetTitle("(100/E0) (dE/dt)");
   h1d->SetStats(kFALSE);  // Eliminate statistics box
   h1d->Draw("HIST");

   // pdg.ascii came from egs4 simulation
   ifstream in;
   in.open("pdg.ascii");

   // Create a new histogramm which egs4.acsii values
   int nb_bins = 40;   
   float x_min = 0;
   float x_max = 20;
   TH1F* h1f = new TH1F("h1f","",nb_bins,x_min,x_max);

   Float_t x, y;
   while (1) {
      in >> x >> y ;
      if (!in.good()) break;
      h1f->Fill(x,y);
   }
   in.close();

   // Draw histogram fill by egs4.acsii values
   h1f->SetLineColor(2);
   h1f->Draw("SAME");

   // Print the histograms legend
   TLegend *legend = new TLegend(0.65,0.55,0.85,0.68);
   legend->AddEntry(h1d,"93ref0","l");
   legend->AddEntry(h1f,"EGS4","L");
   legend->Draw();
}
