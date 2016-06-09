{
   gROOT->Reset();

  /* data 
   * Seltzer, Appl. Radiat. Isot. 42(1991) page 917  
   * L. Ferrer et al. Cancer Bio. Rad. 22-1 (2007)
   */   

   ifstream in;
   in.open("./data/1MeV-DPK.ascii");

   // Create a new histogram
   double x_min = 0;
   double x_max = 1.2;
   
   int nbdata = 0;
   in >> nbdata;
   
   TH1D* h1d = new TH1D("h1d","",nbdata,x_min,x_max);

   Double_t x, y;   
   for ( int i = 0 ; i < nbdata ; i++ ) {
      in >> x >> y ;
      if (!in.good()) break;
      h1d->Fill(x,y);
   }
   in.close();
   
   // Draw histogram
   h1d->SetTitle("Dose point kernel : energy deposition profile, e- 1 MeV");
   h1d->GetXaxis()->SetTitle("d(E/E0)/d(r/r0) along r/r0                             r/r0");
   h1d->GetYaxis()->SetTitle("DPK");
   h1d->SetStats(kFALSE);  // Eliminate statistics box      
   h1d->SetLineColor(2);
   h1d->Draw("L");


   // Draw histogram fill by Geant4 TestEm12 simulation
   TFile f("./1MeV.root");
   TH1D* hroot = (TH1D*) f.Get("8");
   hroot->Draw("SAME HIST");

   // Print the histograms legend
   TLegend *legend = new TLegend(0.7,0.65,0.9,0.78);
   legend->AddEntry(hroot,"G4-9.1-ref-08++ ","l");
   legend->AddEntry(h1d,"Etran simul","l");
   legend->Draw();

}
