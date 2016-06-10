{
   gROOT->Reset();

  /*
   * Mainegra et al. Med. Phys. 32, 685-699 (2005)    
   * Seltzer, Appl. Radiat. Isot. 42 (1991) page 917  
   * L. Ferrer et al. Cancer Bio. Rad. 22-1 (2007)
   */   

   ifstream in;
   in.open("./EGSnrc/100keV-DPK.ascii");
   //in.open("./Etran/100keV-DPK.ascii");
   
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
   h1d->SetTitle("Dose point kernel : energy deposition profile, e- 100 keV");
   h1d->GetXaxis()->SetTitle("d(E/E0)/d(r/r0) along r/r0                             r/r0");
   h1d->GetYaxis()->SetTitle("DPK");
   h1d->SetStats(kFALSE);  // Eliminate statistics box      
   h1d->SetLineColor(kRed);
   h1d->Draw("L");


   // Draw histograms fill by Geant4 TestEm12 simulation
   TFile f1("./100keV_local.root");
   TH1D* h1 = (TH1D*) f1.Get("8");
   h1->SetLineColor(kBlack);      
   h1->Draw("SAME HIST");
/*   
   TFile f2("./10keV.liver.root");
   TH1D* h2 = (TH1D*) f2.Get("8");
   h2->SetLineColor(kBlue);   
   h2->Draw("SAME HIST");
   
   TFile f3("./10keV.penel.root");
   TH1D* h3 = (TH1D*) f3.Get("8");
   h3->SetLineColor(kGreen);   
   h3->Draw("SAME HIST");
*/      
   // Print the histograms legend
   TLegend *legend = new TLegend(0.7,0.6,0.86,0.78);
   legend->AddEntry(h1,"local ","l");
   ///legend->AddEntry(h2,"liver ","l");
   ///legend->AddEntry(h3,"penel","l");      
   legend->AddEntry(h1d,"EGSnrc","l");
   legend->Draw();

}
