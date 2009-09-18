{
   gROOT->Reset();

   // Draw histogram fill by Geant4 TestEm5 simulation
   TFile f("./Si530um.root");
   TH1D* h1d = (TH1D*) f.Get("1");
   h1d->SetTitle("Energy deposited in 530 um of Si by 1 MeV e-");
   h1d->GetXaxis()->SetTitle("Ekine (MeV)");
   h1d->GetYaxis()->SetTitle("nb / MeV");
   h1d->SetStats(kFALSE);  // Eliminate statistics box
   h1d->Draw("HIST");

/*
* e- 1 MeV in Silicon 530 um 
* M.J. Berger et al. NIM 69 (p.181) 1969
* distribution of energy deposition
* (from 110 keV to 1.03 MeV by bin of 10 keV --> 93 bins) 
*/

   ifstream in;
   in.open("530um.ascii");
   
   // First indicate number of data
   int nbdata = 0;
   in >> nbdata;
   
   // Create a new histogram with data.acsii values
   float x_min = 0.110;
   float x_max = 1.030;
   TH1F* h1f = new TH1F("h1f","",nbdata,x_min,x_max);

   Float_t x, y;
   while (1) {
      in >> x >> y ;
      if (!in.good()) break;
      h1f->Fill(x,y);
   }
   in.close();

   // Draw histogram fill by data.acsii values
   h1f->SetLineColor(2);
   h1f->Draw("SAME");

   // Print the histograms legend
   TLegend *legend = new TLegend(0.6,0.6,0.8,0.8);
   legend->AddEntry(h1d,"G4-9.2-ref-04 ","l");
   legend->AddEntry(h1f,"Berger data","L");
   legend->Draw();


}
