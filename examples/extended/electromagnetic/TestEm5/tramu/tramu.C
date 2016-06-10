{
   gROOT->Reset();

   // Draw histograms fill by Geant4 TestEm5 simulation
   TFile f1("./opt3.root");
   TH1D* h1 = (TH1D*) f1.Get("10");
   h1->SetTitle("1 TeV muon in 3 m iron : kinetic energy at exit (GeV)");
   h1->GetXaxis()->SetTitle("Ekine (GeV)");
   h1->GetYaxis()->SetTitle("nb/GeV");
   h1->SetStats(kFALSE);  // Eliminate statistics box
   h1->SetLineColor(kBlue);
   h1->Draw("HIST");
/*   
   TFile f2("./local.root");
   TH1D* h2 = (TH1D*) f2.Get("10");
   h2->SetStats(kFALSE);  // Eliminate statistics box
   h2->SetLineColor(kRed); 
   h2->Draw("SAME HIST");
*/     
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
   h1f->SetLineColor(kRed); 
   h1f->Draw("SAME");

   // Print the histograms legend
   TLegend* legend = new TLegend(0.2,0.55,0.45,0.70);
   legend->AddEntry(h1,"local (Urban90)","l");
   //legend->AddEntry(h2,"local (msc90)","l");   
   legend->AddEntry(h1f,"Mars14 simul ","L");
   legend->Draw();


}
