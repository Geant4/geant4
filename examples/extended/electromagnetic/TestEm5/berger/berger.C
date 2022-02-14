{
   gROOT->Reset();

   // Draw histograms fill by Geant4 TestEm5 simulation
   TFile f1("./local.root");
   TH1D* h1 = (TH1D*) f1.Get("h1");
   h1->SetTitle("Energy deposited in 530 um of Si by 1 MeV e-");
   h1->GetXaxis()->SetTitle("Ekine (MeV)");
   h1->GetYaxis()->SetTitle("nb / MeV");
   h1->SetStats(kFALSE);  // Eliminate statistics box
   h1->SetLineColor(kBlue);
   h1->Draw("HIST");
   
   TFile f2("./opt4.root");
   TH1D* h2 = (TH1D*) f2.Get("h1");
   h2->SetStats(kFALSE);  // Eliminate statistics box
   h2->SetLineColor(kRed);
   h2->Draw("SAME HIST");
/*
* e- 1 MeV in Silicon 530 um 
* M.J. Berger et al. NIM 69 (p.181) 1969
* distribution of energy deposition
* (from 110 keV to 1.03 MeV by bin of 10 keV --> 93 bins) 
*/

  ifstream in;
  in.open("530um.ascii");
   
  TMarker *pt;
  Double_t x, y;   
  // First indicate number of data
  int nbdata = 0;
  in >> nbdata;
  for ( int i = 0 ; i < nbdata ; i++ ) {
      in >> x >> y ;
      if (!in.good()) break;
      pt = new TMarker(x,y,32); // 32 for open triangle-down
      pt->SetMarkerColor(kGreen);
      pt->Draw();
  }   
  in.close();

  // Print the histograms legend
  TLegend *legend = new TLegend(0.6,0.6,0.8,0.8);
  legend->AddEntry(h1,"local","l");
  legend->AddEntry(h2,"opt4 ","l");   
  legend->AddEntry(pt,"Berger data","P");
  legend->Draw();
}
