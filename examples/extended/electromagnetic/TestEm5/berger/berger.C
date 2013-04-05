{
   gROOT->Reset();

   // Draw histograms fill by Geant4 TestEm5 simulation
   TFile f1("./msc93.root");
   TH1D* h1 = (TH1D*) f1.Get("1");
   h1->SetTitle("Energy deposited in 530 um of Si by 1 MeV e-");
   h1->GetXaxis()->SetTitle("Ekine (MeV)");
   h1->GetYaxis()->SetTitle("nb / MeV");
   h1->SetStats(kFALSE);  // Eliminate statistics box
   h1->SetLineColor(1);   // black
   h1->Draw("HIST");
   
   TFile f2("./local.root");
   TH1D* h2 = (TH1D*) f2.Get("1");
   h2->SetStats(kFALSE);  // Eliminate statistics box
   h2->SetLineColor(4);   // blue
   h2->Draw("SAME HIST");
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
   TH1F* h3 = new TH1F("h1f","",nbdata,x_min,x_max);

   Float_t x, y;
   while (1) {
      in >> x >> y ;
      if (!in.good()) break;
      h3->Fill(x,y);
   }
   in.close();

   // Draw histogram fill by data.acsii values
   h3->SetLineColor(2);   // red
   h3->Draw("SAME");

   // Print the histograms legend
   TLegend *legend = new TLegend(0.6,0.6,0.8,0.8);
   legend->AddEntry(h1,"Urban93","l");
   legend->AddEntry(h2,"Urban95","l");   
   legend->AddEntry(h3,"Berger data","L");
   legend->Draw();


}
