
{
   gROOT->Reset();

   // Draw histograms fill by Geant4 TestEm11 simulation
   //
   ///TFile f1("./Al_1033keV_opt3.root");      
   TFile f1("./Ta_1000keV_opt3.root");
   TH1D* h1 = (TH1D*) f1.Get("8");
   h1->SetTitle("Depth dose distribution of 1000 keV e- in Ta");
   h1->GetXaxis()->SetTitle("Edep (Mev.cm2/g) along x/r0                    x/r0");
   h1->GetYaxis()->SetTitle("MeV*cm2/g");
   h1->SetStats(kFALSE);  // Eliminate statistics box
   h1->SetLineColor(kBlack);
   h1->Draw("HIST");
      
/* EGSnrc
* Yann Perrot
*/

   ifstream in;
   ///in.open("./EGSnrc/Al_1033keV_EGSnrc.ascii");
   in.open("./EGSnrc/Ta_1000keV_EGSnrc.ascii");
      
   TMarker *pt;
   Double_t x, y;
   // First indicate number of data
   int nbdata = 0;
   in >> nbdata;
   for ( int i = 0 ; i < nbdata ; i++ ) {
      in >> x >> y ;
      if (!in.good()) break;
      pt = new TMarker(x,y,22); // 22 for triangle TMatker
      pt->SetMarkerColor(kRed);
      pt->Draw();
   }
   in.close();

   // Print the histograms legend
   TLegend* legend = new TLegend(0.6,0.5,0.8,0.70);
   legend->AddEntry(h1,"ref10 ","l");
   legend->AddEntry(pt,"EGSnrc","P");
   legend->Draw();
}
