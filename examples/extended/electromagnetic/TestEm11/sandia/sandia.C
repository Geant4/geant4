
{
   gROOT->Reset();

   // Draw histograms fill by Geant4 TestEm11 simulation
   TFile f1("./Ta_1000keV_opt3.root");
   TH1D* h1 = (TH1D*) f1.Get("8");
   h1->SetTitle("Depth dose distribution of 1000 keV e- in Ta");
   h1->GetXaxis()->SetTitle("Edep (Mev.cm2/g) along x/r0                    x/r0");
   h1->GetYaxis()->SetTitle("MeV*cm2/g");
   h1->SetStats(kFALSE);  // Eliminate statistics box
   h1->SetLineColor(4);   // blue
   h1->Draw("HIST");
/*   
   TFile f2("./Ta.1000keV.opt2.root");
   TH1D* h2 = (TH1D*) f2.Get("8");
   h2->SetStats(kFALSE);  // Eliminate statistics box
   h2->SetLineColor(3);   // green
   h2->Draw("SAME HIST");
*/   
/* data
* G.J.Lockwood et al.
*     Sandia report SAND79-0414.UC-34a, February 1987
* O.Kadri et al. NIM B 258 (2007) 381
*/

   ifstream in;
   in.open("./data/Ta_1000keV.ascii");

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
   TLegend* legend = new TLegend(0.6,0.55,0.8,0.68);
   legend->AddEntry(h1,"ref10-opt3 ","l");
   ////legend->AddEntry(h2,"ref10-opt2 ","l");   
   legend->AddEntry(pt,"Sandia data","P");
   legend->Draw();
}
