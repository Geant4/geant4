// -------------------------------------------------------------------
// -------------------------------------------------------------------

// This macro requires the pdb4dna_output.root file generated from PDB4DNA example

{
  gROOT->Reset();
 
  gStyle->SetOptStat("em");

  TCanvas *c1;
  TPad *pad1, *pad2, *pad3;
  c1 = new TCanvas("c1","PDB DNA outputs",200,10,700,780);
  c1->SetFillColor(0);

  pad1 = new TPad("pad1","pad1",0.02,0.52,0.98,0.98,21);
  pad2 = new TPad("pad2","pad2",0.02,0.02,0.48,0.48,21);
  pad3 = new TPad("pad3","pad3",0.52,0.02,0.98,0.48,21);

  pad1->SetFillColor(0);
  pad1->Draw();
  pad2->SetFillColor(0);
  pad2->Draw();
  pad3->SetFillColor(0);
  pad3->Draw();


  TFile f("pdb4dna_output.root");

  // Draw histograms

  TH1D* hist1 = (TH1D*)f.Get("1");
  pad1->cd();
  hist1->Draw("HIST");

  TH1D* hist2 = (TH1D*)f.Get("2");
  pad2->cd();
  hist2->Draw("HIST");

  TH1D* hist3 = (TH1D*)f.Get("3");
  pad3->cd();
  hist3->Draw("HIST");

  c1->Modified();
  c1->Update();


  // Read stats to get global quantities

  double* pdbStats=new double[4];

  hist1->GetStats(pdbStats);
  cout << "-> Edep in the target : " << pdbStats[2]/1E6 << " MeV" << endl;

  hist2->GetStats(pdbStats);
  cout << "-> Number of SSB : " << pdbStats[2] << endl;

  hist3->GetStats(pdbStats);
  cout << "-> Number of DSB : " << pdbStats[2] << endl;
}
