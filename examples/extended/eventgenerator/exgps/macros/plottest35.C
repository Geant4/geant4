{

  TFile *input_file_1 = new TFile("test35a.root");
  TFile *input_file_2 = new TFile("test35b.root");

  TCanvas *c1 = new TCanvas("c1", "test35", 200, 10, 700, 500);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetLogx();
  c1->SetLogy();

  // histogram for energy spectra
  int n = 41;
  float bin[41];
  
  for (int i = 0; i < n; i++) {
    bin[i] =pow(10,(-2+0.1*i));
  }
  //
  TH1 *h_1 = new TH1D("unbiased","Source spectrum",40,bin);
  TH1 *h_2 = new TH1D("biased","Source spectrum",40,bin);
  input_file_1->cd();
  input_file_1->ls();
  // get the tuple t1
  double energy, weight;
  TTree *t1 = (TTree *) input_file_1->Get("MyTuple");
  t1->SetBranchAddress("Energy", &energy);
  t1->SetBranchAddress("Weight", &weight);
  cout <<t1->GetEntries() << endl;
  for (int i = 0; i < t1->GetEntries(); i++) {
    t1.GetEntry(i);
    //    cout << energy << " " << weight << endl;
    h_1->Fill(energy,weight);
  }
  input_file_2->cd();
  TTree *t2 = (TTree *) input_file_2->Get("MyTuple");
  t2->SetBranchAddress("Energy", &energy);
  t2->SetBranchAddress("Weight", &weight);
  cout <<t2->GetEntries() << endl;
  for (int i = 0; i < t1->GetEntries(); i++) {
    t2.GetEntry(i);
    h_2->Fill(energy,weight);
  }
  //  h_2->SetFillColor(kRed);
  h_2->SetLineStyle(kDashed);
  h_2->SetLineColor(kBlue);
  h_2->Draw();
  h_1->Draw("same") ;
  c1->Update();
  c1->Print("./test35.png");
  
  input_file_1->Close();
  input_file_2->Close();
}
