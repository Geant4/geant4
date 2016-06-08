{
  TFile f("Histo.root");
  f->ls();

  hist1->SetFillColor(kRed);
  hist1->SetMarkerStyle(21);
  hist1->SetMarkerColor(kRed);

  hist2->SetFillColor(kBlue);
  hist2->SetMarkerStyle(21);
  hist2->SetMarkerColor(kBlue);

  hist3->SetFillColor(50);
  hist3->SetMarkerStyle(21);
  hist3->SetMarkerColor(50);

  hist4->SetFillColor(kGreen);
  hist4->SetMarkerStyle(21);
  hist4->SetMarkerColor(kGreen);

  // hist5->SetFillColor(kYellow);
  // hist5->SetMarkerStyle(21);
  // hist5->SetMarkerColor(kYellow);

  // hist6->SetFillColor(7);
  // hist6->SetMarkerStyle(21);
  // hist6->SetMarkerColor(7);

  // TCanvas *c1 = new TCanvas("c1","stacked hists",10,10,750,600);
  TCanvas *c1 = new TCanvas("c1","stacked hists",10,10,1000,800);
  c1->SetFillColor(41);
  c1->Divide(2,2);

  // in top left pad, draw the stack with defaults
  c1->cd(1);
  hist1->Draw();

  // in top right pad, draw the stack in non-stack mode and errors option
  c1->cd(2);
  hist2->Draw();

  c1->cd(3);
  hist3->Draw();

  c1->cd(4);
  hist4->Draw();

  // c1->cd(5);
  // hist5->Draw();
}
