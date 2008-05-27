// all files, all energies Plot.C

{
TString file1 = "test45.root";
cout << "file1= " << file1 << endl;
TFile ff(file1);

c1.cd(1);
h7->SetLineColor(2);
h7->Draw("HISTO SAME");
h1->SetLineColor(2);
h1->Draw("HISTO SAME");
leg[0]->AddEntry(h7, "QGSP_BERT", "l"); 

if(npl>1) {
  c1.cd(2);
  h8->SetLineColor(2);
  h8->Draw("HISTO SAME");
  h2->SetLineColor(2);
  h2->Draw("HISTO SAME");
}
if(npl>2) {
  c1.cd(3);
  h9->SetLineColor(2);
  h9->Draw("HISTO SAME");
  h3->SetLineColor(2);
  h3->Draw("HISTO SAME");
}
if(npl>3) {
  c1.cd(4);
  h10->SetLineColor(2);
  h10->Draw("HISTO SAME");
  h4->SetLineColor(2);
  h4->Draw("HISTO SAME");
}
if(npl>4) {
  c1.cd(5);
  h11->SetLineColor(2);
  h11->Draw("HISTO SAME");
  h5->SetLineColor(2);
  h5->Draw("HISTO SAME");
}
if(npl>5) {
  c1.cd(6);
  h12->SetLineColor(2);
  h12->Draw("HISTO SAME");
  h6->SetLineColor(2);
  h6->Draw("HISTO SAME");
}
}
