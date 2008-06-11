// all files, all energies Plot.C

{
TString file1 = fl[j] + ".root";
cout << "file1= " << file1 << endl;
if(j == 0) TFile ff1(file1);
if(j == 1) TFile ff2(file1);

c1.cd(1);
h1->SetLineColor(coll[j]);
h1->Draw("HISTO SAME");
leg[0]->AddEntry(h1, nam[j], "l"); 
cout << "legend color is " << coll[j] << endl;

if(npl>1) {
  c1.cd(2);
  h2->SetLineColor(coll[j]);
  h2->Draw("HISTO SAME");
cout << "h2 is OK"  << endl;
}
if(npl>2) {
  c1.cd(3);
  h3->SetLineColor(coll[j]);
  h3->Draw("HISTO SAME");
cout << "h3 is OK"  << endl;
}
if(npl>3) {
  c1.cd(4);
  h4->SetLineColor(coll[j]);
  h4->Draw("HISTO SAME");
cout << "h4 is OK"  << endl;
}
if(npl>4) {
  c1.cd(5);
  h5->SetLineColor(coll[j]);
  h5->Draw("HISTO SAME");
}
if(npl>5) {
  c1.cd(6);
  h6->SetLineColor(coll[j]);
  h6->Draw("HISTO SAME");
}
}
