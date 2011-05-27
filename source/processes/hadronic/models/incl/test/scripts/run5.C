void run5()
{
  gROOT->ProcessLine(".x scripts/rootlogon.C");
  gROOT->SetStyle("clearRetro");

  TFile *cpp_f = new TFile("tmp/run5.root");
  TFile *fort_f = new TFile("tmp/run5ref.root");
  TTree *cpp = (TTree*) cpp_f->Get("h101");
  TTree *fort = (TTree*) fort_f->Get("h101");
  cpp->SetLineColor(kRed);

  TCanvas *c1 = new TCanvas();
  c1->SetLogy();
  cpp->Draw("Avv", "Avv");
  fort->Draw("Avv", "Avv", "same");
  cpp->GetHistogram()->SetTitle("n(425 MeV) + 63Cu");
  Double_t xsize = cpp->GetHistogram()->GetXaxis()->GetTitleSize();
  cpp->GetHistogram()->GetXaxis()->SetTitleSize(1.1*xsize);
  cpp->GetHistogram()->GetYaxis()->SetTitleSize(1.1*xsize);
  cpp->GetHistogram()->GetXaxis()->SetTitle("Mass number of the residual nuclei");
  cpp->GetHistogram()->GetYaxis()->SetTitle("Number of particles");

  cpp->Draw("Avv", "Avv > 30");
  fort->Draw("Avv", "Avv > 30", "same");
  cpp->GetHistogram()->SetTitle("n(425 MeV) + 63Cu");
  Double_t xsize = cpp->GetHistogram()->GetXaxis()->GetTitleSize();
  cpp->GetHistogram()->GetXaxis()->SetTitleSize(1.1*xsize);
  cpp->GetHistogram()->GetYaxis()->SetTitleSize(1.1*xsize);
  cpp->GetHistogram()->GetXaxis()->SetTitle("Mass number of the residual nuclei");
  cpp->GetHistogram()->GetYaxis()->SetTitle("Number of particles");

  // Legend
  Double_t legendUpperX = 0.12;
  Double_t legendUpperY = 0.85;
  Double_t legendLowerX = 0.38;
  Double_t legendLowerY = 0.73;
  TLegend *legend = new TLegend(legendUpperX, legendUpperY, legendLowerX, legendLowerY);
  legend->AddEntry(cpp, "C++", "l");
  legend->SetShadowColor(kWhite);
  legend->SetBorderSize(0);
  legend->AddEntry(fort, "FORTRAN", "l");
  legend->Draw();  
  c1->SaveAs("test5.ps(");
  c1->SaveAs("n425Cu.eps");

  cpp->Draw("Zvv");
  fort->Draw("Zvv", "", "same");
  c1->SaveAs("test5.ps");

  cpp->Draw("Massini");
  fort->Draw("Massini", "", "same");
  c1->SaveAs("test5.ps");

  cpp->Draw("Enerj", "Avv == 1 && Zvv == 0");
  fort->Draw("Enerj", "Avv == 1 && Zvv == 0", "same");
  c1->SaveAs("test5.ps");  
  
  c1->SetLogy(0);
  cpp->Draw("Exini");
  fort->Draw("Exini", "", "same");
  c1->SaveAs("test5.ps");

  c1->SetLogy(0);
  cpp->Draw("Enerj", "Avv == 4 && Zvv == 2");
  fort->Draw("Enerj", "Avv == 4 && Zvv == 2", "same");
  c1->SaveAs("test5.ps");

  cpp->Draw("Exini:Massini");
  fort->Draw("Exini:Massini", "", "same");
  c1->SaveAs("test5.ps)");

  INCL();
}

void INCL()
{
  TFile *cf = new TFile("tmp/run5.root");
  TFile *ff = new TFile("tmp/run5ref.root");

  TTree *ct = (TTree *) cf->Get("h101");
  ct->SetLineColor(kRed);
  TTree *ft = (TTree *) ff->Get("h101");

  TCanvas *inclResult = new TCanvas();
  inclResult->Divide(2,2);

  inclResult->cd(1);
  ft->Draw("Massini");
  ct->Draw("Massini", "", "same");

  inclResult->cd(2);
  ft->Draw("Mzini");
  ct->Draw("Mzini", "", "same");

  inclResult->cd(3);
  ft->Draw("Exini");
  ct->Draw("Exini", "", "same");

  inclResult->cd(4);
  ft->Draw("Jremn");
  ct->Draw("Jremn", "", "same");

  TCanvas *bimpactPlot = new TCanvas();
  ft->Draw("Bimpact");
  ct->Draw("Bimpact", "", "same");

  inclResult->SaveAs("run5Remnants.eps");
  bimpactPlot->SaveAs("run5ImpactParameter.eps");
}
