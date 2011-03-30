{
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  TCanvas c("c","Plot Hadr02",-5, 5, 600, 400);
  leg = new TLegend(.75, .75, 1., 1.);
  //  TH1F* hh = c->DrawFrame(0, 0, 30, .19, "");
  //  hh->GetXaxis()->SetTitle();
  //  hh->GetYaxis()->SetTitle(tit[m]);
  //  hh->Draw("AXIS");
  TFile f1("dpmjet_S32_QGSP_BIC.root");
  h17->SetLineColor(2);
  h17->Draw("HIST");
  leg->AddEntry(h17, "DPMJET", "L");
  TFile f2("ftfp_S32_QGSP_BIC.root");
  h17->SetLineColor(3);
  h17->Draw("HIST SAME");
  leg->AddEntry(h17, "FTFP", "L");
  leg->Draw("SAME");
  c.Update();
  c.Print("Apic_A.gif");
}
