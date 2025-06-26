// *********************************************************************
// To execute this macro under ROOT after your simulation ended,
//   1 - launch ROOT (usually type 'root' at your machine's prompt)
//   2 - type '.X plot.C' at the ROOT session prompt
// *********************************************************************

void SetLeafAddress(TNtuple* ntuple, const char* name, void* address);

void plotElastic()
{
  gROOT->Reset();

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetGridColor(kGray+1);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetTitleOffset(1.2, "Y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLineWidth(2);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(42);
  gStyle->SetOptStat(000000);

  TCanvas* c1 = new TCanvas ("c1","",20,20,1500,500);
  c1->Divide(3,1);

  TFile* f = new TFile("dna.root");

  TNtuple* ntuple2;
  ntuple2 = (TNtuple*)f->Get("step");
  bool rowWise2 = true;
  TBranch* eventBranch2 = ntuple2->FindBranch("row_wise_branch");
  if ( ! eventBranch2 ) rowWise2 = false;

  // Cosine plot

  c1->cd(1);
  gPad->SetLogy();

  ntuple2->SetFillStyle(1001);
  ntuple2->SetFillColor(2);
  ntuple2->Draw
   ("(cosTheta)","parentID==0&&trackID==1&&stepID==1","");

  TH1* hist = (TH1*)gPad->GetPrimitive("htemp");
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->SetLabelSize(0.03);
  hist->GetXaxis()->SetTitleSize(0.03);
  hist->GetYaxis()->SetTitleSize(0.03);
  hist->GetXaxis()->SetTitleOffset(1.4);
  hist->GetYaxis()->SetTitleOffset(1.8);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->SetTitle("cos(#theta)");
  hist->GetYaxis()->SetTitle("");
  gPad->SetLogy();
  gPad->SetTicks(1, 1);
  gPad->Modified();

  // Solid angle (Omega) plot
  // (as in dsigma/dOmega vs Omega diff. cross section plots)

  c1->cd(2);
  gPad->SetLogy();
  ntuple2->Draw
   ("2*3.1415926535*(1-cosTheta)","parentID==0&&trackID==1&&stepID==1","");

  TH1* hist2 = (TH1*)gPad->GetPrimitive("htemp");
  hist2->GetXaxis()->SetLabelSize(0.03);
  hist2->GetYaxis()->SetLabelSize(0.03);
  hist2->GetXaxis()->SetTitleSize(0.03);
  hist2->GetYaxis()->SetTitleSize(0.03);
  hist2->GetXaxis()->SetTitleOffset(1.4);
  hist2->GetYaxis()->SetTitleOffset(1.8);
  hist2->GetXaxis()->CenterTitle();
  hist2->GetYaxis()->CenterTitle();
  hist2->GetXaxis()->SetTitle("Solid angle (sr)");
  hist2->GetYaxis()->SetTitle("");
  gPad->SetLogy();
  gPad->SetTicks(1, 1);
  gPad->Modified();

  // Angle plot

  c1->cd(3);
  gPad->SetLogy();
  ntuple2->Draw
   ("acos(cosTheta)*180/3.14159","parentID==0&&trackID==1&&stepID==1","");

  TH1* hist3 = (TH1*)gPad->GetPrimitive("htemp");
  hist3->GetXaxis()->SetLabelSize(0.03);
  hist3->GetYaxis()->SetLabelSize(0.03);
  hist3->GetXaxis()->SetTitleSize(0.03);
  hist3->GetYaxis()->SetTitleSize(0.03);
  hist3->GetXaxis()->SetTitleOffset(1.4);
  hist3->GetYaxis()->SetTitleOffset(1.8);
  hist3->GetXaxis()->CenterTitle();
  hist3->GetYaxis()->CenterTitle();
  hist3->GetYaxis()->SetTitle("");
  hist3->GetXaxis()->SetTitle("Angle (deg)");
  gPad->SetLogy();
  gPad->SetTicks(1, 1);
  gPad->Modified();
}

void SetLeafAddress(TNtuple* ntuple, const char* name, void* address) {
  TLeaf* leaf = ntuple->FindLeaf(name);
  if ( ! leaf ) {
    std::cerr << "Error in <SetLeafAddress>: unknown leaf --> " << name << std::endl;
    return;
  }
  leaf->SetAddress(address);
}
