{
fin = filp[ipart] + filt[itarg] + ".root";

TFile* fff = new TFile(fin);
TH1F*  hhh;

cout << "Opened file " << fin << "  ixs= " << ixs << endl;
if(ixs==0) hhh = (TH1F*)fff->Get("h1");
if(ixs==1) hhh = (TH1F*)fff->Get("h2");
if(ixs==2) hhh = (TH1F*)fff->Get("h3");
if(ixs==3) hhh = (TH1F*)fff->Get("h4");
if(ixs==4) hhh = (TH1F*)fff->Get("h5");
if(ixs==5) hhh = (TH1F*)fff->Get("h6");
if(ixs==6) hhh = (TH1F*)fff->Get("h7");

hhh->SetLineColor(col[itarg]);
hhh->SetLineWidth(2);
leg[ixs]->AddEntry(hhh, targ[itarg], "l");
hhh->Draw("C HIST SAME");
c1->Update();
}
