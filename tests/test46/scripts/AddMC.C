{
fin = dir[idir] + "/" + fil[iener] + ".root";

TFile* fff = new TFile(fin);
TH1F*  hhh;

cout <<"Opened file "<<fin<<"  iplot= "<<iplot<<" idir= "<<idir<<endl;
if(iplot==0) hhh = (TH1F*)fff->Get("h2");
if(iplot==1) hhh = (TH1F*)fff->Get("h7");
if(iplot==2) hhh = (TH1F*)fff->Get("h9");
if(iplot==3) hhh = (TH1F*)fff->Get("h11");

hhh->SetLineColor(col[idir]);
hhh->SetLineWidth(2);
leg[iplot]->AddEntry(hhh, hed[idir], "l");
//hhh->Draw("C HIST SAME");
hhh->Draw("HIST SAME");
//c1->Update();
}
