{
fin = dir[idir] + "/" + fil[iener] + ".root";

TFile* fff = new TFile(fin);
TH1F*  hhh;

cout <<"Opened file "<<fin<<"  iplot= "<<iplot<<" idir= "<<idir<<endl;
if(iplot==3) hhh = (TH1F*)fff->Get("h12");
if(iplot==4) hhh = (TH1F*)fff->Get("h10");
if(iplot==5) hhh = (TH1F*)fff->Get("h11");

hhh->SetLineColor(col[idir]);
//hhh->SetLineWidth(2);
leg[0]->AddEntry(hhh, dir[idir], "l");
//hhh->Draw("C HIST SAME");
hhh->Draw("HIST SAME 9");
cout << "AddMC1 done" << endl;
}
