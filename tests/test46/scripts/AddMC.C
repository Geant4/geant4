{
fin = dir[idir] + "/" + fil[iener] + ".root";

TFile* fff = new TFile(fin);
TH1F*  hhh;

cout <<"Opened file "<<fin<<"  iplot= "<<iplot<<" idir= "<<idir<<endl;
if(iplot+iplot0==0) hhh = (TH1F*)fff->Get("h2");
if(iplot+iplot0==1) hhh = (TH1F*)fff->Get("h7");
if(iplot+iplot0==2) hhh = (TH1F*)fff->Get("h9");
if(iplot+iplot0==3) hhh = (TH1F*)fff->Get("h12");
if(iplot+iplot0==6) hhh = (TH1F*)fff->Get("h4");
if(iplot+iplot0==7) hhh = (TH1F*)fff->Get("h5");

hhh->SetLineColor(col[idir]);
//hhh->SetLineWidth(2);
leg[iplot]->AddEntry(hhh, dir1[idir], "l");
//leg[iplot]->AddEntry(hhh, dir[idir], "l");
//hhh->Draw("C HIST SAME");
hhh->Draw("HIST SAME");
//c1->Update();
//delete fff;
 cout << "AddMC done iplot= " << iplot << " iplot0= " << iplot0 << endl;
}
