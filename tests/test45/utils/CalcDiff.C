{
TString file1 = fl[jj] + ".root";
cout << "file1= " << file1 << endl;
gROOT->ProcessLine(".x  $TEST45/utils/Calc0.C");
TFile ff(file1);
TH1F *h;
Int_t j1 = 11; 
Int_t j2 = 39; 
Double_t z1, z2, z3, z4;

for (i = 0; i < npl; i++) {
  c1.cd(i+1);  
  for (j=0; j<bin; j++) {    
    if (x11[j] < 5.0) {    
      if(i==0) h = (TH1F*)ff.Get("h1");
      if(i==1) h = (TH1F*)ff.Get("h2");
      if(i==2) h = (TH1F*)ff.Get("h3");
      if(i==3) h = (TH1F*)ff.Get("h4");
      if(i==4) h = (TH1F*)ff.Get("h5");
      if(i==5) h = (TH1F*)ff.Get("h6");

      z1 = h->GetBinContent(j+j1);
      z2 = h->GetBinError(j+j1);
    } else {
      if(i==0) h = (TH1F*)ff.Get("h7");
      if(i==1) h = (TH1F*)ff.Get("h8");
      if(i==2) h = (TH1F*)ff.Get("h9");
      if(i==3) h = (TH1F*)ff.Get("h10");
      if(i==4) h = (TH1F*)ff.Get("h11");
      if(i==5) h = (TH1F*)ff.Get("h12");

      z1 = h->GetBinContent(j-j2);
      z2 = h->GetBinError(j-j2);
    }

    if(i==0) {z3 = y11[j]; z4 = erry11[j];}
    if(i==1) {z3 = y22[j]; z4 = erry22[j];}
    if(i==2) {z3 = y33[j]; z4 = erry33[j];}
    if(i==3) {z3 = y44[j]; z4 = erry44[j];}
    if(i==4) {z3 = y55[j]; z4 = erry55[j];}
    if(i==5) {z3 = y66[j]; z4 = erry66[j];}
   
    if(z3 == 0.0) {yy[j] = 0.0; erry[j] = 0.0;}
    else if(z1 == 0.0) {yy[j] = 0.0; erry[j] = z2/z3;}
    else {
      yy[j] = z1/z3;
      z4    = z4/z3;
      z3    = z2/z1;
      erry[j] = yy[j]*sqrt(z3*z3 + z4*z4);
    }
  cout << "i= "<< i  << " j= " << j << " x= " << x11[j] << " y= " << yy[j] << " dx= " << errx11[j] << " dy= " << erry[j] <<endl;    
  }
  
  gr[i] = new TGraphErrors(bin,x11,yy,errx11,erry);
  gr[i]->SetMarkerColor(coll[jj]);
  gr[i]->Draw("P SAME");  
  if (i == 5)  {leg[0]->AddEntry(gr[5], nam1[jj], "p"); }
  cout  << "## i= " << i << "  j= " << j << endl;
}
}
