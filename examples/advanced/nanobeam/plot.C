{
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  gStyle->SetOptStat(1111);
  gStyle->SetPalette(1);

  c1 = new TCanvas ("c1","",200,10,800,800);
  c1.Divide(2,2);

  Float_t x,y,z,theta,phi,nrj;
  Int_t ncols;

  FILE *fp1 = fopen("./results/profile.txt","r");
  if (fp1==0) goto jump1;
  
  TNtuple *ntuple1 = new TNtuple("ntuple1","file1","x:y:z");
  while (1) 
  {
    ncols = fscanf(fp1,"%f %f %f",&x, &y, &z);
      if (ncols < 0) break;    
      ntuple1->Fill(x/1000,y/1000,z/1000);
  }
  fclose(fp1);

  c1.cd(1);
  ntuple1->SetMarkerSize(.2);
  ntuple1->SetMarkerColor(2);
  ntuple1->SetMarkerStyle(20);
  ntuple1->Draw("x:y:z","");
  htemp->GetXaxis()->SetLabelSize(0.025);
  htemp->GetYaxis()->SetLabelSize(0.025);
  htemp->GetZaxis()->SetLabelSize(0.025);
  htemp->GetXaxis()->SetTitleSize(0.035);
  htemp->GetYaxis()->SetTitleSize(0.035);
  htemp->GetZaxis()->SetTitleSize(0.035);
  htemp->GetXaxis()->SetTitleOffset(1.4);
  htemp->GetYaxis()->SetTitleOffset(1.4);
  htemp->GetZaxis()->SetTitleOffset(1.4);
  htemp->GetXaxis()->SetTitle("Z (m)");
  htemp->GetYaxis()->SetTitle("X (mm)");
  htemp->GetZaxis()->SetTitle("Y (mm)");
  htemp->SetTitle("Nanobeam profile");
 
  jump1:
  
  FILE *fp1 = fopen("./results/x.txt","r");
  if (fp1==0) goto jump2;
  FILE *fp2 = fopen("./results/y.txt","r");
  FILE *fp3 = fopen("./results/theta.txt","r");
  FILE *fp4 = fopen("./results/phi.txt","r");
  TNtuple *ntuple2 = new TNtuple("ntuple2","file2","x:y:theta:phi");
  while (1) 
  {
    ncols = fscanf(fp1,"%f",&x);
      if (ncols < 0) break;    
    ncols = fscanf(fp2,"%f",&y);
    ncols = fscanf(fp3,"%f",&theta);
    ncols = fscanf(fp4,"%f",&phi);
    ntuple2->Fill(x,y,theta,phi);
  }
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);

  c1.cd(2);
  gStyle->SetPalette(1);
  ntuple2->SetMarkerColor(2);
  ntuple2->SetMarkerStyle(20);
  gPad->SetLogz();
  ntuple2->Draw("y:x","","hcolz");
  htemp->GetXaxis()->SetLabelSize(0.025);
  htemp->GetYaxis()->SetLabelSize(0.025);
  htemp->GetXaxis()->SetTitleSize(0.035);
  htemp->GetYaxis()->SetTitleSize(0.035);
  htemp->GetXaxis()->SetTitleOffset(1.4);
  htemp->GetYaxis()->SetTitleOffset(1.4);
  htemp->GetXaxis()->SetTitle("X (micrometer)");
  htemp->GetYaxis()->SetTitle("Y (micrometer)");
  htemp->SetTitle("Nanobeam image");
 
  c1.cd(3);
  gStyle->SetPalette(1);
  ntuple2->SetMarkerColor(4);
  ntuple2->SetMarkerStyle(20);
  gPad->SetLogz();
  ntuple2->Draw("theta:x","","hcolz");
  htemp->GetXaxis()->SetLabelSize(0.025);
  htemp->GetYaxis()->SetLabelSize(0.025);
  htemp->GetXaxis()->SetTitleSize(0.035);
  htemp->GetYaxis()->SetTitleSize(0.035);
  htemp->GetXaxis()->SetTitleOffset(1.4);
  htemp->GetYaxis()->SetTitleOffset(1.4);
  htemp->GetXaxis()->SetTitle("X (micrometer)");
  htemp->GetYaxis()->SetTitle("THETA (mrad)");
  htemp->SetTitle("Nanobeam image emittance (X,THETA)");

  c1.cd(4);
  gStyle->SetPalette(1);
  ntuple2->SetMarkerColor(4);
  ntuple2->SetMarkerStyle(20);
  gPad->SetLogz();
  ntuple2->Draw("phi:y","","hcolz");
  htemp->GetXaxis()->SetLabelSize(0.025);
  htemp->GetYaxis()->SetLabelSize(0.025);
  htemp->GetXaxis()->SetTitleSize(0.035);
  htemp->GetYaxis()->SetTitleSize(0.035);
  htemp->GetXaxis()->SetTitleOffset(1.4);
  htemp->GetYaxis()->SetTitleOffset(1.4);
  htemp->GetXaxis()->SetTitle("Y (micrometer)");
  htemp->GetYaxis()->SetTitle("PHI (mrad)");
  htemp->SetTitle("Nanobeam image emittance (Y,PHI)");
   
jump2:
  FILE *fp3 = fopen("./results/grid.txt","r");
  if (fp3==0) goto jump3;
  
  TNtuple *ntuple3 = new TNtuple("ntuple3","file3","x:y:nrj");
  while (1) 
  {
    ncols = fscanf(fp3,"%f %f %f",&x, &y, &nrj);
      if (ncols < 0) break;    
      ntuple3->Fill(x/1000,y/1000,nrj);
  }
  fclose(fp3);

  c1.cd(1);
  gStyle->SetPalette(1);
  ntuple3->SetMarkerSize(.2);
  ntuple3->SetMarkerColor(2);
  ntuple3->SetMarkerStyle(20);
  gPad->SetLogz();
  ntuple3->Draw("y:x","","hcolz");
  htemp->GetXaxis()->SetLabelSize(0.025);
  htemp->GetYaxis()->SetLabelSize(0.025);
  htemp->GetXaxis()->SetTitleSize(0.035);
  htemp->GetYaxis()->SetTitleSize(0.035);
  htemp->GetXaxis()->SetTitleOffset(1.4);
  htemp->GetYaxis()->SetTitleOffset(1.4);
  htemp->GetXaxis()->SetTitle("X (mm)");
  htemp->GetYaxis()->SetTitle("Y (mm)");
  htemp->SetTitle("Grid shadow");
 
jump3:
}
