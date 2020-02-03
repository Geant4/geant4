{
// Read reference data in granero.txt
FILE *fg1=fopen("granero.txt", "r");
Int_t n_points_granero =13;
Float_t x1[n_points_granero], y1[n_points_granero], ratio_liv[n_points_granero];
Float_t x, y;
Int_t ncols_granero;
Int_t nlines1 =0;

while(1)
 {
  ncols_granero = fscanf(fg1,"%f %f",&x, &y);
  if (ncols_granero<0) break;
 // std::cout << "x " << x << std::endl;
  x1[nlines1]=x;
  y1[nlines1]=y;
  nlines1++;
}

fclose(fg1);

// Read the results of the brachytherapy advanced example
// FlexiSorceMacro.mac - livermore
FILE *fg2=fopen("geant4_dose_Flexi_livermore.txt", "r");
Int_t n_points_geant4 =398;
Float_t x2[n_points_geant4], y2[n_points_geant4];
Int_t ncols_geant4_liv;
Int_t nlines2 =0;

while(1)
 {
  ncols_geant4_liv = fscanf(fg2,"%f %f",&x, &y);
  if (ncols_geant4_liv<0) break;
 // std::cout << "x " << x << std::endl;
  x2[nlines2]=x;
  y2[nlines2]=y;
  
  for (int i=0; i<n_points_granero; i++)
   {
   if (x1[i]==x2[nlines2])
   { 
   ratio_liv[i]= y2[nlines2]/y1[i];
  // std::cout << "granero: " <<  x1[i] << "," << y1[i] 
   //          << ", livermore: "<< x2[nlines2] << "," << y2[nlines2]
   //          << " ratio:" << ratio_liv[i] << std::endl;
  }
  }
  nlines2++;
}

fclose(fg2);

// Read the results of the brachytherapy advanced example
//  penelope
FILE *fg3=fopen("geant4_dose_Flexi_penelope.txt", "r");
Float_t x3[n_points_geant4], y3[n_points_geant4], ratio_pen[n_points_granero];
Int_t ncols_geant4_penelope;
Int_t nlines3 =0;

while(1)
 {
  ncols_geant4_penelope = fscanf(fg3,"%f %f",&x, &y);
  if (ncols_geant4_penelope<0) break;
 // std::cout << "x " << x << std::endl;
  x3[nlines3]=x;
  y3[nlines3]=y;
  for (int i=0; i<n_points_granero; i++)
   {
   if (x1[i]==x3[nlines3])
   { 
   ratio_pen[i]= y3[nlines3]/y1[i];
 //  std::cout << "granero: " <<  x1[i] << "," << y1[i] 
  //           << ", penelope: "<< x3[nlines3] << "," << y3[nlines3]
   //          << " ratio:" << ratio_pen[i] << std::endl;
  }
  }
  nlines3++;
}

fclose(fg3);

// FlexiSorceMacro.mac - opt 0
FILE *fg4=fopen("geant4_dose_Flexi_opt0.txt", "r");
Float_t x4[n_points_geant4], y4[n_points_geant4], ratio_opt0[n_points_granero];
Int_t ncols_geant4_opt0;
Int_t nlines4 =0;

while(1)
 {
  ncols_geant4_opt0 = fscanf(fg4,"%f %f",&x, &y);
  if (ncols_geant4_opt0<0) break;
 // std::cout << "x " << x << std::endl;
  x4[nlines4]=x;
  y4[nlines4]=y;
  for (int i=0; i<n_points_granero; i++)
   {
   if (x1[i]==x4[nlines4])
   { 
   ratio_opt0[i]= y4[nlines4]/y1[i];
  // std::cout << "granero: " <<  x1[i] << "," << y1[i] 
   //          << ", opt0: "<< x4[nlines4] << "," << y4[nlines4]
    //         << " ratio:" << ratio_opt0[i] << std::endl;
  }
  }
  nlines4++;
}

fclose(fg4);

// FlexiSorceMacro.mac - opt 3
FILE *fg5=fopen("geant4_dose_Flexi_opt3.txt", "r");
Float_t x5[n_points_geant4], y5[n_points_geant4], ratio_opt3[n_points_granero];
Int_t ncols_geant4_opt3;
Int_t nlines5 =0;

while(1)
 {
  ncols_geant4_opt3 = fscanf(fg5,"%f %f",&x, &y);
  if (ncols_geant4_opt3<0) break;
 // std::cout << "x " << x << std::endl;
  x5[nlines5]=x;
  y5[nlines5]=y;

 for (int i=0; i<n_points_granero; i++)
   {
   if (x1[i]==x5[nlines5])
   { 
   ratio_opt3[i]= y5[nlines5]/y1[i];
 /*  std::cout << "granero: " <<  x1[i] << "," << y1[i] 
             << ", opt3: "<< x5[nlines5] << "," << y5[nlines5]
            << " ratio:" << ratio_opt3[i] << std::endl;*/
  }
  }
  nlines5++;
}

fclose(fg5);

// FlexiSorceMacro.mac - opt 4
FILE *fg6=fopen("geant4_dose_Flexi_opt4.txt", "r");
Float_t x6[n_points_geant4], y6[n_points_geant4], ratio_opt4[n_points_granero];
Int_t ncols_geant4_opt4;
Int_t nlines6 =0;

while(1)
 {
  ncols_geant4_opt4 = fscanf(fg6,"%f %f",&x, &y);
  if (ncols_geant4_opt4<0) break;
 // std::cout << "x " << x << std::endl;
  x6[nlines6]=x;
  y6[nlines6]=y;
  for (int i=0; i<n_points_granero; i++)
   {
   if (x1[i]==x6[nlines6])
   { 
   ratio_opt4[i]= y6[nlines6]/y1[i];
/*   std::cout << "granero: " <<  x1[i] << "," << y1[i] 
              << ", opt4: "<< x6[nlines6] << "," << y6[nlines6]
            << " ratio:" << ratio_opt4[i] << std::endl;*/
  }
  }
  nlines6++;
}

fclose(fg6);

TGraph *gr1 = new TGraph (nlines1, x1, y1);
TGraph *gr2 = new TGraph (nlines2, x2, y2);
TGraph *gr3 = new TGraph (nlines3, x3, y3);
TGraph *gr4 = new TGraph (nlines4, x4, y4);
TGraph *gr5 = new TGraph (nlines5, x5, y5);
TGraph *gr6 = new TGraph (nlines6, x6, y6);

std::cout<< "Livermore" << std::endl;

for (Int_t j=0; j < nlines1; j++)
{
 if ((ratio_liv[j] > 1.03) || (ratio_liv[j] < 0.97)) std::cout<< "Difference above 3% :"<< x1[j] << ", " << ratio_liv[j] << std::endl;
} 
std::cout<< "penelope" << std::endl;
for (Int_t j=0; j < nlines1; j++)
 {
  if ((ratio_pen[j] > 1.03) || (ratio_pen[j] < 0.97)) std::cout<< "Difference above 3% :" << x1[j] << ", " << ratio_pen[j] << std::endl;
} 

std::cout<< "opt0" << std::endl;
for (Int_t j=0; j < nlines1; j++)
 {
 if ((ratio_opt0[j] > 1.03) || (ratio_opt0[j] < 0.97)) std::cout<< "Difference above 3% :" << x1[j] << ", " << ratio_opt0[j] << std::endl;
} 

std::cout<< "opt3" << std::endl;
for (Int_t j=0; j < nlines1; j++)
 {
 if ((ratio_opt3[j] > 1.03) || (ratio_opt3[j] < 0.97)) std::cout<< "Difference above 3% :" << x1[j] << ", " << ratio_opt3[j] << std::endl;
} 
std::cout<< "opt4" << std::endl;
for (Int_t j=0; j < nlines1; j++)
 {
 if ((ratio_opt4[j] > 1.03) || (ratio_opt4[j] < 0.97)) std::cout<< "Difference above 3% :" << x1[j] << ", " << ratio_opt4[j] << std::endl;
} 

TGraph *gr1_ratio = new TGraph (nlines1, x1, ratio_liv);
TGraph *gr2_ratio = new TGraph (nlines1, x1, ratio_pen);
TGraph *gr3_ratio = new TGraph (nlines1, x1, ratio_opt0);
TGraph *gr4_ratio = new TGraph (nlines1, x1, ratio_opt3);
TGraph *gr5_ratio = new TGraph (nlines1, x1, ratio_opt4);

// draw the graph with axis, continuous line, and put
// a * at each point
gr1->SetTitle("Dose rate distribution");
gr1-> GetXaxis()->SetTitle("Distance from the centre (cm)");
gr1->GetYaxis()->SetTitle("Normalised dose rate distribution");
gr1->SetLineWidth(1);
gr1->SetMarkerColor(1);
gr1->SetMarkerStyle(20);
gr1->Draw("AP");

gr2->SetLineWidth(0.3);
gr2->SetMarkerColor(2);
gr2->SetMarkerStyle(21);
gr2->SetMarkerSize(0.2);
gr2->SetLineColor(2);
gr2->Draw("CP");

gr3->SetLineWidth(0.3);
gr3->SetMarkerColor(3);
gr3->SetMarkerStyle(21);
gr3->SetMarkerSize(0.2);
gr3->SetLineColor(3);
gr3->Draw("CP");

gr4->SetLineWidth(0.3);
gr4->SetMarkerColor(4);
gr4->SetMarkerStyle(21);
gr4->SetMarkerSize(0.2);
gr4->SetLineColor(4);
gr4->Draw("CP");

gr5->SetLineWidth(0.3);
gr5->SetMarkerColor(6);
gr5->SetMarkerStyle(21);
gr5->SetMarkerSize(0.2);
gr5->SetLineColor(6);
gr5->Draw("CP");

gr6->SetLineWidth(0.3);
gr6->SetMarkerColor(9);
gr6->SetMarkerStyle(21);
gr6->SetMarkerSize(0.2);
gr6->SetLineColor(9);
gr6->Draw("CP");

TLegend *leg = new TLegend(0.3, 0.5, 0.6, 0.8);
leg->SetFillColor(0);
leg->AddEntry(gr1, "Reference data", "lp");
leg->AddEntry(gr2, "Geant4 - Livermore", "lp");
leg->AddEntry(gr3, "Geant4 - Penelope", "lp");
leg->AddEntry(gr4, "Geant4 - opt0", "lp");
leg->AddEntry(gr5, "Geant4 - opt3", "lp");
leg->AddEntry(gr6, "Geant4 - opt4", "lp");
leg->Draw();

c1 -> Print("Flexi_dose_rate_distribution.jpg");

// ratio plot
gr1_ratio->SetTitle("Dose rate distribution - RATIO");
gr1_ratio-> GetXaxis()->SetTitle("Distance from the centre (cm)");
gr1_ratio->GetYaxis()->SetTitle("Ratio");
gr1_ratio->SetLineWidth(1);
gr2_ratio->SetMarkerSize(0.8);
gr1_ratio->SetMarkerColor(1);
gr1_ratio->SetMarkerStyle(20);
gr1_ratio->Draw("AP");

gr2_ratio->SetLineWidth(0.3);
gr2_ratio->SetMarkerColor(2);
gr2_ratio->SetMarkerStyle(20);
gr2_ratio->SetMarkerSize(0.8);
gr2_ratio->SetLineColor(2);
gr2_ratio->Draw("P");

gr3_ratio->SetLineWidth(0.3);
gr3_ratio->SetMarkerColor(3);
gr3_ratio->SetMarkerStyle(20);
gr3_ratio->SetMarkerSize(0.8);
gr3_ratio->SetLineColor(3);
gr3_ratio->Draw("P");

gr4_ratio->SetLineWidth(0.3);
gr4_ratio->SetMarkerColor(4);
gr4_ratio->SetMarkerStyle(20);
gr4_ratio->SetMarkerSize(0.8);
gr4_ratio->SetLineColor(4);
gr4_ratio->Draw("P");

gr5_ratio->SetLineWidth(0.3);
gr5_ratio->SetMarkerColor(6);
gr5_ratio->SetMarkerStyle(21);
gr5_ratio->SetMarkerSize(1);
gr5_ratio->SetLineColor(6);
gr5_ratio->Draw("P");

TLegend *leg = new TLegend(0.3, 0.5, 0.6, 0.8);
leg->SetFillColor(0);
leg->AddEntry(gr1_ratio, "Geant4 - Livermore", "lp");
leg->AddEntry(gr2_ratio, "Geant4 - Penelope", "lp");
leg->AddEntry(gr3_ratio, "Geant4 - opt0", "lp");
leg->AddEntry(gr4_ratio, "Geant4 - opt3", "lp");
leg->AddEntry(gr5_ratio, "Geant4 - opt4", "lp");
leg->Draw();

c1 -> Print("Flexi_dose_rate_distribution_ratio.jpg");

}
