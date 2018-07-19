// -------------------------------------------------------------------
// $Id: plotDend.C $
// -------------------------------------------------------------------
//
// *********************************************************************
// To execute this macro under ROOT, 
// launch ROOT (usually type 'root' at your machine's prompt)
// This macro needs Dend3DEdep.out file : 
// *********************************************************************
{ 
gROOT->Reset();
gStyle->SetOptStat(0000);

c1 = new TCanvas ("c1","",20,20,1200,600);
c1->Divide(2,1); 

Int_t ncols=0;
Int_t nlines = 0;

FILE * fp = fopen("Dend3DEdep.out","r");
Float_t posX, posY, posZ ;
Float_t distB, distA, EdepR, DoseR;
Float_t distMaxA = -1e-9;
Float_t distMaxB = -1e-9;
Float_t edepMax = -1e-9;
Float_t doseMax = -1e-9;
Float_t edepMin = 1e9;
Float_t doseMin = 1e9;

h1 = new TProfile("Energy", "Energy deposits (keV) in dendritic compartments",1000,-1000,1000,0.001,1000);  
h2 = new TProfile("Dose", "Dose deposits (Gy) in dendritic compartments",1000,-1000,1000,0.001,1000);  
while (1) 
   {
      ncols = fscanf(fp," %f %f %f %f %f %f %f",&posX, &posY, &posZ, &distA, &distB, &EdepR, &DoseR);
      if (ncols < 0) break;
      if (distMaxA < distA ) distMaxA = distA ;
      if (distMaxB < distB ) distMaxB = distB ;
      if (edepMax < EdepR ) edepMax = EdepR ;
      if (doseMax < DoseR ) doseMax = DoseR ;
      if (edepMin > EdepR ) edepMin = EdepR ;
      if (doseMin > DoseR ) doseMin = DoseR ;
      // ....  
      h1->Fill(-distB, EdepR);	// Basal dendrite
      h1->Fill(distA, EdepR);	// Apical dendrite
      h2->Fill(-distB, DoseR);
      h2->Fill(distA, DoseR);
      nlines++;	
   }
fclose(fp);
cout << " Max and Min Energy deposits (keV) ==  " << edepMax << " ; "<< edepMin<<endl;
cout << " Max and Min Dose deposits (Gy) ==  " << doseMax << " ; "<< doseMin<<endl;
cout << " Maximum Basal Distance (um) == " << distMaxB << " "<<endl;
cout << " Maximum Apical Distance (um) == " << distMaxA << " "<<endl;

c1->cd(1); 
h1->Draw("P"); 
//gPad->SetLogy();
h1->SetMarkerSize(2);
h1->SetMarkerColor(4);
h1->SetMarkerStyle(27); 
h1->SetFillStyle(3005);
//h1->GetYaxis()->SetTitle("Energy deposits in basal and apical dendrite (keV)");
h1->GetXaxis()->SetTitle("Distance from Soma (um)");
h1->GetYaxis()->SetRangeUser(edepMin, edepMax+3.);
if (distMaxB > 0.)
{
    TLatex text(-distMaxB,edepMax-1.,"Basal");
    text.DrawClone();
}
if (distMaxA > 0.)
{
    TLatex text(distMaxA/3.,edepMax-2.,"Apical");
    text.DrawClone();
}
h1->GetXaxis()->SetRangeUser(-distMaxB-10., distMaxA+10.);
//h1->Fit("gaus");

c1->cd(2);
h2->Draw("P");	
//gPad->SetLogy();
h2->SetMarkerSize(2);
h2->SetMarkerColor(kRed);
h2->SetMarkerStyle(27); 
h2->SetFillStyle(3005);
//h2->GetYaxis()->SetTitle("Dose deposits in basal and apical dendrite (Gy)");
h2->GetXaxis()->SetTitle("Distance from Soma (um)");
h2->GetYaxis()->SetRangeUser(doseMin, doseMax+0.3);
if (distMaxB > 0.)
{
    TLatex text(-distMaxB,doseMax-0.1,"Basal");
    text.DrawClone();
}
if (distMaxA > 0.)
{
    TLatex text(distMaxA/3.,doseMax-0.2,"Apical");
    text.DrawClone();
}
h2->GetXaxis()->SetRangeUser(-distMaxB-10., distMaxA+10.);

}

