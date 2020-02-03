{
// Read reference data in dolan.txt
FILE * fg1=fopen("dolan.txt", "r");
Int_t n_points_dolan = 8;
Float_t x1[n_points_dolan], y1[n_points_dolan];
Float_t x, y;
Int_t ncols_dolan;
Int_t nlines1 = 0;

while(1)
{
	ncols_dolan = fscanf(fg1, "%f %f", &x, &y);
	if (ncols_dolan<0) break;
	x1[nlines1]= x;
	y1[nlines1] = y;
	nlines1++;
}

fclose(fg1);

// Read the results of the brachyadvanced example
// OncuraSourceMacro.mac with 20 B events

FILE *fg2=fopen("geant4_6711_dose.txt", "r");
Int_t n_points_geant4 = 398;
Float_t x2[n_points_geant4], y2[n_points_geant4];
Int_t ncols_geant4;
Int_t nlines2 = 0;

while(1)
{
	ncols_geant4 = fscanf(fg2, "%f %f", &x, &y);
	if (ncols_geant4<0) break;
	x2[nlines2] = x;
	y2[nlines2] = y;
	nlines2++;
}

fclose(fg2);

TGraph *gr1 = new TGraph (nlines1, x1, y1);
TGraph *gr2 = new TGraph (nlines2, x2, y2);

TCanvas *c1 = new TCanvas("c1", "Graph Draw Options", 200, 10, 600, 400);

gPad->SetLogy();

// Draw the graph with axis, continuous line, and put a '*' at each point

gr1->SetTitle("Dose rate distribution");
gr1->GetXaxis()->SetTitle("Distance from the centre (cm)");
gr1->GetYaxis()->SetTitle("Normalised dose rate distribution");
gr1->SetLineWidth(1);
gr1->SetMarkerColor(1);
gr1->SetMarkerStyle(20);
gr1->Draw("AP");

gr2->SetLineWidth(1);
gr2->SetMarkerColor(2);
gr2->SetMarkerStyle(21);
gr2->SetMarkerSize(0.5);
gr2->SetLineColor(2);
gr2->Draw("CP");

TLegend *leg = new TLegend(0.3, 0.5, 0.6, 0.8);
leg->SetFillColor(0);
leg->AddEntry(gr1, "Reference data", "lp");
leg->AddEntry(gr2, "Geant4 - 20B Events", "lp");
leg->Draw();

}


