// Function for the computation of the channeling efficiency
Double_t ComputeEfficiency(TH1D *h1,Double_t *fPar){
    
    Double_t vXMin = fPar[1] - 3. * fPar[2];
    Double_t vXMax = fPar[1] + 3. * fPar[2];
    
    Int_t vBinMin = h1->FindBin(vXMin);
    Int_t vBinMed = h1->FindBin(fPar[1]);
    Int_t vBinMax = h1->FindBin(vXMax);

    Double_t vEfficiency = h1->Integral(vBinMin+1,vBinMax-1);
    vEfficiency += h1->GetBinContent(vBinMin)*(h1->GetBinLowEdge(vBinMin+1)-vXMin)/h1->GetBinWidth(vBinMin);
    vEfficiency += h1->GetBinContent(vBinMax)*(-h1->GetBinLowEdge(vBinMax)+vXMax)/h1->GetBinWidth(vBinMax);
    vEfficiency /= h1->Integral(0,-1);
    vEfficiency *= 100.;
    
    return vEfficiency;
}

// Function for the computation of channeling efficiency at various incoming angle
Int_t AnalyseChannelingEfficiency(TTree *fTree,Float_t fChannelingMinimum = 35., Float_t fChannelingMaximum = 70.){
    //**//Channeling Gaussian Fit Function
    TF1 *vChanneling = new TF1("vChanneling","gaus",fChannelingMinimum,fChannelingMaximum);
    vChanneling->SetParNames("Const","Mean","Sigma");
    vChanneling->SetLineColor(4);
    vChanneling->SetLineStyle(2);

    TH2D *hChannelingPlot = new TH2D("hChannelingPlot","Deflection Angle vs. Incoming Angle;Horizontal Incoming Angle [#murad];Horizontal Deflection Angle [#murad]",21,-10.5,10.5,256,-127.5,128.5);
    
    TH1F *hChannelingEfficiency = new TH1F("hChannelingEfficiency","G4Channeling;Horizontal Incoming Angle [#murad];Efficiency [%]",21,-10.5,10.5);

    fTree->Draw("(angXout-angXin):angXin>>hChannelingPlot");
    
    Double_t vNormalizationToAmorphous = 0.965; // Normalization for channeling efficiency, see PRSTAB 11, 063501 (2008)
    
    for(int i=2;i<=21;i++){
        TH1D* h1 = hChannelingPlot->ProjectionY("h1",i,i);
        h1->Fit(vChanneling,"QR");
        Double_t *vChannelingParameters;
        vChannelingParameters = vChanneling->GetParameters();
        hChannelingEfficiency->SetBinContent(i,ComputeEfficiency(h1,vChannelingParameters)/vNormalizationToAmorphous);
        h1->Delete();
    }
    hChannelingEfficiency->SetLineColor(3);
    hChannelingEfficiency->SetLineStyle(4);
    hChannelingEfficiency->SetMarkerColor(3);
    hChannelingEfficiency->SetFillStyle(0);
    hChannelingEfficiency->SetMarkerStyle(20);
    hChannelingEfficiency->Draw("PL");

    TGraph* gRoughExperimentalData = new TGraph(11);
    gRoughExperimentalData->SetPoint(	0	,	-10	,	20	);
    gRoughExperimentalData->SetPoint(	1	,	-8	,	38	);
    gRoughExperimentalData->SetPoint(	2	,	-6	,	56	);
    gRoughExperimentalData->SetPoint(	3	,	-4	,	72	);
    gRoughExperimentalData->SetPoint(	4	,	-2	,	80	);
    gRoughExperimentalData->SetPoint(	5	,	0	,	84	);
    gRoughExperimentalData->SetPoint(	6	,	2	,	82	);
    gRoughExperimentalData->SetPoint(	7	,	4	,	78	);
    gRoughExperimentalData->SetPoint(	8	,	6	,	66	);
    gRoughExperimentalData->SetPoint(	9	,	8	,	52	);
    gRoughExperimentalData->SetPoint(	10	,	10	,	37	);

    gRoughExperimentalData->SetLineColor(4);
    gRoughExperimentalData->SetLineStyle(3);
    gRoughExperimentalData->SetFillStyle(0);
    gRoughExperimentalData->SetFillColor(0);
    gRoughExperimentalData->SetMarkerColor(4);
    gRoughExperimentalData->SetMarkerStyle(21);
    gRoughExperimentalData->SetTitle("Phys. Lett. B 680, 129");

    gRoughExperimentalData->Draw("sameCP");

    TLegend *aLegend = new TLegend(0.30,0.15,0.55,0.3);
    aLegend->AddEntry(hChannelingEfficiency);
    aLegend->AddEntry(gRoughExperimentalData);
    aLegend->SetFillStyle(0);
    aLegend->SetLineColor(0);
    aLegend->Draw();

    return 0;
}
