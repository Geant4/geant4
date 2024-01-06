void reco()
{
    const Int_t nCores = 8;
    
    TChain *fHits = new TChain("Hits");
    
    for(Int_t i = 0; i < nCores; i++)
    {
        fHits->AddFile(Form("../build/output0_t%d.root", i));
    }
    
    Int_t entries = fHits->GetEntries();
    
    cout << "Number of hits: " << entries << endl;
    
    Double_t globalTime;
    Double_t x;
    Double_t y;
    Double_t z;
    
    fHits->SetBranchAddress("fGlobalTime", &globalTime);
    fHits->SetBranchAddress("fX", &x);
    fHits->SetBranchAddress("fY", &y);
    fHits->SetBranchAddress("fZ", &z);
    
    Double_t oldGlobalTime = 0.;
    TVector3 oldPos(0., 0., 0.);
    
    TH2F *hPos = new TH2F("hPos", "Recostructed Position in XZ Plane;X [mm]; Z[mm]", 100, -500, 500, 100, -500, 500);
    
    for(Int_t i = 0; i < entries; i++)
    {
        fHits->GetEntry(i);
        
        TVector3 pos(x, y, z);
        
        TVector3 midPoint = 0.5 * (pos + oldPos);
        
        TVector3 LOR = (pos - oldPos).Unit();
        
        //cout << "LOR: " << flush;
        //LOR.Print();
        
        Double_t timeDifference = globalTime - oldGlobalTime;
        
        //cout << "Time difference: " << timeDifference << endl;
        
        Double_t dist = 0.5 * (timeDifference * 299.792458);
        
        TVector3 annihilationPoint = midPoint + LOR * dist;
        
        hPos->Fill(annihilationPoint.X(), annihilationPoint.Z());
        
        //cout << "Annihilation Point: " << flush;
        //annihilationPoint.Print();
        //cout << globalTime - oldGlobalTime << endl;
        //midPoint.Print();
        
        oldPos = pos;
        oldGlobalTime = globalTime;
    }
    
    TCanvas *c1 = new TCanvas();
    hPos->Draw("colz");
    hPos->Print("pos.pdf");
}
