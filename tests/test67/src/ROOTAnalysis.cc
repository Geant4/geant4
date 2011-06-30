//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: ROOTAnalysis.cc,v 1.1 2009/03/21 18:37:27 vnivanch Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#ifdef G4_USE_ROOT

#include "ROOTAnalysis.hh"
#include "G4Version.hh"

ROOTAnalysis* ROOTAnalysis::instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
ROOTAnalysis::ROOTAnalysis() : 
  fFile(0),fHistograms(0),fGraph(0),fGraph2(0),fTree(0)
{
  fPhysicsListName = new TObjString("undefined");
  fVersionName = new TObjString((TString)G4Version);
  fPE = 0;
  fIon = 0;
  fComp = 0;
  fEfficiency = 0;
  fEfficiencyErr = 0;
  fEfficiencyFull = 0;
  fEfficiencyFullErr = 0;
  fPrimaryEnergy = 0;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ROOTAnalysis::~ROOTAnalysis()
{
  //loop over all histograms
  std::map<G4int,TH1D*>::iterator i;
  if (fHistograms)
    {
      for (i=fHistograms->begin(); i != fHistograms->end(); i++)        
	delete i->second;        
      delete fHistograms;
    }
  if (fTree)
    delete fTree;
  if (fGraph)
    delete fGraph;
  if (fGraph2)
    delete fGraph2;
  if (fVersionName)
    delete fVersionName;
  if (fPhysicsListName)
    delete fPhysicsListName;
  if (fFile) 
    delete fFile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ROOTAnalysis* ROOTAnalysis::getInstance()
{
  if (instance == 0) instance = new ROOTAnalysis();
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ROOTAnalysis::BookNewHistogram(G4int runID, G4double primaryEnergy)
{
  if (!fFile)
    {
      //Use the version name for the ROOT file
      TString filename;
      //remove blank spaces
      for (size_t i = 0; i < G4Version.length(); i++)
        if (G4Version[i] != ' ' && G4Version[i] != '$') filename += G4Version[i];
      //check if there is "Name" in front, if so remove 5 chars      
      if (filename[0] == 'N')
	{
	  G4String vs = (G4String) filename;
	  filename = (TString) vs.substr(5,vs.size()-5);
	}
      filename += "_";
      filename += fPhysicsListName->GetString();     
      filename += ".root";
      fFile = new TFile(filename,"RECREATE");
    }
  if (!fHistograms)
    fHistograms = new std::map<G4int,TH1D*>;
  if (!fGraph)
    fGraph = new TGraphErrors();
  if (!fGraph2)
    fGraph2 = new TGraphErrors();
  if (!fTree)
    {
      fTree = new TTree("tree","Global results");
      fTree->Branch("energy",&fPrimaryEnergy,"energy/D");
      fTree->Branch("efficiency",&fEfficiency,"efficiency/D");
      fTree->Branch("efferror",&fEfficiencyErr,"efferror/D");
      fTree->Branch("effTotal",&fEfficiencyFull,"effTotal/D");
      fTree->Branch("effTotalerror",&fEfficiencyFullErr,"effTotalerror/D");
      fTree->Branch("Compton",&fComp,"Compton/D");
      fTree->Branch("Ionisation",&fIon,"Ionisation/D");
      fTree->Branch("PhotoEl",&fPE,"PhotoEl/D");
    }

  if (!(fHistograms->count(runID)))
    {
      TString id = "h";
      id += runID;
      TString title = "Histo ";
      title += primaryEnergy/keV; 
      title += " keV";
      Int_t nbins = (Int_t) ((primaryEnergy+10*keV)/(0.5*keV));
      TH1D* histo = new  TH1D(id,title,nbins,
			      0,(primaryEnergy+10*keV)/keV);
      fHistograms->insert(std::make_pair(runID,histo));
					
    }
  //Zero counters
  fPE = 0;
  fIon = 0;
  fComp = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ROOTAnalysis::AddEventEnergy(G4int runID,G4double val)
{
  if (!(fHistograms->count(runID)))
    {
      G4cout << "Problem! Histogram for run " << runID << " not booked" << G4endl;
      G4Exception();
      return;
    }
  TH1D* theHisto = (fHistograms->find(runID))->second;
  theHisto->Fill(val/keV);
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ROOTAnalysis::AddSecondary(G4String name)
{
  if (name == "compt")
    fComp += 1.0;
  else if (name == "eIoni")
    fIon += 1.0;
  else if (name == "phot")
    fPE += 1.0;
  return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ROOTAnalysis::CloseFile()
{
  if (!fFile) //file not created at all: e.g. for a vis-only execution
    return;
  if (!fFile->IsOpen())
    {
      G4cout << "Problema in ROOTAnalysis::CloseFile()" << G4endl;
      G4Exception();
      return;
    }
  fFile->cd();
  //loop over all histograms
  std::map<G4int,TH1D*>::iterator i;
  if (fHistograms)
    {
      for (i=fHistograms->begin(); i != fHistograms->end(); i++)
        {
	  TH1D* hi = i->second;
          hi->Write(hi->GetName());
        }
    }
  if (fVersionName)
   fVersionName->Write("Version");
  if (fPhysicsListName)
   fPhysicsListName->Write("List");
  if (fGraph)  
   fGraph->Write("PeakEff");
  if (fGraph2)
   fGraph2->Write("FullEff");
  if (fTree)
   fTree->Write(fTree->GetName());
  fFile->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ROOTAnalysis::SetListName(G4String name)
{
  if (!fPhysicsListName)
    {
      G4cout << "Problem in ROOTAnalysis::SetListName()" << G4endl;
      G4Exception();
      return;
    }
  fPhysicsListName->SetString((TString) name);
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ROOTAnalysis::EndOfRun(G4int runID,G4double energy,G4int nEventsTotal,
			    G4int nEventsGood,G4int nEventsFull)
{
  if (!fGraph || !fGraph2 || !nEventsTotal)
    {
      G4cout << "Problem in ROOTAnalysis::AddGraphPoint()" << G4endl;
      G4Exception();
      return;
    }
  fEfficiency = ((Double_t) nEventsGood)/((Double_t) nEventsTotal);
  fPrimaryEnergy = energy/keV;
  fEfficiencyErr = std::sqrt(nEventsGood)/nEventsTotal;
  fEfficiency *= 100.; //so it is in %
  fEfficiencyErr *= 100.;
  fGraph->SetPoint(runID,fPrimaryEnergy,fEfficiency);
  fGraph->SetPointError(runID,0.,fEfficiencyErr);

  fEfficiencyFull = ((Double_t) nEventsFull)/((Double_t) nEventsTotal);
  fEfficiencyFullErr = std::sqrt(nEventsFull)/nEventsTotal;  
  fEfficiencyFull *= 100.; //so it is in %
  fEfficiencyFullErr *= 100.;
  fGraph2->SetPoint(runID,fPrimaryEnergy,fEfficiencyFull);
  fGraph2->SetPointError(runID,0.,fEfficiencyFullErr);

  //Here write ntuple with fluorescence and efficiency
  fPE /= nEventsTotal;
  fIon /= nEventsTotal;
  fComp /= nEventsTotal;
  fTree->Fill(); 

  return;
}

#endif
