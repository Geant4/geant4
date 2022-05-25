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
/// \file electromagnetic/TestEm6/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4SystemOfUnits.hh"  
#include "G4PhysicalConstants.hh"  
#include "G4EmCalculator.hh" 
#include "G4ParticleTable.hh" 
#include "G4ParticleDefinition.hh" 
#include "G4Positron.hh" 
#include "G4AnnihiToMuPair.hh"  
#include "G4eeToHadrons.hh"  
#include "G4eeToHadronsModel.hh"  
#include "G4eBremsstrahlung.hh"  

#include "RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4MuonMinus.hh"

#include "Randomize.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det)
 : G4UserRunAction(),fDetector(det),fProcCounter(0),fAnalysis(0),fMat(0)
{
  fMinE = 40*GeV;
  fMaxE = 10000*GeV;
  fnBin = 10000;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //get material
  //
  fMat = fDetector->GetMaterial();
  G4cout << "###RunAction::BeginOfRunAction  Material:" 
         << fMat->GetName() << G4endl;
 
  fProcCounter = new ProcessesCount;

  fAnalysis = G4AnalysisManager::Instance();
  fAnalysis->SetDefaultFileType("root");

  // Open an output file
  //
  std::stringstream tmp;
  tmp << "testem6_" << aRun->GetRunID();
  G4String fileName = tmp.str();
  fAnalysis->OpenFile(fileName);    
  fAnalysis->SetVerboseLevel(2);
  fAnalysis->SetActivation(true);
    
  // Creating histograms 1,2,3,4,5,6
  // 
  fAnalysis->SetFirstHistoId(1);  
  fAnalysis->CreateH1("h1","1/(1+(theta+[g]+)**2)",100, 0 ,1.);
  fAnalysis->CreateH1("h2","log10(theta+ [g]+)",   100,-3.,1.);
  fAnalysis->CreateH1("h3","log10(theta- [g]-)",   100,-3.,1.);
  fAnalysis->CreateH1("h4","log10(theta+ [g]+ -theta- [g]-)", 100,-3.,1.);
  fAnalysis->CreateH1("h5","xPlus" ,100,0.,1.);
  fAnalysis->CreateH1("h6","xMinus",100,0.,1.);
  
  //creating histogram 7,8,9,10,11 (CrossSectionPerAtom)
  //
  G4double minBin = std::log10(fMinE/GeV);
  G4double maxBin = std::log10(fMaxE/GeV);
  fAnalysis->CreateH1("h7","CrossSectionPerAtom of AnnihiToMuMu (microbarn)",
                      fnBin,minBin,maxBin);
  fAnalysis->CreateH1("h8",
    "CrossSectionPerAtom of AnnihiToTwoGamma (microbarn)",fnBin,minBin,maxBin);
  fAnalysis->CreateH1("h9","CrossSectionPerAtom of AnnihiToHadrons (microbarn)",
                      fnBin,minBin,maxBin);
  fAnalysis->CreateH1("h10",
    "Theoretical CrossSectionPerAtom of AnnihiToTwoGamma (microbarn)",
                      fnBin,minBin,maxBin);
  fAnalysis->CreateH1("h11",
    "Theoretical CrossSectionPerAtom of AnnihiToMuMu (microbarn)",
                      fnBin,minBin,maxBin);

  //creating histogram 12,13,14,15,16(CrossSectionPerVolume)
  //
  fAnalysis->CreateH1("h12","CrossSectionPerVol of Bremsstraulung (1/mm) ",
                      fnBin,minBin,maxBin);
  fAnalysis->CreateH1("h13","CrossSectionPerVol of Ionization (1/mm)",
                      fnBin,minBin,maxBin);
  fAnalysis->CreateH1("h14","CrossSectionPerVol of AnnihiToMuMu (1/mm)",
                      fnBin,minBin,maxBin);
  fAnalysis->CreateH1("h15","CrossSectionPerVol of AnnihiToTwoGamma (1/mm)",
                      fnBin,minBin,maxBin);
  fAnalysis->CreateH1("h16","CrossSectionPerVol of AnnihiToHadrons (1/mm)",
                      fnBin,minBin,maxBin);
  
  //creating histogram 17 (R ratio)
  fAnalysis->CreateH1("h17","R : eeToHadr/eeToMu",fnBin,minBin,maxBin);
    
  G4cout << "\n----> Histogram file is opened in " << fileName << G4endl;   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CountProcesses(G4String procName)
{
  //does the process  already encounted ?
  //
  size_t nbProc = fProcCounter->size();
  size_t i = 0;
  while ((i<nbProc)&&((*fProcCounter)[i]->GetName()!=procName)) i++;
  if (i == nbProc) fProcCounter->push_back( new OneProcessCount(procName));
  
  (*fProcCounter)[i]->Count();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  G4cout << "### RunAction::EndOfRunAction" << G4endl;
  //total number of process calls
  //
  G4cout << "\n Number of process calls --->";
  for (size_t i=0; i< fProcCounter->size(); ++i) {
    G4String procName = (*fProcCounter)[i]->GetName();
    if (procName != "Transportation") {
      G4int count = (*fProcCounter)[i]->GetCounter(); 
      G4cout << "\t" << procName << " : " << count;
    }
  }
  G4cout << G4endl;
  
  //instance EmCalculator
  //
  G4EmCalculator emCal;
  emCal.SetVerbose(0);

  //get positron
  //
  G4String positronName = "e+";
  G4ParticleDefinition* positron = 
    G4ParticleTable::GetParticleTable()->FindParticle(positronName);

  //process name
  //
  G4String annihilName      = "annihil";
  G4String annihiToMuName   = "AnnihiToMuPair";
  G4String annihiToHadrName = "ee2hadr";
  G4String BremName         = "eBrem";
  G4String IoniName         = "eIoni";
  
  //get AnnihiToMuPair
  //
  G4AnnihiToMuPair* annihiToMu = 
    reinterpret_cast<G4AnnihiToMuPair*>(emCal.FindProcess(positron,
                                                          annihiToMuName));

  //parameters for ComputeCrossSection
  //
  G4double atomicZ = 1.;
  G4double atomicA = 2.;  

  //set parameters for theory
  //
  const G4ParticleDefinition* muon = G4MuonMinus::MuonMinus();
  G4double Mu = muon->GetPDGMass();
  G4double Me = electron_mass_c2;
  G4double Re = classic_electr_radius;
  G4double Ru = Re*Me/Mu;
  G4double Eth = 2*Mu*Mu/Me-Me;
  G4PhysicsLogVector v(fMinE, fMaxE, fnBin, false);
  
  //Compute CrossSections and Fill histgrams
  //
  for(G4int i=0; i<=fnBin; ++i) {
    G4double energy = v.Energy(i);
    G4double x = std::log10(energy/GeV);

    //CrossSectionPerAtom
    //
    G4double crs_annihiToMu = 
      annihiToMu->ComputeCrossSectionPerAtom(energy,atomicZ);
    //G4cout << "crs_annihiToMu(mkb)=" << crs_annihiToMu/microbarn << G4endl;
    fAnalysis->FillH1(7,x,crs_annihiToMu/microbarn);
  
    G4double crs_annihil = 
      emCal.ComputeCrossSectionPerAtom(energy,positron,annihilName,
                                       atomicZ,atomicA);
    fAnalysis->FillH1(8,x,crs_annihil/microbarn);

    G4double crs_annihiToHadr = 
      emCal.ComputeCrossSectionPerAtom(energy,positron,annihiToHadrName,
                                       atomicZ,atomicA);
    fAnalysis->FillH1(9,x,crs_annihiToHadr/microbarn);

    //CrossSectionPerVolume
    //
    G4double crsVol_Brem = 
      emCal.ComputeCrossSectionPerVolume(energy,positron,BremName,fMat,100*keV);
    fAnalysis->FillH1(12,x,crsVol_Brem*mm);

    G4double crsVol_Ioni = 
      emCal.ComputeCrossSectionPerVolume(energy,positron,IoniName,fMat,100*keV);
    fAnalysis->FillH1(13,x,crsVol_Ioni*mm);
                
    G4double crsVol_annihiToMu = annihiToMu->CrossSectionPerVolume(energy,fMat);
    fAnalysis->FillH1(14,x,crsVol_annihiToMu*mm);
  
    G4double crsVol_annihil = 
      emCal.ComputeCrossSectionPerVolume(energy,positron,annihilName,fMat);
    fAnalysis->FillH1(15,x,crsVol_annihil*mm);

    G4double crsVol_annihiToHadr = 
      emCal.ComputeCrossSectionPerVolume(energy,positron,annihiToHadrName,fMat);
    fAnalysis->FillH1(16,x,crsVol_annihiToHadr*mm);

    //R ratio
    //
    G4double RR = 0.0;
    if(crsVol_annihiToMu > 0.) RR = crsVol_annihiToHadr/crsVol_annihiToMu;
    fAnalysis->FillH1(17,x,RR);

    //Theoretical calculation
    //
    G4double X1 = energy/Me;
    if(X1>1 && i%1000==0){
      G4double crs_annihil_theory = atomicZ*pi*Re*Re*
        ( (X1*X1+4*X1+1)*G4Log(X1+std::sqrt(X1*X1-1))/(X1*X1-1)
         -(X1+3)/std::sqrt(X1*X1-1) )/(X1+1);
      fAnalysis->FillH1(10,x,crs_annihil_theory/microbarn);
    }

    G4double X2 = Eth/energy;
    if(X2<1. && i%1000==0){
      G4double crs_annihiToMu_theory = 
	atomicZ*pi*Ru*Ru/3*X2*(1+X2/2)*std::sqrt(1-X2);
      fAnalysis->FillH1(11,x,crs_annihiToMu_theory/microbarn);
    }

    //if(i%1000==0)G4cout <<"###energy:" << energy << "/crs_ToMuMu:"  
    //        << crs_annihiToMu << "/crs_ToTwoGamma:"<< crs_annihil 
    //        <<"/crs_ToToHadr:"<<crs_annihiToHadr<< G4endl;
  }
 
  fAnalysis->Write();
  fAnalysis->CloseFile();
  fAnalysis->Clear();

  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
