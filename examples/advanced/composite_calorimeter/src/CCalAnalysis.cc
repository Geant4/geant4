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
///////////////////////////////////////////////////////////////////////////////
// File: CCalAnalysis.cc
// Description: CCalAnalysis interfaces all user analysis code
///////////////////////////////////////////////////////////////////////////////


#include "G4RunManager.hh" 

#include "CCalAnalysis.hh"
#include "CCalutils.hh"

CCalAnalysis* CCalAnalysis::instance = 0;
 
CCalAnalysis::CCalAnalysis() :
  energy(0), hcalE(0), ecalE(0), timeHist(0), lateralProfile(0),
  timeProfile(0)
{

#ifdef debug
  fVerbosity = 1;
#else
  fVerbosity = 0;  
#endif

  numberOfTimeSlices = 200;

  // Create analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);
  analysisManager->SetFirstNtupleId(1);

  // Open an output file
  analysisManager->OpenFile("ccal");
  G4cout << "********************************************" << G4endl
	 << "* o/p file ccal"  << G4endl
	 << "********************************************" << G4endl 
	 << G4endl;
  

  // Create a tuple :
  // Create ntuple
  analysisManager->CreateNtuple("ntuple1", "Event info");
  char tupleid[7];
  for (int i=0;i<28;i++)
    {
      sprintf(tupleid,"hcal%d",i);
      analysisManager->CreateNtupleFColumn(tupleid);
    }
  for (int i=0; i<49; i++) 
    {
      sprintf(tupleid,"ecal%d",i);
      analysisManager->CreateNtupleFColumn(tupleid);
    }
  analysisManager->CreateNtupleFColumn("ELAB");
  analysisManager->CreateNtupleFColumn("XPOS");
  analysisManager->CreateNtupleFColumn("YPOS");
  analysisManager->CreateNtupleFColumn("ZPOS");
  analysisManager->CreateNtupleFColumn("EDEP");
  analysisManager->CreateNtupleFColumn("EDEC");
  analysisManager->CreateNtupleFColumn("EHDC");
  analysisManager->FinishNtuple();

  //
  //Create histograms. Save the ID of the first histogram in each block
  //
  char id[5], ntupletag[50];
  //Energy deposit in Hcal layers 
  for (int i = 0; i<28; i++) {
    sprintf(id,"h%d",i+100);
    sprintf(ntupletag, "Energy Deposit in Hcal Layer%d   in GeV",i);
    G4int histoID = analysisManager->CreateH1(id,ntupletag, 100, 0., 1.0);
    if (!i)
      hcalE = histoID;
  }
  // Energy deposits in Ecal towers
  for (int i = 0; i<49; i++) 
    {
      sprintf(id, "h%d",i+200);
      sprintf(ntupletag, "Energy Deposit in Ecal Tower%d   in GeV",i);
      G4int histoID = analysisManager->CreateH1(id,ntupletag, 100, 0., 1.0);
      if (!i)
	ecalE = histoID;
    }
  // Total energy deposit
  energy = analysisManager->CreateH1("h4000", "Total energy deposited in GeV", 
				     100, 0., 100.0);

  // Time slices	  
  for (int i=0; i<numberOfTimeSlices; i++){
    sprintf(id, "h%d",i+300);
    sprintf(ntupletag, "Time slice %d nsec energy profile in GeV",i);
    G4int histoID = analysisManager->CreateH1(id,ntupletag,100, 0., 100.0);
    if (!i)
      timeHist = histoID;
  }

  // Profile of lateral energy deposit in Hcal
  for (int i = 0; i<70; i++) {
    sprintf(id, "h%d",i+500);
    sprintf(ntupletag, "Lateral energy profile at %d cm  in GeV",i);
    G4int histoID = analysisManager->CreateH1(id, ntupletag, 100, 0., 10.0);
    if (!i)
      lateralProfile = histoID;
  }

  // Time profile 
  timeProfile = analysisManager->CreateH1("h901", "Time Profile in Sensitive Detector", 
						200, 0., 200.);
  analysisManager->CreateH1("h902", "Time Profile in Sensitive+Passive", 
			    200, 0., 200.);
}


CCalAnalysis::~CCalAnalysis() {  
  Finish();
  delete G4AnalysisManager::Instance();
}


void CCalAnalysis::Init() {
}                       

void CCalAnalysis::Finish() 
{;}             

CCalAnalysis* CCalAnalysis::getInstance() {
  if (instance == 0) instance = new CCalAnalysis();
  return instance;
}


// This function fill the 1d histogram of the energies in HCal layers
void CCalAnalysis::InsertEnergyHcal(float* v) 
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  G4double totalFilledEnergyHcal = 0.0;
  for (int i=0; i<28; i++) {
    G4double x = v[i];
    man->FillH1(hcalE+i,x);
    if (fVerbosity)
      {
	G4cout << "Fill Hcal histo " << i << " with " << x << G4endl;
	totalFilledEnergyHcal += x;
      }
  }  

  if (fVerbosity)
    G4cout << 
      "CCalAnalysis::InsertEnergyHcal: Total filled Energy Hcal histo " 
	   << totalFilledEnergyHcal << G4endl;
}


// This function fill the 1d histogram of the energies in ECal layers
void CCalAnalysis::InsertEnergyEcal(float* v) 
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  
  G4double totalFilledEnergyEcal = 0.0;

  for (G4int i=0; i<49; i++) {    
    G4double x = v[i];
    man->FillH1(ecalE+i,x);
    if (fVerbosity)
      {
	G4cout << "Fill Ecal histo " << i << " with " << x << G4endl;
	totalFilledEnergyEcal += x;
      }
  }
  if (fVerbosity)
    G4cout << 
      "CCalAnalysis::InsertEnergyEcal: Total filled Energy Ecal histo " 
	   << totalFilledEnergyEcal << G4endl;
}


// This function fill the 1d histogram of the lateral profile
void CCalAnalysis::InsertLateralProfile(float* v) 
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  G4double totalFilledProfileHcal = 0.0;

  for (G4int i=0; i<70; i++) {    
      G4double x = v[i];
      man->FillH1(lateralProfile+1,x);
      if (fVerbosity)
	{
	  G4cout << "Fill Profile Hcal histo " << i << " with " << x << G4endl;
	  totalFilledProfileHcal += x;
	}    
  }
  if (fVerbosity)
    G4cout << "CCalAnalysis::InsertLateralProfile: Total filled Profile Hcal"
	   << " histo " << totalFilledProfileHcal << G4endl;
}


// This function fill the 1d histogram of the energy 
void CCalAnalysis::InsertEnergy(float v) 
{
 
  G4double x = v;
  G4AnalysisManager::Instance()->FillH1(energy,x);
  if (fVerbosity)
    G4cout << "CCalAnalysis::InsertEnergy: Fill Total energy Hcal histo with " 
	   << x << G4endl;
}


// This function fill the 1d histograms of time profiles at stepping action
void CCalAnalysis::InsertTime(float* v) 
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  
  G4double totalFilledTimeProfile = 0.0;
  for (G4int j=0; j<numberOfTimeSlices; j++) 
    {
      G4double x = v[j];
      man->FillH1(timeHist+j,x);
      if (fVerbosity)
	{
	  G4cout << "Fill Time slice histo " << j << " with " << x << G4endl;
	  totalFilledTimeProfile += x;
	}    

      G4double t = j + 0.5;
      man->FillH1(timeProfile+1,t,x);
      if (fVerbosity)
	G4cout << "Fill Time profile histo 1 with " << t << " " << x << G4endl;
    }
  if (fVerbosity)
    G4cout << "CCalAnalysis::InsertTime: Total filled Time profile histo " 
	   << totalFilledTimeProfile << G4endl;
}


// This function fill the 1d histograms of time profiles in SD
void CCalAnalysis::InsertTimeProfile(int hit, double time, double edep) 
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillH1(timeProfile,time,edep);

  if (fVerbosity)
    G4cout << "CCalAnalysis:: Fill Time Profile with Hit " << hit
	   << " Edeposit " << edep << " Gev at " << time << " ns" << G4endl;
}



void CCalAnalysis::setNtuple(float* HCalE, float* ECalE, float elab, 
			     float x, float y, float z, float edep, 
			     float edec, float edhc) 
{
  G4AnalysisManager* man = G4AnalysisManager::Instance();  
  G4int counter=0;
  for (int i=0; i<28; i++) 
    {
      man->FillNtupleFColumn(counter,HCalE[i]);
      counter++;
    }
  for (int i=0; i<49; i++) 
    {
      man->FillNtupleFColumn(counter,ECalE[i]);
      counter++;
    }
  man->FillNtupleFColumn(counter,elab);
  man->FillNtupleFColumn(counter+1,x);
  man->FillNtupleFColumn(counter+2,y);
  man->FillNtupleFColumn(counter+3,z);
  man->FillNtupleFColumn(counter+4,edep);
  man->FillNtupleFColumn(counter+5,edec);
  man->FillNtupleFColumn(counter+6,edhc);

  man->AddNtupleRow();  

  if (fVerbosity)
    G4cout << "CCalAnalysis:: Fill Ntuple " << G4endl;
}


/* 
   This member reset the histograms and it is called at the begin
   of each run; here we put the inizialization so that the histograms have 
   always the right dimensions depending from the detector geometry
*/
void CCalAnalysis::BeginOfRun(G4int )  
{ 
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  for (int i=0; i<numberOfTimeSlices; i++){
    man->GetH1(timeHist+1)->reset();
  }
}

//============================================================================

//  This member is called at the end of each run 
void CCalAnalysis::EndOfRun(G4int )  
{
  // Save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
}


// This member is called at the end of every event 
void CCalAnalysis::EndOfEvent(G4int flag) {
  // The plotter is updated only if there is some
  // hits in the event
  if (!flag) return;
}


