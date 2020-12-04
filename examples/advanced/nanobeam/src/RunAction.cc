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
// Please cite the following paper if you use this software
// Nucl.Instrum.Meth.B260:20-27, 2007

// #define MATRIX_BOUND_CHECK

#include "RunAction.hh"
#include "Analysis.hh"
#include "G4AutoLock.hh"

namespace 
{
  G4Mutex aMutex = G4MUTEX_INITIALIZER;
}

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* pri)
:fDetector(det),fPrimary(pri)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginOfRunAction(const G4Run* /*aRun*/)
{
   // Vector initialization 
   
   fRow=0;
     
   fXVector = CLHEP::HepVector(32);
   fYVector = CLHEP::HepVector(32);
   fThetaVector = CLHEP::HepVector(32);
   fPhiVector = CLHEP::HepVector(32);
   
   // Histograms
  
   // Get/create analysis manager
   G4cout << "##### Create analysis manager " << "  " << this << G4endl;
  
   G4AnalysisManager* man = G4AnalysisManager::Instance();
  
   G4cout << "Using " << man->GetType() << " analysis manager" << G4endl;

   // Open an output file
   man->OpenFile("nanobeam");
   man->SetFirstHistoId(1);
   man->SetFirstNtupleId(1);
  
   // Create 1st ntuple (id = 1)
   man->CreateNtuple("ntuple0", "BeamProfile");
   man->CreateNtupleDColumn("xIn");
   man->CreateNtupleDColumn("yIn");
   man->CreateNtupleDColumn("zIn");
   man->FinishNtuple();
   G4cout << "Ntuple-1 created" << G4endl;

   // Create 2nd htuple (id = 2)
   man->CreateNtuple("ntuple1","Grid");
   man->CreateNtupleDColumn("xIn");
   man->CreateNtupleDColumn("yIn");
   man->CreateNtupleDColumn("e");
   man->FinishNtuple();
   G4cout << "Ntuple-2 created" << G4endl;
 
   // Create 3rd ntuple (id = 3)   
   man->CreateNtuple("ntuple2","Coef");
   man->CreateNtupleDColumn("xIn");
   man->CreateNtupleDColumn("yIn");
   man->CreateNtupleDColumn("thetaIn");
   man->CreateNtupleDColumn("phiIn");
   man->FinishNtuple();
   G4cout << "Ntuple-3 created" << G4endl;

   return;

}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run* /*aRun*/)
{

if (fDetector->GetCoef()==1)
{
        // 17/12/2013 - thanks to A. Dotti
	// CLHEP Matrix inversion (as for CLHEP version 2.1.4.1)  
        // is not thread safe, we need to protect this code.
        // Since this is not performance critical we simply lock 
        // all this part.
	G4AutoLock l(&aMutex);
        //
	
	CLHEP::HepMatrix m;
	
	// VECTOR READING
	// VECTOR READING

	m = CLHEP::HepMatrix(32,32);
	m = fPrimary->GetMatrix();

	G4cout << G4endl;
	G4cout << "===> NANOBEAM LINE INTRINSIC ABERRATION COEFFICIENTS (units of micrometer and mrad) :" << G4endl;
	G4cout << G4endl;

	int inv;

	m.invert(inv);
	CLHEP::HepVector tmp(32,0);
	tmp=m*fXVector;
	CLHEP::HepVector b;
	b=tmp.sub(2,2);   G4cout << "<x|theta>=" << b << G4endl;
	b=tmp.sub(8,8);   G4cout << "<x|theta*delta>=" << b << G4endl;
	b=tmp.sub(10,10); G4cout << "<x|theta^3>=" << b << G4endl;
	b=tmp.sub(12,12); G4cout << "<x|theta*phi^2>=" << b << G4endl;
	m.invert(inv);

	m.invert(inv);
	tmp = m*fThetaVector;
	m.invert(inv);
	b=tmp.sub(2,2); G4cout << "<x|x>=" << b << G4endl;

	m.invert(inv);
	tmp=m*fYVector;
	b=tmp.sub(3,3);   G4cout << "<y|phi>=" << b << G4endl;
	b=tmp.sub(9,9);   G4cout << "<y|phi*delta>=" << b << G4endl;
	b=tmp.sub(11,11); G4cout << "<y|theta^2*phi>=" << b << G4endl;
	b=tmp.sub(13,13); G4cout << "<y|phi^3>=" << b << G4endl;
	m.invert(inv);

	m.invert(inv);
	tmp = m*fPhiVector;
	m.invert(inv);
	b=tmp.sub(3,3); G4cout << "<y|y>=" << b << G4endl;

}

 // Save histograms
 
 G4AnalysisManager* man = G4AnalysisManager::Instance();
 man->Write();
 man->CloseFile();
 
 // Complete clean-up
 
 delete G4AnalysisManager::Instance();

}
