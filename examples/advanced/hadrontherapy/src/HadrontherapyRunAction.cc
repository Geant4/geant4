//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************

#include "HadrontherapyRunAction.hh"
#include "HadrontherapyEventAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4Timer.hh"
#include "HadrontherapyRunAction.hh"
#ifdef G4ANALYSIS_USE
#include "HadrontherapyAnalysisManager.hh"
#endif

HadrontherapyRunAction::HadrontherapyRunAction(G4String &SDNAME)
{
  extern G4String sensitiveDetectorName;
  sensitiveDetectorName = SDNAME;
  detector = new HadrontherapyDetectorConstruction(sensitiveDetectorName);
}

HadrontherapyRunAction::~HadrontherapyRunAction()
{ 
  delete detector; 
}

void HadrontherapyRunAction::BeginOfRunAction(const G4Run*aRun)
{ 
#ifdef G4ANALYSIS_USE
  HadrontherapyAnalysisManager* analysis = HadrontherapyAnalysisManager::getInstance();
  analysis->book();
#endif  	
  //inizializza la matrice di dose
  extern G4double matrix[40][40][40];

  G4int indexI; 
  G4int indexJ; 
  G4int indexK;

  for (indexI = 0; indexI < 40; indexI++){
    for (indexJ = 0; indexJ < 40; indexJ++){
      for (indexK = 0; indexK < 40; indexK++){
	matrix[indexI][indexJ][indexK] = 0;
      }}};
  
 G4RunManager::GetRunManager()->SetRandomNumberStore(true);


}

void HadrontherapyRunAction::EndOfRunAction(const G4Run* aRun)
{

  extern G4double matrix[40][40][40];
  G4int indexI; 
  G4int indexJ; 
  G4int indexK;
  G4int count;
  count = 0;

#ifdef G4ANALYSIS_USE
  HadrontherapyAnalysisManager* analysis = HadrontherapyAnalysisManager::getInstance();  
  //Scrittura della matrice di dose su file
  for (indexI = 0; indexI < 40; indexI++){
    for (indexJ = 0; indexJ < 40; indexJ++){
      for (indexK = 0; indexK < 40; indexK++){
	
	std::ofstream pmtfile("matrice.out", std::ios::app);
	if(pmtfile.is_open())
	  
	  {
           if (matrix[indexI][indexJ][indexK] != 0)
	     {
	      pmtfile << matrix[indexI][indexJ][indexK]/MeV << G4endl; 
            
              analysis -> energyDeposit3D(count, indexI, indexJ, indexK,
                                         matrix[indexI][indexJ][indexK]/MeV);	              count ++;
	     }
	  }
      }}};

  analysis -> finish();
#endif
  }




