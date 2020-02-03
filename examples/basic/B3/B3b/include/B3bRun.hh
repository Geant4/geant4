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
//
/// \file B3bRun.hh
/// \brief Definition of the B3bRun class

#ifndef B3bRun_h
#define B3bRun_h 1

#include "G4Run.hh"
#include "globals.hh"
#include "G4StatAnalysis.hh"

/// Run class
///
/// In RecordEvent() there is collected information event per event 
/// from Hits Collections, and accumulated statistic for the run 

class B3bRun : public G4Run
{
  public:
    B3bRun();
    virtual ~B3bRun();

    virtual void RecordEvent(const G4Event*);
    virtual void Merge(const G4Run*);
    
  public:
    G4int GetNbGoodEvents() const { return fGoodEvents; }
    G4double GetSumDose()   const { return fSumDose; }    
    G4StatAnalysis GetStatDose() const { return fStatDose; }

  private:
    G4int fCollID_cryst;
    G4int fCollID_patient;   
    G4int fPrintModulo;
    G4int fGoodEvents;
    G4double fSumDose;
    G4StatAnalysis fStatDose;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
