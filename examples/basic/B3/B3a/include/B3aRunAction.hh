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
<<<<<<< HEAD
// $Id: B3aRunAction.hh 94031 2015-11-05 11:54:38Z ihrivnac $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
/// \file B3aRunAction.hh
/// \brief Definition of the B3aRunAction class

#ifndef B3aRunAction_h
#define B3aRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Parameter.hh"
#include "globals.hh"

/// Run action class

class B3aRunAction : public G4UserRunAction
{
  public:
    B3aRunAction();
    virtual ~B3aRunAction();
    
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void CountEvent()           { fGoodEvents += 1; };
    void SumDose(G4double dose) { fSumDose += dose; };  

private:
    G4Parameter<G4int>    fGoodEvents;
    G4Parameter<G4double> fSumDose;  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

