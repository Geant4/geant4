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
// $Id$
// 
/// \file B4bEventAction.hh
/// \brief Definition of the B4bEventAction class

#ifndef B4bEventAction_h
#define B4bEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4GenericMessenger;

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy 
/// deposit and track lengths of charged particles in Absober and Gap layers 
/// stored in B4bRunData object.
///
/// The data member fPrintModulo defines the frequency of printing.
/// the accumulated quantities. Its value can be changed via a command
/// defined using G4GenericMessenger class:
/// - /B4/event/setPrintModulo value

class B4bEventAction : public G4UserEventAction
{
  public:
    B4bEventAction();
    virtual ~B4bEventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
                     
    void SetPrintModulo(G4int value);
    
  private:
    // methods
    void PrintEventStatistics(G4double absoEdep, G4double absoTrackLength,
                              G4double gapEdep, G4double gapTrackLength) const;
    
    // data members
    G4GenericMessenger*  fMessenger;
                     
    G4int  fPrintModulo;
};

// inline functions

inline void B4bEventAction::SetPrintModulo(G4int value) {
  fPrintModulo = value;
}
                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
