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
// $Id: B1EventInformation.hh 66536 2012-12-19 14:32:36Z ihrivnac $
//
/// \file B1EventInformation.hh
/// \brief Definition of the B1EventInformation class

#ifndef B1EventInformation_h
#define B1EventInformation_h 1

#include "G4VUserEventInformation.hh"
#include "globals.hh"

/// Event information class
///
/// It holds data member fEnergySum for accumulating 
/// the event energy deposit for an event.
/// These data are then used in the run action to compute the dose.

class B1EventInformation : public G4VUserEventInformation
{
  public:
    B1EventInformation();
    virtual ~B1EventInformation();
    
    virtual void Print() const;
    
    G4double GetEnergySum() const { return fEnergySum; }
    void AddEDep(G4double eDep) { fEnergySum += eDep; }
     
  private:
    G4double  fEnergySum;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
