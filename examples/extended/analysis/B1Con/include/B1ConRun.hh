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
/// \file B1ConRun.hh
/// \brief Definition of the B1ConRun class

#ifndef B1ConRun_h
#define B1ConRun_h 1

#include "B1Run.hh"
#include "globals.hh"

class G4Event;

/// Run class which extends B1Run
///

class B1ConRun : public B1Run
{
  public:
    B1ConRun();
    virtual ~B1ConRun();

    // method from the base class
    virtual void Merge(const G4Run*);
    virtual void AddEdep (G4double edep); 

    // get methods
    G4int GetNumberOfEvent() const { return (G4int)fEdepEventVector.size(); }
    G4double GetEdepPerEvent(G4int i) const { return fEdepEventVector[i]; }

  private:
    std::vector<G4double> fEdepEventVector;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

