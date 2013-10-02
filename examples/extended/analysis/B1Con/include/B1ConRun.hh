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
// $Id: B1ConRun.hh 66536 2012-12-19 14:32:36Z ihrivnac $
//
/// \file B1ConRun.hh
/// \brief Definition of the B1ConRun class

#ifndef B1ConRun_h
#define B1ConRun_h 1

#include "G4Run.hh"
#include "globals.hh"

class G4Event;

/// Run class
///

class B1ConRun : public G4Run
{
  public:
    B1ConRun();
    virtual ~B1ConRun();

    virtual void RecordEvent(const G4Event*);
    virtual void Merge(const G4Run*);

  public:
    // get methods
    G4double GetEdepRun()  const { return fEdepRun; }
    G4double GetEdep2Run() const { return fEdep2Run; }
    G4int GetNumberOfEvent() const { return (G4int)fEdepEventVector.size(); };
    G4double GetEdepPerEvent(G4int i) const { return fEdepEventVector[i]; };

  private:
    G4double  fEdepRun;
    G4double  fEdep2Run;
    std::vector<G4double> fEdepEventVector;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

