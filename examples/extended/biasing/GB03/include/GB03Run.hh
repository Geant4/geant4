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
/// \file GB03Run.hh
/// \brief Definition of the GB03Run class

#ifndef GB03Run_h
#define GB03Run_h 1

#include "G4Event.hh"
#include "G4Run.hh"

#include <map>
#include <tuple>
#include <utility>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ParticlesCount
{
  public:
    ParticlesCount() = default;
    ~ParticlesCount() = default;

    std::tuple<G4double, G4int, G4double, G4double> GetCounter(const G4String& name)
    {
      return fPartCntr[name];
    };

    void Count(const G4String& name, G4double w, G4double E);

    void Print() const;

    void Merge(const ParticlesCount&);
    void Reset() { fPartCntr.clear(); }

  private:
    std::map<G4String, std::tuple<G4double, G4int, G4double, G4double>> fPartCntr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class GB03Run : public G4Run
{
  public:
    // constructor and destructor.
    GB03Run() = default;
    ~GB03Run() override;

    void CountParticle(const G4String&, G4double, G4double);
    const ParticlesCount* GetPartCounter() const { return &fPartCounter; }

    // Merge run data
    void Merge(const G4Run*) override;

  private:
    ParticlesCount fPartCounter;
};

#endif
