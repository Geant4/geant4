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
/// \file medical/electronScattering2/include/ElectronRun.hh
/// \brief Definition of the ElectronRun class

#ifndef ELECTRONRUN_HH
#define ELECTRONRUN_HH

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4THitsMap.hh"
#include <map>

class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ElectronRun : public G4Run {
    
public:
    ElectronRun();
    virtual ~ElectronRun();
    
    virtual void RecordEvent(const G4Event*);
    virtual void Merge(const G4Run*);
    void DumpData(G4String&) const;

private:
    void Print(const std::vector<G4String>& title,
               const std::map< G4int, std::vector<G4double> >&out,
               G4String&) const;
    
    std::map<G4int, G4THitsMap<G4double>* > fMap;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
