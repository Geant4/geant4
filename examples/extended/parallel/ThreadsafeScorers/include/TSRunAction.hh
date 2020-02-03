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
/// \file parallel/ThreadsafeScorers/include/TSRunAction.hh
/// \brief Definition of the TSRunAction class
//
//
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#ifndef tsrunaction_hh
#define tsrunaction_hh 1

#include "globals.hh"
#include "G4UserRunAction.hh"
#include <vector>
#include <map>
#include <tuple>

class G4Run;
class G4Timer;
class TSDetectorConstruction;

class TSRunAction : public G4UserRunAction
{
public:
    typedef std::tuple<G4double, G4double, G4double>  Compare_t;
    typedef std::map<G4int, Compare_t>                IDcompare_t;
    typedef std::map<G4String, IDcompare_t>           TypeCompare_t;

public:
    TSRunAction();
    virtual ~TSRunAction();

public:
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    virtual G4Run* GenerateRun();

private:
    TSDetectorConstruction* fDetector;
    G4String fName;
    TypeCompare_t fTypeCompare;

};

#endif
