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
//
// Created by ngoc hoang tran on 03/08/2023.
//

#ifndef G4ChemEquilibrium_hh
#define G4ChemEquilibrium_hh 1
#include "globals.hh"
#include "G4MoleculeTable.hh"
#include "G4UnitsTable.hh"

class G4DNAMolecularReactionData;
class G4ChemEquilibrium //for each Equilibrium process
{
public:
  using MolType    = const G4MolecularConfiguration*;
  using Reaction    = const G4DNAMolecularReactionData*;
  explicit G4ChemEquilibrium(const G4int& type, const G4double& time);
  ~G4ChemEquilibrium() = default;
  void Initialize();
  inline G4bool IsStatusChanged()
  {
    if(fStatus == fAddEquilibrium){
      return false;
    }else
    {
      fStatus = fAddEquilibrium;
      if(fVerbose > 0)
      {
        PrintInfo();
      }
      return true;
    }
  }

  inline void Reset()
  {
    fStatus = false;
    fAddEquilibrium = false;
    fEquilibriumTime = 0;
    fGlobalTime = 0;
  }

  inline void SetVerbose(const G4int& verbose)
  {
    fVerbose = verbose;
  }

  inline void SetGlobalTime(const G4double& time)
  {
    fGlobalTime = time;

    if(fGlobalTime - fEquilibriumTime > fEquilibriumDuration && fAddEquilibrium)
    {
      fAddEquilibrium = false;
      if(fVerbose) {
        G4cout << "SetEquilibrium : off " << fRectionType
               << "  fGlobalTime : " << G4BestUnit(fGlobalTime, "Time")
               << "  fEquilibriumTime8 : " << G4BestUnit(fEquilibriumTime, "Time")
               << " fAddEquilibrium : " << fAddEquilibrium << G4endl;
      }
    }
  }
  void SetEquilibrium(Reaction pReaction);

  inline G4bool GetEquilibriumStatus() const
  {
    return fAddEquilibrium;
  }

  void PrintInfo() const;
private:
  G4bool fStatus = false;
  G4bool fAddEquilibrium = false;
  G4double fEquilibriumTime = 0;
  const G4double fEquilibriumDuration = 0;
  G4int fRectionType = 0;
  MolType fReactant1 = nullptr;
  MolType fReactant2 = nullptr;
  MolType fReactantB1 = nullptr;
  MolType fReactantB2 = nullptr;
  G4double fGlobalTime = 0;
  G4int fVerbose = 1;
};

#endif  //
