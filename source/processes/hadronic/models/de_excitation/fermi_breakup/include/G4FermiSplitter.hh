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
// G4FermiBreakUpAN alternative de-excitation model
// by A. Novikov (January 2025)
//

#ifndef G4FERMISPLITTER_HH
#define G4FERMISPLITTER_HH

#include "G4FermiDataTypes.hh"
#include "G4FermiVFragment.hh"

class G4FermiSplitter
{
  public:
    static G4double DecayWeight(const G4FermiFragmentVector& split, G4FermiAtomicMass atomicMass,
                                G4double totalEnergy);

    static G4double SplitFactor(const G4FermiFragmentVector& split, G4FermiAtomicMass atomicMass);

    static G4double KineticFactor(const G4FermiFragmentVector& split, G4double totalEnergy);

    static void GenerateSplits(G4FermiNucleiData nucleiData,
                               std::vector<G4FermiFragmentVector>& splits);

    static std::vector<G4FermiFragmentVector> GenerateSplits(G4FermiNucleiData nucleiData);
};

#endif  // G4FERMISPLITTER_HH
