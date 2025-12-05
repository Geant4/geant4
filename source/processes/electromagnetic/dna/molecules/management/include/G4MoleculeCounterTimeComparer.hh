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
#ifndef G4MOLECULECOUNTERTIMECOMPARER_HH
#define G4MOLECULECOUNTERTIMECOMPARER_HH 1

#include <CLHEP/Units/PhysicalConstants.h>
#include "G4Types.hh"

#include <map>
#include <vector>

class G4MoleculeCounterTimeComparer
{
  public:
    enum TimeComparerType
    {
      FixedPrecision,
      VariablePrecision
    };

  public:
    G4MoleculeCounterTimeComparer();
    virtual ~G4MoleculeCounterTimeComparer() = default;

    G4MoleculeCounterTimeComparer(const G4MoleculeCounterTimeComparer&);
    G4MoleculeCounterTimeComparer& operator=(const G4MoleculeCounterTimeComparer&);

    void SetFixedPrecision(G4double);
    void SetVariablePrecision(const std::vector<G4double>&, const std::vector<G4double>&);
    void SetVariablePrecision(const std::map<G4double, G4double>&);

    G4double GetPrecisionAtTime(G4double) const;

    bool operator()(const G4double&, const G4double&) const;

  private:
    TimeComparerType fType;
    G4double fPrecision{1 * CLHEP::picosecond};
    std::map<G4double, G4double> fVariablePrecision{};

  public:  // Factory
    static G4MoleculeCounterTimeComparer CreateWithFixedPrecision(G4double);
    static G4MoleculeCounterTimeComparer CreateWithVariablePrecision(const std::map<G4double, G4double>&);
};

#endif
