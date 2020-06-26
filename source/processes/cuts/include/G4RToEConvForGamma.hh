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
// G4RToEConvForGamma
//
// Class description:
//
// This class is a Range to Energy Converter for gamma.

// Author: H.Kurashige, 05 October 2002 - First implementation
// --------------------------------------------------------------------
#ifndef G4RToEConvForGamma_hh
#define G4RToEConvForGamma_hh 1

#include <vector>

#include "globals.hh"
#include "G4ios.hh"
#include "G4VRangeToEnergyConverter.hh"

class G4RToEConvForGamma : public G4VRangeToEnergyConverter
{
  public:

    G4RToEConvForGamma();
      // Constructor

    virtual ~G4RToEConvForGamma();
      // Destructor

  protected:

    using G4CrossSectionTable = G4LossTable;

    virtual G4double ComputeLoss( G4double AtomicNumber,
                                  G4double KineticEnergy );
  
    virtual void BuildRangeVector( const G4Material* aMaterial,
                                   G4RangeVector* rangeVector );
      // The Range Table

    void BuildAbsorptionLengthVector( const G4Material* aMaterial,
                                      G4RangeVector* rangeVector );
 
    G4double ComputeCrossSection( G4double AtomicNumber,
                                  G4double KineticEnergy );

    G4double Z = -1.0;  
    G4double s200keV = 0.0, s1keV = 0.0;
    G4double tmin = 0.0, tlow = 0.0; 
    G4double smin = 0.0, slow = 0.0;
    G4double cmin = 0.0, clow = 0.0, chigh = 0.0;
};

// ------------------
// Inline methods
// ------------------

inline 
G4double G4RToEConvForGamma::ComputeLoss(G4double AtomicNumber,
                                         G4double KineticEnergy) 
{
  return ComputeCrossSection(AtomicNumber,KineticEnergy);
}

inline 
void G4RToEConvForGamma::BuildRangeVector(const G4Material* aMaterial,
                                          G4RangeVector* rangeVector)
{
  BuildAbsorptionLengthVector(aMaterial, rangeVector);
}

#endif
