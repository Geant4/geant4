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
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//
// Class Description
//  This class is a Range to Energy Converter for gamma.
//
// ------------------------------------------------------------
//   First Implementation          5 Oct. 2002  H.Kurahige
// ------------------------------------------------------------

#ifndef G4RToEConvForGamma_h
#define G4RToEConvForGamma_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "G4VRangeToEnergyConverter.hh"


class G4RToEConvForGamma : public G4VRangeToEnergyConverter
{
  public: // with description
  //  constructor
  G4RToEConvForGamma();

  public:
  //  destructor
  virtual ~G4RToEConvForGamma();

  public: // with description 
  // calculate energy cut from given range cut for the material
  // virtual G4double convert(G4double rangeCut, const G4Material* material); 

  protected:
    virtual G4double ComputeLoss( G4double AtomicNumber,
                                  G4double KineticEnergy
				  ) ;
  
  //-------------- Range Table ------------------------------------------
    virtual void BuildRangeVector( const G4Material* aMaterial,
				   G4RangeVector* rangeVector);

    typedef G4LossTable G4CrossSectionTable;
    void BuildAbsorptionLengthVector( const G4Material* aMaterial,
				      G4RangeVector* rangeVector);
 
    G4double ComputeCrossSection( G4double AtomicNumber,
				  G4double KineticEnergy
				  );

    G4double Z;  
    G4double s200keV, s1keV;
    G4double tmin, tlow; 
    G4double smin, slow;
    G4double cmin, clow, chigh;

};

inline 
 G4double G4RToEConvForGamma::ComputeLoss(G4double AtomicNumber,
					  G4double KineticEnergy) 
{
  return ComputeCrossSection(AtomicNumber,KineticEnergy);
}

inline 
 void G4RToEConvForGamma::BuildRangeVector(
                                const G4Material* aMaterial,
                                G4RangeVector* rangeVector)
{
  BuildAbsorptionLengthVector(aMaterial, rangeVector);
}



#endif









