//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4RToEConvForGamma.hh,v 1.2 2002-12-16 11:15:43 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "g4std/vector"

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
				  ) const;
  
  //-------------- Range Table ------------------------------------------
    virtual void BuildRangeVector( const G4Material* aMaterial,
				   G4double       maxEnergy,
				   G4double       aMass,
				   G4RangeVector* rangeVector);

    typedef G4LossTable G4CrossSectionTable;
    void BuildAbsorptionLengthVector( const G4Material* aMaterial,
				      G4double       maxEnergy,
				      G4double       aMass,
				      G4RangeVector* rangeVector);
 
    G4double ComputeCrossSection( G4double AtomicNumber,
				  G4double KineticEnergy
				  ) const;


};

inline 
 G4double G4RToEConvForGamma::ComputeLoss(G4double AtomicNumber,
					  G4double KineticEnergy) const
{
  return ComputeCrossSection(AtomicNumber,KineticEnergy);
}

inline 
 void G4RToEConvForGamma::BuildRangeVector(
                                const G4Material* aMaterial,
                                G4double       maxEnergy,     
                                G4double       aMass,
                                G4RangeVector* rangeVector)
{
  BuildAbsorptionLengthVector(aMaterial, maxEnergy, aMass, rangeVector);
}



#endif









