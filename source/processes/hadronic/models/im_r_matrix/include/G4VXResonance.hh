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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:     G4VXResonance
//
//      Author:        Maria Grazia Pia (MariaGrazia.Pia@genova.infn.it)
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
//      Kinetic Model, resonance cross section class
//
// -------------------------------------------------------------------

#ifndef G4VXRESONANCE_HH
#define G4VXRESONANCE_HH

#include "globals.hh"
#include "G4VCrossSectionSource.hh"
#include "G4Clebsch.hh"

class G4ParticleDefinition;
class G4KineticTrack;
class G4ParticleDefinition;

class G4VXResonance : public G4VCrossSectionSource
{

public:

  G4VXResonance();

  virtual ~G4VXResonance();

  G4bool operator==(const G4VXResonance &right) const;
  G4bool operator!=(const G4VXResonance &right) const;


protected:

  G4VXResonance(const G4VXResonance &right);
  G4VXResonance& operator=(const G4VXResonance &right);
  
  G4double DegeneracyFactor(const G4KineticTrack& trk1, 
			    const G4KineticTrack& trk2,
			    G4double iSpinOut1,
			    G4double iSpinOut2) const;
  
  G4double DetailedBalance(const G4KineticTrack& trk1, 
			   const G4KineticTrack& trk2, 
			   G4int isoOut1, G4int isoOut2,
			   G4double iSpinOut1, G4double iSpinOut2,
			   G4double mOut1, G4double mOut2) const;

  G4double IsospinCorrection(const G4KineticTrack& trk1, 
			     const G4KineticTrack& trk2,
			     G4int isoOut1, G4int isoOut2,
			     G4double iSpinOut1, G4double iSpinOut2) const;

private:  

  G4Clebsch clebsch;

};

#endif


















