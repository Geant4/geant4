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


















