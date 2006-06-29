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

#ifndef G4XRESONANCE_HH
#define G4XRESONANCE_HH

#include "globals.hh"
#include "G4PhysicsVector.hh"
#include "G4VCrossSectionSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4VXResonanceTable.hh"
#include "G4VXResonance.hh"

class G4KineticTrack;

class G4XResonance : public G4VXResonance

{

public:

  G4XResonance(const G4ParticleDefinition* in1, 
	       const G4ParticleDefinition* in2,
	       G4int iIsospinOut1, G4double iSpinOut1, G4double massOut1,
	       G4int iIsospinOut2, G4double iSpinOut2, G4double massOut2,
	       G4String subType1, G4String subType2,
	       const G4VXResonanceTable& sigmaTable);

  virtual ~G4XResonance();

  G4bool operator==(const G4XResonance &right) const;
  G4bool operator!=(const G4XResonance &right) const;

  virtual G4double CrossSection(const G4KineticTrack& trk1, const G4KineticTrack& trk2) const;

  virtual const G4CrossSectionVector* GetComponents() const { return 0; }

  virtual G4String Name() const;


protected:


private:  

  G4XResonance(const G4XResonance &right);
  G4XResonance& operator=(const G4XResonance &right);
  
  G4int isoOut1;
  G4double iSpinOut1;
  G4double mOut1;
  
  G4int isoOut2;
  G4double iSpinOut2;
  G4double mOut2;

  // Owned pointer
  G4PhysicsVector* table;

  G4String name;

};

#endif


















