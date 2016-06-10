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
// $Id: G4ComponentSAIDTotalXS.hh 67988 2013-03-13 10:52:45Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4ComponentSAIDTotalXS
//
// Authors:  G.Folger, V.Ivanchenko, D.Wright
//
// Modifications:
//
 
//
// Class Description
// Total, elastic and inelastic hadron/nucleon cross sections
// from SAID database, G4SAIDXSDATA environment variable
// should be defined
// Class Description - End

#ifndef G4ComponentSAIDTotalXS_h
#define G4ComponentSAIDTotalXS_h 1

#include "G4VComponentCrossSection.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

enum G4SAIDCrossSectionType 
{ 
  saidUnknown = 0, 
  saidPP = 1, 
  saidNP = 2, 
  saidPIPP = 3,
  saidPINP = 4,
  saidPINP_PI0N = 5, 
  saidPINP_ETAN = 6,
  saidGP_PI0P = 7,
  saidGP_PIPN = 8,
  saidGN_PINP = 9,
  saidGN_PI0N = 10,
  saidGP_ETAP = 11,
  saidGP_ETAPP = 12,
  numberOfSaidXS = 13
};

class G4PhysicsVector;

class G4ComponentSAIDTotalXS : public G4VComponentCrossSection
{
public: //with description

  G4ComponentSAIDTotalXS();

  virtual ~G4ComponentSAIDTotalXS();

  virtual
  G4double GetTotalElementCrossSection(const G4ParticleDefinition*,
				       G4double kinEnergy, 
				       G4int /*Z*/, G4double /*N*/);

  virtual
  G4double GetTotalIsotopeCrossSection(const G4ParticleDefinition*,
				       G4double kinEnergy,
				       G4int /*Z*/, G4int /*N*/);

  virtual
  G4double GetInelasticElementCrossSection(const G4ParticleDefinition*,
					   G4double kinEnergy, 
					   G4int /*Z*/, G4double /*N*/);

  virtual
  G4double GetInelasticIsotopeCrossSection(const G4ParticleDefinition*,
					   G4double kinEnergy, 
					   G4int /*Z*/, G4int /*N*/);

  virtual
  G4double GetElasticElementCrossSection(const G4ParticleDefinition*,
					 G4double kinEnergy, 
					 G4int /*Z*/, G4double /*N*/);

  virtual
  G4double GetElasticIsotopeCrossSection(const G4ParticleDefinition*,
					 G4double kinEnergy, 
					 G4int /*Z*/, G4int /*N*/);

  G4double GetChargeExchangeCrossSection(const G4ParticleDefinition* prim,
					 const G4ParticleDefinition* sec,
					 G4double kinEnergy, 
					 G4int /*Z*/, G4int /*N*/);

  virtual
  void Description() const;

private:

  G4SAIDCrossSectionType GetType(const G4ParticleDefinition* prim,
				 const G4ParticleDefinition* sec,
				 G4int Z, G4int N);

  void Initialise(G4SAIDCrossSectionType tp);

  void ReadData(G4int index, G4PhysicsVector*,
		const G4String&, const G4String&);

  void PrintWarning(const G4ParticleDefinition* prim,
		    const G4ParticleDefinition* sec,
		    G4int /*Z*/, G4int /*N*/,
		    const G4String&, const G4String&);

  G4ComponentSAIDTotalXS & operator=(const G4ComponentSAIDTotalXS &right);
  G4ComponentSAIDTotalXS(const G4ComponentSAIDTotalXS&);

  static const G4String fnames[numberOfSaidXS];
  G4PhysicsVector* elastdata[numberOfSaidXS];
  G4PhysicsVector* inelastdata[numberOfSaidXS];

};

#endif
