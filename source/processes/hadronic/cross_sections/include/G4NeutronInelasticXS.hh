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
// $Id: G4NeutronInelasticXS.hh,v 1.3 2009/11/19 11:52:27 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4NeutronInelasticXS
//
// Author  Ivantchenko, Geant4, 3-AUG-09
//
// Modifications:
//
 

#ifndef G4NeutronInelasticXS_h
#define G4NeutronInelasticXS_h 1

#include "G4VCrossSectionDataSet.hh"
#include "globals.hh"

class G4DynamicParticle;
class G4ParticleDefinition;
class G4Element;
class G4PhysicsVector;
class G4GlauberGribovCrossSection;

class G4NeutronInelasticXS : public G4VCrossSectionDataSet
{
public: // With Description

  G4int Z;
  G4NeutronInelasticXS();

  virtual ~G4NeutronInelasticXS();

  // The following methods need to be implemented for each new data set.
  virtual
  G4bool IsApplicable(const G4DynamicParticle*, const G4Element*);

  virtual
  G4bool IsZAApplicable(const G4DynamicParticle*, 
			G4double /*Z*/, G4double /*A*/);

  virtual
  G4double GetCrossSection(const G4DynamicParticle*, 
			   const G4Element*, 
	 		   G4double aTemperature = 0.);

  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual
  void DumpPhysicsTable(const G4ParticleDefinition&);

public: // Without Description

  inline void SetVerboseLevel(G4int value)
  {
    verboseLevel = value;
  }
  inline G4int GetVerboseLevel()
  {
    return verboseLevel;
  }

private: // Without Description

  void Initialise(G4int Z, G4DynamicParticle* dp = 0, const char* = 0);

  G4NeutronInelasticXS & operator=(const G4NeutronInelasticXS &right);
  G4NeutronInelasticXS(const G4NeutronInelasticXS&);
  
  G4GlauberGribovCrossSection* ggXsection;

  G4PhysicsVector* data[93];
  G4double         coeff[93];

  G4bool  isInitialized;

};

#endif
