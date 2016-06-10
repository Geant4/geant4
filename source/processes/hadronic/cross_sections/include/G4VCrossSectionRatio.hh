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
// $Id: G4VCrossSectionRatio.hh 66241 2012-12-13 18:34:42Z gunter $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4VCrossSectionRatio
//
// Author:  V.Ivanchenko 30.10.2014
//
// Modifications:
//
 
//
// Class Description
// This is a base class for hadronic cross section ratio
// quasi-elastic/inelastic or diffractive/elastic
// Class Description - End

#ifndef G4VCrossSectionRatio_h
#define G4VCrossSectionRatio_h 1

#include "globals.hh"

class G4ParticleDefinition;

class G4VCrossSectionRatio
{
public: 

  G4VCrossSectionRatio(const G4String& nam = "", G4int verb = 0);

  virtual ~G4VCrossSectionRatio();

  virtual
  G4double ComputeRatio(const G4ParticleDefinition*,
			G4double kinEnergy, 
			G4int /*Z*/, G4int /*A*/) = 0;

  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual
  void DumpPhysicsTable(const G4ParticleDefinition&);

  virtual 
  void Description() const;

  inline void SetVerboseLevel(G4int value);

  inline G4int GetVerboseLevel() const;

  inline const G4String& GetName() const;

private:

  G4VCrossSectionRatio & operator=(const G4VCrossSectionRatio &right);
  G4VCrossSectionRatio(const G4VCrossSectionRatio&);

  G4int verboseLevel;
  const G4String name;
};

inline void G4VCrossSectionRatio::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline G4int G4VCrossSectionRatio::GetVerboseLevel() const
{
  return verboseLevel;
}

inline const G4String& G4VCrossSectionRatio::GetName() const
{
  return name;
}

#endif
