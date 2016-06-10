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
// $Id: G4VComponentCrossSection.hh 66241 2012-12-13 18:34:42Z gunter $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4VComponentCrossSection
//
// Authors:  G.Folger, V.Ivanchenko, D.Wright
//
// Modifications:
//
 
//
// Class Description
// This is a base class for hadronic cross section data source
// Class Description - End

#ifndef G4VComponentCrossSection_h
#define G4VComponentCrossSection_h 1

#include "G4ParticleDefinition.hh"
#include "G4Element.hh"
#include "globals.hh"

class G4VComponentCrossSection
{
public: //with description

  G4VComponentCrossSection(const G4String& nam = "");

  virtual ~G4VComponentCrossSection();

  inline
  G4double GetTotalElementCrossSection(const G4ParticleDefinition*,
				       G4double kinEnergy,
				       const G4Element*);

  virtual
  G4double GetTotalElementCrossSection(const G4ParticleDefinition*,
				       G4double kinEnergy, 
				       G4int /*Z*/, G4double /*N*/) = 0;

  virtual
  G4double GetTotalIsotopeCrossSection(const G4ParticleDefinition*,
				       G4double kinEnergy,
				       G4int /*Z*/, G4int /*N*/) = 0;

  inline
  G4double GetInelasticElementCrossSection(const G4ParticleDefinition*,
					   G4double kinEnergy, 
					   const G4Element*);

  virtual
  G4double GetInelasticElementCrossSection(const G4ParticleDefinition*,
					   G4double kinEnergy, 
					   G4int /*Z*/, G4double /*N*/) = 0;

  virtual
  G4double GetInelasticIsotopeCrossSection(const G4ParticleDefinition*,
					   G4double kinEnergy, 
					   G4int /*Z*/, G4int /*N*/) = 0;

  inline
  G4double GetElasticElementCrossSection(const G4ParticleDefinition*,
					 G4double kinEnergy, 
					 const G4Element*);

  virtual
  G4double GetElasticElementCrossSection(const G4ParticleDefinition*,
					 G4double kinEnergy, 
					 G4int /*Z*/, G4double /*N*/) = 0;

  virtual
  G4double GetElasticIsotopeCrossSection(const G4ParticleDefinition*,
					 G4double kinEnergy, 
					 G4int /*Z*/, G4int /*N*/) = 0;

  virtual
  G4double ComputeQuasiElasticRatio(const G4ParticleDefinition*,
				    G4double kinEnergy, 
				    G4int /*Z*/, G4int /*N*/);

  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual
  void DumpPhysicsTable(const G4ParticleDefinition&);

  virtual
  void Description() const;

  inline void SetVerboseLevel(G4int value);

  inline G4int GetVerboseLevel() const;

  inline G4double GetMinKinEnergy() const;

  inline void SetMinKinEnergy(G4double value);

  inline G4double GetMaxKinEnergy() const;

  inline void SetMaxKinEnergy(G4double value);

  inline const G4String& GetName() const;

private:

  G4VComponentCrossSection & operator=(const G4VComponentCrossSection &right);
  G4VComponentCrossSection(const G4VComponentCrossSection&);

  G4int verboseLevel;

  G4double minKinEnergy;
  G4double maxKinEnergy;

  const G4String name;
};

inline G4double 
G4VComponentCrossSection::GetTotalElementCrossSection(
         const G4ParticleDefinition* p,
	 G4double kinEnergy, 
	 const G4Element* elm)
{
  return GetTotalElementCrossSection(p,kinEnergy,
				     (G4int)elm->GetZ(),elm->GetN());
}

inline G4double 
G4VComponentCrossSection::GetInelasticElementCrossSection(
         const G4ParticleDefinition* p,
	 G4double kinEnergy, 
	 const G4Element* elm)
{
  return GetInelasticElementCrossSection(p,kinEnergy,
					 (G4int)elm->GetZ(),elm->GetN());
}

inline G4double 
G4VComponentCrossSection::GetElasticElementCrossSection(
         const G4ParticleDefinition* p,
	 G4double kinEnergy, 
	 const G4Element* elm)
{
  return GetElasticElementCrossSection(p,kinEnergy,
				       (G4int)elm->GetZ(),elm->GetN());
}

inline void G4VComponentCrossSection::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline G4int G4VComponentCrossSection::GetVerboseLevel() const
{
  return verboseLevel;
}

inline void G4VComponentCrossSection::SetMinKinEnergy(G4double value)
{
  minKinEnergy = value;
}

inline G4double G4VComponentCrossSection::GetMinKinEnergy() const
{
  return minKinEnergy;
}

inline void G4VComponentCrossSection::SetMaxKinEnergy(G4double value)
{
  maxKinEnergy = value;
}

inline G4double G4VComponentCrossSection::GetMaxKinEnergy() const
{
  return maxKinEnergy;
}

inline const G4String& G4VComponentCrossSection::GetName() const
{
  return name;
}

#endif
