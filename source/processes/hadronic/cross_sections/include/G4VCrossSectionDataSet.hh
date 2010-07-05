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
// $Id: G4VCrossSectionDataSet.hh,v 1.14 2010-07-05 13:39:11 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4VCrossSectionDataSet
//
// Author  F.W. Jones, TRIUMF, 20-JAN-97
//
// Modifications:
// 23.01.2009 V.Ivanchenko move constructor and destructor to source
// 05.07.2010 V.Ivanchenko added name, min and max energy limit and
//            corresponding access methods
//
 
//
// Class Description
// This is a base class for hadronic cross section data sets.  Users may 
// derive specialized cross section classes and register them with the
// appropriate process, or use provided data sets.
// Class Description - End

#ifndef G4VCrossSectionDataSet_h
#define G4VCrossSectionDataSet_h 1

#include "G4DynamicParticle.hh"
#include "G4Element.hh"


class G4VCrossSectionDataSet
{
public: //with description

  G4VCrossSectionDataSet(const G4String& nam = "");

  virtual ~G4VCrossSectionDataSet();

  // The following methods need to be implemented for each new data set.
  virtual
  G4bool IsApplicable(const G4DynamicParticle*, const G4Element*) = 0;

  virtual
  G4bool IsZAApplicable(const G4DynamicParticle*, G4double /*Z*/, G4double /*A*/);

  virtual
  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int /*Z*/, G4int /*N*/);

  virtual
  G4double GetCrossSection(const G4DynamicParticle*, 
			   const G4Element*, 
			   G4double aTemperature = 0.) = 0;

  virtual
  G4double GetIsoCrossSection(const G4DynamicParticle*, const G4Isotope*, 
			      G4double aTemperature = 0.);

  virtual
  G4double GetIsoZACrossSection(const G4DynamicParticle*, G4double /*Z*/,
				G4double /*A*/, G4double aTemperature = 0.);

  virtual
  G4double GetZandACrossSection(const G4DynamicParticle*, G4int /*Z*/,
				G4int /*A*/, G4double aTemperature = 0.);

  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&) = 0;

  virtual
  void DumpPhysicsTable(const G4ParticleDefinition&) = 0;

public: // Without Description

  inline void SetVerboseLevel(G4int value);

  inline G4double GetMinKinEnergy() const;

  inline void SetMinKinEnergy(G4double value);

  inline G4double GetMaxKinEnergy() const;

  inline void SetMaxKinEnergy(G4double value);

  inline const G4String& GetName() const;

protected:

  G4int verboseLevel;

private:

  G4VCrossSectionDataSet & operator=(const G4VCrossSectionDataSet &right);
  G4VCrossSectionDataSet(const G4VCrossSectionDataSet&);

  G4double minKinEnergy;
  G4double maxKinEnergy;

  const G4String name;
};

inline void G4VCrossSectionDataSet::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline void G4VCrossSectionDataSet::SetMinKinEnergy(G4double value)
{
  minKinEnergy = value;
}

inline G4double G4VCrossSectionDataSet::GetMinKinEnergy() const
{
  return minKinEnergy;
}

inline void G4VCrossSectionDataSet::SetMaxKinEnergy(G4double value)
{
  maxKinEnergy = value;
}

inline G4double G4VCrossSectionDataSet::GetMaxKinEnergy() const
{
  return maxKinEnergy;
}

inline const G4String& G4VCrossSectionDataSet::GetName() const
{
  return name;
}

#endif
