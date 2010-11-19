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
// $Id: G4CrossSectionElastic.hh,v 1.2 2010-11-19 11:12:11 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4CrossSectionElastic
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 19.11.2010
// Modifications:
//
// Class Description:
//
// Wrapper for elastic cross section build from a component
//
// -------------------------------------------------------------------
//

#ifndef G4CrossSectionElastic_h
#define G4CrossSectionElastic_h 1

#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"

class G4ParticleDefinition;
class G4VComponentCrossSection;

class G4CrossSectionElastic : public G4VCrossSectionDataSet
{
public:

  G4CrossSectionElastic(const G4ParticleDefinition*, G4VComponentCrossSection*,
			G4int zmin = 0, G4int zmax = 256, 
			G4double Emin = 0.0, G4double Emax = DBL_MAX);

  virtual ~G4CrossSectionElastic();
   
  virtual
  G4bool IsApplicable(const G4DynamicParticle*, const G4Element*);

  virtual
  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A);

  virtual
  G4double GetCrossSection(const G4DynamicParticle*, 
			   const G4Element*, G4double aTemperature = 0.);

  virtual
  G4double GetZandACrossSection(const G4DynamicParticle*, G4int /*Z*/,
                                G4int /*A*/, G4double aTemperature = 0.);

  inline const G4ParticleDefinition* GetParticle() const { return particle; }

private:

  G4CrossSectionElastic & operator=(const G4CrossSectionElastic &right);
  G4CrossSectionElastic(const G4CrossSectionElastic&);

  const G4ParticleDefinition* particle;
  G4VComponentCrossSection* component;
  G4int Zmin;
  G4int Zmax;

};

#endif
