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
// $Id: G4CrossSectionElastic.hh 66241 2012-12-13 18:34:42Z gunter $
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
#include <iostream>

class G4ParticleDefinition;
class G4VComponentCrossSection;
class G4DynamicParticle;
class G4Element;
class G4Material;
class G4NistManager;

class G4CrossSectionElastic : public G4VCrossSectionDataSet
{
public:

  G4CrossSectionElastic(G4VComponentCrossSection*,
			G4int zmin = 0, G4int zmax = 256, 
			G4double Emin = 0.0, G4double Emax = DBL_MAX);

  virtual ~G4CrossSectionElastic();
   
  virtual
  G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z,
			     const G4Material* mat = 0);

  virtual
  G4double GetElementCrossSection(const G4DynamicParticle*, 
				  G4int Z, 
				  const G4Material* mat = 0);

  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual
  void DumpPhysicsTable(const G4ParticleDefinition&);

  virtual void CrossSectionDescription(std::ostream&) const;

private:

  G4CrossSectionElastic & operator=(const G4CrossSectionElastic &right);
  G4CrossSectionElastic(const G4CrossSectionElastic&);

  G4NistManager* nist;
  G4VComponentCrossSection* component;
  G4int Zmin;
  G4int Zmax;

};

#endif
