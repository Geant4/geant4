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
// $Id: G4NeutronElasticXS.hh 93682 2015-10-28 10:09:49Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4NeutronElasticXS
//
// Author  Ivantchenko, Geant4, 3-AUG-09
//
// Modifications:
//
 
// Class Description:
// This is a base class for neutron elastic hadronic cross section based on
// data files from G4NEUTRONXSDATA data set 
// Class Description - End

#ifndef G4NeutronElasticXS_h
#define G4NeutronElasticXS_h 1

#include "G4VCrossSectionDataSet.hh"
#include "globals.hh"
#include <vector>

const G4int MAXZEL = 93;

class G4DynamicParticle;
class G4ParticleDefinition;
class G4Element;
class G4PhysicsVector;
class G4ComponentGGHadronNucleusXsc;
class G4HadronNucleonXsc;

class G4NeutronElasticXS : public G4VCrossSectionDataSet
{
public: // With Description

  G4NeutronElasticXS();

  virtual ~G4NeutronElasticXS();
    
  static const char* Default_Name() {return "G4NeutronElasticXS";}

  virtual
  G4bool IsElementApplicable(const G4DynamicParticle*, 
			     G4int Z, const G4Material*);

  virtual
  G4double GetElementCrossSection(const G4DynamicParticle*, 
				  G4int Z, const G4Material* mat=0); 

  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual void CrossSectionDescription(std::ostream&) const;

private: 

  void Initialise(G4int Z, G4DynamicParticle* dp = 0, const char* = 0);

  G4NeutronElasticXS & operator=(const G4NeutronElasticXS &right);
  G4NeutronElasticXS(const G4NeutronElasticXS&);
  
  G4ComponentGGHadronNucleusXsc* ggXsection;
  G4HadronNucleonXsc* fNucleon;

  const G4ParticleDefinition* proton;

  static std::vector<G4PhysicsVector*>* data;
  static G4double  coeff[MAXZEL];

  G4bool  isMaster;

};

#endif
