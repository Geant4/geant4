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
// $Id: G4NeutronInelasticXS.hh 93682 2015-10-28 10:09:49Z gcosmo $
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

// Class Description:
// This is a base class for neutron inelastic hadronic cross section based on
// data files from G4NEUTRONXSDATA data set 
// Class Description - End
 
#ifndef G4NeutronInelasticXS_h
#define G4NeutronInelasticXS_h 1

#include "G4VCrossSectionDataSet.hh"
#include "globals.hh"
#include "G4ElementData.hh"
#include <vector>

const G4int MAXZINEL = 93;

class G4DynamicParticle;
class G4ParticleDefinition;
class G4Element;
class G4PhysicsVector;
class G4ComponentGGHadronNucleusXsc;
class G4HadronNucleonXsc;

class G4NeutronInelasticXS : public G4VCrossSectionDataSet
{
public: 

  G4NeutronInelasticXS();

  virtual ~G4NeutronInelasticXS();

  static const char* Default_Name() {return "G4NeutronInelasticXS";}
    
  virtual
  G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z,
			     const G4Material*);

  virtual
  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,
			 const G4Element*, const G4Material*);

  virtual
  G4double GetElementCrossSection(const G4DynamicParticle*, 
				  G4int Z, const G4Material* mat=0);

  virtual
  G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
                              const G4Isotope* iso,
                              const G4Element* elm,
                              const G4Material* mat);

  virtual G4Isotope* SelectIsotope(const G4Element*, G4double kinEnergy);

  virtual
  void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual void CrossSectionDescription(std::ostream&) const;

private: 

  void Initialise(G4int Z, G4DynamicParticle* dp = 0, const char* = 0);

  G4PhysicsVector* RetrieveVector(std::ostringstream& in, G4bool warn);

  G4double IsoCrossSection(G4double ekin, G4int Z, G4int A);

  G4NeutronInelasticXS & operator=(const G4NeutronInelasticXS &right);
  G4NeutronInelasticXS(const G4NeutronInelasticXS&);
  
  G4ComponentGGHadronNucleusXsc* ggXsection;
  G4HadronNucleonXsc* fNucleon;

  const G4ParticleDefinition* proton;

  G4bool  isMaster;

  static G4ElementData* data;
  std::vector<G4double> temp;

  static G4double  coeff[MAXZINEL];

  static const G4int amin[MAXZINEL];
  static const G4int amax[MAXZINEL];
};

#endif
