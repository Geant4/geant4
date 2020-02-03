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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4NeutronInelasticXS
//
// Author  Ivantchenko, Geant4, 3-AUG-09
//

// Class Description:
// This is a base class for neutron inelastic hadronic cross section based on
// data files from G4PARTICLEXSDATA data set 
// Class Description - End
 
#ifndef G4NeutronInelasticXS_h
#define G4NeutronInelasticXS_h 1

#include "G4VCrossSectionDataSet.hh"
#include "globals.hh"
#include "G4ElementData.hh"
#include "G4Threading.hh"
#include <vector>

const G4int MAXZINEL = 93;

class G4DynamicParticle;
class G4ParticleDefinition;
class G4Element;
class G4PhysicsVector;
class G4ComponentGGHadronNucleusXsc;
class G4NistManager;

class G4NeutronInelasticXS : public G4VCrossSectionDataSet
{
public: 

  explicit G4NeutronInelasticXS();

  ~G4NeutronInelasticXS() final;

  static const char* Default_Name() {return "G4NeutronInelasticXS";}
    
  G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z,
			     const G4Material*) final;

  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,
			 const G4Element*, const G4Material*) final;

  G4double GetElementCrossSection(const G4DynamicParticle*, 
				  G4int Z, const G4Material*) final;

  G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
                              const G4Isotope* iso,
                              const G4Element* elm,
                              const G4Material* mat) final;

  const G4Isotope* SelectIsotope(const G4Element*, 
                                 G4double kinEnergy, G4double logE) final;

  void BuildPhysicsTable(const G4ParticleDefinition&) final;

  void CrossSectionDescription(std::ostream&) const final;

private: 

  void Initialise(G4int Z);

  void InitialiseOnFly(G4int Z);

  const G4String& FindDirectoryPath();

  const G4PhysicsVector* GetPhysicsVector(G4int Z);

  G4PhysicsVector* RetrieveVector(std::ostringstream& in, G4bool warn);

  G4double IsoCrossSection(G4double ekin, G4double logekin, G4int Z, G4int A);

  G4NeutronInelasticXS & operator=(const G4NeutronInelasticXS &right);
  G4NeutronInelasticXS(const G4NeutronInelasticXS&);
  
  G4ComponentGGHadronNucleusXsc* ggXsection;
  G4NistManager* nist;

  const G4ParticleDefinition* neutron;

  std::vector<G4double> temp;

  G4bool  isMaster;

  static G4ElementData* data;
  static G4double   coeff[MAXZINEL];
  static G4double    aeff[MAXZINEL];
  static const G4int amin[MAXZINEL];
  static const G4int amax[MAXZINEL];
  static G4String gDataDirectory;

#ifdef G4MULTITHREADED
  static G4Mutex neutronInelasticXSMutex;
#endif
};

#endif
