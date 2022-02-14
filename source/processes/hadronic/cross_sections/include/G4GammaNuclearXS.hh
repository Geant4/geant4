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
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:    G4GammaNuclearXS
//
// Authors  V.Ivantchenko, Geant4, 20 October 2020
//          B.Kutsenko, BINP/NSU, 10 August 2021
//
// Modification 


// Class Description:
// This is a base class for gamma-nuclear cross section based on
// data files from IAEA Evaluated Photonuclear Data Library (IAEA/PD-2019) 
// https://www-nds.iaea.org/photonuclear/
// Database review - https://www.sciencedirect.com/science/article/pii/S0090375219300699?via%3Dihub
// Class Description - End

#ifndef G4GammaNuclearXS_h
#define G4GammaNuclearXS_h 1

#include "G4VCrossSectionDataSet.hh"
#include "globals.hh"
#include "G4ElementData.hh"
#include "G4PhysicsVector.hh"
#include "G4Threading.hh"
#include "G4IsotopeList.hh"
#include <vector>

class G4DynamicParticle;
class G4ParticleDefinition;
class G4Element;
class G4VComponentCrossSection;

class G4GammaNuclearXS final : public G4VCrossSectionDataSet
{
public: 

  explicit G4GammaNuclearXS();

  ~G4GammaNuclearXS() final;
    
  static const char* Default_Name() {return "GammaNuclearXS";}

  G4bool IsElementApplicable(const G4DynamicParticle*, 
			     G4int Z, const G4Material*) final;

  G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A,
			 const G4Element*, const G4Material* mat) final;

  G4double GetElementCrossSection(const G4DynamicParticle*, 
			          G4int Z, 
                                  const G4Material* mat = nullptr) final; 

  G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A,
                              const G4Isotope* iso = nullptr,
                              const G4Element* elm = nullptr,
                              const G4Material* mat = nullptr) final;

  const G4Isotope* SelectIsotope(const G4Element*, 
                                 G4double kinEnergy, G4double logE) final;

  void BuildPhysicsTable(const G4ParticleDefinition&) final;

  G4double IsoCrossSection(G4double ekin, G4int Z, G4int A);

  G4double ElementCrossSection(G4double ekin, G4int Z);

  void CrossSectionDescription(std::ostream&) const final;
      
  G4GammaNuclearXS & operator=(const G4GammaNuclearXS &right) = delete;
  G4GammaNuclearXS(const G4GammaNuclearXS&) = delete;
  
private: 

  void Initialise(G4int Z);

  void InitialiseOnFly(G4int Z);

  const G4String& FindDirectoryPath();

  inline G4PhysicsVector* GetPhysicsVector(G4int Z);

  G4PhysicsVector* RetrieveVector(std::ostringstream& in, G4bool isElement, G4int Z);
  
  G4VCrossSectionDataSet* ggXsection = nullptr;

  std::vector<G4double> temp;
  const G4ParticleDefinition* gamma;
  
  G4bool isMaster = false;

  static const G4int MAXZGAMMAXS = 95;
  static G4ElementData* data;
  // Upper limit of the linear transition between IAEA database and CHIPS model
  static const G4int rTransitionBound = 150.*CLHEP::MeV; 
  // The list of elements with non-linear parametrisation for better precision 
  const G4int freeVectorException[11] = {4, 6, 7, 8, 27, 39, 45, 65, 67, 69, 73}; 
  // CHIPS photonuclear model had a problem with high energy parametrisation 
  // for isotopes of H and He, coefficient is needed to connect isotope cross 
  // section with element cross-section on high energies.
  static G4double coeff[3][3]; 
  static G4double xs150[MAXZGAMMAXS];
  static G4String gDataDirectory;
  
#ifdef G4MULTITHREADED
  static G4Mutex gNuclearXSMutex;
#endif
};

inline
G4PhysicsVector* G4GammaNuclearXS::GetPhysicsVector(G4int Z)
{
  G4PhysicsVector* pv = data->GetElementData(Z);
  if(pv == nullptr) { 
    InitialiseOnFly(Z);
    pv = data->GetElementData(Z);
  }
  return pv;
}

#endif
