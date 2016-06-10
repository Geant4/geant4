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
// $Id: G4LivermorePolarizedRayleighModel.hh 93810 2015-11-02 11:27:56Z gcosmo $
//
// Author: Sebastien Incerti
//         30 October 2008
//         on base of G4LowEnergyPolarizedRayleigh developed by R. Capra
//

#ifndef G4LivermorePolarizedRayleighModel_h
#define G4LivermorePolarizedRayleighModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4ProductionCutsTable.hh"

//#include "G4CrossSectionHandler.hh"
//#include "G4LogLogInterpolation.hh"
//#include "G4CompositeEMDataSet.hh"

#include "G4VEMDataSet.hh"
#include "G4CompositeEMDataSet.hh"



class G4LivermorePolarizedRayleighModel : public G4VEmModel
{

public:

  G4LivermorePolarizedRayleighModel(const G4ParticleDefinition* p = 0, 
		     const G4String& nam = "LivermorePolarizedRayleigh");

  virtual ~G4LivermorePolarizedRayleighModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual void InitialiseLocal(const G4ParticleDefinition*, 
			       G4VEmModel* masterModel);
  
  virtual void InitialiseForElement(const G4ParticleDefinition*, G4int Z);
  
  virtual G4double ComputeCrossSectionPerAtom(
					      const G4ParticleDefinition*,
					      G4double kinEnergy, 
					      G4double Z, 
					      G4double A=0, 
					      G4double cut=0,
					      G4double emax=DBL_MAX);
  
  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

private:

  G4ParticleChangeForGamma* fParticleChange;
  G4int verboseLevel;
  G4bool isInitialised;

  G4double lowEnergyLimit;  
  //G4double highEnergyLimit; 
  
  G4LivermorePolarizedRayleighModel & operator=(const  G4LivermorePolarizedRayleighModel &right);
  G4LivermorePolarizedRayleighModel(const  G4LivermorePolarizedRayleighModel&);

  // cross section 

  void ReadData(size_t Z, const char* path = 0);
  
  static G4VEMDataSet* formFactorData;
  static G4int maxZ;
  static G4LPhysicsFreeVector* dataCS[101];

  // Polarization 
  
  //   Generates \f$cos \left ( \theta\right )\f$ of the scattered photon
  //   incomingPhotonEnergy The energy of the incoming photon
  //   zAtom Atomic number
  //   \f$cos \left ( \theta\right )\f$
  G4double GenerateCosTheta(G4double incomingPhotonEnergy, G4int zAtom) const;
   
  //   Generates \f$\phi\f$ of the scattered photon
  //   cosTheta \f$cos \left ( \theta\right )\f$ of the scattered photon
  //    \f$\phi\f$
  G4double GeneratePhi(G4double cosTheta) const;

  //   Generates the polarization direction \f$\beta\f$ in the plane x, y relative to the x direction
  //   \f$\beta\f$
  G4double GeneratePolarizationAngle(void) const;
 
  G4ThreeVector GetPhotonPolarization(const G4DynamicParticle&  photon);


};


#endif
