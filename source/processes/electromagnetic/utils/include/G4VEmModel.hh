//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4VEmModel.hh,v 1.17 2004/03/01 14:16:57 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4VEmModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.01.2002
//
// Modifications:
//
// 23-12-02 V.Ivanchenko change interface before move to cut per region
// 24-01-03 Cut per region (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 25-02-03 Add sample theta and displacement (V.Ivanchenko)
// 23-07-03 Replace G4Material by G4MaterialCutCouple in dE/dx and CrossSection
//          calculation (V.Ivanchenko)
// 01-03-04 L.Urban signature changed in SampleCosineTheta 
//

//
// Class Description:
//
// Abstract interface to energy loss models

// -------------------------------------------------------------------
//

#ifndef G4VEmModel_h
#define G4VEmModel_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4DataVector.hh"

class G4PhysicsTable;
class G4Region;

class G4VEmModel
{

public:

  G4VEmModel(const G4String& nam): name(nam) {};

  virtual ~G4VEmModel() {};

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&) = 0;

  virtual G4double HighEnergyLimit(const G4ParticleDefinition*) = 0;

  virtual G4double LowEnergyLimit(const G4ParticleDefinition*) = 0;

  virtual void SetHighEnergyLimit(G4double) = 0;

  virtual void SetLowEnergyLimit(G4double) = 0;

  virtual G4double MinEnergyCut(const G4ParticleDefinition*,
                                const G4MaterialCutsCouple*) = 0;

  virtual G4bool IsInCharge(const G4ParticleDefinition*) = 0;

  virtual G4double ComputeDEDX(const G4MaterialCutsCouple*,
                               const G4ParticleDefinition*,
                                     G4double kineticEnergy,
                                     G4double cutEnergy) = 0;

  virtual G4double CrossSection(const G4MaterialCutsCouple*,
                                const G4ParticleDefinition*,
                                      G4double kineticEnergy,
                                      G4double cutEnergy,
                                      G4double maxEnergy) = 0;

  virtual G4DynamicParticle* SampleSecondary(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double tmax) = 0;

  virtual std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double tmax) = 0;

  virtual G4double MaxSecondaryEnergy(
				const G4DynamicParticle* dynParticle) = 0;

  const G4String& GetName() const {return name;};

  // Methods for msc simulation
  virtual G4double GeomPathLength(G4PhysicsTable*,
                            const G4MaterialCutsCouple*,
		            const G4ParticleDefinition*,
		                  G4double&,
			          G4double,
			          G4double,
    			          G4double truePathLength) {return truePathLength;};
  // G4double parameters: kinEnergy, lambda, range,
  // G4PhysicsTable: theLambdaTable

  virtual G4double TrueStepLength(G4double geomStepLength) {return geomStepLength;};

  virtual G4double SampleCosineTheta(G4double,G4double,G4double ) {return 1.0;};
  //  trueStepLength + Tkin + lambda   

  virtual G4double SampleDisplacement() {return 0.0;};

  virtual void DefineForRegion(const G4Region*) {};

protected:

  virtual G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
    				            G4double kineticEnergy) = 0;

private:

  //  hide assignment operator
  G4VEmModel & operator=(const  G4VEmModel &right);
  G4VEmModel(const  G4VEmModel&);

  const G4String  name;
};

#endif

