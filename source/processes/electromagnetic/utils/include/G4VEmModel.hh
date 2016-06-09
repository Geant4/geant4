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
// $Id: G4VEmModel.hh,v 1.35 2005/07/25 18:13:35 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01-patch-01 $
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
// 23-04-04 L.urban signature of SampleCosineTheta changed back 
// 17-11-04 Add method CrossSectionPerAtom (V.Ivanchenko)
// 14-03-05 Reduce number of pure virtual methods and make inline part separate (V.Ivanchenko)
// 24-03-05 Remove IsInCharge and add G4VParticleChange in the constructor (V.Ivanchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 15-04-05 optimize internal interface for msc (V.Ivanchenko)
// 08-05-05 A -> N (V.Ivanchenko)
// 25-07-05 Move constructor and destructor to the body (V.Ivanchenko)
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
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4DataVector.hh"
#include "G4VEmFluctuationModel.hh"
#include "Randomize.hh"

class G4PhysicsTable;
class G4Region;
class G4VParticleChange;

class G4VEmModel
{

public:

  G4VEmModel(const G4String& nam);

  virtual ~G4VEmModel();

  //------------------------------------------------------------------------
  // Virtual methods to be implemented for the concrete model
  //------------------------------------------------------------------------

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&) = 0;


  virtual std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin = 0.0,
                                      G4double tmax = DBL_MAX) = 0;

  //------------------------------------------------------------------------
  // Methods with standard implementation; may be overwritten if needed 
  //------------------------------------------------------------------------

  virtual G4double MinEnergyCut(const G4ParticleDefinition*,
                                const G4MaterialCutsCouple*);

  virtual G4double ComputeDEDX(const G4MaterialCutsCouple*,
			       const G4ParticleDefinition*,
			       G4double kineticEnergy,
			       G4double cutEnergy = DBL_MAX);

  virtual G4double CrossSection(const G4MaterialCutsCouple*,
				const G4ParticleDefinition*,
				G4double kineticEnergy,
				G4double cutEnergy = 0.0,
				G4double maxEnergy = DBL_MAX);

  virtual G4double ComputeDEDXPerVolume(
                               const G4Material*,
			       const G4ParticleDefinition*,
                                     G4double kineticEnergy,
                                     G4double cutEnergy = DBL_MAX);

  virtual G4double CrossSectionPerVolume(
                               const G4Material*,
			       const G4ParticleDefinition*,
                                     G4double kineticEnergy,
			             G4double cutEnergy = 0.0,
                                     G4double maxEnergy = DBL_MAX);

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double Z, 
                                      G4double A, 
                                      G4double cutEnergy = 0.0,
                                      G4double maxEnergy = DBL_MAX);

  // Methods for msc simulation
  virtual G4double GeomPathLength(G4PhysicsTable* theLambdaTable,
                            const G4MaterialCutsCouple*,
		            const G4ParticleDefinition*,
		                  G4double& kinEnergy,
			          G4double  lambda,
			          G4double  range,
    			          G4double  truePathLength);

  virtual G4double TrueStepLength(G4double geomStepLength);

  virtual void DefineForRegion(const G4Region*);

protected:

  virtual G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
				      G4double kineticEnergy);

  //------------------------------------------------------------------------
  // Generic methods common to all models
  //------------------------------------------------------------------------

public:

  void SetParticleChange(G4VParticleChange*, G4VEmFluctuationModel*);

  G4VEmFluctuationModel* GetModelOfFluctuations();

  G4double HighEnergyLimit();

  G4double LowEnergyLimit();

  void SetHighEnergyLimit(G4double);

  void SetLowEnergyLimit(G4double);

  G4double MaxSecondaryKinEnergy(const G4DynamicParticle* dynParticle);

  const G4Element* SelectRandomAtom(const G4Material*,
				    const G4ParticleDefinition*,
				    G4double kineticEnergy,
				    G4double cutEnergy = 0.0,
				    G4double maxEnergy = DBL_MAX);

  const G4String& GetName() const;

private:

  //  hide assignment operator
  G4VEmModel & operator=(const  G4VEmModel &right);
  G4VEmModel(const  G4VEmModel&);

  G4double        lowLimit;
  G4double        highLimit;
  G4double        xsec[40];

  G4VEmFluctuationModel* fluc;

  const G4String  name;

protected:

  G4VParticleChange*  pParticleChange;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmModel::HighEnergyLimit()
{
  return highLimit;
}

inline G4double G4VEmModel::LowEnergyLimit()
{
  return lowLimit;
}

inline void G4VEmModel::SetHighEnergyLimit(G4double val)
{
  highLimit = val;
}

inline void G4VEmModel::SetLowEnergyLimit(G4double val)
{
  lowLimit = val;
}

inline void G4VEmModel::SetParticleChange(G4VParticleChange* p,  
                                          G4VEmFluctuationModel* f = 0)
{
  if(p && pParticleChange != p) pParticleChange = p;
  fluc = f;
}

inline G4VEmFluctuationModel* G4VEmModel::GetModelOfFluctuations()
{
  return fluc;
}

inline G4double G4VEmModel::MinEnergyCut(const G4ParticleDefinition*,
                                         const G4MaterialCutsCouple*)
{
  return 0.0;
}

inline G4double G4VEmModel::ComputeDEDXPerVolume(
                                        const G4Material*,
					const G4ParticleDefinition*,
                                              G4double,
                                              G4double)
{
  return 0.0;
}

inline G4double G4VEmModel::ComputeDEDX(const G4MaterialCutsCouple* c,
                                        const G4ParticleDefinition* p,
                                              G4double kinEnergy,
                                              G4double cutEnergy)
{
  return ComputeDEDXPerVolume(c->GetMaterial(),p,kinEnergy,cutEnergy);
}

inline G4double G4VEmModel::CrossSection(const G4MaterialCutsCouple* c,
                                         const G4ParticleDefinition* p,
                                               G4double kinEnergy,
                                               G4double cutEnergy,
                                               G4double maxEnergy)
{
  return CrossSectionPerVolume(c->GetMaterial(),p,kinEnergy,cutEnergy,maxEnergy);
}

inline G4double G4VEmModel::ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double, G4double, G4double, G4double, G4double)
{
  return 0.0;
}

inline const G4Element* G4VEmModel::SelectRandomAtom(
				    const G4Material* material,
				    const G4ParticleDefinition* pd,
				    G4double kinEnergy,
				    G4double tcut,
				    G4double tmax)
{
  const G4ElementVector* theElementVector = material->GetElementVector();
  const G4Element* elm = (*theElementVector)[0];
  G4int nelm = material->GetNumberOfElements() - 1;
  if(nelm > 0) {
    G4double x = G4UniformRand()*CrossSectionPerVolume(material,pd,kinEnergy,tcut,tmax);
    G4int i = -1;
    do {i++;} while (x > xsec[i] && i < nelm);
    elm = (*theElementVector)[i];
  }
  return elm;
}

inline G4double G4VEmModel::MaxSecondaryKinEnergy(const G4DynamicParticle* dynParticle)
{
  return MaxSecondaryEnergy(dynParticle->GetDefinition(), dynParticle->GetKineticEnergy());
}

inline G4double G4VEmModel::MaxSecondaryEnergy(const G4ParticleDefinition*,
    				                     G4double kineticEnergy)
{
  return kineticEnergy;
}

inline const G4String& G4VEmModel::GetName() const 
{
  return name;
}

// Methods for msc simulation
inline G4double G4VEmModel::GeomPathLength(G4PhysicsTable*,
                                const G4MaterialCutsCouple*,
	  	                const G4ParticleDefinition*,
		                      G4double&,
			              G4double,
			              G4double,
    			              G4double truePathLength) 
{
  return truePathLength;
}

inline G4double G4VEmModel::TrueStepLength(G4double geomStepLength) 
{
  return geomStepLength;
}

inline void G4VEmModel::DefineForRegion(const G4Region*) 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

