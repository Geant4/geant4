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
// $Id: G4VEmModel.hh,v 1.26 2005-03-18 12:48:31 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4DataVector.hh"

class G4PhysicsTable;
class G4Region;

class G4VEmModel
{

public:

  G4VEmModel(const G4String& nam): name(nam) {};

  virtual ~G4VEmModel() {};

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&) = 0;

  virtual G4double HighEnergyLimit(const G4ParticleDefinition*);

  virtual G4double LowEnergyLimit(const G4ParticleDefinition*);

  virtual void SetHighEnergyLimit(G4double);

  virtual void SetLowEnergyLimit(G4double);

  virtual G4double MinEnergyCut(const G4ParticleDefinition*,
                                const G4MaterialCutsCouple*);

  virtual G4bool IsInCharge(const G4ParticleDefinition*) = 0;

  virtual G4double ComputeDEDXPerVolume(
                               const G4Material*,
			       const G4ParticleDefinition*,
                                     G4double kineticEnergy,
                                     G4double cutEnergy = DBL_MAX);

  virtual G4double ComputeDEDX(const G4MaterialCutsCouple*,
                               const G4ParticleDefinition*,
                                     G4double kineticEnergy,
                                     G4double cutEnergy = DBL_MAX);

  virtual G4double CrossSectionPerVolume(
                               const G4Material*,
			       const G4ParticleDefinition*,
                                     G4double kineticEnergy,
			             G4double cutEnergy = 0.0,
                                     G4double maxEnergy = DBL_MAX);

  virtual G4double CrossSection(const G4MaterialCutsCouple*,
                                const G4ParticleDefinition*,
                                      G4double kineticEnergy,
                                      G4double cutEnergy = 0.0,
                                      G4double maxEnergy = DBL_MAX);

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
                                      G4double& kinEnergy, 
                                      G4double& Z, 
                                      G4double& A, 
                                      G4double  cutEnergy = 0.0,
                                      G4double  maxEnergy = DBL_MAX);

  virtual G4DynamicParticle* SampleSecondary(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin = 0.0,
                                      G4double tmax = DBL_MAX) = 0;

  virtual std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin = 0.0,
                                      G4double tmax = DBL_MAX) = 0;

  virtual void CurrentAtomAndShell(G4int& Z, G4int& shell);

  virtual void SetDeexcitation(G4bool fluo, G4bool auger);

  G4double MaxSecondaryKinEnergy(const G4DynamicParticle* dynParticle);

  const G4String& GetName() const;

  // Methods for msc simulation
  virtual G4double GeomPathLength(G4PhysicsTable* theLambdaTable,
                            const G4MaterialCutsCouple*,
		            const G4ParticleDefinition*,
		                  G4double& kinEnergy,
			          G4double  lambda,
			          G4double  range,
    			          G4double  truePathLength);

  virtual G4double TrueStepLength(G4double geomStepLength);

  virtual G4double SampleCosineTheta(G4double trueStepLength,
                                     G4double kinEnergy);

  virtual G4double SampleDisplacement();

  virtual void DefineForRegion(const G4Region*);

  virtual void SetDynamicParticle(const G4DynamicParticle*);

protected:

  virtual G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
    				               G4double kineticEnergy);

private:

  //  hide assignment operator
  G4VEmModel & operator=(const  G4VEmModel &right);
  G4VEmModel(const  G4VEmModel&);

  const G4String  name;
  G4double        lowLimit;
  G4double        highLimit;

protected:

  G4int           currentZ;
  G4int           currentShell;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4VEmModel::HighEnergyLimit(const G4ParticleDefinition*)
{
  return highLimit;
}

inline G4double G4VEmModel::LowEnergyLimit(const G4ParticleDefinition*)
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

inline G4double G4VEmModel::CrossSectionPerVolume(
                                        const G4Material* material,
					const G4ParticleDefinition* p,
					      G4double ekin,
					      G4double emin,
                                              G4double emax)
{
  G4double cross = 0.0;
  const G4ElementVector* theElementVector = material->GetElementVector() ;
  const G4double* theAtomNumDensityVector = material->GetVecNbOfAtomsPerVolume();
  size_t nelm = material->GetNumberOfElements();
  for (size_t i=0; i<nelm; i++) {
    const G4Element* elm = (*theElementVector)[i];
    G4double Z = elm->GetZ(); 
    G4double A = elm->GetA();
    cross += theAtomNumDensityVector[i]*ComputeCrossSectionPerAtom(p, ekin, Z, A, emin, emax); 
  }
  return cross;
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
                                      G4double&, G4double&, G4double&, G4double, G4double)
{
  return 0.0;
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

inline void G4VEmModel::CurrentAtomAndShell(G4int& Z, G4int& shell)
{
  Z = currentZ;
  shell = currentShell;
}

inline void G4VEmModel::SetDeexcitation(G4bool, G4bool)
{}

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
};

inline G4double G4VEmModel::SampleCosineTheta(G4double,G4double ) 
{
  return 1.0;
}

inline G4double G4VEmModel::SampleDisplacement() 
{
  return 0.0;
}

inline void G4VEmModel::DefineForRegion(const G4Region*) 
{}

inline void G4VEmModel::SetDynamicParticle(const G4DynamicParticle*) 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

