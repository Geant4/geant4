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
// $Id: G4WentzelVIModel.hh,v 1.29 2010-05-27 14:22:05 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4WentzelVIModel
//
// Author:        V.Ivanchenko 
//
// Creation date: 09.04.2008 from G4MuMscModel
//
// Modifications:
// 27-05-2010 V.Ivanchenko added G4WentzelOKandVIxSection class to
//              compute cross sections and sample scattering angle
//
// Class Description:
//
// Implementation of the model of multiple scattering based on
// G.Wentzel, Z. Phys. 40 (1927) 590.
// H.W.Lewis, Phys Rev 78 (1950) 526.
// J.M. Fernandez-Varea et al., NIM B73 (1993) 447.
// L.Urban, CERN-OPEN-2006-077.

// -------------------------------------------------------------------
//

#ifndef G4WentzelVIModel_h
#define G4WentzelVIModel_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4VMscModel.hh"
#include "G4PhysicsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4WentzelOKandVIxSection.hh"

class G4ParticleDefinition;
class G4LossTableManager;
class G4Pow;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4WentzelVIModel : public G4VMscModel
{

public:

  G4WentzelVIModel(const G4String& nam = "WentzelVIUni");

  virtual ~G4WentzelVIModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
					      G4double KineticEnergy,
					      G4double AtomicNumber,
					      G4double AtomicWeight=0., 
					      G4double cut = DBL_MAX,
					      G4double emax= DBL_MAX);

  virtual void SampleScattering(const G4DynamicParticle*, G4double safety);

  virtual G4double ComputeTruePathLengthLimit(const G4Track& track,
					      G4PhysicsTable* theLambdaTable,
					      G4double currentMinimalStep);

  virtual G4double ComputeGeomPathLength(G4double truePathLength);

  virtual G4double ComputeTrueStepLength(G4double geomStepLength);

private:

  G4double ComputeXSectionPerVolume();

  inline G4double GetLambda(G4double kinEnergy);

  inline void SetupParticle(const G4ParticleDefinition*);

  inline void DefineMaterial(const G4MaterialCutsCouple*);

  //  hide assignment operator
  G4WentzelVIModel & operator=(const  G4WentzelVIModel &right);
  G4WentzelVIModel(const  G4WentzelVIModel&);

  G4LossTableManager*       theManager;
  G4ParticleChangeForMSC*   fParticleChange;
  G4WentzelOKandVIxSection* wokvi;
  G4Pow*                    fG4pow;

  G4PhysicsTable*           theLambdaTable;
  const G4DataVector*       currentCuts;

  G4double tlimitminfix;
  G4double invsqrt12;

  // cache kinematics
  G4double preKinEnergy;
  G4double tPathLength;
  G4double zPathLength;
  G4double lambdaeff;
  G4double currentRange; 

  // data for single scattering mode
  G4double xtsec;
  std::vector<G4double> xsecn;
  std::vector<G4double> prob;
  G4int    nelments;

  G4double numlimit;

  // cache material
  G4int    currentMaterialIndex;
  const G4MaterialCutsCouple* currentCouple;
  const G4Material* currentMaterial;

  // single scattering parameters
  G4double cosThetaMin;
  G4double cosThetaMax;
  G4double cosTetMaxNuc;

  // projectile
  const G4ParticleDefinition* particle;
  G4double lowEnergyLimit;

  // flags
  G4bool   isInitialized;
  G4bool   inside;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4WentzelVIModel::DefineMaterial(const G4MaterialCutsCouple* cup) 
{ 
  if(cup != currentCouple) {
    currentCouple = cup;
    currentMaterial = cup->GetMaterial();
    currentMaterialIndex = currentCouple->GetIndex(); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4double G4WentzelVIModel::GetLambda(G4double e)
{
  G4double x;
  if(theLambdaTable) {
    x = ((*theLambdaTable)[currentMaterialIndex])->Value(e);
  } else {
    x = CrossSection(currentCouple,particle,e,
		     (*currentCuts)[currentMaterialIndex]);
  }
  if(x > DBL_MIN) { x = 1./x; }
  else            { x = DBL_MAX; }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4WentzelVIModel::SetupParticle(const G4ParticleDefinition* p)
{
  // Initialise mass and charge
  if(p != particle) {
    particle = p;
    wokvi->SetupParticle(p);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

