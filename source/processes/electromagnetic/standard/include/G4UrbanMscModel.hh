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
// $Id: G4UrbanMscModel.hh,v 1.5 2006-03-07 16:57:46 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4UrbanMscModel
//
// Author:        Laszlo Urban
//
// Creation date: 01.03.2001
//
// Modifications:
//
// 27-03-03 Move model part from G4MultipleScattering (V.Ivanchenko)
// 27-03-03 Rename (V.Ivanchenko)
//
// 05-08-03 angle distribution has been modified (L.Urban)
// 26-11-03 new data member currentRange (L.Urban)
// 01-03-04 changes in data members + signature changed in SampleCosineTheta
// 11-03-04 changes in data members (L.Urban)
// 23-04-04 changes in data members and in signature of SampleCosineTheta
//          (L.Urban)
// 17-08-04 name of data member facxsi changed to factail (L.Urban)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 15-04-05 optimize internal interface - add SampleSecondaries method (V.Ivanchenko)
// 11-08-05 computation of lateral correlation added (L.Urban)
// 02-10-05 nuclear size correction computation removed, the correction
//          included in the (theoretical) tabulated values (L.Urban)
// 16-02-06 data members b and xsi have been removed (L.Urban)
// 17-02-06 Save table of transport cross sections not mfp (V.Ivanchenko)
// 07-03-06 Create G4UrbanMscModel and move there step limit calculation (V.Ivanchenko)
//

//
// Class Description:
//
// Implementation of the model of multiple scattering based on
// H.W.Lewis Phys Rev 78 (1950) 526 and L.Urban model

// -------------------------------------------------------------------
//

#ifndef G4UrbanMscModel_h
#define G4UrbanMscModel_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4VEmModel.hh"
#include "G4PhysicsTable.hh"

class G4ParticleChangeForMSC;
class G4Navigator;
class G4LossTableManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4UrbanMscModel : public G4VEmModel
{

public:

  G4UrbanMscModel(G4double facrange, G4double dtrl, G4double tkinlimit, 
		  G4double facgeom, G4double factail, 
		  G4bool samplez, G4bool stepAlg, 
		  const G4String& nam = "UrbanMscUni");

  virtual ~G4UrbanMscModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double ComputeCrossSectionPerAtom( 
                             const G4ParticleDefinition* particle,
                                   G4double KineticEnergy,
                                   G4double AtomicNumber,
                                   G4double AtomicWeight=0., 
				   G4double cut =0.,
				   G4double emax=DBL_MAX);

  virtual std::vector<G4DynamicParticle*>* SampleSecondaries(
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
                                   G4double length,
                                   G4double safety);

  virtual G4double ComputeTruePathLengthLimit(
                             const G4Track& track,
			           G4PhysicsTable* theLambdaTable,
			           G4double currentMinimalStep);

  virtual G4double ComputeGeomPathLength(G4double truePathLength);

  virtual G4double ComputeTrueStepLength(G4double geomStepLength);

  void SetLateralDisplasmentFlag(G4bool val);

private:

  G4double SampleCosineTheta(G4double trueStepLength, G4double KineticEnergy);

  G4double SampleDisplacement();

  G4double LatCorrelation();

  G4double GetLambda(G4double kinEnergy);

  G4double GeomLimit(const G4Track& track);

  void SetParticle(const G4ParticleDefinition* p);

  //  hide assignment operator
  G4UrbanMscModel & operator=(const  G4UrbanMscModel &right);
  G4UrbanMscModel(const  G4UrbanMscModel&);

  const G4ParticleDefinition* particle;
  G4ParticleChangeForMSC*     fParticleChange;
  G4Navigator*                navigator;
  G4PhysicsTable*             theLambdaTable;
  const G4MaterialCutsCouple* couple;
  G4LossTableManager*         theManager;

  G4double mass;
  G4double charge;

  G4double taubig;
  G4double tausmall;
  G4double taulim;
  G4double currentTau;
  G4double dtrl;
  G4double factail ;

  G4double Tkinlimit;
  G4double Tlimit;
  G4double facrange;
  G4double tlimit;
  G4double tlimitmin;
  G4double geombig;
  G4double geommin;
  G4double facgeom;
  G4double safety;
  G4double facsafety;

  G4double lambda0;
  G4double lambdaeff;
  G4double tPathLength;
  G4double zPathLength;
  G4double geomLength;
  G4double par1,par2,par3 ;

  G4double stepmin ;

  G4double currentKinEnergy;
  G4double currentRange; 
  G4double currentRadLength;

  G4int    currentMaterialIndex;

  G4bool   samplez;
  G4bool   latDisplasment;
  G4bool   steppingAlgorithm;
  G4bool   isInitialized;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4UrbanMscModel::SetLateralDisplasmentFlag(G4bool val) 
{ 
  latDisplasment = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4double G4UrbanMscModel::GetLambda(G4double e)
{
  G4double x;
  if(theLambdaTable) {
    G4bool b;
    x = ((*theLambdaTable)[currentMaterialIndex])->GetValue(e, b);
  } else {
    x = CrossSection(couple,particle,e);
  }
  if(x > DBL_MIN) x = 1./x;
  else            x = DBL_MAX;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4UrbanMscModel::SetParticle(const G4ParticleDefinition* p)
{
  if (p != particle) {
    particle = p;
    mass = p->GetPDGMass();
    charge = p->GetPDGCharge()/eplus;
    Tlimit = Tkinlimit*electron_mass_c2/mass;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

