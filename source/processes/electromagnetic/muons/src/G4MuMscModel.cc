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
// $Id: G4MuMscModel.cc,v 1.4 2007-10-29 18:50:27 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4MuMscModel
//
// Author:      Laszlo Mu
//
// Creation date: 03.03.2001
//
// Modifications:
//
// 27-03-03 Move model part from G4MultipleScattering80 (V.Ivanchenko)
//

// Class Description:
//
// Implementation of the model of multiple scattering based on
// H.W.Lewis Phys Rev 78 (1950) 526 and others

// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MuMscModel.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4TransportationManager.hh"
#include "G4SafetyHelper.hh"
#include "G4eCoulombScatteringModel.hh"

#include "G4Poisson.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4MuMscModel::G4MuMscModel(G4double frange, 
			   G4double thetaMax, 
			   G4double tMax,  
			   const G4String& nam)
  : G4eCoulombScatteringModel(0.0,thetaMax,false,tMax,nam),
    theLambdaTable(0),
    theLambda1Table(0),
    //    theLambda2Table(0),
    dtrl(0.05),
    facrange(frange),
    thetaLimit(thetaMax),
    numlimit(0.2),
    couple(0),
    isInitialized(false),
    buildTables(true),
    inside(false),
    newrun(true)
{
  invsqrt12 = 1./sqrt(12.);
  tlimitminfix = 1.e-6*mm;
  xSection = currentRange = 0.0;
  theManager = G4LossTableManager::Instance(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MuMscModel::~G4MuMscModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuMscModel::Initialise(const G4ParticleDefinition* p,
			      const G4DataVector& v)
{

  // set values of some data members
  if(!isInitialized) {
    SetupParticle(p);
    if(p->GetParticleName() == "GenericIon") buildTables = false;

    if (pParticleChange)
      fParticleChange = reinterpret_cast<G4ParticleChangeForMSC*>(pParticleChange);
    else
      fParticleChange = new G4ParticleChangeForMSC();

    safetyHelper = G4TransportationManager::GetTransportationManager()
      ->GetSafetyHelper();
    safetyHelper->InitialiseHelper();
  }
  G4eCoulombScatteringModel::Initialise(p, v);
  newrun = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuMscModel::BuildTables()
{
  newrun = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuMscModel::ComputeCrossSectionPerAtom( 
                             const G4ParticleDefinition* p,
			     G4double kinEnergy,
			     G4double Z, G4double A,
			     G4double cutEnergy, G4double)
{
  if(p == particle && kinEnergy == tkin && Z == targetZ &&
     cutEnergy == ecut) return xSection;
  ecut = cutEnergy;
  xSection = 0.0;
  SetupParticle(p);
  G4double ekin = std::max(keV, kinEnergy);
  SetupTarget(Z, A, ekin);

  G4double tmax = tkin;
  if(p == theElectron) tmax *= 0.5;
  else if(p != thePositron) {
    G4double ratio = electron_mass_c2/mass;
    tmax = 2.0*mom2/
      (electron_mass_c2*(1.0 + ratio*(tkin/mass + 1.0) + ratio*ratio));
  }
  G4double t = std::min(cutEnergy, tmax);
  G4double mom21 = t*(t + 2.0*electron_mass_c2);
  t = tkin - t;
  G4double mom22 = t*(t + 2.0*mass);
  cosTetMaxElec = (mom2 + mom22 - mom21)*0.5/sqrt(mom2*mom22);
  if(cosTetMaxElec < cosTetMaxNuc) cosTetMaxElec = cosTetMaxNuc;

  if(cosTetMaxElec < 1.0) {
    G4double x2 = screenZ/(1.0 - cosTetMaxElec + screenZ);
    xSection += coeff*Z*chargeSquare*invbeta2*(x2 - 1.0 - log(x2))/mom2;
  }
  //  G4cout << "cut= " << ecut << " e= " << tkin << " croosE= " 
  //  << xSection/barn << G4endl;

  // limit integral because of nuclear size effect
  G4double costm = cosTetMaxNuc;
  if(formfactA > 2.01) {
    G4double ctet = sqrt(1.0 - 2.0/formfactA);
    if(ctet > costm) costm = ctet; 
  }

  if(costm < 1.0) {
    G4double x1 = 1.0 - costm; 
    G4double x2 = screenZ/(x1 + screenZ);
    G4double x3 = 1.0/(1.0 + formfactA*x1);
    xSection += 
      coeff*Z*Z*chargeSquare*invbeta2*(log(x3/x2) - 2.0 + x2 + x3)/mom2; 
  }
  //  G4cout << " croosE= " << xSection/barn << " screenZ= " 
  //	 << screenZ << " formF= " << formfactA << G4endl;

  return xSection; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuMscModel::ComputeLambda1(const G4MaterialCutsCouple*,
				      const G4ParticleDefinition*,
				      G4double kinEnergy)
{
  G4double x = 0.0;
  if(kinEnergy < 0.0) x = 1.0;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuMscModel::ComputeTruePathLengthLimit(
                             const G4Track& track,
			     G4PhysicsTable* theTable,
			     G4double currentMinimalStep)
{
  G4double tlimit = currentMinimalStep;
  const G4DynamicParticle* dp = track.GetDynamicParticle();

  // initialisation for 1st step  
  if(track.GetCurrentStepNumber() == 1) {
    inside = false;
    SetupParticle(dp->GetDefinition());
    theLambdaTable = theTable;
  }

  // initialisation for each step  
  G4double kinEnergy = dp->GetKineticEnergy();
  DefineMaterial(track.GetMaterialCutsCouple());
  lambda0 = GetLambda(kinEnergy);
  currentRange = 
    theManager->GetRangeFromRestricteDEDX(particle,kinEnergy,couple);

  // extra check for abnormal situation
  // this check needed to run MSC with eIoni and eBrem inactivated
  if(tlimit > currentRange) tlimit = currentRange;

  // stop here if small range particle
  if(inside) return tlimit;   

  // pre step
  G4StepPoint* sp = track.GetStep()->GetPreStepPoint();
  G4StepStatus stepStatus = sp->GetStepStatus();
  G4double presafety = sp->GetSafety();

  // compute presafety again if presafety <= 0 and no boundary
  // i.e. when it is needed for optimization purposes
  if(stepStatus != fGeomBoundary && presafety < tlimitminfix) 
    presafety = safetyHelper->ComputeSafety(sp->GetPosition()); 

  //  G4cout << "G4MuMscModel::ComputeTruePathLengthLimit tlimit= " 
  //	 <<tlimit<<" safety= " << presafety
  //	 << " range= " <<currentRange<<G4endl;

  // far from geometry boundary
  if(currentRange < presafety) {
    inside = true;  
    
    // limit mean scattering angle
  } else {
    tlimit = std::min(facrange*lambda0, tlimit);
  }

  //  G4cout << "tlimit= " << tlimit  
  //	 << " currentMinimalStep= " << currentMinimalStep << G4endl;
  return tlimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuMscModel::ComputeGeomPathLength(G4double truelength)
{
  tPathLength  = truelength;
  zPathLength  = tPathLength;

  G4double tau = tPathLength/lambda0;
  lambdaeff    = lambda0;
  //G4cout << "ComputeGeomPathLength: tLength= " << tPathLength
  //	 << " lambda0= " << lambda0 << " tau= " << tau << G4endl; 
  // small step
  if(tau < numlimit) {
    par1 = -1. ;  
    par2 = par3 = 0. ;  
    zPathLength *= (1.0 - 0.5*tau + tau*tau/6.0);

    // medium step
  } else if(tPathLength < currentRange*dtrl) {
    zPathLength = lambda0*(1.0 - exp(-tau));

  } else if(tkin < mass) {

    par1 = 1./currentRange;
    par2 = 1./(par1*lambda0);
    par3 = 1.+ par2;
    lambdaeff = 1.0/(par1*par3);
    G4double x = tPathLength/currentRange;
    G4double x1;
    if(x < numlimit) x1 = x*(1.0  - 0.5*x + x*x/3.0);
    else             x1 = log(1.0 - x); 
    zPathLength = lambdaeff*(1.-exp(par3*x1));

  } else {

    G4double T1 = theManager->GetEnergy(particle,currentRange-tPathLength,couple);
    G4double lambda1 = GetLambda(T1);

    par1 = (lambda0-lambda1)/(lambda0*tPathLength) ;
    par2 = 1./(par1*lambda0) ;
    par3 = 1.+ par2 ;
    lambdaeff = 1.0/(par1*par3);
    zPathLength = lambdaeff*(1.-exp(par3*log(lambda1/lambda0)));
  }

  //  if(zPathLength > lambda0) zPathLength = lambda0;
  if(zPathLength > tPathLength) zPathLength = tPathLength;

  return zPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuMscModel::ComputeTrueStepLength(G4double geomStepLength)
{
  // step defined other than transportation 
  if(geomStepLength == zPathLength) return tPathLength;

  tPathLength  = geomStepLength;
  zPathLength  = geomStepLength;
  G4double tau = geomStepLength/lambda0;
  if(tau < numlimit) {
    tPathLength *= (1.0 + 0.5*tau - tau*tau/3.0); 
  
  } else if(par1 <  0.) {
    tPathLength = -lambda0*log(1.0 - tau); 

  } else {
    G4double x = par1*par3*geomStepLength;
    if(x < numlimit) 
      tPathLength = (1.- exp(- x*(1.- 0.5*x + x*x/3.0)/par3))/par1 ;
    else if (x < 1.0)  
      tPathLength = (1.-exp(log(1.- x)/par3))/par1;
    else 
      tPathLength = currentRange;
  }
  if(tPathLength < geomStepLength) tPathLength = geomStepLength;

  return tPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuMscModel::SampleScattering(const G4DynamicParticle* dynParticle,
				    G4double safety)
{
  if(dynParticle->GetKineticEnergy() == 0.0) return;
  G4double rms  = sqrt(tPathLength/lambdaeff);
  G4double dirx, diry, dirz, tetx, tety;
  //  G4cout << "G4MuMscModel::SampleSecondaries: tstep(mm)= " << truestep/mm
  //	 << " lambdaeff= " << lambdaeff
  //	 << " rms= " << rms << G4endl;
  do {
    tetx = G4RandGauss::shoot(0.0,rms);
    tety = G4RandGauss::shoot(0.0,rms);
    dirx = sin(tetx);
    diry = sin(tety);
    dirz = 1.0 - dirx*dirx - diry*diry;
  } while(dirz < 0.0);

  G4ThreeVector oldDirection = dynParticle->GetMomentumDirection();
  G4ThreeVector newDirection(dirx,diry,sqrt(dirz));
  newDirection.rotateUz(oldDirection);
  fParticleChange->ProposeMomentumDirection(newDirection);

  if (latDisplasment && safety > tlimitminfix) {

    G4double rx = zPathLength*(0.5*tetx + invsqrt12*G4RandGauss::shoot(0.0,rms));
    G4double ry = zPathLength*(0.5*tety + invsqrt12*G4RandGauss::shoot(0.0,rms));
    G4double r  = sqrt(rx*rx + ry*ry);
/*
    G4cout << "G4MuMscModel::SampleSecondaries: e(MeV)= " << kineticEnergy
	   << " sinTheta= " << sth << " r(mm)= " << r
           << " trueStep(mm)= " << truestep 
           << " geomStep(mm)= " << zPathLength
           << G4endl;
*/
    
    G4ThreeVector latDirection(rx,ry,0.0);
    latDirection.rotateUz(oldDirection);

    G4ThreeVector Position = *(fParticleChange->GetProposedPosition());
    G4double fac = 1.;
    if(r >  safety) {
      //  ******* so safety is computed at boundary too ************
      G4double newsafety = safetyHelper->ComputeSafety(Position);
      if(r > newsafety)
	fac = newsafety/r ;
    }  

    if(fac > 0.) {
      // compute new endpoint of the Step
      G4ThreeVector newPosition = Position+fac*r*latDirection;

      // definitely not on boundary
      if(1. == fac) {
	safetyHelper->ReLocateWithinVolume(newPosition);
	    
      } else {
	// check safety after displacement
	G4double postsafety = safetyHelper->ComputeSafety(newPosition);

	// displacement to boundary
	if(postsafety <= 0.0) {
	  safetyHelper->Locate(newPosition, newDirection);

	  // not on the boundary
	} else { 
	      safetyHelper->ReLocateWithinVolume(newPosition);
	}
      }
      fParticleChange->ProposePosition(newPosition);
    } 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuMscModel::SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                     const G4MaterialCutsCouple*,
				     const G4DynamicParticle*,
				     G4double,
				     G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
