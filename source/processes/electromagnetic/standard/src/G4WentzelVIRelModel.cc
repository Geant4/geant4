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
// $Id: G4WentzelVIRelModel.cc 88979 2015-03-17 10:10:21Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4WentzelVIRelModel
//
// Author:      V.Ivanchenko 
//
// Creation date: 08.06.2012 from G4WentzelVIRelModel
//
// Modifications:
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4WentzelVIRelModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4LossTableManager.hh"
#include "G4Pow.hh"
#include "G4NistManager.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4WentzelVIRelModel::G4WentzelVIRelModel(G4bool combined) :
  G4VMscModel("WentzelVIRel"),
  numlimit(0.1),
  currentCouple(0),
  cosThetaMin(1.0),
  isCombined(combined),
  inside(false),
  singleScatteringMode(false)
{
  invsqrt12 = 1./sqrt(12.);
  tlimitminfix = 1.e-6*mm;
  lowEnergyLimit = 1.0*eV;
  particle = 0;
  nelments = 5;
  xsecn.resize(nelments);
  prob.resize(nelments);
  theManager = G4LossTableManager::Instance();
  fNistManager = G4NistManager::Instance();
  fG4pow = G4Pow::GetInstance();
  wokvi = new G4WentzelVIRelXSection(combined);

  preKinEnergy = tPathLength = zPathLength = lambdaeff = currentRange = 
    xtsec = 0;
  currentMaterialIndex = 0;
  cosThetaMax = cosTetMaxNuc = 1.0;

  fParticleChange = 0;
  currentCuts = 0;
  currentMaterial = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4WentzelVIRelModel::~G4WentzelVIRelModel()
{
  delete wokvi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WentzelVIRelModel::Initialise(const G4ParticleDefinition* p,
				     const G4DataVector& cuts)
{
  // reset parameters
  SetupParticle(p);
  currentRange = 0.0;

  if(isCombined) {
    G4double tet = PolarAngleLimit();
    if(tet >= pi)      { cosThetaMax = -1.0; }
    else if(tet > 0.0) { cosThetaMax = cos(tet); }
  }

  wokvi->Initialise(p, cosThetaMax);
  /*  
  G4cout << "G4WentzelVIRelModel: " << particle->GetParticleName()
         << "  1-cos(ThetaLimit)= " << 1 - cosThetaMax 
	 << G4endl;
  */
  currentCuts = &cuts;

  // set values of some data members
  fParticleChange = GetParticleChangeForMSC(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIRelModel::ComputeCrossSectionPerAtom( 
                             const G4ParticleDefinition* p,
			     G4double kinEnergy,
			     G4double Z, G4double,
			     G4double cutEnergy, G4double)
{
  G4double cross = 0.0;
  if(p != particle) { SetupParticle(p); }
  if(kinEnergy < lowEnergyLimit) { return cross; }
  if(!CurrentCouple()) {
    G4Exception("G4WentzelVIRelModel::ComputeCrossSectionPerAtom", "em0011",
		FatalException, " G4MaterialCutsCouple is not defined");
    return 0.0;
  }
  DefineMaterial(CurrentCouple());
  cosTetMaxNuc = wokvi->SetupKinematic(kinEnergy, currentMaterial);
  if(cosTetMaxNuc < 1.0) {
    G4double cost = wokvi->SetupTarget(G4lrint(Z), cutEnergy);
    cross = wokvi->ComputeTransportCrossSectionPerAtom(cost);
    /*
    if(p->GetParticleName() == "e-")      
    G4cout << "G4WentzelVIRelModel::CS: Z= " << G4int(Z) 
           << " e(MeV)= " << kinEnergy 
	   << " 1-cosN= " << 1 - cosTetMaxNuc << " cross(bn)= " << cross/barn
	   << " " << particle->GetParticleName() << G4endl;
    */
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WentzelVIRelModel::StartTracking(G4Track* track)
{
  SetupParticle(track->GetDynamicParticle()->GetDefinition());
  inside = false;
  G4VEmModel::StartTracking(track);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIRelModel::ComputeTruePathLengthLimit(
                             const G4Track& track,
			     G4double& currentMinimalStep)
{
  G4double tlimit = currentMinimalStep;
  const G4DynamicParticle* dp = track.GetDynamicParticle();
  G4StepPoint* sp = track.GetStep()->GetPreStepPoint();
  G4StepStatus stepStatus = sp->GetStepStatus();
  singleScatteringMode = false;
  //G4cout << "G4WentzelVIRelModel::ComputeTruePathLengthLimit stepStatus= " 
  //	 << stepStatus << G4endl;


  // initialisation for each step, lambda may be computed from scratch
  preKinEnergy  = dp->GetKineticEnergy();
  DefineMaterial(track.GetMaterialCutsCouple());
  lambdaeff = GetTransportMeanFreePath(particle,preKinEnergy);
  currentRange = GetRange(particle,preKinEnergy,currentCouple);
  cosTetMaxNuc = wokvi->SetupKinematic(preKinEnergy, currentMaterial);

  // extra check for abnormal situation
  // this check needed to run MSC with eIoni and eBrem inactivated
  if(tlimit > currentRange) { tlimit = currentRange; }

  // stop here if small range particle
  if(inside || tlimit < tlimitminfix) { 
    return ConvertTrueToGeom(tlimit, currentMinimalStep); 
  }

  // pre step
  G4double presafety = sp->GetSafety();
  // far from geometry boundary
  if(currentRange < presafety) {
    inside = true;  
    return ConvertTrueToGeom(tlimit, currentMinimalStep);
  }

  // compute presafety again if presafety <= 0 and no boundary
  // i.e. when it is needed for optimization purposes
  if(stepStatus != fGeomBoundary && presafety < tlimitminfix) {
    presafety = ComputeSafety(sp->GetPosition(), tlimit); 
    if(currentRange < presafety) {
      inside = true;  
      return ConvertTrueToGeom(tlimit, currentMinimalStep);
    }
  }
  /*  
  G4cout << "e(MeV)= " << preKinEnergy/MeV
	 << "  " << particle->GetParticleName() 
	 << " CurLimit(mm)= " << tlimit/mm <<" safety(mm)= " << presafety/mm
	 << " R(mm)= " <<currentRange/mm
	 << " L0(mm^-1)= " << lambdaeff*mm 
	 <<G4endl;
  */

  // natural limit for high energy
  G4double rlimit = std::max(facrange*currentRange, 
			     0.7*(1.0 - cosTetMaxNuc)*lambdaeff);

  // low-energy e-
  if(cosThetaMax > cosTetMaxNuc) {
    rlimit = std::min(rlimit, facsafety*presafety);
  }
   
  // cut correction
  G4double rcut = currentCouple->GetProductionCuts()->GetProductionCut(1);
  //G4cout << "rcut= " << rcut << " rlimit= " << rlimit 
  //   << " presafety= " << presafety 
  // << " 1-cosThetaMax= " <<1-cosThetaMax 
  // << " 1-cosTetMaxNuc= " << 1-cosTetMaxNuc
  // << G4endl;
  if(rcut > rlimit) { rlimit = std::min(rlimit, rcut*sqrt(rlimit/rcut)); }

  if(rlimit < tlimit) { tlimit = rlimit; }

  tlimit = std::max(tlimit, tlimitminfix);

  // step limit in infinite media
  tlimit = std::min(tlimit, 50*currentMaterial->GetRadlen()/facgeom);

  //compute geomlimit and force few steps within a volume
  if(steppingAlgorithm == fUseDistanceToBoundary && stepStatus == fGeomBoundary)
    {
      G4double geomlimit = ComputeGeomLimit(track, presafety, currentRange);
      tlimit = std::min(tlimit, geomlimit/facgeom);
    } 

  /*  
  G4cout << particle->GetParticleName() << " e= " << preKinEnergy
	 << " L0= " << lambdaeff << " R= " << currentRange
	 << "tlimit= " << tlimit  
  	 << " currentMinimalStep= " << currentMinimalStep << G4endl;
  */
  return ConvertTrueToGeom(tlimit, currentMinimalStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIRelModel::ComputeGeomPathLength(G4double truelength)
{
  tPathLength  = truelength;
  zPathLength  = tPathLength;

  if(lambdaeff > 0.0 && lambdaeff != DBL_MAX) {
    G4double tau = tPathLength/lambdaeff;
    //G4cout << "ComputeGeomPathLength: tLength= " << tPathLength
    //	 << " Leff= " << lambdaeff << " tau= " << tau << G4endl; 
    // small step
    if(tau < numlimit) {
      zPathLength *= (1.0 - 0.5*tau + tau*tau/6.0);

      // medium step
    } else {
      G4double e1 = 0.0;
      if(currentRange > tPathLength) {
	e1 = GetEnergy(particle,currentRange-tPathLength,currentCouple);
      }
      e1 = 0.5*(e1 + preKinEnergy);
      cosTetMaxNuc = wokvi->SetupKinematic(e1, currentMaterial);
      lambdaeff = GetTransportMeanFreePath(particle,e1);
      zPathLength = lambdaeff*(1.0 - G4Exp(-tPathLength/lambdaeff));
    }
  } else { lambdaeff = DBL_MAX; }
  //G4cout<<"Comp.geom: zLength= "<<zPathLength
  // <<" tLength= "<<tPathLength<<G4endl;
  return zPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIRelModel::ComputeTrueStepLength(G4double geomStepLength)
{
  // initialisation of single scattering x-section
  xtsec = 0.0;
  cosThetaMin = cosTetMaxNuc;

  //G4cout << "Step= " << geomStepLength << "  Lambda= " <<  lambdaeff 
  //	 << " 1-cosThetaMaxNuc= " << 1 - cosTetMaxNuc << G4endl;
  // pathalogical case
  if(lambdaeff == DBL_MAX) { 
    singleScatteringMode = true;
    zPathLength  = geomStepLength;
    tPathLength  = geomStepLength;
    cosThetaMin = 1.0;

    // normal case
  } else {

    // small step use only single scattering
    static const G4double singleScatLimit = 1.0e-7;
    if(geomStepLength < lambdaeff*singleScatLimit*(1.0 - cosTetMaxNuc)) {
      singleScatteringMode = true;
      cosThetaMin = 1.0;
      zPathLength  = geomStepLength;
      tPathLength  = geomStepLength;

      // step defined by transportation 
    } else if(geomStepLength != zPathLength) { 

      // step defined by transportation 
      zPathLength  = geomStepLength;
      G4double tau = geomStepLength/lambdaeff;
      tPathLength  = zPathLength*(1.0 + 0.5*tau + tau*tau/3.0); 

      // energy correction for a big step
      if(tau > numlimit) {
	G4double e1 = 0.0;
	if(currentRange > tPathLength) {
	  e1 = GetEnergy(particle,currentRange-tPathLength,currentCouple);
	}
	e1 = 0.5*(e1 + preKinEnergy);
	cosTetMaxNuc = wokvi->SetupKinematic(e1, currentMaterial);
	lambdaeff = GetTransportMeanFreePath(particle,e1);
	tau = zPathLength/lambdaeff;
      
	if(tau < 0.999999) { tPathLength = -lambdaeff*G4Log(1.0 - tau); } 
	else               { tPathLength = currentRange; }
      }
    }
  }

  // check of step length
  // define threshold angle between single and multiple scattering 
  if(!singleScatteringMode) { cosThetaMin = 1.0 - 1.5*tPathLength/lambdaeff; }

  // recompute transport cross section - do not change energy
  // anymore - cannot be applied for big steps
  if(cosThetaMin > cosTetMaxNuc) {

    // new computation
    G4double cross = ComputeXSectionPerVolume();
    //G4cout << "%%%% cross= " << cross << "  xtsec= " << xtsec << G4endl;
    if(cross <= 0.0) {
      singleScatteringMode = true;
      tPathLength = zPathLength; 
      lambdaeff = DBL_MAX;
    } else if(xtsec > 0.0) {

      lambdaeff = 1./cross; 
      G4double tau = zPathLength*cross;
      if(tau < numlimit) { 
	tPathLength = zPathLength*(1.0 + 0.5*tau + tau*tau/3.0); 
      } 
      else if(tau < 0.999999) { tPathLength = -lambdaeff*G4Log(1.0 - tau); } 
      else                    { tPathLength = currentRange; }

      if(tPathLength > currentRange) { tPathLength = currentRange; }
    } 
  }

  /*      
  G4cout <<"Comp.true: zLength= "<<zPathLength<<" tLength= "<<tPathLength
	 <<" Leff(mm)= "<<lambdaeff/mm<<" sig0(1/mm)= " << xtsec <<G4endl;
  G4cout << particle->GetParticleName() << " 1-cosThetaMin= " << 1-cosThetaMin
	 << " 1-cosTetMaxNuc= " << 1-cosTetMaxNuc 
	 << " e(MeV)= " << preKinEnergy/MeV << "  "  << singleScatteringMode 
	 << G4endl;
  */
  return tPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector& 
G4WentzelVIRelModel::SampleScattering(const G4ThreeVector& oldDirection,
				      G4double safety)
{
  fDisplacement.set(0.0,0.0,0.0);
  //G4cout << "!##! G4WentzelVIRelModel::SampleScattering for " 
  //	 << particle->GetParticleName() << G4endl;

  // ignore scattering for zero step length and energy below the limit
  if(preKinEnergy < lowEnergyLimit || tPathLength <= 0.0) 
    { return fDisplacement; }
  
  G4double invlambda = 0.0;
  if(lambdaeff < DBL_MAX) { invlambda = 0.5/lambdaeff; }

  // use average kinetic energy over the step
  G4double cut = (*currentCuts)[currentMaterialIndex];
  /*
  G4cout <<"SampleScat: E0(MeV)= "<< preKinEnergy/MeV
  	 << " Leff= " << lambdaeff <<" sig0(1/mm)= " << xtsec 
 	 << " xmsc= " <<  tPathLength*invlambda 
	 << " safety= " << safety << G4endl;
  */

  // step limit due msc
  G4double x0 = tPathLength;
  //  const G4double thinlimit = 0.5; 
  static const G4double thinlimit = 0.1; 
  G4int nMscSteps = 1;
  // large scattering angle case - two step approach
  if(tPathLength*invlambda > thinlimit && safety > tlimitminfix) { 
    x0 *= 0.5; 
    nMscSteps = 2; 
  } 

  // step limit due to single scattering
  G4double x1 = 2*tPathLength;
  if(0.0 < xtsec) { x1 = -G4Log(G4UniformRand())/xtsec; }

  const G4ElementVector* theElementVector = 
    currentMaterial->GetElementVector();
  G4int nelm = currentMaterial->GetNumberOfElements();

  // geometry
  G4double sint, cost, phi;
  G4ThreeVector temp(0.0,0.0,1.0);

  // current position and direction relative to the end point
  // because of magnetic field geometry is computed relatively to the 
  // end point of the step 
  G4ThreeVector dir(0.0,0.0,1.0);
  fDisplacement.set(0.0,0.0,-zPathLength);
  G4double mscfac = zPathLength/tPathLength;

  // start a loop 
  G4double x2 = x0;
  G4double step, z;
  G4bool singleScat;
  /*
    G4cout << "Start of the loop x1(mm)= " << x1 << "  x2(mm)= " << x2 
  	 << " 1-cost1= " << 1 - cosThetaMin << "  " << singleScatteringMode 
  	 << " xtsec= " << xtsec << G4endl;
  */
  do {

    // single scattering case
    if(x1 < x2) { 
      step = x1;
      singleScat = true;
    } else {
      step = x2;
      singleScat = false;
    }

    // new position
    fDisplacement += step*mscfac*dir;

    if(singleScat) {

      // select element
      G4int i = 0;
      if(nelm > 1) {
	G4double qsec = G4UniformRand()*xtsec;
	for (; i<nelm; ++i) { if(xsecn[i] >= qsec) { break; } }
      }
      G4double cosTetM = 
	wokvi->SetupTarget(G4lrint((*theElementVector)[i]->GetZ()), cut);
      //G4cout << "!!! " << cosThetaMin << "  " << cosTetM 
      // << "  " << prob[i] << G4endl;
      temp = wokvi->SampleSingleScattering(cosThetaMin, cosTetM, prob[i]);

      // direction is changed
      temp.rotateUz(dir);
      dir = temp;

      // new proposed step length
      x1  = -G4Log(G4UniformRand())/xtsec; 
      x2 -= step; 
      if(x2 <= 0.0) { --nMscSteps; }

    // multiple scattering
    } else { 
      --nMscSteps;
      x1 -= step;
      x2  = x0;

      // for multiple scattering x0 should be used as a step size
      if(!singleScatteringMode) {
        G4double z0 = x0*invlambda;

	// correction to keep first moment
 
	// sample z in interval 0 - 1
        if(z0 > 5.0) { z = G4UniformRand(); }
	else {
	  G4double zzz = 0.0;
	  if(z0 > 0.01) { zzz = G4Exp(-1.0/z0); }
	  z = -z0*G4Log(1.0 - (1.0 - zzz)*G4UniformRand());
	  //  /(1.0 - (1.0/z0 + 1.0)*zzz); 
	}

	cost = 1.0 - 2.0*z/*factCM*/;
	if(cost > 1.0)       { cost = 1.0; }
	else if(cost < -1.0) { cost =-1.0; }
	sint = sqrt((1.0 - cost)*(1.0 + cost));
	phi  = twopi*G4UniformRand();
	G4double vx1 = sint*cos(phi);
	G4double vy1 = sint*sin(phi);

	// lateral displacement  
	if (latDisplasment && x0 > tlimitminfix) {
	  G4double rms = invsqrt12*sqrt(2*z0);
          G4double r   = x0*mscfac;
	  G4double dx  = r*(0.5*vx1 + rms*G4RandGauss::shoot(0.0,1.0));
	  G4double dy  = r*(0.5*vy1 + rms*G4RandGauss::shoot(0.0,1.0));
	  G4double dz;
	  G4double d   = r*r - dx*dx - dy*dy;
	  if(d >= 0.0)  { dz = sqrt(d) - r; }
	  else          { dx = dy = dz = 0.0; }

	  // change position
	  temp.set(dx,dy,dz);
	  temp.rotateUz(dir); 
	  fDisplacement += temp;
	}
	// change direction
	temp.set(vx1,vy1,cost);
	temp.rotateUz(dir);
	dir = temp;
      }
    }
  } while (0 < nMscSteps);
    
  dir.rotateUz(oldDirection);

  //G4cout << "G4WentzelVIRelModel sampling of scattering is done" << G4endl;
  // end of sampling -------------------------------

  fParticleChange->ProposeMomentumDirection(dir);

  // lateral displacement  
  fDisplacement.rotateUz(oldDirection);

  /*            
	 G4cout << " r(mm)= " << fDisplacement.mag() 
		<< " safety= " << safety
		<< " trueStep(mm)= " << tPathLength
		<< " geomStep(mm)= " << zPathLength
		<< " x= " << fDisplacement.x() 
		<< " y= " << fDisplacement.y() 
		<< " z= " << fDisplacement.z()
		<< G4endl;
  */

  //G4cout<< "G4WentzelVIRelModel::SampleScattering end NewDir= " 
  // << dir<< G4endl;
  return fDisplacement;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIRelModel::ComputeXSectionPerVolume()
{
  // prepare recomputation of x-sections
  const G4ElementVector* theElementVector = currentMaterial->GetElementVector();
  const G4double* theAtomNumDensityVector = 
    currentMaterial->GetVecNbOfAtomsPerVolume();
  G4int nelm = currentMaterial->GetNumberOfElements();
  if(nelm > nelments) {
    nelments = nelm;
    xsecn.resize(nelm);
    prob.resize(nelm);
  }
  G4double cut = (*currentCuts)[currentMaterialIndex];
  //  cosTetMaxNuc = wokvi->GetCosThetaNuc();

  // check consistency
  xtsec = 0.0;
  if(cosTetMaxNuc > cosThetaMin) { return 0.0; }

  // loop over elements
  G4double xs = 0.0;
  for (G4int i=0; i<nelm; ++i) {
    G4double costm = 
      wokvi->SetupTarget(G4lrint((*theElementVector)[i]->GetZ()), cut);
    G4double density = theAtomNumDensityVector[i];

    G4double esec = 0.0;
    if(costm < cosThetaMin) {  

      // recompute the transport x-section
      if(1.0 > cosThetaMin) {
	xs += density*wokvi->ComputeTransportCrossSectionPerAtom(cosThetaMin);
      }
      // recompute the total x-section
      G4double nucsec = wokvi->ComputeNuclearCrossSection(cosThetaMin, costm);
      esec = wokvi->ComputeElectronCrossSection(cosThetaMin, costm);
      nucsec += esec;
      if(nucsec > 0.0) { esec /= nucsec; }
      xtsec += nucsec*density;
    }
    xsecn[i] = xtsec;
    prob[i]  = esec;
    //G4cout << i << "  xs= " << xs << " xtsec= " << xtsec 
    //       << " 1-cosThetaMin= " << 1-cosThetaMin 
    //	   << " 1-cosTetMaxNuc2= " <<1-cosTetMaxNuc2<< G4endl;
  }
  
  //G4cout << "ComputeXS result:  xsec(1/mm)= " << xs 
  //	 << " txsec(1/mm)= " << xtsec <<G4endl; 
  return xs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
