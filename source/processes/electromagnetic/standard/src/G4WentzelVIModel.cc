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
// $Id: G4WentzelVIModel.cc,v 1.1 2008-04-13 17:20:06 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4WentzelVIModel
//
// Author:      V.Ivanchenko 
//
// Creation date: 09.04.2008 from G4MuMscModel
//
// Modifications:
//
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

#include "G4WentzelVIModel.hh"
#include "Randomize.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4TransportationManager.hh"
#include "G4SafetyHelper.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4WentzelVIModel::G4WentzelVIModel(G4double thetaMax, 
				   G4double tMax,  
				   const G4String& nam) :
  G4VMscModel(nam),
  theLambdaTable(0),
  theLambda2Table(0),
  numlimit(0.2),
  lowBinEnergy(keV),
  highBinEnergy(PeV),
  nbins(60),
  nwarnings(0),
  nwarnlimit(50),
  currentCouple(0),
  cosThetaMin(1.0),
  cosThetaMax(cos(thetaMax)),
  q2Limit(tMax),
  alpha2(fine_structure_const*fine_structure_const),
  isInitialized(false),
  newrun(true),
  inside(false)
{
  invsqrt12 = 1./sqrt(12.);
  tlimitminfix = 1.e-6*mm;
  theManager = G4LossTableManager::Instance(); 
  fNistManager = G4NistManager::Instance();
  theElectron = G4Electron::Electron();
  thePositron = G4Positron::Positron();
  theProton   = G4Proton::Proton();
  a0 = alpha2*electron_mass_c2*electron_mass_c2/(0.885*0.885);
  G4double p0 = electron_mass_c2*classic_electr_radius;
  coeff  = twopi*p0*p0;
  constn = 6.937e-6/(MeV*MeV);
  tkin = targetZ = targetA = mom2 = DBL_MIN;
  ecut = etag = DBL_MAX;
  particle = 0;
  for(size_t j=0; j<100; j++) {
    FF[j]    = 0.0;
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4WentzelVIModel::~G4WentzelVIModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WentzelVIModel::Initialise(const G4ParticleDefinition* p,
			      const G4DataVector& cuts)
{
  SetupParticle(p);
  newrun = true;
  tkin = targetZ = targetA = mom2 = DBL_MIN;
  ecut = etag = DBL_MAX;
  xSection = currentRange = 0.0;
  // set values of some data members
  if(!isInitialized) {
    isInitialized = true;


    if (pParticleChange)
      fParticleChange = reinterpret_cast<G4ParticleChangeForMSC*>(pParticleChange);
    else
      fParticleChange = new G4ParticleChangeForMSC();

    safetyHelper = G4TransportationManager::GetTransportationManager()
      ->GetSafetyHelper();
    safetyHelper->InitialiseHelper();
  }
  currentCuts = &cuts;
  cosThetaLimit = cosThetaMax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIModel::ComputeCrossSectionPerAtom( 
                             const G4ParticleDefinition* p,
			     G4double kinEnergy,
			     G4double Z, G4double A,
			     G4double cutEnergy, G4double)
{
  if(p == particle && kinEnergy == tkin && Z == targetZ &&
     cutEnergy == ecut) return xSection;
  xSection = 0.0;
  SetupParticle(p);
  G4double ekin = std::max(keV, kinEnergy);
  SetupKinematic(ekin, cutEnergy);
  SetupTarget(Z, A, ekin);
  /*
  G4cout << "CS: e= " << tkin << " cosEl= " << cosTetMaxElec 
	 << " cosN= " << cosTetMaxNuc
	 << " cosH= " << cosTetMaxHad << G4endl;
  */
  G4double x, y, x1, x2, x3, x4;

  // scattering off electrons
  if(cosTetMaxElec < 1.0) {
    x = (1.0 - cosTetMaxElec)/screenZ;
    if(x < numlimit) y = 0.5*x*x*(1.0 - 1.3333333*x + 1.5*x*x); 
    else             y = log(1.0 + x) - x/(1.0 + x);
    if(y < 0.0) {
      nwarnings++;
      if(nwarnings < nwarnlimit /*&& y < -1.e-10*/) {
	G4cout << "Electron scattering <0 for L1 " << y 
	       << " e(MeV)= " << tkin << " p(MeV/c)= " << sqrt(mom2) 
	       << " Z= " << Z << "  " 
	       << particle->GetParticleName() << G4endl;
	G4cout << " z= " << 1.0-cosTetMaxElec << " screenZ= " << screenZ 
	       << " x= " << x << G4endl;
      }
      y = 0.0;
    }
    xSection += y/Z;
  }
  /*
  G4cout << "G4WentzelVIModel:XS per A " << " Z= " << Z << " e(MeV)= " << kinEnergy/MeV 
	 << " cut(MeV)= " << ecut/MeV  
  	 << " zmaxE= " << (1.0 - cosTetMaxElec)/screenZ 
	 << " zmaxN= " << (1.0 - cosTetMsxNuc)/screenZ << G4endl;
  */

  // scattering off nucleus
  G4double costm = std::max(cosTetMaxNuc,cosTetMaxHad);
  if(costm < 1.0) {
    x  = 1.0 - costm;
    x1 = screenZ*formfactA;
    x2 = 1.0/(1.0 - x1); 
    x3 = x/screenZ;
    x4 = formfactA*x;
    if(x3 < numlimit) {
      y = 0.5*x3*x3*x2*x2*x2*(1.0 - 1.333333*x3 + 1.5*x3*x3 - 1.5*x1);
    } else {
      y = ((1.0 + x1)*x2*log((1. + x3)/(1. + x4)) 
	   - x3/(1. + x3) - x4/(1. + x4))*x2*x2; 
    }
    if(y < 0.0) {
      nwarnings++;
      if(nwarnings < nwarnlimit /*&& y < -1.e-10*/) { 
	G4cout << "Nuclear scattering <0 for L1 " << y 
	       << " e(MeV)= " << tkin << " Z= " << Z << "  " 
	       << particle->GetParticleName() << G4endl;
	G4cout << " formfactA= " << formfactA << " screenZ= " << screenZ 
	       << " x= " << " x1= " << x1 << " x2= " << x2 
	       << " x3= " << x3 << " x4= " << x4 <<G4endl;
      }
      y = 0.0;
    }
    xSection += y; 
  }
  xSection *= (coeff*Z*Z*chargeSquare*invbeta2/mom2); 
  //  G4cout << "   XStotal= " << xSection/barn << " screenZ= " << screenZ 
  //	 << " formF= " << formfactA << " for " << p->GetParticleName() << G4endl;
  return xSection; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIModel::ComputeTruePathLengthLimit(
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

  // initialisation for each step, lambda may be computed from scratch
  cosThetaLimit = cosThetaMax;
  preKinEnergy  = dp->GetKineticEnergy();
  DefineMaterial(track.GetMaterialCutsCouple());
  lambda0 = GetLambda(preKinEnergy);
  currentRange = 
    theManager->GetRangeFromRestricteDEDX(particle,preKinEnergy,currentCouple);

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

  //  G4cout << "G4WentzelVIModel::ComputeTruePathLengthLimit tlimit= " 
  //	 <<tlimit<<" safety= " << presafety
  //	 << " range= " <<currentRange<<G4endl;

  // far from geometry boundary
  if(currentRange < presafety) {
    inside = true;  
    
    // limit mean scattering angle
  } else {
    G4double rlimit = facrange*lambda0;
    G4double rcut = currentCouple->GetProductionCuts()->GetProductionCut(1);
    if(rcut > rlimit) rlimit = std::pow(2.0*rcut*rcut*lambda0,0.33333333);
    if(rlimit < tlimit) tlimit = rlimit;
  }
  /*
  G4cout << particle->GetParticleName() << " e= " << preKinEnergy
	 << " L0= " << lambda0 << " R= " << currentRange
	 << "tlimit= " << tlimit  
  	 << " currentMinimalStep= " << currentMinimalStep << G4endl;
  */
  return tlimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIModel::ComputeGeomPathLength(G4double truelength)
{
  tPathLength  = truelength;
  zPathLength  = tPathLength;
  lambdaeff    = lambda0;

  if(lambda0 > 0.0) {
    G4double tau = tPathLength/lambda0;
    //G4cout << "ComputeGeomPathLength: tLength= " << tPathLength
    //	 << " lambda0= " << lambda0 << " tau= " << tau << G4endl; 
    // small step
    if(tau < numlimit) {
      zPathLength *= (1.0 - 0.5*tau + tau*tau/6.0);

      // medium step
    } else {
      zPathLength = lambda0*(1.0 - exp(-tPathLength/lambda0));
      /*
      G4double e1 = 0.0;
      if(currentRange > tPathLength) {
	e1 = theManager->GetEnergy(particle,
				   currentRange-tPathLength,
				   currentCouple);
      }
      lambdaeff = GetLambda(0.5*(e1 + preKinEnergy));
      zPathLength = lambdaeff*(1.0 - exp(-tPathLength/lambdaeff));
      */
    }
  }
  return zPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIModel::ComputeTrueStepLength(G4double geomStepLength)
{
  // step defined other than transportation 
  if(geomStepLength == zPathLength) return tPathLength;

  // step defined by transportation 
  tPathLength  = geomStepLength;
  zPathLength  = geomStepLength;
  G4double tau = zPathLength/lambdaeff;
  if(tau < numlimit) {
    tPathLength *= (1.0 + 0.5*tau + tau*tau/3.0); 

  } else {
    G4double e1 = 0.0;
    if(currentRange > zPathLength) {
      e1 = theManager->GetEnergy(particle,
				 currentRange-zPathLength,
				 currentCouple);
    }
    lambdaeff = GetLambda(0.5*(e1 + preKinEnergy));
    tau = zPathLength/lambdaeff;

    if(tau < 0.999999) tPathLength = -lambdaeff*log(1.0 - tau); 
    else               tPathLength = currentRange;

    if(tPathLength < zPathLength || tPathLength > currentRange) {
      tPathLength = zPathLength;
    }
  }
  return tPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WentzelVIModel::SampleScattering(const G4DynamicParticle* dynParticle,
					G4double safety)
{
  G4double kinEnergy = dynParticle->GetKineticEnergy();
  if(kinEnergy <= DBL_MIN || tPathLength <= DBL_MIN) return;
  
  G4double x1 = 0.5*tPathLength/lambdaeff;

  // define threshold angle as 2 sigma of central value
  cosThetaMin = 1.0 - 2.0*x1;

  // prepare kinematics
  cosThetaLimit = std::max(cosThetaMin,cosThetaMax);
  G4double e  = 0.5*(preKinEnergy + kinEnergy);
  const G4ParticleDefinition* part = dynParticle->GetDefinition();
  G4double cut = (*currentCuts)[currentMaterialIndex];
  SetupParticle(part);
  G4double ekin = std::max(keV, e);
  SetupKinematic(ekin, cut);
  /*
  G4cout << "SampleScat: cosTetMaxNuc= " << cosTetMaxNuc 
	 << " cosThetaMin= " << cosThetaMin 
	 << " cosThetaMax= " << cosThetaMax 
	 << " cosTetMaxElec= " << cosTetMaxElec
	 << " L0= " << lambda0 << " Leff= " << lambdaeff 
	 << G4endl;
  */
  G4double xsec = 0.0;
  if(cosThetaMin < cosTetMaxNuc) {
    cosThetaMin = cosTetMaxNuc;

    // recompute transport cross section
  } else {
    //G4cout << " cosTetMaxElec= " << cosTetMaxElec << " x1= " << x1 << G4endl;
    G4double x = CrossSection(currentCouple,particle,ekin);
    if(x > DBL_MIN) x1 = 0.5*tPathLength*x;
    else            x1 = 0.0;
    // compute cross section for the tail distribution 
    //G4cout << " cosTetMaxElec= " << cosTetMaxElec << G4endl;
    xsec = ComputeXSectionPerVolume();
  }
  // result of central part sampling
  G4double z; 
    
  //G4cout<<part->GetParticleName()<<" e= "<<tkin<<"  x1= "<<x1<<G4endl;
  
  do {
    z = -x1*log(G4UniformRand());
  } while (z > 1.0); 

  // cost is sampled ------------------------------
  G4double cost = 1.0 - 2.0*z;
  if(cost < -1.0) cost = -1.0;
  else if(cost > 1.0) cost = 1.0;
  G4double sint = sqrt((1.0 - cost)*(1.0 + cost));

  G4double phi  = twopi*G4UniformRand();

  G4double dirx = sint*cos(phi);
  G4double diry = sint*sin(phi);
  /*
  G4cout << "G4WentzelVIModel: step(mm)= " << tPathLength/mm
	 << " e= " << e << " Epre= " << preKinEnergy
	 << " sint= " << sint << " cost= " << cost
	 << " xsec= " << xsec << G4endl;
  */
  G4ThreeVector oldDirection = dynParticle->GetMomentumDirection();
  G4ThreeVector newDirection(dirx,diry,cost);
  G4ThreeVector newDirection1(0.0,0.0,1.0);
  G4ThreeVector pos(0.0,0.0,-zPathLength);
  G4ThreeVector dir(0.0,0.0,1.0);
  G4bool isscat = false;

  // sample scattering for large angle -------------
  if(xsec <= DBL_MIN) {
    pos += dir*zPathLength;
  } else {
    G4double t = tPathLength;
    do{
      G4double x  = -log(G4UniformRand())/xsec;      
      pos += dir*(zPathLength*std::min(x,t)/tPathLength);
      t -= x;
      if(t > 0.0) {
	G4double zz1 = 1.0;
	G4double qsec = G4UniformRand()*xsec;

	G4double costm = cosTetMaxNuc;
	G4double formf = 0.0;
	size_t nelm = currentMaterial->GetNumberOfElements();
	const G4ElementVector* theElementVector = 
	  currentMaterial->GetElementVector();

	// scattering off nucleus
	if(qsec >= xsece2) {
	  qsec -= xsece2;
	  for (size_t i=0; i<nelm; i++) {
	    if(xsecn[i] > qsec || nelm-1 == i) {
	      const G4Element* elm = (*theElementVector)[i];
	      SetupTarget(elm->GetZ(), elm->GetN(), tkin);
	      formf = formfactA;
              if(cosTetMaxHad > costm) costm = cosTetMaxHad;
	      //G4cout << "Reserford nuc Z= " << elm->GetZ() << G4endl;
	      break;
	    }
	  }
	  // scattering off electrons
	} else {
	  costm = cosTetMaxElec;
	  for (size_t i=0; i<nelm; i++) {
	    if(xsece[i] > qsec || nelm-1 == i) {
	      const G4Element* elm = (*theElementVector)[i];
	      SetupTarget(elm->GetZ(), elm->GetN(), tkin);
	      //G4cout << "Reserford elec Z= " << elm->GetZ() << G4endl;
	      break;
	    }
	  }
	}
	if(cosThetaMin > costm) {

	  G4double w1 = 1. - cosThetaMin + screenZ;
	  G4double w2 = 1. - costm + screenZ;
	  G4double w3 = cosThetaMin - costm;
	  G4double grej, zz; 
	  do {
	    zz = w1*w2/(w1 + G4UniformRand()*w3) - screenZ;
	    grej = 1.0/(1.0 + formf*zz);
	  } while ( G4UniformRand() > grej*grej );  
	  if(zz < 0.0) zz = 0.0;
	  else if(zz > 2.0) zz = 2.0;
	  zz1 = 1.0 - zz;
	}
        if(zz1 < 1.0) {
	  isscat = true;
	  //G4cout << "Reserford zz1= " << zz1 << G4endl;
	  sint = sqrt((1.0 - zz1)*(1.0 + zz1));
	  //G4cout << "sint= " << sint << G4endl;
	  phi  = twopi*G4UniformRand();

	  G4double vx1 = sint*cos(phi);
	  G4double vy1 = sint*sin(phi);
        
	  newDirection1.set(vx1,vy1,zz1);
	  newDirection1.rotateUz(dir);
	  dir = newDirection1;
	}
      }
    } while (t > 0.0); 
  }
  if(isscat) newDirection.rotateUz(dir);
  newDirection.rotateUz(oldDirection);

  dir.rotateUz(oldDirection);

  // and of sampling -------------------------------

  fParticleChange->ProposeMomentumDirection(newDirection);

  if (latDisplasment && safety > tlimitminfix) {
    G4double rms = sqrt(2.0*x1);
    G4double dx = zPathLength*(0.5*dirx + invsqrt12*G4RandGauss::shoot(0.0,rms));
    G4double dy = zPathLength*(0.5*diry + invsqrt12*G4RandGauss::shoot(0.0,rms));
    G4double dz;
    G4double d = (dx*dx + dy*dy)/(zPathLength*zPathLength);
    if(d < 0.2)       dz = -zPathLength*d*(1.0 + 0.25*d);
    else if(d >= 1.0) dz = -zPathLength;
    else              dz = -zPathLength*(1.0 - sqrt(1.0 - d));

    G4ThreeVector dv = G4ThreeVector(dx,dy,dz);
    if(isscat) dv.rotateUz(dir);
    pos += dv;
    pos.rotateUz(oldDirection);

    G4double r2 = pos.mag2();
    G4double r  = sqrt(r2);
    /*
    G4cout << " r(mm)= " << r << " safety= " << safety
           << " trueStep(mm)= " << tPathLength
           << " geomStep(mm)= " << zPathLength
           << G4endl;
    */
    //    G4ThreeVector displacement(rx,ry,rz);
    // displacement.rotateUz(oldDirection);

    G4ThreeVector Position = *(fParticleChange->GetProposedPosition());
    G4double fac= 1.;
    if(r > safety) {
      //  ******* so safety is computed at boundary too ************
      G4double newsafety = safetyHelper->ComputeSafety(Position);
      if(r > newsafety) fac = newsafety/r ;
    }  

    if(fac > 0.) {
      // compute new endpoint of the Step
      //      G4ThreeVector newPosition = Position + fac*displacement;
      G4ThreeVector newPosition = Position + fac*pos;

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

G4double G4WentzelVIModel::ComputeXSectionPerVolume()
{
  const G4ElementVector* theElementVector = 
    currentMaterial->GetElementVector();
  const G4double* theAtomNumDensityVector = 
    currentMaterial->GetVecNbOfAtomsPerVolume();
  size_t nelm = currentMaterial->GetNumberOfElements();

  xsece2 = 0.0;
  xsecn2 = 0.0;

  G4double fac = coeff*chargeSquare*invbeta2/mom2;

  for (size_t i=0; i<nelm; i++) {
    const G4Element* elm = (*theElementVector)[i];
    G4double Z = elm->GetZ();
    G4double A = elm->GetN();
    SetupTarget(Z, A, tkin);
    G4double den = fac*theAtomNumDensityVector[i]*Z;

    G4double x1 = 1.0 - cosThetaMin + screenZ;
    /*    
    G4cout<<"ComputeXS: Z= "<<Z<<" den= "<<den<<" cosEl= "<< cosTetMaxElec 
          << " cosThtaMin= " << cosThetaMin << G4endl;
    G4cout<<"cosTetMaxNuc= "<<cosTetMaxNuc<<" cosTetMaxHad= "<<cosTetMaxHad<<G4endl;
    */
    // scattering off electrons
    if(cosTetMaxElec < cosThetaMin) {

      // Reserford part
      G4double z1 = 1.0 - cosTetMaxElec + screenZ;
      G4double z2 = (cosThetaMin - cosTetMaxElec)/x1; 
      xsece2  += den*z2/z1;
    }
    den *= Z;

    // scattering off nucleaus
    G4double costm = std::max(cosTetMaxNuc,cosTetMaxHad);
    if(costm < cosThetaMin) {

      // Reserford part
      G4double s  = screenZ*formfactA;
      G4double w  = 1.0 + 2.0*s;
      G4double z1 = 1.0 - costm + screenZ;
      G4double d  = (1.0 - s)/formfactA;
      G4double x4 = x1 + d;
      G4double z4 = z1 + d;

      xsecn2  += den*w*((cosThetaMin - costm)*(1.0/(x1*z1) + 1.0/(x4*z4)) 
		       - 2.0*log(z1*x4/(x1*z4))/d);
    }
    xsece[i] = xsece2;
    xsecn[i] = xsecn2;
    //G4cout << i << "  xsece2= " << xsece2 << "  xsecn2= " << xsecn2 << G4endl;
  }
  G4double xsec = xsece2 + xsecn2;
  
  //G4cout << "ComputeXS result:  xsec= " << xsec << G4endl; 

  return xsec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WentzelVIModel::ComputeMaxElectronScattering(G4double cutEnergy)
{
  ecut = cutEnergy;
  G4double tmax = tkin;
  if(particle == theElectron) tmax *= 0.5;
  else if(particle != thePositron) {
    G4double ratio = electron_mass_c2/mass;
    G4double tau = tkin/mass;
    tmax = 2.0*electron_mass_c2*tau*(tau + 2.)/
      (1.0 + 2.0*ratio*(tau + 1.0) + ratio*ratio); 
  }
  cosTetMaxElec = cosTetMaxNuc;
  G4double t = std::min(cutEnergy, tmax);
  G4double mom21 = t*(t + 2.0*electron_mass_c2);
  G4double t1 = tkin - t;
  if(t1 > 0.0) {
    G4double mom22 = t1*(t1 + 2.0*mass);
    G4double ctm = (mom2 + mom22 - mom21)*0.5/sqrt(mom2*mom22);
    if(ctm > cosTetMaxNuc && ctm < 1.0) cosTetMaxElec = ctm;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WentzelVIModel::SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                     const G4MaterialCutsCouple*,
				     const G4DynamicParticle*,
				     G4double,
				     G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
