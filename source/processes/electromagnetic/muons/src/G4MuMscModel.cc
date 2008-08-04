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
// $Id: G4MuMscModel.cc,v 1.26 2008-08-04 09:07:23 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4MuMscModel
//
// Author:      V.Ivanchenko on base of Laszlo Urban model
//
// Creation date: 26.10.2007
//
// Modifications:
//
// 22-02-08 Introduce sampling of the tail function (V.Ivanchenko)
//

// Class Description:
//
// Implementation of the model of multiple scattering based on
// H.W.Lewis Phys Rev 78 (1950) 526; 
// J.M. Fernandez-Varea et al., NIM B73 (1993) 447;
// G.Wentzel, Z. Phys. 40 (1927) 590.

// -------------------------------------------------------------------
//


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MuMscModel.hh"
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

G4MuMscModel::G4MuMscModel(const G4String& nam) :
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
  q2Limit(TeV*TeV),
  alpha2(fine_structure_const*fine_structure_const),
  isInitialized(false),
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
  tkin = targetZ = mom2 = DBL_MIN;
  ecut = etag = DBL_MAX;
  particle = 0;
  for(size_t j=0; j<100; j++) {
    FF[j]    = 0.0;
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MuMscModel::~G4MuMscModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuMscModel::Initialise(const G4ParticleDefinition* p,
			      const G4DataVector& cuts)
{
  SetupParticle(p);
  tkin = targetZ = mom2 = DBL_MIN;
  ecut = etag = DBL_MAX;
  xSection = currentRange = 0.0;
  cosThetaMax = cos(PolarAngleLimit());

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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuMscModel::ComputeCrossSectionPerAtom( 
                             const G4ParticleDefinition* p,
			     G4double kinEnergy,
			     G4double Z, G4double,
			     G4double cutEnergy, G4double)
{
  if(p == particle && kinEnergy == tkin && Z == targetZ &&
     cutEnergy == ecut) return xSection;
  xSection = 0.0;
  SetupParticle(p);
  G4double ekin = std::max(keV, kinEnergy);
  SetupKinematic(ekin, cutEnergy);
  SetupTarget(Z, ekin);

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
  G4cout << "G4MuMscModel:XS per A " << " Z= " << Z << " e(MeV)= " << kinEnergy/MeV 
	 << " cut(MeV)= " << ecut/MeV  
  	 << " zmaxE= " << (1.0 - cosTetMaxElec)/screenZ 
	 << " zmaxN= " << (1.0 - cosTetMsxNuc)/screenZ << G4endl;
  */

  // scattering off nucleus
  if(cosTetMaxNuc < 1.0) {
    x  = 1.0 - cosTetMaxNuc;
    x1 = screenZ*formfactA;
    x2 = 1.0/(1.0 - x1); 
    x3 = x/screenZ;
    x4 = formfactA*x;
    if(x3 < numlimit && x1 < numlimit) {
      y = 0.5*x3*x3*x2*x2*x2*(1.0 - 1.333333*x3 + 1.5*x3*x3 
			      - 1.5*x1 + 3.0*x1*x1 + 2.666666*x3*x1);
    } else {
      y  = ((1.0 + x1)*x2*log((1. + x3)/(1. + x4)) 
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
  preKinEnergy = dp->GetKineticEnergy();
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

  //  G4cout << "G4MuMscModel::ComputeTruePathLengthLimit tlimit= " 
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

G4double G4MuMscModel::ComputeGeomPathLength(G4double truelength)
{
  tPathLength  = truelength;
  zPathLength  = tPathLength;
  lambdaeff    = lambda0;

  if(lambda0 > 0.0) {
    G4double tau = tPathLength/lambda0;
    //G4cout << "ComputeGeomPathLength: tLength= " << tPathLength
    //	 << " lambda0= " << lambda0 << " tau= " << tau << G4endl; 
    // small step
    //    par1 = -1.;  
    //    par2 = par3 = 0.;  
    if(tau < numlimit) {
      zPathLength *= (1.0 - 0.5*tau + tau*tau/6.0);

      // medium step
    } else {
      G4double e1 = theManager->GetEnergy(particle,
					  currentRange-tPathLength,
					  currentCouple);
      lambdaeff = GetLambda(0.5*(e1 + preKinEnergy));
      zPathLength = lambdaeff*(1.0 - exp(-tPathLength/lambdaeff));
    }
  }
  /*
  } else if(tPathLength < currentRange*dtrl) {
    zPathLength = lambda0*(1.0 - exp(-tau));

  } else if(tkin < mass) {

    par1 = 1./currentRange;
    par2 = 1./(par1*lambda0);
    par3 = 1.+ par2;
    lambdaeff = 1.0/(par1*par3);
    G4double x = tPathLength/currentRange;
    G4double x1;
    if(x < numlimit) x1 = x*(1.0  + 0.5*x + x*x/3.0);
    else             x1 = log(1.0 - x); 

    zPathLength = lambdaeff*(1.-exp(par3*x1));

  } else {

    G4double T1 = theManager->GetEnergy(particle,
					currentRange-tPathLength,
					currentCouple);
    G4double lambda1 = GetLambda(T1);

    par1 = (lambda0-lambda1)/(lambda0*tPathLength) ;
    par2 = 1./(par1*lambda0) ;
    par3 = 1.+ par2 ;
    lambdaeff = 1.0/(par1*par3);
    zPathLength = lambdaeff*(1.-exp(par3*log(lambda1/lambda0)));
  }

  //  if(zPathLength > lambda0) zPathLength = lambda0;
  if(zPathLength > tPathLength) zPathLength = tPathLength;
  */
  return zPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MuMscModel::ComputeTrueStepLength(G4double geomStepLength)
{
  // step defined other than transportation 
  if(geomStepLength == zPathLength) return tPathLength;

  tPathLength  = geomStepLength;
  zPathLength  = geomStepLength;
  G4double tau = zPathLength/lambdaeff;
  if(tau < numlimit) {
    tPathLength *= (1.0 + 0.5*tau + tau*tau/3.0); 

  } else {
    G4double e1 = theManager->GetEnergy(particle,
					currentRange-zPathLength,
					currentCouple);
    lambdaeff = GetLambda(0.5*(e1 + preKinEnergy));
    tau = zPathLength/lambdaeff;

    if(tau < 0.999999) tPathLength = -lambdaeff*log(1.0 - tau); 
    else               tPathLength = currentRange;

    if(tPathLength < zPathLength) tPathLength = zPathLength;
  }
  if(tPathLength > currentRange) tPathLength = currentRange;
  /*
  } else {
    G4double x = par1*par3*geomStepLength;
    if(x > 1.0) tPathLength = currentRange;
    else { 
      G4double x1;
      if(x < numlimit) x1 = x*(1.+ 0.5*x + x*x/3.0)/par3;
      else if(x < 1.0) x1 = log(1.- x)/par3;

      if(x1 < numlimit)tPathLength = x1*(1.0 - 0.5*x1 + x1*x1/6.0)/par1;
      else             tPathLength = (1.- exp(- x1))/par1;
    }      
  }

  if(tPathLength < zPathLength) tPathLength = zPathLength;
  */
  return tPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuMscModel::SampleScattering(const G4DynamicParticle* dynParticle,
				    G4double safety)
{
  G4double kinEnergy = dynParticle->GetKineticEnergy();
  if(kinEnergy <= DBL_MIN || tPathLength <= DBL_MIN) return;
  
  G4double x1 = 0.5*tPathLength/lambdaeff;

  G4double e  = 0.5*(preKinEnergy + kinEnergy);

  const G4ParticleDefinition* part = dynParticle->GetDefinition();
  G4double cut = (*currentCuts)[currentMaterialIndex];
  SetupParticle(part);
  G4double ekin = std::max(keV, e);
  SetupKinematic(ekin, cut);

  // result of sampling
  G4double z;

  // Gaussian part for combined algorithm ---------

  // define threshold angle as 2 sigma of central value
  cosThetaMin = 1.0 - 4.0*x1;
  /*
  G4cout << "cosTmin= " << cosThetaMin << " cosTmax= " 
	 << cosThetaMax << G4endl;
  */
  
  if(cosThetaMin < cosTetMaxNuc) cosThetaMin = cosTetMaxNuc;

  // The compute cross section for the tail distribution and 
  // correction to the width of the central part 

  G4double xsec = ComputeXSectionPerVolume();
  if(zcorr < x1) x1 -= zcorr;

  if(x1 > 0.1) x1 /= (1.0 - exp(-1.0/x1));
 
  /*  
  G4cout << part->GetParticleName() << " e= " << tkin << "  x1= " 
  	 << x1 << "  zcorr= " << zcorr << G4endl;
  */
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
  G4cout << "G4MuMscModel: step(mm)= " << tPathLength/mm
	 << " e= " << e << " Epre= " << preKinEnergy
  	 << " lambdaeff= " << lambdaeff
	 << " sint= " << sint << " cost= " << cost
	 << " xsec= " << xsec
  	 << " x1= " << x1 << G4endl;
  */
  G4ThreeVector oldDirection = dynParticle->GetMomentumDirection();
  G4ThreeVector newDirection(dirx,diry,cost);
  newDirection.rotateUz(oldDirection);
  G4double rx = 0.0;
  G4double ry = 0.0;
  G4double rz = 0.0;

  //G4double xsec = -1.0;

  // sample scattering for large angle -------------
  if(xsec > DBL_MIN) {
    G4double t = tPathLength;
    do{
      t += log(G4UniformRand())/xsec;
      if(t > 0.0) {
	G4double zz1 = 1.0;
	G4double qsec = G4UniformRand()*xsec;

	// uniform part
	if(qsec < xsece1) {
	  zz1 -= G4UniformRand()*(1.0 - cosThetaMin);
	  //G4cout << "Uniform zz1= " << zz1 << G4endl;

	  // Reserford part
	} else {
	  qsec -= (xsece1 + xsece2);
	  G4double costm = cosTetMaxNuc;
	  G4double formf = 0.0;
	  size_t nelm = currentMaterial->GetNumberOfElements();
	  const G4ElementVector* theElementVector = 
	    currentMaterial->GetElementVector();

	  // scattering off nucleus
	  if(qsec >= 0.0) {
	    for (size_t i=0; i<nelm; i++) {
	      if(xsecn[i] > qsec || nelm-1 == i) {
		const G4Element* elm = (*theElementVector)[i];
		SetupTarget(elm->GetZ(), tkin);
		formf = formfactA;
		//G4cout << "Reserford nuc Z= " << elm->GetZ() << G4endl;
		break;
	      }
	    }
	    // scattering off electrons
	  } else {
	    qsec += xsece2;
	    costm = cosTetMaxElec;
	    for (size_t i=0; i<nelm; i++) {
	      if(xsece[i] > qsec || nelm-1 == i) {
		const G4Element* elm = (*theElementVector)[i];
		SetupTarget(elm->GetZ(), tkin);
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
	  //G4cout << "Reserford zz1= " << zz1 << G4endl;
	}
	sint = sqrt((1.0 - zz1)*(1.0 + zz1));
	//G4cout << "sint= " << sint << G4endl;
	phi  = twopi*G4UniformRand();

	G4double vx1 = sint*cos(phi);
	G4double vy1 = sint*sin(phi);
	G4double zr = zPathLength*t/tPathLength;
	rx += vx1*zr; 
	ry += vy1*zr;
	rz -= (1.0 - zz1)*zr;
	G4ThreeVector newDirection1(vx1,vy1,zz1);
	newDirection1.rotateUz(newDirection);
	newDirection = newDirection1;
      }
    } while (t > 0.0);
  }

  // and of sampling -------------------------------

  fParticleChange->ProposeMomentumDirection(newDirection);

  if (latDisplasment && safety > tlimitminfix) {
    G4double rms = sqrt(2.0*x1);
    G4double dx = zPathLength*(0.5*dirx + invsqrt12*G4RandGauss::shoot(0.0,rms));
    G4double dy = zPathLength*(0.5*diry + invsqrt12*G4RandGauss::shoot(0.0,rms));
    rx += dx;
    ry += dy;
    G4double dz  = (dx*dx + dy*dy)/(zPathLength*zPathLength);
    if(dz < 0.2)      rz -= zPathLength*dz*(1.0 + 0.25*dz);
    else if(dz > 1.0) rz  = -zPathLength;
    else              rz -= zPathLength*(1.0 - sqrt(1.0 - dz));

    G4double r2 = rx*rx + ry*ry + rz*rz;
    G4double r  = sqrt(r2);
    /*
    G4cout << " r(mm)= " << r << " safety= " << safety
           << " trueStep(mm)= " << tPathLength
           << " geomStep(mm)= " << zPathLength
           << G4endl;
    */
    G4ThreeVector displacement(rx,ry,rz);
    displacement.rotateUz(oldDirection);

    G4ThreeVector Position = *(fParticleChange->GetProposedPosition());
    G4double fac= 1.;
    if(r > safety) {
      //  ******* so safety is computed at boundary too ************
      G4double newsafety = safetyHelper->ComputeSafety(Position);
      if(r > newsafety) fac = newsafety/r ;
    }  

    if(fac > 0.) {
      // compute new endpoint of the Step
      G4ThreeVector newPosition = Position + fac*displacement;

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

G4double G4MuMscModel::ComputeXSectionPerVolume()
{
  const G4ElementVector* theElementVector = 
    currentMaterial->GetElementVector();
  const G4double* theAtomNumDensityVector = 
    currentMaterial->GetVecNbOfAtomsPerVolume();
  size_t nelm = currentMaterial->GetNumberOfElements();

  xsece1 = 0.0;
  xsece2 = 0.0;
  xsecn2 = 0.0;
  zcorr  = 0.0;

  G4double fac = coeff*chargeSquare*invbeta2/mom2;

  for (size_t i=0; i<nelm; i++) {
    const G4Element* elm = (*theElementVector)[i];
    G4double Z = elm->GetZ();
    SetupTarget(Z, tkin);
    G4double den = fac*theAtomNumDensityVector[i]*Z;

    G4double x  = 1.0 - cosThetaMin;
    G4double x1 = x + screenZ;
    G4double x2 = 1.0/(x1*x1);
    G4double x3 = 1.0 + x*formfactA;
    /*
    G4cout << "x= " << x << " den= " << den << " cosE= " << cosTetMaxElec << G4endl;
    G4cout << "cosThtaMin= " << cosThetaMin << G4endl;
    G4cout << "cosTetMaxNuc= " << cosTetMaxNuc << " q2Limit= " << q2Limit << G4endl;
    */
    // scattering off electrons
    if(cosTetMaxElec < cosThetaMin) {

      // flat part
      G4double s = den*x2*x;
      xsece1 += s;
      zcorr  += 0.5*x*s;

      // Reserford part
      G4double z1 = 1.0 - cosTetMaxElec + screenZ;
      G4double z2 = (cosThetaMin - cosTetMaxElec)/x1; 
      if(z2 < 0.2) s = z2*(x - 0.5*z2*(x - screenZ))/x1;
      else         s = log(1.0 + z2)  - screenZ*z2/z1;
      xsece2  += den*z2/z1;
      zcorr   += den*s;
    }
    den *= Z;

    //G4cout << "Z= " << Z<< " cosL= " << cosTetMaxNuc << " cosMin= " << cosThetaMin << G4endl;
    // scattering off nucleaus
    if(cosTetMaxNuc < cosThetaMin) {

      // flat part
      G4double s = den*x2*x/(x3*x3);
      xsece1 += s;
      zcorr  += 0.5*x*s;

      // Reserford part
      s  = screenZ*formfactA;
      G4double w  = 1.0 + 2.0*s;
      G4double z1 = 1.0 - cosTetMaxNuc + screenZ;
      G4double d  = (1.0 - s)/formfactA;
      G4double x4 = x1 + d;
      G4double z4 = z1 + d;
      G4double t1 = 1.0/(x1*z1);
      G4double t4 = 1.0/(x4*z4);
      G4double w1 = cosThetaMin - cosTetMaxNuc;
      G4double w2 = log(z1*x4/(x1*z4));

      den *= w;     
      xsecn2  += den*(w1*(t1 + t4) - 2.0*w2/d);
      zcorr   += den*(w*w2 - w1*(screenZ*t1 + t4/formfactA));
    }
    xsece[i] = xsece2;
    xsecn[i] = xsecn2;
    //    G4cout << i << "  xsece2= " << xsece2 << "  xsecn2= " << xsecn2 << G4endl;
  }
  G4double xsec = xsece1 + xsece2 + xsecn2;
  /*
    G4cout << "xsece1= " << xsece1 << "  xsece2= " << xsece2 
    << "  xsecn2= " << xsecn2 
	 << " zsec= " << zcorr*0.5*tPathLength << G4endl;
  */
  zcorr *= 0.5*tPathLength;

  return xsec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuMscModel::ComputeMaxElectronScattering(G4double cutEnergy)
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
    if(ctm > cosTetMaxNuc)  cosTetMaxElec = ctm;
    if(cosTetMaxElec > 1.0) cosTetMaxElec = 1.0;
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
