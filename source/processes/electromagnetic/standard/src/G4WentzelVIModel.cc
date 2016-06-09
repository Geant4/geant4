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
// $Id: G4WentzelVIModel.cc,v 1.16 2008/11/19 11:47:50 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
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

G4WentzelVIModel::G4WentzelVIModel(const G4String& nam) :
  G4VMscModel(nam),
  theLambdaTable(0),
  theLambda2Table(0),
  numlimit(0.2),
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
  nelments = 5;
  xsecn.resize(nelments);
  prob.resize(nelments);
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
  // reset parameters
  SetupParticle(p);
  tkin = targetZ = mom2 = DBL_MIN;
  ecut = etag = DBL_MAX;
  currentRange = 0.0;
  cosThetaMax = cos(PolarAngleLimit());
  currentCuts = &cuts;

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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIModel::ComputeCrossSectionPerAtom( 
                             const G4ParticleDefinition* p,
			     G4double kinEnergy,
			     G4double Z, G4double,
			     G4double cutEnergy, G4double)
{
  SetupParticle(p);
  G4double ekin = std::max(lowEnergyLimit, kinEnergy);
  SetupKinematic(ekin, cutEnergy);
  SetupTarget(Z, ekin);
  G4double xsec = ComputeTransportXSectionPerVolume();
  /*  
  G4cout << "CS: e= " << tkin << " cosEl= " << cosTetMaxElec2 
	 << " cosN= " << cosTetMaxNuc2 << " xsec(bn)= " << xsec/barn
	 << " " << particle->GetParticleName() << G4endl;
  */
  return xsec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIModel::ComputeTransportXSectionPerVolume()
{
  G4double xSection = 0.0;
  G4double x, y, x1, x2, x3, x4;

  // scattering off electrons
  if(cosTetMaxElec2 < 1.0) {
    x = (1.0 - cosTetMaxElec2)/screenZ;
    if(x < numlimit) y = 0.5*x*x*(1.0 - 1.3333333*x + 1.5*x*x); 
    else             y = log(1.0 + x) - x/(1.0 + x);
    if(y < 0.0) {
      nwarnings++;
      if(nwarnings < nwarnlimit /*&& y < -1.e-10*/) {
	G4cout << "Electron scattering <0 for L1 " << y 
	       << " e(MeV)= " << tkin << " p(MeV/c)= " << sqrt(mom2) 
	       << " Z= " << targetZ << "  " 
	       << particle->GetParticleName() << G4endl;
	G4cout << " z= " << 1.0-cosTetMaxElec2 << " screenZ= " << screenZ 
	       << " x= " << x << G4endl;
      }
      y = 0.0;
    }
    xSection += y/targetZ;
  }
  /*
  G4cout << "G4WentzelVIModel:XS per A " << " Z= " << Z << " e(MeV)= " << kinEnergy/MeV 
	 << " cut(MeV)= " << ecut/MeV  
  	 << " zmaxE= " << (1.0 - cosTetMaxElec)/screenZ 
	 << " zmaxN= " << (1.0 - cosTetMsxNuc)/screenZ << G4endl;
  */

  // scattering off nucleus
  if(cosTetMaxNuc2 < 1.0) {
    x  = 1.0 - cosTetMaxNuc2;
    x1 = screenZ*formfactA;
    x2 = 1.0/(1.0 - x1); 
    x3 = x/screenZ;
    x4 = formfactA*x;
    // low-energy limit
    if(x3 < numlimit && x1 < numlimit) {
      y = 0.5*x3*x3*x2*x2*x2*(1.0 - 1.333333*x3 + 1.5*x3*x3 - 1.5*x1
      			      + 3.0*x1*x1 + 2.666666*x3*x1);
      // high energy limit
    } else if(1.0 < x1) {
      x4 = x1*(1.0 + x3);
      y  = x3*(1.0 + 0.5*x3 - (2.0 - x1)*(1.0 + x3 + x3*x3/3.0)/x4)/(x4*x4);
      // middle energy 
    } else {
      y = ((1.0 + x1)*x2*log((1. + x3)/(1. + x4)) 
	   - x3/(1. + x3) - x4/(1. + x4))*x2*x2; 
    }
    if(y < 0.0) {
      nwarnings++;
      if(nwarnings < nwarnlimit /*&& y < -1.e-10*/) { 
	G4cout << "Nuclear scattering <0 for L1 " << y 
	       << " e(MeV)= " << tkin << " Z= " << targetZ << "  " 
	       << particle->GetParticleName() << G4endl;
	G4cout << " formfactA= " << formfactA << " screenZ= " << screenZ 
	       << " x= " << " x1= " << x1 << " x2= " << x2 
	       << " x3= " << x3 << " x4= " << x4 <<G4endl;
      }
      y = 0.0;
    }
    xSection += y; 
  }
  xSection *= (coeff*targetZ*targetZ*chargeSquare*invbeta2/mom2); 
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
  G4StepPoint* sp = track.GetStep()->GetPreStepPoint();
  G4StepStatus stepStatus = sp->GetStepStatus();

  // initialisation for 1st step  
  if(stepStatus == fUndefined) {
    inside = false;
    SetupParticle(dp->GetDefinition());
    theLambdaTable = theTable;
  }

  // initialisation for each step, lambda may be computed from scratch
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
  G4double presafety = sp->GetSafety();

  // compute presafety again if presafety <= 0 and no boundary
  // i.e. when it is needed for optimization purposes
  if(stepStatus != fGeomBoundary && presafety < tlimitminfix) 
    presafety = safetyHelper->ComputeSafety(sp->GetPosition()); 
  /*
  G4cout << "G4WentzelVIModel::ComputeTruePathLengthLimit tlimit= " 
 	 <<tlimit<<" safety= " << presafety
  	 << " range= " <<currentRange<<G4endl;
  */
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
      //      zPathLength = lambda0*(1.0 - exp(-tPathLength/lambda0));
      G4double e1 = 0.0;
      if(currentRange > tPathLength) {
	e1 = theManager->GetEnergy(particle,
				   currentRange-tPathLength,
				   currentCouple);
      }
      lambdaeff = GetLambda(0.5*(e1 + preKinEnergy));
      zPathLength = lambdaeff*(1.0 - exp(-tPathLength/lambdaeff));
    }
  }
  //G4cout<<"Comp.geom: zLength= "<<zPathLength<<" tLength= "<<tPathLength<<G4endl;
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
  tPathLength *= (1.0 + 0.5*tau + tau*tau/3.0); 

  if(tau > numlimit) {
    G4double e1 = 0.0;
    if(currentRange > tPathLength) {
      e1 = theManager->GetEnergy(particle,
				 currentRange-tPathLength,
				 currentCouple);
    }
    lambdaeff = GetLambda(0.5*(e1 + preKinEnergy));
    tau = zPathLength/lambdaeff;

    if(tau < 0.999999) tPathLength = -lambdaeff*log(1.0 - tau); 
    else               tPathLength = currentRange;

    if(tPathLength < zPathLength) tPathLength = zPathLength;
  }
  if(tPathLength > currentRange) tPathLength = currentRange;
  //G4cout<<"Comp.true: zLength= "<<zPathLength<<" tLength= "<<tPathLength<<G4endl;
  return tPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WentzelVIModel::SampleScattering(const G4DynamicParticle* dynParticle,
					G4double safety)
{
  //G4cout << "!##! G4WentzelVIModel::SampleScattering for " 
  //	 << particle->GetParticleName() << G4endl;
  G4double kinEnergy = dynParticle->GetKineticEnergy();
  if(kinEnergy <= DBL_MIN || tPathLength <= DBL_MIN) return;

  G4double ekin = preKinEnergy;
  if(ekin - kinEnergy > ekin*dtrl) {
    ekin = 0.5*(preKinEnergy + kinEnergy);
    lambdaeff = GetLambda(ekin);
  }  
  
  G4double x1 = 0.5*tPathLength/lambdaeff;
  G4double cut= (*currentCuts)[currentMaterialIndex];
  /*  
  G4cout <<"SampleScat: E0(MeV)= "<< preKinEnergy<<" Eeff(MeV)= "<<ekin/MeV
	 << " L0= " << lambda0 << " Leff= " << lambdaeff 
	 << " x1= " << x1 << " safety= " << safety << G4endl;
  */

  G4double xsec = 0.0;
  G4bool largeAng = false;

  // large scattering angle case
  if(x1 > 0.5) {
    x1 *= 0.5;
    largeAng = true;

    // normal case
  } else {

    // define threshold angle as 2 sigma of central value
    cosThetaMin = 1.0 - 3.0*x1;

    // for low-energy e-,e+ no limit
    ekin = std::max(ekin, lowEnergyLimit);
    SetupKinematic(ekin, cut);
  
    // recompute transport cross section
    if(cosThetaMin > cosTetMaxNuc) {

      xsec = ComputeXSectionPerVolume();

      if(xtsec > DBL_MIN) x1 = 0.5*tPathLength*xtsec;
      else                x1 = 0.0;

      /*      
	G4cout << "cosTetMaxNuc= " << cosTetMaxNuc 
	<< " cosThetaMin= " << cosThetaMin 
	<< " cosThetaMax= " << cosThetaMax 
	<< " cosTetMaxElec2= " << cosTetMaxElec2 << G4endl;
	G4cout << "Recomputed xsec(1/mm)= " << xsec << " x1= " << x1 << G4endl;
      */
    }
  }

  // result of central part sampling 
  G4double z; 
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
  
  //G4cout << "G4WentzelVIModel: step(mm)= " << tPathLength/mm
  //	 << " sint= " << sint << " cost= " << cost<< G4endl;
  
  G4ThreeVector oldDirection = dynParticle->GetMomentumDirection();
  G4ThreeVector newDirection(dirx,diry,cost);
  G4ThreeVector temp(0.0,0.0,1.0);
  G4ThreeVector pos(0.0,0.0,-zPathLength);
  G4ThreeVector dir(0.0,0.0,1.0);
  G4bool isscat = false;

  // sample MSC scattering for large angle
  // extra central scattering for holf step
  if(largeAng) {
    isscat = true;
    pos.setZ(-0.5*zPathLength);
    do {
      z = -x1*log(G4UniformRand());
    } while (z > 1.0); 
    cost = 1.0 - 2.0*z;
    if(std::abs(cost) > 1.0) cost = 1.0;

    sint = sqrt((1.0 - cost)*(1.0 + cost));
    phi  = twopi*G4UniformRand();

    // position and direction for secondary scattering
    dir.set(sint*cos(phi),sint*sin(phi),cost);
    pos += 0.5*dir*zPathLength;
    x1 *= 2.0;
  }

  // sample Reserford scattering for large angle
  if(xsec > DBL_MIN) {
    G4double t = tPathLength;
    G4int nelm = currentMaterial->GetNumberOfElements();
    const G4ElementVector* theElementVector = 
      currentMaterial->GetElementVector();
    do{
      G4double x  = -log(G4UniformRand())/xsec;      
      pos += dir*(zPathLength*std::min(x,t)/tPathLength);
      t -= x;
      if(t > 0.0) {
	G4double zz1 = 1.0;
	G4double qsec = G4UniformRand()*xsec;

	// scattering off nucleus
        G4int i = 0;
	if(nelm > 1) {
	  for (; i<nelm; i++) {if(xsecn[i] >= qsec) break;}
	  if(i >= nelm) i = nelm - 1;
	}
	SetupTarget((*theElementVector)[i]->GetZ(), tkin);
        G4double formf = formfactA;
        G4double costm = cosTetMaxNuc2;
        if(prob[i] > 0.0) {
	  if(G4UniformRand() <= prob[i]) {
	    formf = 0.0;
	    costm = cosTetMaxElec2;
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
	  //G4cout << "Reserford zz1= " << zz1 << " t= " << t << G4endl;
	  sint = sqrt((1.0 - zz1)*(1.0 + zz1));
	  //G4cout << "sint= " << sint << G4endl;
	  phi  = twopi*G4UniformRand();
	  G4double vx1 = sint*cos(phi);
	  G4double vy1 = sint*sin(phi);
	  temp.set(vx1,vy1,zz1);
	  temp.rotateUz(dir);
	  dir = temp;
	}
      }
    } while (t > 0.0); 
  }
  if(isscat) newDirection.rotateUz(dir);
  newDirection.rotateUz(oldDirection);

  //G4cout << "G4WentzelVIModel sampling of scattering is done" << G4endl;
  // end of sampling -------------------------------

  fParticleChange->ProposeMomentumDirection(newDirection);

  if (latDisplasment && safety > tlimitminfix) {
    G4double rms = invsqrt12*sqrt(2.0*x1);
    G4double dx = zPathLength*(0.5*dirx + rms*G4RandGauss::shoot(0.0,1.0));
    G4double dy = zPathLength*(0.5*diry + rms*G4RandGauss::shoot(0.0,1.0));
    G4double dz;
    G4double d = (dx*dx + dy*dy)/(zPathLength*zPathLength);
    if(d < numlimit)  dz = -0.5*zPathLength*d*(1.0 + 0.25*d);
    else if(d < 1.0)  dz = -zPathLength*(1.0 - sqrt(1.0 - d));
    else {
      dx = dy = dz = 0.0;
    }

    temp.set(dx,dy,dz);
    if(isscat) temp.rotateUz(dir);
    pos += temp;
   
    pos.rotateUz(oldDirection);

    G4double r = pos.mag();

    /*    
    G4cout << " r(mm)= " << r << " safety= " << safety
           << " trueStep(mm)= " << tPathLength
           << " geomStep(mm)= " << zPathLength
           << G4endl;
    */

    if(r > tlimitminfix) {
      G4ThreeVector Position = *(fParticleChange->GetProposedPosition());
      G4double fac= 1.;
      if(r >= safety) {
	//  ******* so safety is computed at boundary too ************
	G4double newsafety = 
	  safetyHelper->ComputeSafety(Position) - tlimitminfix;
        if(newsafety <= 0.0) fac = 0.0;
	else if(r > newsafety) fac = newsafety/r ;
        //G4cout << "NewSafety= " << newsafety << " fac= " << fac 
	// << " r= " << r << " sint= " << sint << " pos " << Position << G4endl;
      }  

      if(fac > 0.) {
	// compute new endpoint of the Step
	G4ThreeVector newPosition = Position + fac*pos;

	// check safety after displacement
	G4double postsafety = safetyHelper->ComputeSafety(newPosition);

	// displacement to boundary
	if(postsafety <= 0.0) {
	  safetyHelper->Locate(newPosition, newDirection);

	  // not on the boundary
	} else { 
	  safetyHelper->ReLocateWithinVolume(newPosition);
	  // if(fac < 1.0) G4cout << "NewPosition " << newPosition << G4endl;
	}
     
	fParticleChange->ProposePosition(newPosition);
      } 
    }
  }
  //G4cout << "G4WentzelVIModel::SampleScattering end" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4WentzelVIModel::ComputeXSectionPerVolume()
{
  const G4ElementVector* theElementVector = 
    currentMaterial->GetElementVector();
  const G4double* theAtomNumDensityVector = 
    currentMaterial->GetVecNbOfAtomsPerVolume();
  G4int nelm = currentMaterial->GetNumberOfElements();
  if(nelm > nelments) {
    nelments = nelm;
    xsecn.resize(nelments);
    prob.resize(nelments);
  }

  xtsec = 0.0;
  G4double xs = 0.0;

  G4double fac = coeff*chargeSquare*invbeta2/mom2;

  for (G4int i=0; i<nelm; i++) {
    SetupTarget((*theElementVector)[i]->GetZ(), tkin);
    G4double density = theAtomNumDensityVector[i];
    G4double cosnm = cosTetMaxNuc2;
    G4double cosem = cosTetMaxElec2;

    // recompute the angular limit 
    cosTetMaxNuc2  = std::max(cosnm,cosThetaMin); 
    cosTetMaxElec2 = std::max(cosem,cosThetaMin); 
    xtsec += ComputeTransportXSectionPerVolume()*density;
    // return limit back
    cosTetMaxElec2 = cosem;
    cosTetMaxNuc2  = cosnm;

    G4double esec = 0.0;
    G4double nsec = 0.0;
    G4double x1 = 1.0 - cosThetaMin + screenZ;
    G4double f  = fac*targetZ*density; 

    // scattering off electrons
    if(cosThetaMin > cosem) {
      esec = f*(cosThetaMin - cosem)/(x1*(1.0 - cosem + screenZ));
    }

    // scattering off nucleaus
    if(cosThetaMin > cosnm) {

      // Reserford part
      G4double s  = screenZ*formfactA;
      G4double z1 = 1.0 - cosnm + screenZ;
      G4double d  = (1.0 - s)/formfactA;

      // check numerical limit
      if(d < numlimit*x1) {
	G4double x2 = x1*x1;
	G4double z2 = z1*z1;
	nsec = (1.0/(x1*x2) - 1.0/(z1*z2) - d*1.5*(1.0/(x2*x2) - 1.0/(z2*z2)))/
	  (3.0*formfactA*formfactA);
      } else {
	G4double x2 = x1 + d;
	G4double z2 = z1 + d;
	nsec = (1.0 + 2.0*s)*((cosThetaMin - cosnm)*(1.0/(x1*z1) + 1.0/(x2*z2)) -
			   2.0*log(z1*x2/(z2*x1))/d);
      }
      nsec *= f*targetZ;
    }
    nsec += esec;
    if(nsec > 0.0) esec /= nsec;
    xs += nsec;
    xsecn[i] = xs;
    prob[i]  = esec;
    //G4cout << i << "  xs= " << xs << " cosThetaMin= " << cosThetaMin 
    //	   << " costm= " << costm << G4endl;
  }
  
  //G4cout << "ComputeXS result:  xsec(1/mm)= " << xs 
  //<< " txsec(1/mm)= " << xtsec <<G4endl; 
  return xs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
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
    
    //G4cout << "x= " << x << " den= " << den << " cosE= " << cosTetMaxElec << G4endl;
    //G4cout << "cosThtaMin= " << cosThetaMin << G4endl;
    //G4cout << "cosTetMaxNuc= " << cosTetMaxNuc << " q2Limit= " << q2Limit << G4endl;
    
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
 
    //G4cout << "xsece1= " << xsece1 << "  xsece2= " << xsece2 
    //<< "  xsecn2= " << xsecn2 
	// << " zsec= " << zcorr*0.5*tPathLength << G4endl;
  zcorr *= 0.5*tPathLength;

  return xsec;
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WentzelVIModel::ComputeMaxElectronScattering(G4double cutEnergy)
{
  ecut = cutEnergy;
  G4double tmax = tkin;
  cosTetMaxElec = 1.0;
  if(mass > MeV) {
    G4double ratio = electron_mass_c2/mass;
    G4double tau = tkin/mass;
    tmax = 2.0*electron_mass_c2*tau*(tau + 2.)/
      (1.0 + 2.0*ratio*(tau + 1.0) + ratio*ratio); 
    cosTetMaxElec = 1.0 - std::min(cutEnergy, tmax)*electron_mass_c2/mom2;
  } else {

    if(particle == theElectron) tmax *= 0.5;
    G4double t = std::min(cutEnergy, tmax);
    G4double mom21 = t*(t + 2.0*electron_mass_c2);
    G4double t1 = tkin - t;
    //G4cout <<"tkin=" <<tkin<<" tmax= "<<tmax<<" t= " 
    //<<t<< " t1= "<<t1<<" cut= "<<ecut<<G4endl;
    if(t1 > 0.0) {
      G4double mom22 = t1*(t1 + 2.0*mass);
      G4double ctm = (mom2 + mom22 - mom21)*0.5/sqrt(mom2*mom22);
      if(ctm < 1.0) cosTetMaxElec = ctm;
    }
  }
  if(cosTetMaxElec < cosTetMaxNuc) cosTetMaxElec = cosTetMaxNuc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4WentzelVIModel::SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                     const G4MaterialCutsCouple*,
				     const G4DynamicParticle*,
				     G4double,
				     G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
