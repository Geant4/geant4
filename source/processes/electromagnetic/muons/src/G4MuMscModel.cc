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
// $Id: G4MuMscModel.cc,v 1.16 2008-02-22 14:38:35 vnivanch Exp $
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
#include "G4Electron.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4TransportationManager.hh"
#include "G4SafetyHelper.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PhysicsLogVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4MuMscModel::G4MuMscModel(G4double frange, 
			   G4double thetaMax, 
			   G4double tMax,  
			   const G4String& nam)
  : G4eCoulombScatteringModel(0.0,thetaMax,false,tMax,nam),
    theLambdaTable(0),
    theLambda2Table(0),
    dtrl(0.05),
    facrange(frange),
    numlimit(0.2),
    lowBinEnergy(keV),
    highBinEnergy(PeV),
    nbins(60),
    nwarnings(0),
    nwarnlimit(50),
    currentCouple(0),
    isInitialized(false),
    buildTables(true),
    newrun(true),
    inside(false)
{
  invsqrt12 = 1./sqrt(12.);
  tlimitminfix = 1.e-6*mm;
  theManager = G4LossTableManager::Instance(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MuMscModel::~G4MuMscModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuMscModel::Initialise(const G4ParticleDefinition* p,
			      const G4DataVector& cuts)
{
  SetupParticle(p);
  newrun = true;
  xSection = currentRange = targetZ = ecut = tkin = 0.0;
  // set values of some data members
  if(!isInitialized) {
    isInitialized = true;
    if(p->GetParticleType() == "nucleus") buildTables = false;

    if (pParticleChange)
      fParticleChange = reinterpret_cast<G4ParticleChangeForMSC*>(pParticleChange);
    else
      fParticleChange = new G4ParticleChangeForMSC();

    safetyHelper = G4TransportationManager::GetTransportationManager()
      ->GetSafetyHelper();
    safetyHelper->InitialiseHelper();
  }
  G4eCoulombScatteringModel::Initialise(p, cuts);
  currentCuts = &cuts;
  if(buildTables)
    theLambda2Table = G4PhysicsTableHelper::PreparePhysicsTable(theLambda2Table);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuMscModel::BuildTables()
{
  //G4cout << "G4MuMscModel::BuildTables flags newrun= " << newrun 
  //	 << "  buildTables= " << buildTables << G4endl;
  newrun = false;
  if(!buildTables) return;

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  G4double e, s, cut;

  for(size_t i=0; i<numOfCouples; i++) {

    if (theLambda2Table->GetFlag(i)) {

      // create physics vector and fill it
      DefineMaterial(theCoupleTable->GetMaterialCutsCouple(i));
      cut = (*currentCuts)[currentMaterialIndex];
      G4PhysicsVector* aVector =
	new G4PhysicsLogVector(lowBinEnergy, highBinEnergy, nbins);
      for(G4int j=0; j<nbins; j++) {
        e = aVector->GetLowEdgeEnergy(j);
        s = ComputeLambda2(e, cut);
        //G4cout << j << "  " << currentCouple->GetMaterial()->GetName() 
	//       << "  e(MeV)= " << e << " cut(MeV)= " << cut
	//       << " L2= " << s << G4endl;  
	aVector->PutValue(j, s);
      }
      
      G4PhysicsTableHelper::SetPhysicsVector(theLambda2Table, i, aVector);
    }
  }
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
    if(ctm > cosTetMaxNuc && ctm <= 1.0) cosTetMaxElec = ctm;
  }

  G4double x, y, x1, x2, x3, x4;
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
	 << " zmaxN= " << (1.0 - cosTetMaxNuc)/screenZ << G4endl;
  */
  if(cosTetMaxNuc < 1.0) {
    x  = 1.0 - cosTetMaxNuc;
    x1 = screenZ*formfactA;
    x2 = 1.0/(1.0 - x1); 
    x3 = x/screenZ;
    x4 = formfactA*x;
    if(x3 < numlimit) 
      y = 0.5*x3*x3*x2*x2*x2*(1.0 - 1.333333*x3 + 1.5*x3*x3 
			      - 1.5*x1 + 3.0*x1*x1 + 2.666666*x3*x1);
    else {
      y  = ((1.0 + x1)*x2*log((1. + x3)/(1. + x4)) - x3/(1. + x3) - x4/(1. + x4))*x2*x2; 
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

G4double G4MuMscModel::ComputeLambda2(G4double kinEnergy,
				      G4double cutEnergy)
{
  G4double res = 0.0;
  G4double ekin = std::max(keV, kinEnergy);

  const G4Material* mat = currentCouple->GetMaterial();
  const G4ElementVector* theElementVector = mat->GetElementVector();
  const G4double* theAtomNumDensityVector = mat->GetVecNbOfAtomsPerVolume();
  size_t nelm = mat->GetNumberOfElements();

  G4double tmax = ekin;
  if(particle == theElectron) tmax *= 0.5;
  else if(particle != thePositron) {
    G4double ratio = electron_mass_c2/mass;
    G4double tau = ekin/mass;
    tmax = 2.0*electron_mass_c2*tau*(tau + 2.)/
      (1.0 + 2.0*ratio*(tau + 1.0) + ratio*ratio); 
  }
  cosTetMaxElec = cosTetMaxNuc;
  G4double t = std::min(cutEnergy, tmax);
  G4double mom0  = ekin*(ekin + mass*2.0);
  G4double mom21 = t*(t + 2.0*electron_mass_c2);
  G4double t1 = ekin - t;
  if(t1 > 0.0) {
    G4double mom22 = t1*(t1 + 2.0*mass);
    G4double ctm = (mom0 + mom22 - mom21)*0.5/sqrt(mom0*mom22);
    if(ctm > cosTetMaxNuc && ctm <= 1.0) cosTetMaxElec = ctm;
  }

  G4double x, x1, x2, x3, x4, x5, y;

  for (size_t i=0; i<nelm; i++) {
    const G4Element* elm = (*theElementVector)[i];
    G4double Z = elm->GetZ();
    SetupTarget(Z, elm->GetN(), ekin);
    G4double s = 0.0;  
    if(cosTetMaxElec < 1.0) {
      x = 1.0 - cosTetMaxElec;
      x1= x/screenZ;
      if(x1 < numlimit) {
	y = x*x1*x1*(1.0 - 1.5*x1 + 1.8*x1*x1)/3.0;
      } else {
	y = screenZ*(x1*(2.0 + x1)/(1.0 + x1) - 2.0*log(1.0 + x1));
      }
      if(y < 0.0) {
        nwarnings++;
        if(nwarnings < nwarnlimit) {
	  G4cout << "Electron scattering <0 for L2 " << y 
		 << " e(MeV)= " << tkin << " p(MeV/c)= " << sqrt(mom2)
		 << " Z= " << Z << "  " 
		 << particle->GetParticleName() << G4endl;
          G4cout << " z= " << " screenZ= " << screenZ 
		 << " x= " << x << " x1= " << x1 << G4endl;
	}
        y = 0.0;
      }
      s = y/Z;
    }
    //  G4cout << "cut= " << cut << " e= " << tkin << " croosE= " 
    //  << xSection/barn << G4endl;

    // limit main integral because of nuclear size effect

    if(cosTetMaxNuc < 1.0) {
      // Compute second moment
      x   = 1.0 - cosTetMaxNuc;
      x1  = screenZ*formfactA;
      x2  = 1.0/(1.0 - x1); 
      x3  = x/screenZ;
      x4  = 1.0 + x*formfactA;
      x5  = 1.0 + x3;
      if(x3 < numlimit) {
        G4double x32 = x3*x3;
	y = x*x32*x2*x2*x2*(1.0 - 1.5*x3 - 3.0*x1 + 
			      3.0*x3*x1 + 3.0*x1*x1 + 1.8*x32)/3.0;
      } else {
	y = screenZ*(x3/x5 + x3/x4 - 2.0*x2*log(x5/x4))*x2*x2; 
      }
      if(y < 0.0) {
        nwarnings++;
        if(nwarnings < nwarnlimit) { 
	  G4cout << "Nuclear scattering <0 for L2 " << y 
		 << " e(MeV)= " << tkin << " Z= " << Z << "  " 
		 << particle->GetParticleName() << G4endl;
          G4cout << " formfactA= " << formfactA << " screenZ= " << screenZ 
		 << " x= " << x << " x1= " << x1 << " x2= " << x2 
		 << " x3-1= " << x3-1.0 << " x4-1.0= " << x4-1.0 <<G4endl;
	}
        y = 0.0;
      }
      s += y;
    }
    res += s*Z*Z*theAtomNumDensityVector[i]; 
  }
  res  *= coeff*chargeSquare*invbeta2/mom2;
  //  G4cout << "L2 e= " << ekin << " sig2= " << res
  //	 << " " << particle->GetParticleName() << " in " << mat->GetName() << G4endl;
  return res; 
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
    if(newrun) BuildTables();
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
  //  G4double t1 = zPathLength/lambdaeff;

  G4double e  = 0.5*(preKinEnergy + kinEnergy);
  /*
  G4double t2 = zPathLength*GetLambda2(e); 
  G4double x2;  
  if(t1 < numlimit && t2 < numlimit) {
    x2 = 0.25*t2*(1.0 + t1*(2.0*t1/t2 + 0.75*t2/t1 - 3.0));
  } else {
    x2 = x1 - (1.0 - exp(-3.0*(t1 - 0.5*t2)))/6.0;
  }
  */
  // result of sampling
  G4double z;

  // Gaussian distribution ---------------------
  /*
  G4double x  = G4UniformRand();
  G4double z1 = 1.0/x1;
  G4double z2;
  if(z1 < numlimit) z2 = x + (1.0 - x)*(1.0 - z1 + z1*z1*0.5);
  else              z2 = x + (1.0 - x)*exp(-z1);
  z = -x1*log(z2);
  */
  // Penelope algorithm ------------------------
  /*
  G4double y1  = 1.0 - x1;
  G4double x3  = (x2 - x1*x1)/(x1*y1);
  //G4cout << "t1= "<<t1<<" t2= "<<t2<<" x1= "<<x1<<" x2= "<<x2<<" x3= "<<x3<<G4endl;

  if(x3 <= 0.0 || x3 >= 0.33) {
    nwarnings++;
    if(nwarnings < nwarnlimit) 
      G4cout << "G4MuMscModel::SampleScattering: ePre(MeV)= " << preKinEnergy/MeV
	     << " ePost(MeV)= " << kinEnergy/MeV
	     << " <x>= " << x1 << " x2= " << x2
	     << " x3= " << x3
	     << G4endl;
    x3 = std::min(1.0/y1,0.16666);
  }
  G4double x4 = 0.25*(3.0*x3 + sqrt(x3*(x3 + 8.0)))/(1.0 - x3);
  //  G4double z1;  
  G4double x  = G4UniformRand(); 
  if(x <= y1) {
    z = x1*pow(x/y1, x4);
  } else {
    z = 1.0 - y1*pow((1.0 - x)/x1, x4);
  }
  //  G4cout << " x4= " << x4 << " z= " << z << G4endl;
  */
  /*  
  // Home made algorithm --------------------------
  G4double x3 = 2.0*x1*x1/x2;
  G4double x4 = 0.0;
  if(x3 > 1.0) x4 = x3 - 1.0 + sqrt(x3*(x3 - 1.0));
  G4double z1 = x1/(1.0 + x4);
  G4double grej = std::max(1.0, 1.0 + x4*(1.0/z1 - 1.0));
  //G4cout << "z1= "<<z1<<" grej= "<<grej<<" x3= "<<x3<<" x1= "<<x1<<" x2= "<<x2<<G4endl;
  
  G4double x;
  do {
    z = -z1*log(G4UniformRand());
    x = G4UniformRand(); 
  } while (z > 1.0 || grej*x > 1.0 + x4*(z/z1 - 1.0));
  */
  // Gaussian part for combined algorithm ---------

  // define threshold angle as 2 sigma of central value
  G4double cosTetSave = cosThetaMin;
  cosThetaMin = 1.0 - 4.0*x1;
  //  G4cout << "cosTmin= " << cosThetaMin << " cosTmax= " << cosThetaMax << G4endl;
 
  if(cosThetaMin < cosThetaMax) cosThetaMin = cosThetaMax;

  // The compute cross section for the tail distribution and 
  // correction to the width of the central part 
  const G4Material* mat = currentCouple->GetMaterial();
  const G4ParticleDefinition* part = dynParticle->GetDefinition();
  G4double cut = (*currentCuts)[currentMaterialIndex];
  G4double xsec = ComputeXSectionPerVolume(mat,part,e,cut);
  if(zcorr < x1) x1 -= zcorr;
  //G4cout << part->GetParticleName() << " e= " << e << "  x1= " 
  //	 << x1 << "  zcorr= " << zcorr << G4endl;

  G4double x  = G4UniformRand();
  G4double z1 = 1.0/x1;
  G4double z2;
  if(z1 < numlimit) z2 = x + (1.0 - x)*(1.0 - z1 + z1*z1*0.5);
  else              z2 = x + (1.0 - x)*exp(-z1);

  z = -x1*log(z2);

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
  	 << " x1= " << x1 << G4endl;
  */
  G4ThreeVector oldDirection = dynParticle->GetMomentumDirection();
  G4ThreeVector newDirection(dirx,diry,cost);
  newDirection.rotateUz(oldDirection);
  G4double rx = 0.0;
  G4double ry = 0.0;

  // sample scattering for large angle -------------

  if(xsec > DBL_MIN) {
    G4double x = -log(G4UniformRand())/xsec;
    if(x < tPathLength) {
      x = tPathLength - x;
      G4double zz1 = 1.0;
      if(G4UniformRand()*xsec < xsec1) {
        zz1 -= G4UniformRand()*(1.0 - cosThetaMin);
      } else {
        zz1 = SampleCosineTheta(mat,part,e,cut);
      }
      if(std::abs(zz1) < 1.0) {
	sint = sqrt((1.0 - zz1)*(1.0 + zz1));

	phi  = twopi*G4UniformRand();

	G4double dirx1 = sint*cos(phi);
	G4double diry1 = sint*sin(phi);
	rx = dirx1*x; 
	ry = diry1*x; 
	G4ThreeVector newDirection1(dirx1,diry1,zz1);
	newDirection1.rotateUz(newDirection);
	newDirection = newDirection1;
      }      
    }
  }
  cosThetaMin = cosTetSave;
  // and of sampling -------------------------------

  fParticleChange->ProposeMomentumDirection(newDirection);

  if (latDisplasment && safety > tlimitminfix) {
    G4double rms= sqrt(2.0*x1);
    rx += zPathLength*(0.5*dirx + invsqrt12*G4RandGauss::shoot(0.0,rms));
    ry += zPathLength*(0.5*diry + invsqrt12*G4RandGauss::shoot(0.0,rms));
    G4double r2 = rx*rx + ry*ry;
    G4double r  = sqrt(r2);
    // protection against situation when geometry length > true step length
    if(zPathLength*zPathLength + r2 > tPathLength*tPathLength) {
      G4double r0 = sqrt((tPathLength - zPathLength)*(tPathLength + zPathLength));
      rx *= r0/r;
      ry *= r0/r;
      r   = r0;
    }
    /*
    G4cout << " r(mm)= " << r << " safety= " << safety
           << " trueStep(mm)= " << tPathLength
           << " geomStep(mm)= " << zPathLength
           << G4endl;
    */
    G4ThreeVector latDirection(rx,ry,0.0);
    latDirection.rotateUz(oldDirection);

    G4ThreeVector Position = *(fParticleChange->GetProposedPosition());
    G4double fac= 1.;
    if(r > safety) {
      //  ******* so safety is computed at boundary too ************
      G4double newsafety = safetyHelper->ComputeSafety(Position);
      if(r > newsafety) fac = newsafety/r ;
    }  

    if(fac > 0.) {
      // compute new endpoint of the Step
      G4ThreeVector newPosition = Position + fac*latDirection;

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

G4double G4MuMscModel::ComputeXSectionPerVolume(const G4Material* mat,
						const G4ParticleDefinition* p,
						G4double kinEnergy, 
						G4double cut)
{
  G4double ekin = std::max(keV, kinEnergy);

  const G4ElementVector* theElementVector = mat->GetElementVector();
  const G4double* theAtomNumDensityVector = mat->GetVecNbOfAtomsPerVolume();
  size_t nelm = mat->GetNumberOfElements();

  xsec1 = 0.0;
  xsec2 = 0.0;
  zcorr = 0.0;
  G4double s, x, x1, x2, x3, z, z1, z2, z3, z4;

  for (size_t i=0; i<nelm; i++) {
    const G4Element* elm = (*theElementVector)[i];
    G4double Z = elm->GetZ();
    G4double A = elm->GetN();
    SetupTarget(Z, A, ekin);
    // initialize and compute cross section on nucleus 
    xsec2 += theAtomNumDensityVector[i]*
      G4eCoulombScatteringModel::ComputeCrossSectionPerAtom(p,ekin,Z,A,cut,ekin);
    x  = 1.0 - cosThetaMin;
    x1 = x + screenZ;
    x2 = 1.0/(x1*x1);
    x3 = 1.0 + x*formfactA;
    // scattering off nucleaus
    s  = x2*x/(x3*x3);
    z4 = s;
    //G4cout << "#0 s= " << s << G4endl;

    // scattering off electrons
    if(cosTetMaxElec <= cosThetaMin) {
      s += x2*x/Z;
      z2 = x/x1;
      z  = 0.5*z2*z2;
      z2 = 1.0 - std::max(cosTetMaxElec,cosThetaMax) + screenZ;
      z3 = (z2 - x1)/x1;
      if(z3 > 0.0) {
	if(z3 < 0.2) z += z3*(x - 0.5*z3*(x - screenZ));
        else         z += log(1.0 + z3)  - screenZ*z3/(x1*(1.0 + z3));
      } 
    } else {
      z1 = 1.0 - cosTetMaxElec;
      z2 = z1 + screenZ;
      s += z1/(z2*z2*Z);
      z3 = z1/z2;
      z  = 0.5*z3*z3;
    }
    z /= Z;
    //G4cout << "#1 s= " << s << " z= " << z << G4endl;

    // scattering off nucleus
    xsec1 += s*Z*Z*theAtomNumDensityVector[i]; 

    z += 0.5*x*z4;

    //
    z2 = 1.0 - cosThetaMax + screenZ;
    z3 = (z2 - x1)/x1;

    if(z3 > 0.0) {
      z1 = 2.0*screenZ*formfactA;
      z4 = z3*formfactA*x1/x3;
      if(z3 < 0.2) 
        z += (1.0 + z1)*z3*(x/x1 + z1 - z3*((x - screenZ)/x1 + 0.5*z1));      
      else 
        z += (1.0 + z1)*(log(1.0 + z3) - z3*screenZ/(x1*(1.0 + z3)));

      if(z4 < 0.2) z -= (1.0 + 2.0*z1)*z4*(1.0 - 0.5*z4);
      else         z -= (1.0 + 2.0*z1)*log(1.0 + z4);
      z -= (1.0 + z1)*z3*formfactA*x1/((1.0 + z4)*x3*x3);
    }

    //G4cout << "z1= "<<z1 << " z2= " << z2<< " x1= " << x1 << " z4= " << z4 << G4endl;

    // common factor
    zcorr += z*Z*Z*theAtomNumDensityVector[i];
    //G4cout << "#2 zcorr= " << zcorr << " z= " << z << G4endl;
  }
  G4double fac = coeff*chargeSquare*invbeta2/mom2;
  xsec1 *= fac;
  G4double xsec = xsec1 + xsec2;
  //G4cout << "xsec1= " << xsec1 << "  xsec2= " << xsec2 << " zsec= " << zcorr*0.5*fac << G4endl;

  zcorr *= 0.5*fac*tPathLength;

  return xsec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuMscModel::SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                     const G4MaterialCutsCouple*,
				     const G4DynamicParticle*,
				     G4double,
				     G4double)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
