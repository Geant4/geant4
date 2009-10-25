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
// $Id: G4WeMoSoftMscModel.cc,v 1.3 2009-10-25 16:46:22 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4WeMoSoftMscModel
//
// Author:      V.Grichine based on G4WentzelVIModel 
//
// Creation date: 31.07.2009 from G4WentzelVIModel
//
// Modifications:
//
//
// Class Description:
//
//
///////////////////////////////////////////////////////////////////////////////////////


#include "G4WeMoSoftMscModel.hh"
#include "Randomize.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4ElementVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

//////////////////////////////////////////////////////////////////////////////////

G4double G4WeMoSoftMscModel::fScreenRSquare[] = {0.0};
G4double G4WeMoSoftMscModel::fFormFactor[]    = {0.0};

using namespace std;

G4WeMoSoftMscModel::G4WeMoSoftMscModel(const G4String& nam) :
  G4VMscModel(nam),
  theLambdaTable(0),
  theLambda2Table(0),
  fNumLimit(0.2),
  fnBins(60),
  fnWarnings(0),
  fnWarnLimit(50),
  fCurrentCouple(0),
  fCosThetaMin(1.0),
  fq2Limit(TeV*TeV),
  fAlpha2(fine_structure_const*fine_structure_const),
  fInitialized(false),
  fInside(false)
{
  fInvSqrt12      = 1./sqrt(12.);
  ftLimitMinFix   = 1.e-6*mm;
  theManager      = G4LossTableManager::Instance(); 
  fNistManager    = G4NistManager::Instance();
  theElectron     = G4Electron::Electron();
  thePositron     = G4Positron::Positron();
  theProton       = G4Proton::Proton();
  fLowEnergyLimit = 0.1*keV;
  G4double p0     = electron_mass_c2*classic_electr_radius;
  fCoeff          = twopi*p0*p0;

  fTkin = fTargetZ = fMom2 = DBL_MIN;

  feCut = feTag = DBL_MAX;

  fParticle  = 0;
  fnElements = 5;

  fXsc.resize(fnElements);
  fProb.resize(fnElements);

  // Thomas-Fermi screening radii
  // Formfactors from A.V. Butkevich et al., NIM A 488 (2002) 282

  if( 0.0 == fScreenRSquare[0] ) 
  {
    G4double a0 = electron_mass_c2/0.88534; 
    G4double constn = 6.937e-6/(MeV*MeV);

    fScreenRSquare[0] = fAlpha2*a0*a0;

    for( G4int j = 1; j < 100; j++ ) 
    {
      G4double x        = a0*fNistManager->GetZ13(j);
      fScreenRSquare[j] = fAlpha2*x*x;
      x                 = fNistManager->GetA27(j); 
      fFormFactor[j]    = constn*x*x;
    } 
  }
}

///////////////////////////////////////////////////////////////////////////

G4WeMoSoftMscModel::~G4WeMoSoftMscModel()
{}

///////////////////////////////////////////////////////////////////////////
//
//

void G4WeMoSoftMscModel::Initialise(const G4ParticleDefinition* p,
				  const G4DataVector& cuts)
{
  // reset parameters

  SetupParticle(p);

  fTkin = fTargetZ = fMom2 = DBL_MIN;
  feCut = feTag = DBL_MAX;

  fCurrentRange = 0.0;
  fCosThetaMax  = cos(PolarAngleLimit());
  fCurrentCuts  = &cuts;

  // set values of some data members

  if(!fInitialized) 
  {
    fInitialized   = true;
    fParticleChange = GetParticleChangeForMSC();
    InitialiseSafetyHelper();
  }
}

///////////////////////////////////////////////////////////////////////////////
//
//

G4double G4WeMoSoftMscModel::ComputeCrossSectionPerAtom( 
                             const G4ParticleDefinition* p,
			     G4double kinEnergy,
			     G4double Z, G4double,
			     G4double cutEnergy, G4double)
{
  SetupParticle(p);
  G4double ekin = std::max(fLowEnergyLimit, kinEnergy);
  SetupKinematic(ekin, cutEnergy);
  SetupTarget(Z, ekin);
  G4double xsec = ComputeTransportXSectionPerAtom();
  /*  
  G4cout << "CS: e= " << fTkin << " cosEl= " << fCosTetMaxElec2 
	 << " cosN= " << fCosTetMaxNuc2 << " xsec(bn)= " << xsec/barn
	 << " " << fParticle->GetParticleName() << G4endl;
  */
  return xsec;
}

////////////////////////////////////////////////////////////////////////////////
//
//

G4double G4WeMoSoftMscModel::ComputeTransportXSectionPerAtom()
{
  G4double xSection = 0.0;
  G4double x, y, x1, x2, x3, x4;

  // scattering off electrons
  if( fCosTetMaxElec2 < 1.0 ) 
  {
    x = ( 1.0 - fCosTetMaxElec2 )/fScreenZ;

    if( x < fNumLimit ) y = 0.5*x*x*(1.0 - 1.3333333*x + 1.5*x*x); 
    else               y = log(1.0 + x) - x/(1.0 + x);

    if( y < 0.0 ) 
    {
      fnWarnings++;

      if(fnWarnings < fnWarnLimit /*&& y < -1.e-10*/) 
      {
	G4cout << "Electron scattering <0 for L1 " << y 
	       << " e(MeV)= " << fTkin << " p(MeV/c)= " << sqrt(fMom2) 
	       << " Z= " << fTargetZ << "  " 
	       << fParticle->GetParticleName() << G4endl;
	G4cout << " z= " << 1.0-fCosTetMaxElec2 << " fScreenZ= " << fScreenZ 
	       << " x= " << x << G4endl;
      }
      y = 0.0;
    }
    xSection = y;
  }
  /*  
  G4cout << "G4WeMoSoftMsc:XS per A " << " Z= " << fTargetZ 
	 << " e(MeV)= " << fTkin/MeV << " XSel= " << xSection
	 << " cut(MeV)= " << feCut/MeV  
  	 << " zmaxE= " << (1.0 - fCosTetMaxElec)/fScreenZ 
	 << " zmaxN= " << (1.0 - fCosTetMaxNuc2)/fScreenZ 
         << " costm= " << fCosTetMaxNuc2 << G4endl;
  */
  // scattering off nucleus

  if(fCosTetMaxNuc2 < 1.0) 
  {
    x  = 1.0 - fCosTetMaxNuc2;
    x1 = fScreenZ*fFormFactA;
    x2 = 1.0 - x1; 
    x3 = x/fScreenZ;
    x4 = fFormFactA*x;

    // low-energy limit

    if(x3 < fNumLimit && x1 < fNumLimit) 
    {
      y = 0.5*x3*x3*(1.0 - 1.3333333*x3 + 1.5*x3*x3 - 1.5*x1
		     + 3.0*x1*x1 + 2.666666*x3*x1)/(x2*x2*x2);      
    } 
    else if(x2 <= 0.0)  // high energy limit
    {
      x4 = x1*(1.0 + x3);
      y  = x3*(1.0 + 0.5*x3 - (2.0 - x1)*(1.0 + x3 + x3*x3/3.0)/x4)/(x4*x4);     
    } 
    else  // middle energy
    {
      y = ((1.0 + x1)*x2*log((1. + x3)/(1. + x4)) 
	   - x3/(1. + x3) - x4/(1. + x4))/(x2*x2); 
    }
    //G4cout << "y= " << y << " x1= " <<x1<<"  x2= " <<x2
    //	   <<"  x3= "<<x3<<"  x4= " << x4<<G4endl;

    if(y < 0.0) 
    {
      fnWarnings++;

      if(fnWarnings < fnWarnLimit /*&& y < -1.e-10*/) 
      { 
	G4cout << "Nuclear scattering <0 for L1 " << y 
	       << " e(MeV)= " << fTkin << " Z= " << fTargetZ << "  " 
	       << fParticle->GetParticleName() << G4endl;
	G4cout << " fFormFactA= " << fFormFactA << " fScreenZ= " << fScreenZ 
	       << " x= " << " x1= " << x1 << " x2= " << x2 
	       << " x3= " << x3 << " x4= " << x4 <<G4endl;
      }
      y = 0.0;
    }
    xSection += y*fTargetZ; 
  }
  xSection *= fKinFactor;

  /*
  G4cout << "Z= " << fTargetZ << " XStot= " << xSection/barn 
	 << " fScreenZ= " << fScreenZ << " formF= " << fFormFactA 
	 << " for " << fParticle->GetParticleName() 
  	 << " m= " << fMass << " 1/v= " << sqrt(fInvBeta2) << " p= " << sqrt(fMom2)
	 << " x= " << x 
	 << G4endl;
  */
  return xSection; 
}

////////////////////////////////////////////////////////////////////////////////
//
//

G4double G4WeMoSoftMscModel::ComputeTruePathLengthLimit(
                             const G4Track& track,
			     G4PhysicsTable* theTable,
			     G4double currentMinimalStep)
{
  G4double tlimit             = currentMinimalStep;
  const G4DynamicParticle* dp = track.GetDynamicParticle();
  G4StepPoint* sp             = track.GetStep()->GetPreStepPoint();
  G4StepStatus stepStatus     = sp->GetStepStatus();

  // initialisation for 1st step 
 
  if(stepStatus == fUndefined) 
  {
    fInside = false;
    SetupParticle(dp->GetDefinition());
    theLambdaTable = theTable;
  }

  // initialisation for each step, lambda may be computed from scratch

  fPreKinEnergy  = dp->GetKineticEnergy();
  DefineMaterial(track.GetMaterialCutsCouple());
  fLambda1       = GetLambda(fPreKinEnergy);
  fCurrentRange  = 
    theManager->GetRangeFromRestricteDEDX(fParticle,fPreKinEnergy,fCurrentCouple);

  // extra check for abnormal situation
  // this check needed to run MSC with eIoni and eBrem inactivated

  if(tlimit > fCurrentRange) tlimit = fCurrentRange;

  // stop here if small range particle

  if(fInside) return tlimit;   

  // pre step

  G4double presafety = sp->GetSafety();

  // compute presafety again if presafety <= 0 and no boundary
  // i.e. when it is needed for optimization purposes

  if(stepStatus != fGeomBoundary && presafety < ftLimitMinFix) 
    presafety = ComputeSafety(sp->GetPosition(), tlimit); 

  /*
  G4cout << "G4WeMoSoftMscModel::ComputeTruePathLengthLimit tlimit= " 
 	 <<tlimit<<" safety= " << presafety
  	 << " range= " <<fCurrentRange<<G4endl;
  */

  // far from geometry boundary

  if(fCurrentRange < presafety) 
  {
    fInside = true;      
  } 
  else // limit mean scattering angle
  {
    G4double rlimit = facrange*fLambda1;
    G4double rcut   = fCurrentCouple->GetProductionCuts()->GetProductionCut(1);

    if(rcut > rlimit) rlimit = std::pow(rcut*rcut*rlimit,0.33333333);

    rlimit = std::min(rlimit, facgeom*fCurrentMaterial->GetRadlen());

    if(rlimit < tlimit) tlimit = rlimit;
  }
  /*
  G4cout << fParticle->GetParticleName() << " e= " << fPreKinEnergy
	 << " L0= " << fLambda1 << " R= " << fCurrentRange
	 << "tlimit= " << tlimit  
  	 << " currentMinimalStep= " << currentMinimalStep << G4endl;
  */
  return tlimit;
}

/////////////////////////////////////////////////////////////////////////////////
//
//

G4double G4WeMoSoftMscModel::ComputeGeomPathLength(G4double truelength)
{
  fTruePathLength  = truelength;
  fGeomPathLength  = fTruePathLength;
  fLambdaEff    = fLambda1;

  if(fLambda1 > 0.0) 
  {
    G4double tau = fTruePathLength/fLambda1;

    //G4cout << "ComputeGeomPathLength: tLength= " << fTruePathLength
    //	 << " fLambda1= " << fLambda1 << " tau= " << tau << G4endl; 
    // small step

    if( tau < fNumLimit ) 
    {
      fGeomPathLength *= (1.0 - 0.5*tau + tau*tau/6.0);

      // medium step
    } 
    else 
    {
      //      fGeomPathLength = fLambda1*(1.0 - exp(-fTruePathLength/fLambda1));

      G4double e1 = 0.0;

      if( fCurrentRange > fTruePathLength ) 
      {
	e1 = theManager->GetEnergy( fParticle,
				    fCurrentRange-fTruePathLength,
				    fCurrentCouple              );
      }
      fLambdaEff = GetLambda(0.5*(e1 + fPreKinEnergy));
      fGeomPathLength = fLambdaEff*(1.0 - exp(-fTruePathLength/fLambdaEff));
    }
  }
  //G4cout<<"Comp.geom: zLength= "<<fGeomPathLength<<" tLength= "<<fTruePathLength<<G4endl;

  return fGeomPathLength;
}

///////////////////////////////////////////////////////////////////////////////
//
//

G4double G4WeMoSoftMscModel::ComputeTrueStepLength(G4double geomStepLength)
{
  // step defined other than transportation 

  if(geomStepLength == fGeomPathLength) return fTruePathLength;

  // step defined by transportation 

  fTruePathLength  = geomStepLength;
  fGeomPathLength  = geomStepLength;
  G4double tau = fGeomPathLength/fLambdaEff;
  fTruePathLength *= (1.0 + 0.5*tau + tau*tau/3.0); 

  if(tau > fNumLimit) 
  {
    G4double e1 = 0.0;

    if( fCurrentRange > fTruePathLength ) 
    {
      e1 = theManager->GetEnergy(fParticle,
				 fCurrentRange-fTruePathLength,
				 fCurrentCouple);
    }
    fLambdaEff = GetLambda(0.5*(e1 + fPreKinEnergy));
    tau       = fGeomPathLength/fLambdaEff;

    if(tau < 0.999999) fTruePathLength = -fLambdaEff*log(1.0 - tau); 
    else               fTruePathLength = fCurrentRange;

    if(fTruePathLength < fGeomPathLength) fTruePathLength = fGeomPathLength;
  }
  if(fTruePathLength > fCurrentRange) fTruePathLength = fCurrentRange;

  //G4cout<<"Comp.true: zLength= "<<fGeomPathLength<<" tLength= "<<fTruePathLength<<G4endl;

  return fTruePathLength;
}

////////////////////////////////////////////////////////////////////////////////////
//
//

void G4WeMoSoftMscModel::SampleScattering(const G4DynamicParticle* dynParticle,
					G4double safety)
{
  //G4cout << "!##! G4WeMoSoftMscModel::SampleScattering for " 
  //	 << fParticle->GetParticleName() << G4endl;
  
  G4double kinEnergy = dynParticle->GetKineticEnergy();

  if(kinEnergy <= DBL_MIN || fTruePathLength <= DBL_MIN) return;

  G4double ekin = fPreKinEnergy;

  if(ekin - kinEnergy > ekin*dtrl) 
  {
    ekin = 0.5*(fPreKinEnergy + kinEnergy);
    fLambdaEff = GetLambda(ekin);
  }  
  
  G4double x1 = 0.5*fTruePathLength/fLambdaEff;
  G4double cut= (*fCurrentCuts)[fCurrentMaterialIndex];
  /*  
  G4cout <<"SampleScat: E0(MeV)= "<< fPreKinEnergy<<" Eeff(MeV)= "<<ekin/MeV
	 << " L0= " << fLambda1 << " Leff= " << fLambdaEff 
	 << " x1= " << x1 << " safety= " << safety << G4endl;
  */

  G4double xsec = 0.0;
  G4bool largeAng = false;

  // large scattering angle case

  if(x1 > 0.5) 
  {
    x1      *= 0.5;
    largeAng = true;    
  } 
  else // normal case
  {
    // define threshold angle between single and multiple scattering 

    fCosThetaMin = 1.0 - 3.0*x1;

    // for low-energy e-,e+ no limit

    ekin = std::max(ekin, fLowEnergyLimit);
    SetupKinematic(ekin, cut);
  
    // recompute transport cross section

    if(fCosThetaMin > fCosTetMaxNuc) 
    {
      xsec = ComputeXSectionPerVolume();

      if(fTransportXsc > DBL_MIN) x1 = 0.5*fTruePathLength*fTransportXsc;
      else                x1 = 0.0;

      /*      
	G4cout << "fCosTetMaxNuc= " << fCosTetMaxNuc 
	<< " fCosThetaMin= " << fCosThetaMin 
	<< " fCosThetaMax= " << fCosThetaMax 
	<< " fCosTetMaxElec2= " << fCosTetMaxElec2 << G4endl;
	G4cout << "Recomputed xsec(1/mm)= " << xsec << " x1= " << x1 << G4endl;
      */
    }
  }

  // result of central part sampling 

  G4double z;
 
  do 
  {
    z = -x1*log(G4UniformRand());
  } 
  while (z > 1.0); 

  // cost is sampled 

  G4double              cost = 1.0 - 2.0*z;
  if( cost < -1.0 )     cost = -1.0;
  else if( cost > 1.0 ) cost = 1.0;

  G4double sint = sqrt((1.0 - cost)*(1.0 + cost));

  G4double phi  = twopi*G4UniformRand();

  G4double dirx = sint*cos(phi);
  G4double diry = sint*sin(phi);
  
  //G4cout << "G4WeMoSoftMscModel: step(mm)= " << fTruePathLength/mm
  //	 << " sint= " << sint << " cost= " << cost<< G4endl;
  
  G4ThreeVector oldDirection = dynParticle->GetMomentumDirection();
  G4ThreeVector newDirection(dirx,diry,cost);
  G4ThreeVector temp(0.0,0.0,1.0);
  G4ThreeVector pos(0.0,0.0,-fGeomPathLength);
  G4ThreeVector dir(0.0,0.0,1.0);
  G4bool isscat = false;

  // sample MSC scattering for large angle
  // extra central scattering for holf step

  if(largeAng) 
  {
    isscat = true;
    pos.setZ(-0.5*fGeomPathLength);

    do 
    {
      z = -x1*log(G4UniformRand());
    } 
    while (z > 1.0); 

    cost = 1.0 - 2.0*z;

    if(std::abs(cost) > 1.0) cost = 1.0;

    sint = sqrt((1.0 - cost)*(1.0 + cost));
    phi  = twopi*G4UniformRand();

    // position and direction for secondary scattering

    dir.set(sint*cos(phi),sint*sin(phi),cost);
    pos += 0.5*dir*fGeomPathLength;
    x1 *= 2.0;
  }

  // sample Rutherford scattering for large angle

  if(xsec > DBL_MIN) 
  {
    G4double t = fTruePathLength;
    G4int nelm = fCurrentMaterial->GetNumberOfElements();
    const G4ElementVector* theElementVector = 
      fCurrentMaterial->GetElementVector();

    do
    {
      G4double x = -log(G4UniformRand())/xsec;      
      pos       += dir*(fGeomPathLength*std::min(x,t)/fTruePathLength);
      t         -= x;

      if( t > 0.0 ) 
      {
	G4double zz1  = 1.0;
	G4double qsec = G4UniformRand()*xsec;

	// scattering off nucleus

        G4int i = 0;

	if(nelm > 1) 
        {
	  for (; i < nelm; i++ ) 
          {
            if( fXsc[i] >= qsec ) break;
          }
	  if(i >= nelm) i = nelm - 1;
	}
	SetupTarget((*theElementVector)[i]->GetZ(), fTkin);
        G4double formf = fFormFactA;
        G4double costm = fCosTetMaxNuc2;

        if(fProb[i] > 0.0) 
        {
	  if(G4UniformRand() <= fProb[i]) 
          {
	    formf = 0.0;
	    costm = fCosTetMaxElec2;
	  }
	}
	if( fCosThetaMin > costm ) 
        {
	  G4double w1 = 1. - fCosThetaMin + fScreenZ;
	  G4double w2 = 1. - costm + fScreenZ;
	  G4double w3 = fCosThetaMin - costm;

	  G4double grej, zz; 

	  do 
          {
	    zz = w1*w2/(w1 + G4UniformRand()*w3) - fScreenZ;
	    grej = 1.0/(1.0 + formf*zz);
	  } 
          while ( G4UniformRand() > grej*grej ); 
 
	  if(zz < 0.0) zz = 0.0;
	  else if(zz > 2.0) zz = 2.0;

	  zz1 = 1.0 - zz;
	}
        if(zz1 < 1.0) 
        {
	  isscat = true;

	  //G4cout << "Rutherford zz1= " << zz1 << " t= " << t << G4endl;

	  sint = sqrt((1.0 - zz1)*(1.0 + zz1));

	  //G4cout << "sint= " << sint << G4endl;

	  phi          = twopi*G4UniformRand();
	  G4double vx1 = sint*cos(phi);
	  G4double vy1 = sint*sin(phi);
	  temp.set(vx1,vy1,zz1);
	  temp.rotateUz(dir);
	  dir = temp;
	}
      }
    } 
    while (t > 0.0); 
  }
  if(isscat) newDirection.rotateUz(dir);

  newDirection.rotateUz(oldDirection);

  //G4cout << "G4WeMoSoftMscModel sampling of scattering is done" << G4endl;
  // end of sampling -------------------------------

  fParticleChange->ProposeMomentumDirection(newDirection);

  if (latDisplasment && safety > ftLimitMinFix) 
  {
    G4double rms = fInvSqrt12*sqrt(2.0*x1);
    G4double dx = fGeomPathLength*(0.5*dirx + rms*G4RandGauss::shoot(0.0,1.0));
    G4double dy = fGeomPathLength*(0.5*diry + rms*G4RandGauss::shoot(0.0,1.0));
    G4double dz;
    G4double d = (dx*dx + dy*dy)/(fGeomPathLength*fGeomPathLength);

    if(      d < fNumLimit)  dz = -0.5*fGeomPathLength*d*(1.0 + 0.25*d);
    else if( d < 1.0     )  dz = -fGeomPathLength*(1.0 - sqrt(1.0 - d));
    else 
    {
                  dx = dy = dz = 0.0;
    }

    temp.set(dx,dy,dz);

    if( isscat ) temp.rotateUz(dir);

    pos += temp;
   
    pos.rotateUz(oldDirection);

    G4double r = pos.mag();

    /*    
    G4cout << " r(mm)= " << r << " safety= " << safety
           << " trueStep(mm)= " << fTruePathLength
           << " geomStep(mm)= " << fGeomPathLength
           << G4endl;
    */

    if(r > ftLimitMinFix) 
    {
      pos /= r;
      ComputeDisplacement(fParticleChange, pos, r, safety);
    }
  }
  //G4cout << "G4WeMoSoftMscModel::SampleScattering end" << G4endl;
}

///////////////////////////////////////////////////////////////////////////////////
//
//

G4double G4WeMoSoftMscModel::ComputeXSectionPerVolume()
{
  const G4ElementVector* theElementVector = 
    fCurrentMaterial->GetElementVector();
  const G4double* theAtomNumDensityVector = 
    fCurrentMaterial->GetVecNbOfAtomsPerVolume();
  G4int nelm = fCurrentMaterial->GetNumberOfElements();

  if(nelm > fnElements) 
  {
    fnElements = nelm;
    fXsc.resize(fnElements);
    fProb.resize(fnElements);
  }

  fTransportXsc = 0.0;
  G4double xs = 0.0;

  for (G4int i=0; i<nelm; i++) 
  {
    SetupTarget((*theElementVector)[i]->GetZ(), fTkin);
    G4double density = theAtomNumDensityVector[i];
    G4double cosnm = fCosTetMaxNuc2;
    G4double cosem = fCosTetMaxElec2;

    // recompute the angular limit
 
    fCosTetMaxNuc2  = std::max(cosnm,fCosThetaMin); 
    fCosTetMaxElec2 = std::max(cosem,fCosThetaMin); 
    fTransportXsc += ComputeTransportXSectionPerAtom()*density;

    // return limit back

    fCosTetMaxElec2 = cosem;
    fCosTetMaxNuc2  = cosnm;

    G4double esec  = 0.0;
    G4double nsec  = 0.0;
    G4double x1    = 1.0 - fCosThetaMin + fScreenZ;
    G4double f     = fKinFactor*density; 

    // scattering off electrons

    if(fCosThetaMin > cosem) 
    {
      esec = f*(fCosThetaMin - cosem)/(x1*(1.0 - cosem + fScreenZ));
    }

    // scattering off nucleaus

    if(fCosThetaMin > cosnm) 
    {

      // Rutherford part

      G4double s  = fScreenZ*fFormFactA;
      G4double z1 = 1.0 - cosnm + fScreenZ;
      G4double s1 = 1.0 - s;
      G4double d  = s1/fFormFactA;

      // check numerical limit

      if(d < fNumLimit*x1) 
      {
	G4double x2 = x1*x1;
	G4double z2 = z1*z1;
	nsec = (1.0/(x1*x2) - 1.0/(z1*z2) - d*1.5*(1.0/(x2*x2) - 1.0/(z2*z2)))/
	  (3.0*fFormFactA*fFormFactA);
      } 
      else 
      {
	G4double x2 = x1 + d;
	G4double z2 = z1 + d;
	nsec = (1.0/x1 - 1.0/z1 + 1.0/x2 - 1.0/z2 - 2.0*log(z1*x2/(z2*x1))/d)/(s1*s1);
      }
      nsec *= f*fTargetZ;
    }
    nsec += esec;

    if(nsec > 0.0) esec /= nsec;

    xs      += nsec;
    fXsc[i] = xs;
    fProb[i]  = esec;

    //G4cout << i << "  xs= " << xs << " fCosThetaMin= " << fCosThetaMin 
    //	   << " costm= " << costm << G4endl;
  }
  
  //G4cout << "ComputeXS result:  xsec(1/mm)= " << xs 
  //<< " txsec(1/mm)= " << fTransportXsc <<G4endl; 
  return xs;
}

//////////////////////////////////////////////////////////////////////////////////
//
//

/*
G4double G4MuMscModel::ComputeXSectionPerVolume()
{
  const G4ElementVector* theElementVector = 
    fCurrentMaterial->GetElementVector();
  const G4double* theAtomNumDensityVector = 
    fCurrentMaterial->GetVecNbOfAtomsPerVolume();
  size_t nelm = fCurrentMaterial->GetNumberOfElements();

  xsece1 = 0.0;
  xsece2 = 0.0;
  fXsc2 = 0.0;
  zcorr  = 0.0;

  G4double fac = fCoeff*fChargeSquare*fInvBeta2/fMom2;

  for ( size_t i=0; i<nelm; i++   ) 
  {
    const G4Element* elm = (*theElementVector)[i];
    G4double Z = elm->GetZ();
    SetupTarget(Z, fTkin);
    G4double den = fac*theAtomNumDensityVector[i]*Z;

    G4double x  = 1.0 - fCosThetaMin;
    G4double x1 = x + fScreenZ;
    G4double x2 = 1.0/(x1*x1);
    G4double x3 = 1.0 + x*fFormFactA;
    
    //G4cout << "x= " << x << " den= " << den << " cosE= " << fCosTetMaxElec << G4endl;
    //G4cout << "cosThtaMin= " << fCosThetaMin << G4endl;
    //G4cout << "fCosTetMaxNuc= " << fCosTetMaxNuc << " fq2Limit= " << fq2Limit << G4endl;
    
    // scattering off electrons

    if(fCosTetMaxElec < fCosThetaMin) 
    {

      // flat part
      G4double s = den*x2*x;
      xsece1 += s;
      zcorr  += 0.5*x*s;

      // Rutherford part
      G4double z1 = 1.0 - fCosTetMaxElec + fScreenZ;
      G4double z2 = (fCosThetaMin - fCosTetMaxElec)/x1; 
      if(z2 < 0.2) s = z2*(x - 0.5*z2*(x - fScreenZ))/x1;
      else         s = log(1.0 + z2)  - fScreenZ*z2/z1;
      xsece2  += den*z2/z1;
      zcorr   += den*s;
    }
    den *= Z;

    //G4cout << "Z= " << Z<< " cosL= " << fCosTetMaxNuc << " cosMin= " << fCosThetaMin << G4endl;
    // scattering off nucleaus

    if(fCosTetMaxNuc < fCosThetaMin) 
    {

      // flat part
      G4double s = den*x2*x/(x3*x3);
      xsece1 += s;
      zcorr  += 0.5*x*s;

      // Rutherford part

      s  = fScreenZ*fFormFactA;
      G4double w  = 1.0 + 2.0*s;
      G4double z1 = 1.0 - fCosTetMaxNuc + fScreenZ;
      G4double d  = (1.0 - s)/fFormFactA;
      G4double x4 = x1 + d;
      G4double z4 = z1 + d;
      G4double t1 = 1.0/(x1*z1);
      G4double t4 = 1.0/(x4*z4);
      G4double w1 = fCosThetaMin - fCosTetMaxNuc;
      G4double w2 = log(z1*x4/(x1*z4));

      den *= w;     
      xsecn2  += den*(w1*(t1 + t4) - 2.0*w2/d);
      zcorr   += den*(w*w2 - w1*(fScreenZ*t1 + t4/fFormFactA));
    }
    xsece[i] = xsece2;
    fXsc[i] = xsecn2;
    //    G4cout << i << "  xsece2= " << xsece2 << "  xsecn2= " << xsecn2 << G4endl;
  }
  G4double xsec = xsece1 + xsece2 + xsecn2;
 
    //G4cout << "xsece1= " << xsece1 << "  xsece2= " << xsece2 
    //<< "  xsecn2= " << xsecn2 
	// << " zsec= " << zcorr*0.5*fTruePathLength << G4endl;
  zcorr *= 0.5*fTruePathLength;

  return xsec;
}
*/

///////////////////////////////////////////////////////////////////////////////////
//
//

void G4WeMoSoftMscModel::ComputeMaxElectronScattering(G4double cutEnergy)
{
  feCut          = cutEnergy;
  G4double tmax = fTkin;
  fCosTetMaxElec = 1.0;

  if( fMass > MeV ) 
  {
    G4double ratio = electron_mass_c2/fMass;
    G4double tau   = fTkin/fMass;
    tmax           = 2.0*electron_mass_c2*tau*(tau + 2.)/
              (1.0 + 2.0*ratio*(tau + 1.0) + ratio*ratio); 

    fCosTetMaxElec  = 1.0 - std::min(cutEnergy, tmax)*electron_mass_c2/fMom2;
  } 
  else 
  {
    if(fParticle == theElectron) tmax *= 0.5;

    G4double t = std::min(cutEnergy, tmax);
    G4double mom21 = t*(t + 2.0*electron_mass_c2);
    G4double t1 = fTkin - t;

    //G4cout <<"fTkin=" <<fTkin<<" tmax= "<<tmax<<" t= " 
    //<<t<< " t1= "<<t1<<" cut= "<<feCut<<G4endl;

    if(t1 > 0.0) 
    {
      G4double mom22 = t1*(t1 + 2.0*fMass);
      G4double ctm = (fMom2 + mom22 - mom21)*0.5/sqrt(fMom2*mom22);
      if(ctm <  1.0) fCosTetMaxElec = ctm;
      if(ctm < -1.0) fCosTetMaxElec = -1.0;
    }
  }
  if( fCosTetMaxElec < fCosTetMaxNuc )  fCosTetMaxElec = fCosTetMaxNuc;
}

//
//
/////////////////////////////////////////////////////////////////////////////////
