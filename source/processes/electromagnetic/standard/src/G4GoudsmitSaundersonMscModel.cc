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
// $Id: G4GoudsmitSaundersonMscModel.cc,v 1.3 2009-03-19 14:17:51 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4GoudsmitSaundersonMscModel
//
// Author:        Omrane Kadri
//
// Creation date: 20.02.2009
//
// Modifications:
// 04.03.2009 V.Ivanchenko cleanup and format according to Geant4 EM style
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//REFERENCES:
//Ref.1:E. Benedito et al.,"Mixed simulation ... cross-sections", NIMB 174 (2001) pp 91-110;
//Ref.2:I. Kawrakow et al.,"On the condensed ... transport",NIMB 142 (1998) pp 253-280;
//Ref.3:I. Kawrakow et al.,"On the representation ... calculations",NIMB 134 (1998) pp 325-336;
//Ref.4:Bielajew et al.,".....", NIMB 173 (2001) 332-343;
//Ref.5:F. Salvat et al.,"ELSEPA--Dirac partial ...molecules", Comp. Phys. Comm. 165 (2005) pp 157-190;
//Ref.6:G4UrbanMscModel G4_v9.1Ref09; 
//Ref.7:G4eCoulombScatteringModel G4_v9.1Ref09.
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4GoudsmitSaundersonTable.hh"

#include "G4ParticleChangeForMSC.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4DynamicParticle.hh"
#include "G4DataInterpolation.hh" 
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4LossTableManager.hh"
#include "G4Track.hh"
#include "G4PhysicsTable.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"

G4double G4GoudsmitSaundersonMscModel::ener[106] = {-1.};
G4double G4GoudsmitSaundersonMscModel::TCSE[103][106];
G4double G4GoudsmitSaundersonMscModel::FTCSE[103][106];
G4double G4GoudsmitSaundersonMscModel::TCSP[103][106];
G4double G4GoudsmitSaundersonMscModel::FTCSP[103][106];

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4GoudsmitSaundersonMscModel::G4GoudsmitSaundersonMscModel(const G4String& nam)
  : G4VMscModel(nam),lowKEnergy(0.1*keV),highKEnergy(GeV),isInitialized(false)
{ 
  fr=0.02,rangeinit=0.,masslimite=0.6*MeV;
  particle=0;tausmall=1.e-16;taulim=1.e-6;tlimit=1.e10*mm;
  tlimitmin=10.e-6*mm;geombig=1.e50*mm;geommin=1.e-3*mm,tgeom=geombig;
  tlimitminfix=1.e-6*mm;stepmin=tlimitminfix;lambdalimit=1.*mm;smallstep=1.e10;
  theManager=G4LossTableManager::Instance();
  inside=false;insideskin=false;
  samplez=false;

  GSTable = new G4GoudsmitSaundersonTable();

  if(ener[0] < 0.0){ 
    G4cout << "### G4GoudsmitSaundersonMscModel loading ELSEPA data" << G4endl;
    LoadELSEPAXSections();
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4GoudsmitSaundersonMscModel::~G4GoudsmitSaundersonMscModel()
{
  delete GSTable;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4GoudsmitSaundersonMscModel::Initialise(const G4ParticleDefinition* p,
					      const G4DataVector&)
{ 
  skindepth=skin*stepmin;
  SetParticle(p);
  if(isInitialized) return;
  if (pParticleChange) fParticleChange = 
			 reinterpret_cast<G4ParticleChangeForMSC*>(pParticleChange);
  else   fParticleChange = new G4ParticleChangeForMSC();

  InitialiseSafetyHelper();
  isInitialized=true;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4GoudsmitSaundersonMscModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition* p,
                          G4double kineticEnergy,G4double Z, G4double, G4double, G4double)
{  
  //Build cross section table : Taken from Ref.7
  G4double cs=0.0;
  G4double kinEnergy = kineticEnergy;
  if(kinEnergy<lowKEnergy) kinEnergy=lowKEnergy;
  if(kinEnergy>highKEnergy)kinEnergy=highKEnergy;

  //value0=Lambda0;value1=Lambda1
  G4double value0,value1;
  CalculateIntegrals(p,Z,kinEnergy,value0,value1);
  
  if(value1 > 0.0) cs = 1./value1;

  return cs;
}  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4GoudsmitSaundersonMscModel::SampleScattering(const G4DynamicParticle* dynParticle,
				       G4double safety)
{    
  G4double kineticEnergy = dynParticle->GetKineticEnergy();
  if((kineticEnergy <= 0.0) || (tPathLength <= tlimitminfix)) return ;

  G4double scrA,llambda0,llambda1;
  G4double cosTheta1,sinTheta1,cosTheta2,sinTheta2;
  G4double phi1,cosPhi1,sinPhi1,phi2,cosPhi2,sinPhi2;
  G4double us,vs,ws,q1,Gamma,Eta,delta,nu,nu0,nu1,nu2,nu_interm,x_coord,y_coord,z_coord;

  ///////////////////////////////////////////
  // Effective energy and path-length from Eq. 4.7.15+16 of Ref.4
  G4double  eloss = theManager->GetEnergy(particle,tPathLength,currentCouple);
  G4double ee       = kineticEnergy - 0.5*eloss;
  G4double ttau     = ee/electron_mass_c2;

  G4double ttau2    = ttau*ttau;
  G4double epsilonpp= eloss/ee;
  G4double temp2  = 0.166666*(4+ttau*(6+ttau*(7+ttau*(4+ttau))))*(epsilonpp/(ttau+1)/(ttau+2))*(epsilonpp/(ttau+1)/(ttau+2));
  G4double cst1=epsilonpp*epsilonpp*(6+10*ttau+5*ttau2)/(24*ttau2+48*ttau+72);

  kineticEnergy *= (1 - cst1);
  tPathLength *= (1 - temp2);
  ///////////////////////////////////////////
  // additivity rule for mixture xsection calculation
  const G4Material* mat = currentCouple->GetMaterial();
  G4int nelm = mat->GetNumberOfElements();
  const G4ElementVector* theElementVector = mat->GetElementVector();
  const G4double* theFraction = mat->GetFractionVector();
  G4double atomPerVolume = mat->GetTotNbOfAtomsPerVolume();
  llambda0 =0.;llambda1=0.;
  for(G4int i=0;i<nelm;i++)
    {
      G4double l0,l1;
      CalculateIntegrals(particle,(*theElementVector)[i]->GetZ(),kineticEnergy,l0,l1);
      llambda0 += (theFraction[i]/l0);
      llambda1 += (theFraction[i]/l1);
    } 
  llambda0 =1./llambda0;
  llambda1 =1./llambda1;
  G4double g1=llambda0/llambda1;
  G4double x1,x0;

  x0=g1/2.;
  do
    {  
      x1 = x0-(x0*((1.+x0)*std::log(1.+1./x0)-1.0)-g1/2.)/( (1.+2.*x0)*std::log(1.+1./x0)-2.0);// x1=x0-f(x0)/f'(x0)
      delta = std::abs( x1 - x0 );    
      x0 = x1;  // new approximation becomes the old approximation for the next iteration
    } while (delta > 1e-10);
  scrA = x1;

  G4double lambdan=0.;
  if((tPathLength>0.)&&(llambda0>0.))lambdan=std::max((atomPerVolume*tPathLength/llambda0),1.0e-12);

  // Ref.2 subsection 4.4 "The best solution found"
  // Sample first substep scattering angle
  SampleCosineTheta(0.5*lambdan,scrA,cosTheta1,sinTheta1);

  phi1  = twopi*G4UniformRand();
  cosPhi1 = cos(phi1);
  sinPhi1 = sin(phi1);

  // Sample second substep scattering angle
  SampleCosineTheta(0.5*lambdan,scrA,cosTheta2,sinTheta2);
  phi2  = twopi*G4UniformRand();
  cosPhi2 = cos(phi2);
  sinPhi2 = sin(phi2);

  // Scattering direction
  us = sinTheta2*(cosTheta1*cosPhi1*cosPhi2 - sinPhi1*sinPhi2) + cosTheta2*sinTheta1*cosPhi1;
  vs = sinTheta2*(cosTheta1*sinPhi1*cosPhi2 + cosPhi1*sinPhi2) + cosTheta2*sinTheta1*sinPhi1;
  ws = cosTheta1*cosTheta2 - sinTheta1*sinTheta2*cosPhi2;  

  G4ThreeVector oldDirection = dynParticle->GetMomentumDirection();
  G4ThreeVector newDirection(us,vs,ws);
  newDirection.rotateUz(oldDirection);
  fParticleChange->ProposeMomentumDirection(newDirection);

  if((safety > tlimitminfix)&&(latDisplasment))
    {  
      // Scattering coordinates
      if(scrA<DBL_MIN)scrA=DBL_MIN;
      if(llambda0<DBL_MIN)llambda0=DBL_MIN;
      q1       = 2.*scrA*((1. + scrA)*log(1. + 1./scrA) - 1.);
      if(q1<DBL_MIN)q1=DBL_MIN;
      Gamma    = 6.*scrA*(1. + scrA)*((1. + 2.*scrA)*log(1. + 1./scrA)* - 2.)/q1;
      Eta      = atomPerVolume*tPathLength/llambda0;
      delta    = 0.90824829 - Eta*(0.102062073-Gamma*0.026374715);

      nu = G4UniformRand(); 
      nu = std::sqrt(nu);
      nu0 = (1.0 - nu)/2.;
      nu1 = nu*delta;
      nu2 = nu*(1.0-delta);
      nu_interm = 1.0 - nu0 - nu1 - nu2;
      x_coord=nu1*sinTheta1*cosPhi1+nu2*sinTheta2*(cosPhi1*cosPhi2-cosTheta1*sinPhi1*sinPhi2)+nu_interm*us;
      y_coord=nu1*sinTheta1*sinPhi1+nu2*sinTheta2*(sinPhi1*cosPhi2+cosTheta1*cosPhi1*sinPhi2)+nu_interm*vs;
      z_coord=nu0+nu1*cosTheta1+nu2*cosTheta2+ nu_interm*ws  ;

      G4double r=sqrt(x_coord*x_coord+y_coord*y_coord+z_coord*z_coord);
      if(r<DBL_MIN){x_coord=0.;y_coord=0.;z_coord=0.;}
      else {
	x_coord /=r; 
	y_coord /=r; 
	z_coord /=r; 
      }

      G4double check= 1.- tPathLength/zPathLength;
      if(check<=0.)    r=0.;
      else if(r>check) r=check;

      r *=tPathLength;
      
      if(r > tlimitminfix) {

        G4ThreeVector latDirection(x_coord,y_coord,z_coord);
        latDirection.rotateUz(oldDirection);
 
	ComputeDisplacement(fParticleChange, latDirection, r, safety);

      }     
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4GoudsmitSaundersonMscModel::SampleCosineTheta(G4double lambdan, G4double scrA,
						     G4double &cost, G4double &sint)
{
  G4double u,Qn1,r1,tet;

  G4double xi=0.;
  G4double epsilon1=G4UniformRand();
  if(epsilon1<(exp(-lambdan))) xi=0.0;// no scattering
  else 
    { 
      Qn1=2.* lambdan *scrA*((1.+scrA)*log(1.+1./scrA)-1.);
      if(epsilon1<((1.+lambdan)*exp(-lambdan)))//just one collision 
	{ 
	  xi = G4UniformRand();
	  xi= 2.*scrA*xi/(1.-xi + scrA); 
	}  
      else    //two or more collisions
	{
	  if((lambdan<1.)||(Qn1<0.001))//plural scatt. or small angle scatt.
	    {
	      G4double xi1,lambdai=0.;
	      G4int i=0;
	      do {xi1=G4UniformRand();
	      lambdai -=std::log(xi1);
	      xi +=2.*scrA*xi1/(1.-xi1 + scrA);
	      i++;
	      }while((lambdai<lambdan)&&(i<30));
	    }
	  else {
	    if(Qn1>0.5)xi=G4UniformRand();//isotropic distribution
	    else{// procedure described by Benedito in Ref.1
	      do{r1=G4UniformRand();
	      u=G4UniformRand();
	      u=GSTable->SampleTheta(lambdan,scrA,u);
	      xi = 2.*u;
	      tet=acos(1.-xi);
	      }while(tet*r1*r1>sin(tet));
	    }    
	  }
	}
    }       	 

  if(!((xi>=0.)&&(xi<=2.)))xi=0.;
  cost=(1. - xi);
  sint=sqrt(xi*(2.-xi));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Cubic spline log-log interpolation of Lambda0 and Lambda1
// Screening parameter calculated according to Eq. 37 of Ref.1 
void G4GoudsmitSaundersonMscModel::CalculateIntegrals(const G4ParticleDefinition* p,G4double Z, 
						      G4double kinEnergy,G4double &Lam0,
						      G4double &Lam1)
{ 
  //////// BEGIN OF: LAMBDA CALCULATION ////////////////////////////////
  G4double TCSEForThisAtom[106],FTCSEForThisAtom[106],TCSPForThisAtom[106],FTCSPForThisAtom[106];
  G4double summ00=0.0;
  G4double summ10=0.0;
  G4double InterpolatedValue=0.0;
  
  
  G4int  iZ = G4int(Z);
  if(iZ > 103) iZ = 103;
  for(G4int i=0;i<106;i++)
    {
      TCSEForThisAtom[i]=TCSE[iZ-1][i];FTCSEForThisAtom[i]=FTCSE[iZ-1][i];
      TCSPForThisAtom[i]=TCSP[iZ-1][i];FTCSPForThisAtom[i]=FTCSP[iZ-1][i];
    }

  G4double kineticE = kinEnergy;
  if(kineticE<lowKEnergy)kineticE=lowKEnergy;
  if(kineticE>highKEnergy)kineticE=highKEnergy;
  kineticE /= eV;
    
  if(p==G4Electron::Electron())        
    {
      MyValue= new G4DataInterpolation(ener,TCSEForThisAtom,106,0.0,0);  
      InterpolatedValue = MyValue ->CubicSplineInterpolation(std::log(kineticE));
      delete  MyValue;
      summ00 = std::exp(InterpolatedValue);  
      MyValue= new G4DataInterpolation(ener,FTCSEForThisAtom,106,0.0,0);  
      InterpolatedValue = MyValue ->CubicSplineInterpolation(std::log(kineticE));
      delete  MyValue;
      summ10 = std::exp(InterpolatedValue);
    }
  if(p==G4Positron::Positron())        
    {
      MyValue= new G4DataInterpolation(ener,TCSPForThisAtom,106,0.0,0);  
      InterpolatedValue = MyValue ->CubicSplineInterpolation(std::log(kineticE));
      delete  MyValue;
      summ00 = std::exp(InterpolatedValue);  
      MyValue= new G4DataInterpolation(ener,FTCSPForThisAtom,106,0.0,0);  
      InterpolatedValue = MyValue ->CubicSplineInterpolation(std::log(kineticE));
      delete  MyValue;
      summ10 = std::exp(InterpolatedValue);
    }

  summ00 *=barn;
  summ10 *=barn;

  Lam0=1./((1.+1./Z)*summ00);
  Lam1=1./((1.+1./Z)*summ10);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//t->g->t step transformations taken from Ref.6 
G4double G4GoudsmitSaundersonMscModel::ComputeTruePathLengthLimit(const G4Track& track,
								  G4PhysicsTable* theTable,
								  G4double currentMinimalStep)
{
  tPathLength = currentMinimalStep;
  G4StepPoint* sp = track.GetStep()->GetPreStepPoint();
  G4StepStatus stepStatus = sp->GetStepStatus();

  const G4DynamicParticle* dp = track.GetDynamicParticle();

  if(stepStatus == fUndefined) {
    inside = false;
    insideskin = false;
    tlimit = geombig;
    SetParticle( dp->GetDefinition() );
  }

  theLambdaTable = theTable;
  currentCouple = track.GetMaterialCutsCouple();
  currentMaterialIndex = currentCouple->GetIndex();
  currentKinEnergy = dp->GetKineticEnergy();
  currentRange = 
    theManager->GetRangeFromRestricteDEDX(particle,currentKinEnergy,currentCouple);

  lambda1 = GetLambda(currentKinEnergy);

  // stop here if small range particle
  if(inside) return tPathLength;            
  
  if(tPathLength > currentRange) tPathLength = currentRange;

  G4double presafety = sp->GetSafety();

  // far from geometry boundary
  if(currentRange < presafety)
    {
      inside = true;
      return tPathLength;  
    }

  // standard  version
  //
  if (steppingAlgorithm == fUseDistanceToBoundary)
    {
      //compute geomlimit and presafety 
      G4double geomlimit = ComputeGeomLimit(track, presafety, tPathLength);
   
      // is far from boundary
      if(currentRange <= presafety)
	{
	  inside = true;
	  return tPathLength;   
	}

      smallstep += 1.;
      insideskin = false;

      if((stepStatus == fGeomBoundary) || (stepStatus == fUndefined))
	{
          rangeinit = currentRange;
          if(stepStatus == fUndefined) smallstep = 1.e10;
          else  smallstep = 1.;

	  // constraint from the geometry 
	  if((geomlimit < geombig) && (geomlimit > geommin))
	    {
	      if(stepStatus == fGeomBoundary)  
	        tgeom = geomlimit/facgeom;
	      else
	        tgeom = 2.*geomlimit/facgeom;
	    }
            else
              tgeom = geombig;

          //define stepmin here (it depends on lambda!)
          //rough estimation of lambda_elastic/lambda_transport
          G4double rat = currentKinEnergy/MeV ;
          rat = 1.e-3/(rat*(10.+rat)) ;
          //stepmin ~ lambda_elastic
          stepmin = rat*lambda1;
          skindepth = skin*stepmin;

          //define tlimitmin
          tlimitmin = 10.*stepmin;
          if(tlimitmin < tlimitminfix) tlimitmin = tlimitminfix;

        }

      //step limit 
      tlimit = facrange*rangeinit;              
      if(tlimit < facsafety*presafety)
        tlimit = facsafety*presafety; 

      //lower limit for tlimit
      if(tlimit < tlimitmin) tlimit = tlimitmin;

      if(tlimit > tgeom) tlimit = tgeom;

      // shortcut
      if((tPathLength < tlimit) && (tPathLength < presafety) &&
         (smallstep >= skin) && (tPathLength < geomlimit-0.999*skindepth))
	return tPathLength;   

      // step reduction near to boundary
      if(smallstep < skin)
	{
	  tlimit = stepmin;
	  insideskin = true;
	}
      else if(geomlimit < geombig)
	{
	  if(geomlimit > skindepth)
	    {
	      if(tlimit > geomlimit-0.999*skindepth)
		tlimit = geomlimit-0.999*skindepth;
	    }
	  else
	    {
	      insideskin = true;
	      if(tlimit > stepmin) tlimit = stepmin;
	    }
	}

      if(tlimit < stepmin) tlimit = stepmin;

      if(tPathLength > tlimit) tPathLength = tlimit  ; 

    }
    // for 'normal' simulation with or without magnetic field 
    //  there no small step/single scattering at boundaries
  else if(steppingAlgorithm == fUseSafety)
    {
      // compute presafety again if presafety <= 0 and no boundary
      // i.e. when it is needed for optimization purposes
      if((stepStatus != fGeomBoundary) && (presafety < tlimitminfix)) 
        presafety = ComputeSafety(sp->GetPosition(),tPathLength); 

      // is far from boundary
      if(currentRange < presafety)
        {
          inside = true;
          return tPathLength;  
        }

      if((stepStatus == fGeomBoundary) || (stepStatus == fUndefined))
      {
        rangeinit = currentRange;
        fr = facrange;
        // 9.1 like stepping for e+/e- only (not for muons,hadrons)
        if(mass < masslimite) 
        {
          if(lambda1 > currentRange)
            rangeinit = lambda1;
          if(lambda1 > lambdalimit)
            fr *= 0.75+0.25*lambda1/lambdalimit;
        }

        //lower limit for tlimit
        G4double rat = currentKinEnergy/MeV ;
        rat = 1.e-3/(rat*(10.+rat)) ;
        tlimitmin = 10.*lambda1*rat;
        if(tlimitmin < tlimitminfix) tlimitmin = tlimitminfix;
      }
      //step limit
      tlimit = fr*rangeinit;               

      if(tlimit < facsafety*presafety)
        tlimit = facsafety*presafety;

      //lower limit for tlimit
      if(tlimit < tlimitmin) tlimit = tlimitmin;

      if(tPathLength > tlimit) tPathLength = tlimit;
    }
  
  // version similar to 7.1 (needed for some experiments)
  else
    {
      if (stepStatus == fGeomBoundary)
	{
	  if (currentRange > lambda1) tlimit = facrange*currentRange;
	  else                        tlimit = facrange*lambda1;

	  if(tlimit < tlimitmin) tlimit = tlimitmin;
	  if(tPathLength > tlimit) tPathLength = tlimit;
	}
    }
  return tPathLength ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4GoudsmitSaundersonMscModel::ComputeGeomPathLength(G4double)
{
  par1 = -1. ;  
  par2 = par3 = 0. ;  

  //  do the true -> geom transformation
  zPathLength = tPathLength;

  // z = t for very small tPathLength
  if(tPathLength < tlimitminfix) return zPathLength;

  // this correction needed to run MSC with eIoni and eBrem inactivated
  // and makes no harm for a normal run
  if(tPathLength > currentRange)
    tPathLength = currentRange ;

  G4double tau   = tPathLength/lambda1 ;

  if ((tau <= tausmall) || insideskin) {
    zPathLength  = tPathLength;
    if(zPathLength > lambda1) zPathLength = lambda1;
    return zPathLength;
  }

  G4double zmean = tPathLength;
  if (tPathLength < currentRange*dtrl) {
    if(tau < taulim) zmean = tPathLength*(1.-0.5*tau) ;
    else             zmean = lambda1*(1.-exp(-tau));
  } else if(currentKinEnergy < mass) {
    par1 = 1./currentRange ;
    par2 = 1./(par1*lambda1) ;
    par3 = 1.+par2 ;
    if(tPathLength < currentRange)
      zmean = (1.-exp(par3*log(1.-tPathLength/currentRange)))/(par1*par3) ;
    else
      zmean = 1./(par1*par3) ;
  } else {
    G4double T1 = theManager->GetEnergy(particle,currentRange-tPathLength,currentCouple);

    lambda11 = GetLambda(T1);

    par1 = (lambda1-lambda11)/(lambda1*tPathLength) ;
    par2 = 1./(par1*lambda1) ;
    par3 = 1.+par2 ;
    zmean = (1.-exp(par3*log(lambda11/lambda1)))/(par1*par3) ;
  }

  zPathLength = zmean ;
  //  sample z
  if(samplez)
  {
    const G4double  ztmax = 0.99, onethird = 1./3. ;
    G4double zt = zmean/tPathLength ;

    if (tPathLength > stepmin && zt < ztmax)              
    {
      G4double u,cz1;
      if(zt >= onethird)
      {
        G4double cz = 0.5*(3.*zt-1.)/(1.-zt) ;
        cz1 = 1.+cz ;
        G4double u0 = cz/cz1 ;
        G4double grej ;
        do {
            u = exp(log(G4UniformRand())/cz1) ;
            grej = exp(cz*log(u/u0))*(1.-u)/(1.-u0) ;
           } while (grej < G4UniformRand()) ;
      }
      else
      {
        cz1 = 1./zt-1.;
        u = 1.-exp(log(G4UniformRand())/cz1) ;
      }
      zPathLength = tPathLength*u ;
    }
  }
  if(zPathLength > lambda1) zPathLength = lambda1;

 return zPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4GoudsmitSaundersonMscModel::ComputeTrueStepLength(G4double geomStepLength)
{
  // step defined other than transportation 
  if(geomStepLength == zPathLength && tPathLength <= currentRange)
    return tPathLength;

  // t = z for very small step
  zPathLength = geomStepLength;
  tPathLength = geomStepLength;
  if(geomStepLength < tlimitminfix) return tPathLength;
  
  // recalculation
  if((geomStepLength > lambda1*tausmall) && !insideskin)
    {
      if(par1 <  0.)
	tPathLength = -lambda1*log(1.-geomStepLength/lambda1) ;
      else 
	{
	  if(par1*par3*geomStepLength < 1.)
	    tPathLength = (1.-exp(log(1.-par1*par3*geomStepLength)/par3))/par1 ;
	  else 
	    tPathLength = currentRange;
	}  
    }
  if(tPathLength < geomStepLength) tPathLength = geomStepLength;

  return tPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4GoudsmitSaundersonMscModel::LoadELSEPAXSections()
{ 
  ///////////////////////////////////////
  //Total & first transport x sections of e-/e+ from ELSEPA code
  G4String filename = "XSECTIONS.dat";

    char* path = getenv("G4LEDATA");
    if (!path)
      {
        G4String excep = "G4GoudsmitSaundersonTable: G4LEDATA environment variable not set";
        G4Exception(excep);
      }

    G4String pathString(path);
    G4String dirFile = pathString + "/msc_GS/" + filename;
    FILE *infile;
    infile = fopen(dirFile,"r"); 
    if (infile == 0)
      {
	G4String excep = "G4GoudsmitSaunderson - data files: " + dirFile + " not found";
	G4Exception(excep);
      }

    // Read parameters from tables and take logarithms
    G4float aRead;
    for(G4int i=0 ; i<106 ;i++){
	  fscanf(infile,"%f\t",&aRead);
          if(aRead > 0.0) aRead = std::log(aRead);
          else  aRead = 0.0;
          ener[i]=aRead;
    }        
    for(G4int j=0;j<103;j++){
      for(G4int i=0;i<106;i++){
	  fscanf(infile,"%f\t",&aRead);
          if(aRead > 0.0) aRead = std::log(aRead);
          else  aRead = 0.0;
	  TCSE[j][i]=aRead;
	}        
     }
    for(G4int j=0;j<103;j++){
      for(G4int i=0;i<106;i++){
	  fscanf(infile,"%f\t",&aRead);
          if(aRead > 0.0) aRead = std::log(aRead);
          else  aRead = 0.0;
	  FTCSE[j][i]=aRead;      
	}        
      }    
    for(G4int j=0;j<103;j++){
      for(G4int i=0;i<106;i++){
	  fscanf(infile,"%f\t",&aRead);
          if(aRead > 0.0) aRead = std::log(aRead);
          else  aRead = 0.0;
	  TCSP[j][i]=aRead;      
	}        
     }
    for(G4int j=0;j<103;j++){
      for(G4int i=0;i<106;i++){
	  fscanf(infile,"%f\t",&aRead);
          if(aRead > 0.0) aRead = std::log(aRead);
          else  aRead = 0.0;
	  FTCSP[j][i]=aRead;      
	}        
     }

    fclose(infile);
   //End loading XSections and Energies

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
