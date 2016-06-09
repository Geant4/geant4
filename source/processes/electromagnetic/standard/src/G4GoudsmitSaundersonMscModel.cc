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
// $Id: G4GoudsmitSaundersonMscModel.cc,v 1.20 2009/12/16 17:50:11 gunter Exp $
// GEANT4 tag $Name: geant4-09-03 $
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
// 15.04.2009 O.Kadri: cleanup: discard no scattering and single scattering theta 
//                     sampling from SampleCosineTheta() which means the splitting 
//                     step into two sub-steps occur only for msc regime
//
// 12.06.2009 O.Kadri: linear log-log extrapolation of lambda0 & lambda1 between 1 GeV - 100 TeV
//                     adding a theta min limit due to screening effect of the atomic nucleus
// 26.08.2009 O.Kadri: Cubic Spline interpolation was replaced with polynomial method
//                     within CalculateIntegrals method
// 05.10.2009 O.Kadri: tuning small angle theta distributions
//                     assuming the case of lambdan<1 as single scattering regime
//                     tuning theta sampling for theta below the screening angle
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
//REFERENCES:
//Ref.1:E. Benedito et al.,"Mixed simulation ... cross-sections", NIMB 174 (2001) pp 91-110;
//Ref.2:I. Kawrakow et al.,"On the condensed ... transport",NIMB 142 (1998) pp 253-280;
//Ref.3:I. Kawrakow et al.,"On the representation ... calculations",NIMB 134 (1998) pp 325-336;
//Ref.4:Bielajew et al.,".....", NIMB 173 (2001) 332-343;
//Ref.5:F. Salvat et al.,"ELSEPA--Dirac partial ...molecules", Comp.Phys.Comm.165 (2005) pp 157-190;
//Ref.6:G4UrbanMscModel G4_v9.1Ref09; 
//Ref.7:G4eCoulombScatteringModel G4_v9.1Ref09.
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4GoudsmitSaundersonTable.hh"

#include "G4ParticleChangeForMSC.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4DynamicParticle.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4LossTableManager.hh"
#include "G4Track.hh"
#include "G4PhysicsTable.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"

using namespace std;

G4double G4GoudsmitSaundersonMscModel::ener[] = {-1.};
G4double G4GoudsmitSaundersonMscModel::TCSE[103][106] ;
G4double G4GoudsmitSaundersonMscModel::FTCSE[103][106] ;
G4double G4GoudsmitSaundersonMscModel::TCSP[103][106] ;
G4double G4GoudsmitSaundersonMscModel::FTCSP[103][106] ;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4GoudsmitSaundersonMscModel::G4GoudsmitSaundersonMscModel(const G4String& nam)
  : G4VMscModel(nam),lowKEnergy(0.1*keV),highKEnergy(100.*TeV),isInitialized(false)
{ 
  fr=0.02,rangeinit=0.,masslimite=0.6*MeV,
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
  fParticleChange = GetParticleChangeForMSC();
  InitialiseSafetyHelper();
  isInitialized=true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4GoudsmitSaundersonMscModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition* p,
                       G4double kineticEnergy,G4double Z, G4double, G4double, G4double)
{  
  //Build cross section table : Taken from Ref.7
  G4double cs=0.0;
  G4double kinEnergy = kineticEnergy;
  if(kinEnergy<lowKEnergy) kinEnergy=lowKEnergy;
  if(kinEnergy>highKEnergy)kinEnergy=highKEnergy;

  G4double value0,value1;
  CalculateIntegrals(p,Z,kinEnergy,value0,value1);
  
  if(value1 > 0.0) cs = 1./value1;

  return cs;
}  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4GoudsmitSaundersonMscModel::SampleScattering(const G4DynamicParticle* dynParticle,
					       G4double safety)
{    
  G4double kineticEnergy = dynParticle->GetKineticEnergy();
  if((kineticEnergy <= 0.0) || (tPathLength <= tlimitminfix)) return ;

  G4double cosTheta1,sinTheta1,cosTheta2,sinTheta2;
  G4double phi1,phi2,cosPhi1=1.0,sinPhi1=0.0,cosPhi2=1.0,sinPhi2=0.0;
  G4double q1,Gamma,Eta,delta,nu,nu0,nu1,nu2;

  ///////////////////////////////////////////
  // Effective energy and path-length from Eq. 4.7.15+16 of Ref.4
  G4double  eloss = theManager->GetEnergy(particle,tPathLength,currentCouple);
  if(eloss>0.5*kineticEnergy)eloss=kineticEnergy-eloss;//exchange possibility between target atomic e- and incident particle
  G4double ee       = kineticEnergy - 0.5*eloss;
  G4double ttau     = ee/electron_mass_c2;
  G4double ttau2    = ttau*ttau;
  G4double epsilonpp= eloss/ee;
  G4double cst1=epsilonpp*epsilonpp*(6+10*ttau+5*ttau2)/(24*ttau2+48*ttau+72);

  kineticEnergy *= (1 - cst1);
  ///////////////////////////////////////////
  // additivity rule for mixture and compound xsection calculation
  const G4Material* mat = currentCouple->GetMaterial();
  G4int nelm = mat->GetNumberOfElements();
  const G4ElementVector* theElementVector = mat->GetElementVector();
  const G4double* theFraction = mat->GetFractionVector();
  G4double atomPerVolume = mat->GetTotNbOfAtomsPerVolume();
  G4double llambda0=0.,llambda1=0.;
  for(G4int i=0;i<nelm;i++)
    {
      G4double l0,l1;
      CalculateIntegrals(particle,(*theElementVector)[i]->GetZ(),kineticEnergy,l0,l1);
      llambda0 += (theFraction[i]/l0);
      llambda1 += (theFraction[i]/l1);
    } 
  if(llambda0>DBL_MIN) llambda0 =1./llambda0;
  if(llambda1>DBL_MIN) llambda1 =1./llambda1;
  G4double g1=0.0;
  if(llambda1>DBL_MIN) g1 = llambda0/llambda1;

  G4double x1,x0;
  x0=g1/2.;
  do
    {  
      x1 = x0-(x0*((1.+x0)*log(1.+1./x0)-1.0)-g1/2.)/( (1.+2.*x0)*log(1.+1./x0)-2.0);// x1=x0-f(x0)/f'(x0)
      delta = std::abs( x1 - x0 );    
      x0 = x1;  // new approximation becomes the old approximation for the next iteration
    } while (delta > 1e-10);
  G4double scrA = x1;

  G4double us=0.0,vs=0.0,ws=1.0,x_coord=0.0,y_coord=0.0,z_coord=1.0;
  G4double lambdan=0.;
  G4bool mscatt=false,noscatt=false;

  if(llambda0>0.)lambdan=atomPerVolume*tPathLength/llambda0;
  if((lambdan<=1.0e-12))return;

  G4double epsilon1=G4UniformRand();
  G4double expn = exp(-lambdan);
  if((epsilon1<expn)||insideskin)// no scattering 
    {noscatt=true;}
  else if((epsilon1<((1.+lambdan)*expn)||(lambdan<1.)))
    {
      mscatt=false;
      ws=G4UniformRand();
      ws= 1.-2.*scrA*ws/(1.-ws + scrA);
      if(acos(ws)<sqrt(scrA))//small angle approximation for theta less than screening angle
      {G4int i=0;
      do{i++;
      ws=1.+0.5*atomPerVolume*tPathLength*log(G4UniformRand())/llambda1;
      }while((fabs(ws)>1.)&&(i<20));//i<20 to avoid time consuming during the run
      if(i==19)ws=cos(sqrt(scrA));
      }
      G4double phi0=twopi*G4UniformRand(); 
      us=sqrt(1.-ws*ws)*cos(phi0);
      vs=sqrt(1.-ws*ws)*sin(phi0);
      G4double rr=G4UniformRand();
      x_coord=(rr*us);
      y_coord=(rr*vs);
      z_coord=((1.-rr)+rr*ws);
    }
  else
    {
      mscatt=true;
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
    }
    
  G4ThreeVector oldDirection = dynParticle->GetMomentumDirection();
  G4ThreeVector newDirection(us,vs,ws);
  newDirection.rotateUz(oldDirection);
  fParticleChange->ProposeMomentumDirection(newDirection);
  
  if((safety > tlimitminfix)&&(latDisplasment))
    {  

      if(mscatt)
	{
	  if(scrA<DBL_MIN)scrA=DBL_MIN;
	  if(llambda1<DBL_MIN)llambda1=DBL_MIN;
	  q1 = 2.*scrA*((1. + scrA)*log(1. + 1./scrA) - 1.);
	  if(q1<DBL_MIN)q1=DBL_MIN;
	  Gamma  = 6.*scrA*(1. + scrA)*((1. + 2.*scrA)*log(1. + 1./scrA) - 2.)/q1;
	  Eta    = atomPerVolume*tPathLength/llambda1;      
	  delta  = 0.90824829 - Eta*(0.102062073-Gamma*0.026374715);

	  nu = G4UniformRand(); 
	  nu = std::sqrt(nu);
	  nu0 = (1.0 - nu)/2.;
	  nu1 = nu*delta;
	  nu2 = nu*(1.0-delta);
	  x_coord=(nu1*sinTheta1*cosPhi1+nu2*sinTheta2*(cosPhi1*cosPhi2-cosTheta1*sinPhi1*sinPhi2)+nu0*us);
	  y_coord=(nu1*sinTheta1*sinPhi1+nu2*sinTheta2*(sinPhi1*cosPhi2+cosTheta1*cosPhi1*sinPhi2)+nu0*vs);
	  z_coord=(nu0+nu1*cosTheta1+nu2*cosTheta2+ nu0*ws)  ;
	}

      // displacement is computed relatively to the end point
      if(!noscatt)z_coord -= 1.0;//for noscatt zcoord z_coord !=0.
      G4double rr = sqrt(x_coord*x_coord+y_coord*y_coord+z_coord*z_coord);
      G4double r  = rr*zPathLength;
      /*
      G4cout << "G4GS::SampleSecondaries: e(MeV)= " << kineticEnergy
	     << " sinTheta= " << sqrt(1.0 - ws*ws) << " r(mm)= " << r
	     << " trueStep(mm)= " << tPathLength
	     << " geomStep(mm)= " << zPathLength
	     << G4endl;
      */
      if(r > tlimitminfix) {

        G4ThreeVector latDirection(x_coord/rr,y_coord/rr,z_coord/rr);
        latDirection.rotateUz(oldDirection);

	ComputeDisplacement(fParticleChange, latDirection, r, safety);
      }     
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4GoudsmitSaundersonMscModel::SampleCosineTheta(G4double lambdan, G4double scrA,
						G4double &cost, G4double &sint)
{
  G4double u,Qn1,r1,tet;
  G4double xi=0.;
  Qn1=2.* lambdan *scrA*((1.+scrA)*log(1.+1./scrA)-1.);

if (Qn1<0.001)xi=-0.5*Qn1*log(G4UniformRand());
else if(Qn1>0.5)xi=2.*G4UniformRand();//isotropic distribution
else
{
      // procedure described by Benedito in Ref.1
      do{
	r1=G4UniformRand();
	u=GSTable->SampleTheta(lambdan,scrA,G4UniformRand());
	xi = 2.*u;
	tet=acos(1.-xi);
      }while(tet*r1*r1>sin(tet));
}

  if(xi<0.)xi=0.;
  if(xi>2.)xi=2.;
  cost=(1. - xi);
  sint=sqrt(xi*(2.-xi));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Polynomial log-log interpolation of Lambda0 and Lambda1 between 100 eV - 1 GeV
// linear log-log extrapolation between 1 GeV - 100 TeV

void 
G4GoudsmitSaundersonMscModel::CalculateIntegrals(const G4ParticleDefinition* p,G4double Z, 
						 G4double kinEnergy,G4double &Lam0,
						 G4double &Lam1)
{ 
  G4double summ00=0.0;
  G4double summ10=0.0;
  G4double x1,x2,y1,y2,acoeff,bcoeff;
  G4double kineticE = kinEnergy;
  if(kineticE<lowKEnergy)kineticE=lowKEnergy;
  if(kineticE>highKEnergy)kineticE=highKEnergy;
  kineticE /= eV;
  G4double logE=log(kineticE);
  
  G4int  iZ = G4int(Z);
  if(iZ > 103) iZ = 103;

  G4int enerInd=0;
  for(G4int i=1;i<106;i++)
  {
  if((logE>=ener[i])&&(logE<ener[i+1])){enerInd=i;break;}
  }

  if(p==G4Electron::Electron())        
    {
    if(kineticE<=1.0e+9)//Interpolation of the form y=ax²+b
      {
	x1=ener[enerInd];
	x2=ener[enerInd+1];       
	y1=TCSE[iZ-1][enerInd];
	y2=TCSE[iZ-1][enerInd+1];
	acoeff=(y2-y1)/(x2*x2-x1*x1);
	bcoeff=y2-acoeff*x2*x2;
	summ00=acoeff*logE*logE+bcoeff;
	summ00 =exp(summ00);
	y1=FTCSE[iZ-1][enerInd];
	y2=FTCSE[iZ-1][enerInd+1];
	acoeff=(y2-y1)/(x2*x2-x1*x1);
	bcoeff=y2-acoeff*x2*x2;
	summ10=acoeff*logE*logE+bcoeff;
	summ10 =exp(summ10);
      }
    else  //Interpolation of the form y=ax+b
      {  
	x1=ener[104];
	x2=ener[105];       
	y1=TCSE[iZ-1][104];
	y2=TCSE[iZ-1][105];
	summ00=(y2-y1)*(logE-x1)/(x2-x1)+y1;
	summ00 =exp(summ00);
	y1=FTCSE[iZ-1][104];
	y2=FTCSE[iZ-1][105];
	summ10=(y2-y1)*(logE-x1)/(x2-x1)+y1;
	summ10 =exp(summ10);
      }
    }
  if(p==G4Positron::Positron())        
    {
    if(kinEnergy<=1.0e+9)
      {
	x1=ener[enerInd];
	x2=ener[enerInd+1];       
	y1=TCSP[iZ-1][enerInd];
	y2=TCSP[iZ-1][enerInd+1];
	acoeff=(y2-y1)/(x2*x2-x1*x1);
	bcoeff=y2-acoeff*x2*x2;
	summ00=acoeff*logE*logE+bcoeff;
	summ00 =exp(summ00);
	y1=FTCSP[iZ-1][enerInd];
	y2=FTCSP[iZ-1][enerInd+1];
	acoeff=(y2-y1)/(x2*x2-x1*x1);
	bcoeff=y2-acoeff*x2*x2;
	summ10=acoeff*logE*logE+bcoeff;
	summ10 =exp(summ10);
      }
    else
      {  
	x1=ener[104];
	x2=ener[105];       
	y1=TCSP[iZ-1][104];
	y2=TCSP[iZ-1][105];
	summ00=(y2-y1)*(logE-x1)/(x2-x1)+y1;
	summ00 =exp(summ00);
	y1=FTCSP[iZ-1][104];
	y2=FTCSP[iZ-1][105];
	summ10=(y2-y1)*(logE-x1)/(x2-x1)+y1;
	summ10 =exp(summ10);
      }
    }
    
  summ00 *=barn;
  summ10 *=barn;

  Lam0=1./((1.+1./Z)*summ00);
  Lam1=1./((1.+1./Z)*summ10);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//t->g->t step transformations taken from Ref.6 

G4double 
G4GoudsmitSaundersonMscModel::ComputeTruePathLengthLimit(const G4Track& track,
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

  //G4cout << "G4GS::StepLimit tPathLength= " 
  //	 <<tPathLength<<" safety= " << presafety
  //       << " range= " <<currentRange<< " lambda= "<<lambda1
  //	 << " Alg: " << steppingAlgorithm <<G4endl;

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

	  //G4cout << "rangeinit= " << rangeinit << " stepmin= " << stepmin
	  //	 << " tlimitmin= " << tlimitmin << " geomlimit= " << geomlimit <<G4endl;
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

        }

      //step limit 
      tlimit = facrange*rangeinit;              
      if(tlimit < facsafety*presafety)
        tlimit = facsafety*presafety; 

      //lower limit for tlimit
      if(tlimit < tlimitmin) tlimit = tlimitmin;

      if(tlimit > tgeom) tlimit = tgeom;

      //G4cout << "tgeom= " << tgeom << " geomlimit= " << geomlimit  
      //      << " tlimit= " << tlimit << " presafety= " << presafety << G4endl;

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

      if(tPathLength > tlimit) tPathLength = tlimit; 

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
  //G4cout << "tPathLength= " << tPathLength  
  // << " currentMinimalStep= " << currentMinimalStep << G4endl;
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
    G4double T1 = theManager->GetEnergy(particle,currentRange-tPathLength,
					currentCouple);

    lambda11 = GetLambda(T1);

    par1 = (lambda1-lambda11)/(lambda1*tPathLength) ;
    par2 = 1./(par1*lambda1) ;
    par3 = 1.+par2 ;
    zmean = (1.-exp(par3*log(lambda11/lambda1)))/(par1*par3) ;
  }

  zPathLength = zmean ;
  //  sample z
  if(samplez) {

    const G4double  ztmax = 0.99, onethird = 1./3. ;
    G4double zt = zmean/tPathLength ;

    if (tPathLength > stepmin && zt < ztmax) {

      G4double u,cz1;
      if(zt >= onethird) {

        G4double cz = 0.5*(3.*zt-1.)/(1.-zt) ;
        cz1 = 1.+cz ;
        G4double u0 = cz/cz1 ;
        G4double grej ;
        do {
	  u = exp(log(G4UniformRand())/cz1) ;
	  grej = exp(cz*log(u/u0))*(1.-u)/(1.-u0) ;
	} while (grej < G4UniformRand()) ;

      } else {
        cz1 = 1./zt-1.;
        u = 1.-exp(log(G4UniformRand())/cz1) ;
      }
      zPathLength = tPathLength*u ;
    }
  }
  if(zPathLength > lambda1) zPathLength = lambda1;
  //G4cout << "zPathLength= " << zPathLength << " lambda1= " << lambda1 << G4endl;

  return zPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double 
G4GoudsmitSaundersonMscModel::ComputeTrueStepLength(G4double geomStepLength)
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
  //G4cout << "tPathLength= " << tPathLength << " step= " << geomStepLength << G4endl;

  return tPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Total & first transport x sections for e-/e+ generated from ELSEPA code

void G4GoudsmitSaundersonMscModel::LoadELSEPAXSections()
{ 
  G4String filename = "XSECTIONS.dat";

  char* path = getenv("G4LEDATA");
  if (!path)
    {
      G4String excep = "G4GoudsmitSaundersonTable: G4LEDATA environment variable not set properly";
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
    if(aRead > 0.0) aRead = log(aRead);
    else  aRead = 0.0;
    ener[i]=aRead;
  }        
  for(G4int j=0;j<103;j++){
    for(G4int i=0;i<106;i++){
      fscanf(infile,"%f\t",&aRead);
      if(aRead > 0.0) aRead = log(aRead);
      else  aRead = 0.0;
      TCSE[j][i]=aRead;
    }        
  }
  for(G4int j=0;j<103;j++){
    for(G4int i=0;i<106;i++){
      fscanf(infile,"%f\t",&aRead);
      if(aRead > 0.0) aRead = log(aRead);
      else  aRead = 0.0;
      FTCSE[j][i]=aRead;      
    }        
  }    
  for(G4int j=0;j<103;j++){
    for(G4int i=0;i<106;i++){
      fscanf(infile,"%f\t",&aRead);
      if(aRead > 0.0) aRead = log(aRead);
      else  aRead = 0.0;
      TCSP[j][i]=aRead;      
    }        
  }
  for(G4int j=0;j<103;j++){
    for(G4int i=0;i<106;i++){
      fscanf(infile,"%f\t",&aRead);
      if(aRead > 0.0) aRead = log(aRead);
      else  aRead = 0.0;
      FTCSP[j][i]=aRead;      
    }        
  }

  fclose(infile);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
