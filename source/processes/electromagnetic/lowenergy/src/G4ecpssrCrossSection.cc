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
//$Id: G4ecpssrCrossSection.cc,v 1.5 2008/12/18 13:01:32 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// Author: Haifa Ben Abdelouahed
//         
//
// History:
// -----------
//  21 Apr 2008   H. Ben Abdelouahed   1st implementation
//  21 Apr 2008   MGP        Major revision according to a design iteration
//
// -------------------------------------------------------------------
// Class description:
// Low Energy Electromagnetic Physics, Cross section, p ionisation, K shell
// Further documentation available from http://www.ge.infn.it/geant4/lowE
// -------------------------------------------------------------------


#include "globals.hh"
#include "G4ecpssrCrossSection.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4NistManager.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include <math.h>

G4ecpssrCrossSection::G4ecpssrCrossSection()
{ }

G4ecpssrCrossSection::~G4ecpssrCrossSection()
{ }

//---------------------------------this "ExpIntFunction" function allows fast evaluation of the n order exponential integral function En(x)------

G4double G4ecpssrCrossSection::ExpIntFunction(G4int n,G4double x)

{
  G4int i;
  G4int ii;
  G4int nm1;
  G4double a;
  G4double b;
  G4double c;
  G4double d;
  G4double del;
  G4double fact;
  G4double h;
  G4double psi;
  G4double ans = 0;
  const G4double euler= 0.5772156649;
  const G4int maxit= 100;
  const G4double fpmin = 1.0e-30;
  const G4double eps = 1.0e-7;
  nm1=n-1;
  if (n<0 || x<0.0 || (x==0.0 && (n==0 || n==1)))
  G4cout << "bad arguments in ExpIntFunction" << G4endl;
  else {
       if (n==0) ans=std::exp(-x)/x;
        else {
           if (x==0.0) ans=1.0/nm1;
              else {
                   if (x > 1.0) {
                         	b=x+n;
                        	c=1.0/fpmin;
                        	d=1.0/b;
				h=d;
				for (i=1;i<=maxit;i++) {
				  a=-i*(nm1+i);
				  b +=2.0;
				  d=1.0/(a*d+b);
				  c=b+a/c;
				  del=c*d;
				  h *=del;
				      if (std::fabs(del-1.0) < eps) {
					ans=h*std::exp(-x);
					return ans;
				      }
				}
		   } else {
		     ans = (nm1!=0 ? 1.0/nm1 : -std::log(x)-euler);
		     fact=1.0;
		     for (i=1;i<=maxit;i++) {
		       fact *=-x/i;
		       if (i !=nm1) del = -fact/(i-nm1);
		       else {
			 psi = -euler;
			 for (ii=1;ii<=nm1;ii++) psi +=1.0/ii;
			 del=fact*(-std::log(x)+psi);
		       }
		       ans += del;
		       if (std::fabs(del) < std::fabs(ans)*eps) return ans;
		     }
		   }
	      }
	}
  }
return ans;
}
//-----------------------------------------------------------------------------------------------------------

 
G4double G4ecpssrCrossSection::CalculateCrossSection(G4int zTarget,G4int zIncident, G4double energyIncident)
 
 //this K-CrossSection calculation method is done according to W.Brandt and G.Lapicki, Phys.Rev.A23(1981)//					                

{

  G4NistManager* massManager = G4NistManager::Instance();   

  G4AtomicTransitionManager*  transitionManager =  G4AtomicTransitionManager::Instance();

  G4double  massIncident; 

 if (zIncident == 1)
    {
    G4Proton* aProtone = G4Proton::Proton();
    
   massIncident = aProtone->GetPDGMass(); 
    }
  else
    {
      if (zIncident == 2)
	{
	  G4Alpha* aAlpha = G4Alpha::Alpha();
	  
	   massIncident  = aAlpha->GetPDGMass(); 
	}
      else
	{ 
	  G4cout << "we can treat only Proton or Alpha incident particles " << G4endl;
	  massIncident =0.;
	}
    }
 
  G4double kBindingEnergy = transitionManager->Shell(zTarget,0)->BindingEnergy();
         
  G4double massTarget = (massManager->GetAtomicMassAmu(zTarget))*amu_c2;
 
  G4double systemMass =((massIncident*massTarget)/(massIncident+massTarget))/electron_mass_c2;//the mass of the system (projectile, target) 

 const G4double zkshell= 0.3;

  G4double screenedzTarget = zTarget-zkshell;                                 // screenedzTarget is the screened nuclear charge of the target

  const G4double rydbergMeV= 13.6e-6;
     
  G4double tetaK = kBindingEnergy/((screenedzTarget*screenedzTarget)*rydbergMeV);  //tetaK denotes the reduced binding energy of the electron

  const G4double bohrPow2Barn=(Bohr_radius*Bohr_radius)/barn ;        
  
  G4double sigma0 = 8.*pi*(zIncident*zIncident)*bohrPow2Barn*std::pow(screenedzTarget,-4.);     //sigma0 is the initial cross section of K shell at stable state

  //---------------------------------------------------------------------------------------------------------------------

 G4double velocity = CalculateVelocity( zTarget, zIncident, energyIncident);   //is the scaled velocity parameter of the system  

  //---------------------------------------------------------------------------------------------------------------------
  
 const G4double kAnalyticalApproximation= 1.5; 
 
  G4double x = kAnalyticalApproximation/velocity;
 
  G4double electrIonizationEnergy;                                         
                                       
  if ( x<0.035)  
    {
      electrIonizationEnergy= 0.75*pi*(std::log(1./(x*x))-1.);  
    }
  else 
    { 
      if ( x<3.) 
	{ 
	  electrIonizationEnergy =std::exp(-2.*x)/(0.031+(0.213*std::pow(x,0.5))+(0.005*x)-(0.069*std::pow(x,3./2.))+(0.324*x*x));
	}
     
      else  
	{
	 electrIonizationEnergy =2.*std::exp(-2.*x)/std::pow(x,1.6); }
    }

  G4double hFunction =(electrIonizationEnergy*2.)/(tetaK*std::pow(velocity,3)); //hFunction represents the correction for polarization effet
    
  G4double gFunction = (1.+(9.*velocity)+(31.*velocity*velocity)+(98.*std::pow(velocity,3.))+(12.*std::pow(velocity,4.))+(25.*std::pow(velocity,5.))
			+(4.2*std::pow(velocity,6.))+(0.515*std::pow(velocity,7.)))/std::pow(1.+velocity,9.); //gFunction represents the correction for binding effet
 
  //-----------------------------------------------------------------------------------------------------------------------------

  G4double sigmaPSS = 1.+(((2.*zIncident)/(screenedzTarget*tetaK))*(gFunction-hFunction)); //describes the perturbed stationnairy state of the affected atomic electon
 
 //----------------------------------------------------------------------------------------------------------------------------
  
  const G4double cNaturalUnit= 1/fine_structure_const;  // it's the speed of light according to Atomic-Unit-System
  
  G4double ykFormula=0.4*(screenedzTarget/cNaturalUnit)*(screenedzTarget/cNaturalUnit)/(velocity/sigmaPSS);
 
  G4double relativityCorrection = std::pow((1.+(1.1*ykFormula*ykFormula)),0.5)+ykFormula;// the relativistic correction parameter

  G4double reducedVelocity = velocity*std::pow(relativityCorrection,0.5);  // presents the reduced collision velocity parameter 

  G4double universalFunction = (std::pow(2.,9.)/45.)*std::pow(reducedVelocity/sigmaPSS,8.)*std::pow((1.+(1.72*(reducedVelocity/sigmaPSS)*(reducedVelocity/sigmaPSS))),-4.);// is the reduced universal cross section

  //----------------------------------------------------------------------------------------------------------------------

  G4double sigmaPSSR = (sigma0/(sigmaPSS*tetaK))*universalFunction; //sigmaPSSR is the straight-line K-shell ionization cross section
  
  //-----------------------------------------------------------------------------------------------------------------------

  G4double pssDeltaK = (4./(systemMass*sigmaPSS*tetaK))*(sigmaPSS/velocity)*(sigmaPSS/velocity);

  G4double energyLoss = std::pow(1-pssDeltaK,0.5); //energyLoss incorporates the straight-line energy-loss 

  G4double energyLossFunction = (std::pow(2.,-9)/8.)*((((9.*energyLoss)-1.)*std::pow(1.+energyLoss,9.))+(((9.*energyLoss)+1.)*std::pow(1.-energyLoss,9.)));//energy loss function 

  //----------------------------------------------------------------------------------------------------------------------------------------------

  G4double coulombDeflection = (4.*pi*zIncident/systemMass)*std::pow(tetaK*sigmaPSS,-2.)*std::pow(velocity/sigmaPSS,-3.)*(zTarget/screenedzTarget); //incorporates Coulomb deflection parameter 
 
  G4double cParameter = 2.*coulombDeflection/(energyLoss*(energyLoss+1.));
  

  G4double coulombDeflectionFunction = 9.*ExpIntFunction(10,cParameter);                         //this function describes Coulomb-deflection effect

  //--------------------------------------------------------------------------------------------------------------------------------------------------
 

  G4double crossSection =  energyLossFunction* coulombDeflectionFunction*sigmaPSSR;  //this ECPSSR cross section is estimated at perturbed-stationnairy-state(PSS)
                                                                                    //and it's reduced by the energy-loss(E),the Coulomb deflection(C),
                                                                                   //and the relativity(R) effects

  //--------------------------------------------------------------------------------------------------------------------------------------------------   

  return crossSection;
}

G4double G4ecpssrCrossSection::CalculateVelocity(G4int zTarget, G4int zIncident,  G4double energyIncident) 
			                     
{  

  G4AtomicTransitionManager*  transitionManager =  G4AtomicTransitionManager::Instance();

  G4double kBindingEnergy = (transitionManager->Shell(zTarget,0)->BindingEnergy())/MeV;
 
 G4double  massIncident; 

 if (zIncident == 1)
    {
    G4Proton* aProtone = G4Proton::Proton();
    
   massIncident = aProtone->GetPDGMass(); 
    }
  else
    {
      if (zIncident == 2)
	{
	  G4Alpha* aAlpha = G4Alpha::Alpha();
	  
	   massIncident  = aAlpha->GetPDGMass(); 
	}
      else
	{ 
	  G4cout << "we can treat only Proton or Alpha incident particles " << G4endl;
	  massIncident =0.;
	}
    }

 const G4double zkshell= 0.3;      
  
 G4double screenedzTarget = zTarget- zkshell;                                  

 const G4double rydbergMeV= 13.6e-6;
 
G4double tetaK = kBindingEnergy/(screenedzTarget*screenedzTarget*rydbergMeV);            
  
G4double velocity =(2./(tetaK*screenedzTarget))*std::pow(((energyIncident*electron_mass_c2)/(massIncident*rydbergMeV)),0.5);

  return velocity;
}

