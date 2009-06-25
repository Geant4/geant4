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
//$Id: G4ecpssrLiCrossSection.cc,v 1.4 2009-06-25 15:52:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Haifa Ben Abdelouahed
//         
//
// History:
// -----------
//  23 Apr 2008   H. Ben Abdelouahed   1st implementation
//  28 Apr 2008   MGP        Major revision according to a design iteration
//  29 Apr 2009   ALF Updated Desing for Integration
//  02 May 2009   ALF + Haifa L1,L2,L3 Extensions
//
// -------------------------------------------------------------------
// Class description:
// Low Energy Electromagnetic Physics, Cross section, p and alpha ionisation, L shell
// Further documentation available from http://www.ge.infn.it/geant4/lowE
// -------------------------------------------------------------------


#include "globals.hh"
#include "G4ecpssrLiCrossSection.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4NistManager.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include <math.h>

G4ecpssrLiCrossSection::G4ecpssrLiCrossSection()
{ 


 // Storing FLi data needed for 0.2 to 3.0  velocities region

    char *path = getenv("G4LEDATA");
 
    if (!path)
    G4Exception("G4ecpssrLCrossSection::CalculateCrossSection: G4LEDDATA environment variable not set");

     std::ostringstream fileName1;
     std::ostringstream fileName2;
    fileName1 << path << "/FL1.dat";
    fileName2 << path << "/FL2.dat";
   
    std::ifstream FL1(fileName1.str().c_str());
    std::ifstream FL2(fileName1.str().c_str());
   
   
    if (!FL1) G4Exception("G4ecpssrLCrossSection::CalculateCrossSection: error opening FL1 data file");
    if (!FL2) G4Exception("G4ecpssrLCrossSection::CalculateCrossSection: error opening FL2 data file");
  
 dummyVec.push_back(0.);

    while(!FL1.eof())
    {
	double x1;
	double y1;
	
	FL1>>x1>>y1;

	//  Mandatory vector initialization
        if (x1 != dummyVec.back()) 
        { 
          dummyVec.push_back(x1); 
          aVecMap[x1].push_back(-1.);
        }
	  
        FL1>>FL1Data[x1][y1];

        if (y1 != aVecMap[x1].back()) aVecMap[x1].push_back(y1);
    }
 while(!FL2.eof())
    {
	double x2;
	double y2;
	
	FL2>>x2>>y2;

	//  Mandatory vector initialization
        if (x2 != dummyVec.back()) 
        { 
          dummyVec.push_back(x2); 
          aVecMap[x2].push_back(-1.);
        }
	  
        FL2>>FL2Data[x2][y2];

        if (y2 != aVecMap[x2].back()) aVecMap[x2].push_back(y2);    
    }


}

G4ecpssrLiCrossSection::~G4ecpssrLiCrossSection()
{ }

//---------------------------------this "ExpIntFunction" function allows fast evaluation of the n order exponential integral function En(x)------

G4double G4ecpssrLiCrossSection::ExpIntFunction(G4int n,G4double x)

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
       if (n==0) ans=exp(-x)/x;
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
				      if (fabs(del-1.0) < eps) {
					ans=h*exp(-x);
					return ans;
				      }
				}
		   } else {
		     ans = (nm1!=0 ? 1.0/nm1 : -log(x)-euler);
		     fact=1.0;
		     for (i=1;i<=maxit;i++) {
		       fact *=-x/i;
		       if (i !=nm1) del = -fact/(i-nm1);
		       else {
			 psi = -euler;
			 for (ii=1;ii<=nm1;ii++) psi +=1.0/ii;
			 del=fact*(-log(x)+psi);
		       }
		       ans += del;
		       if (fabs(del) < fabs(ans)*eps) return ans;
		     }
		   }
	      }
	}
  }
return ans;
}
//-----------------------------------------------------------------------------------------------------------

 
G4double G4ecpssrLiCrossSection::CalculateL1CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)
 
 //this L-CrossSection calculation method is done according to W.Brandt and G.Lapicki, Phys.Rev.A23(1981)//					                

{

  G4NistManager* massManager = G4NistManager::Instance();   

  G4AtomicTransitionManager*  transitionManager =  G4AtomicTransitionManager::Instance();


  G4int zIncident = 0; 
  G4Proton* aProtone = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();

  if (massIncident == aProtone->GetPDGMass() )
  {

   zIncident = (G4int)((aProtone->GetPDGCharge())/eplus); 

   G4cout << "zincident:" << zIncident << G4endl;
  }
  else
    {
      if (massIncident == aAlpha->GetPDGMass())
	{
	  
	  zIncident  =(G4int) ((aAlpha->GetPDGCharge())/eplus); 
	  
	  G4cout << "zincident:" << zIncident << G4endl;
	}
      else
	{ 
	  G4cout << "we can treat only Proton or Alpha incident particles " << G4endl;
	  massIncident =0.;
	}
    }


  G4double l1BindingEnergy = transitionManager->Shell(zTarget,1)->BindingEnergy();
  




         
  G4double massTarget = (massManager->GetAtomicMassAmu(zTarget))*amu_c2;
 
  G4double systemMass =((massIncident*massTarget)/(massIncident+massTarget))/electron_mass_c2;//the mass of the system (projectile, target) 

  const G4double zlshell= 4.15;

  G4double screenedzTarget = zTarget-zlshell;         // screenedzTarget is the screened nuclear charge of the target

  const G4double rydbergMeV= 13.6056923e-6;

  const G4double nl= 2.;        // nl is the quantum number of the L shell 
     
  G4double tetal1 = (l1BindingEnergy*nl*nl)/((screenedzTarget*screenedzTarget)*rydbergMeV);  //tetal1 denotes the reduced L1-shell-binding-energy of the electron






  const G4double bohrPow2Barn=(Bohr_radius*Bohr_radius)/barn ;        
  
  G4double sigma0 = 8.*pi*(zIncident*zIncident)*bohrPow2Barn*pow(screenedzTarget,-4.);     //sigma0 is the initial cross section of L shell at stable state

  //---------------------------------------------------------------------------------------------------------------------

  G4double velocityl1 = CalculateVelocity(1, zTarget, massIncident, energyIncident);   //is the scaled velocity parameter of the system  



  
  //---------------------------------------------------------------------------------------------------------------------
  
  const G4double l1AnalyticalApproximation= 1.5; 



  G4double x1 = nl*l1AnalyticalApproximation/velocityl1;





  //-----------------------------------------x of l1 sub shell--------------------------------------
 
 G4double electrIonizationEnergyl1;
                                       
                                       
  if ( x1<0.035 && x1>= 0.)  
    {
      electrIonizationEnergyl1= 0.75*pi*(log(1./(x1*x1))-1.);  
    }
  else 
    { 
      if ( x1<3.&& x1>=0.035) 
	{ 
	  electrIonizationEnergyl1 =exp(-2.*x1)/(0.031+(0.213*pow(x1,0.5))+(0.005*x1)-(0.069*pow(x1,3./2.))+(0.324*x1*x1));
	}
     
      else  
	{
	  if ( x1<=11.&& x1>=3.) {
	    electrIonizationEnergyl1 =2.*exp(-2.*x1)/pow(x1,1.6);
	  }
	  else {
	    electrIonizationEnergyl1 =0.;
	      }
	}
    }


 
 

  //-------------------------------------------------------- h and g functions for l1 -------------------------------------------------


  G4double hFunctionl1 =(electrIonizationEnergyl1*2.*nl)/(tetal1*pow(velocityl1,3)); //hFunction represents the correction for polarization effet
    
 
  G4double gFunctionl1 = (1.+(9.*velocityl1)+(31.*velocityl1*velocityl1)+(49.*pow(velocityl1,3.))+(162.*pow(velocityl1,4.))+(63.*pow(velocityl1,5.))
			  +(18.*pow(velocityl1,6.))+(1.97*pow(velocityl1,7.)))/pow(1.+velocityl1,9.); //gFunction represents the correction for binding effet




  //-----------------------------------------------------------------------------------------------------------------------------

  G4double sigmaPSS_l1 = 1.+(((2.*zIncident)/(screenedzTarget*tetal1))*(gFunctionl1-hFunctionl1)); //describes the perturbed stationnairy state of the affected atomic electon
 
 

 



 //----------------------------------------------------------------------------------------------------------------------------
  //const G4double cNaturalUnit= 1/fine_structure_const;  // it's the speed of light according to Atomic-Unit-System
  //--------------------------------------------------------------------------------------------------------------

  //G4double yl1Formula=0.4*(screenedzTarget/cNaturalUnit)*(screenedzTarget/cNaturalUnit)/(nl*(velocityl1/sigmaPSS_l1));

 

  //-----------------------------------------------Relativity effect correction L1 -------------------------------------------------------------

  //G4double relativityCorrectionl1 = pow((1.+(1.1*yl1Formula*yl1Formula)),0.5)+yl1Formula;// the relativistic correction parameter
  //G4double reducedVelocityl1 = velocityl1*pow(relativityCorrectionl1,0.5);  // presents the reduced collision velocity parameter 







  //-------------------------------------------------------------------------------------------------------------------
   //------------------------------------------------------------UNIVERSAL FUNCTION ---------------------------------
  //------------------------------------------------------------------------------------------------------------
                                                  // is the reduced universal cross section

  G4double universalFunction_l1 ;

 
 
    //-------------------------------------------------------------------------------------------------------------------
   //------------------------------------------------------------ LIMITS OF ECPSSR MODEL ---------------------------
  //------------------------------------------------------------------------------------------------------------


  G4double L1etaOverTheta2 = (energyIncident*electron_mass_c2)/(massIncident*rydbergMeV*screenedzTarget*screenedzTarget)/(sigmaPSS_l1*tetal1)/(sigmaPSS_l1*tetal1);


     universalFunction_l1 = FunctionFL1((sigmaPSS_l1*tetal1), L1etaOverTheta2);
        



     //-----------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------   PSSR Li CROSS SECTION  -------------
   //-----------------------------------------------------------------------------------------------------------------------

   G4double sigmaPSSR_l1 = (sigma0/(sigmaPSS_l1*tetal1))*universalFunction_l1; //sigmaPSSR is the straight-line L-shell ionization cross section



   //----------------------------------------------------------------------------------------------------------------------
  
  G4double pssDeltal1 = (4./(systemMass*sigmaPSS_l1*tetal1))*(sigmaPSS_l1/velocityl1)*(sigmaPSS_l1/velocityl1);



  G4double energyLossl1 = pow(1-pssDeltal1,0.5); //energyLoss incorporates the straight-line energy-loss 





     //-----------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------ENERGY LOSS CORRECTION-------------
   //-----------------------------------------------------------------------------------------------------------------------

/* 
 G4double energyLossFunction_L1 = (pow(2.,-9)/8.)*((((9.*energyLossl1)-1.)*pow(1.+energyLossl1,9.))+(((9.*energyLossl1)+1.)*pow(1.-energyLossl1,9.))); 
  G4double energyLossFunction_L2 = (pow(2.,-11)/10.)*((((11.*energyLossl2)-1.)*pow(1.+energyLossl2,11.))+(((11.*energyLossl2)+1.)*pow(1.-energyLossl2,11.)));
  G4double energyLossFunction_L3 = (pow(2.,-11)/10.)*((((11.*energyLossl3)-1.)*pow(1.+energyLossl3,11.))+(((11.*energyLossl3)+1.)*pow(1.-energyLossl3,11.)));
*/


     //-----------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------COULOMB DEFLECTION CORRECTION-------------
   //-----------------------------------------------------------------------------------------------------------------------


  G4double coulombDeflectionl1 = (4.*pi*zIncident/systemMass)*pow(tetal1*sigmaPSS_l1,-2.)*pow(velocityl1/sigmaPSS_l1,-3.)*(zTarget/screenedzTarget); //incorporates Coulomb deflection parameter 




  G4double cParameterl1 = 2.*coulombDeflectionl1/(energyLossl1*(energyLossl1+1.));



  G4double coulombDeflectionFunction_l1 = 9.*ExpIntFunction(10,cParameterl1);                         //this function describes Coulomb-deflection effect



      //--------------------------------------------------------------------------------------------------------------------------------------------------
     //-----------------------------------------------------------------------------------------------------------------------
    //---------------------------------------------------------------------------   ECPSSR Li CROSS SECTION  -------------
   //-----------------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------------------------------------------
 
  /*
   G4double crossSection_L1 = energyLossFunction_L1 * coulombDeflectionFunction_l1 * sigmaPSSR_l1;  //this ECPSSR cross section is estimated at perturbed-stationnairy-state(PSS)
   G4double crossSection_L2 = energyLossFunction_L2 * coulombDeflectionFunction_l2 * sigmaPSSR_l2;         //and it's reduced by the energy-loss(E),the Coulomb deflection(C),
   G4double crossSection_L3 = energyLossFunction_L3 * coulombDeflectionFunction_l3 * sigmaPSSR_l3;              //and the relativity(R) effects
  */

   G4double crossSection_L1 =  coulombDeflectionFunction_l1 * sigmaPSSR_l1; 

 
 if (crossSection_L1 >= 0) {
    return crossSection_L1;
  }
  else {return 0;}
}

G4double G4ecpssrLiCrossSection::CalculateL2CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)
 
 //this L-CrossSection calculation method is done according to W.Brandt and G.Lapicki, Phys.Rev.A23(1981)//					                

{


  G4NistManager* massManager = G4NistManager::Instance();   

  G4AtomicTransitionManager*  transitionManager =  G4AtomicTransitionManager::Instance();


  G4int zIncident = 0; 
  G4Proton* aProtone = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();

  if (massIncident == aProtone->GetPDGMass() )
  {

   zIncident =(G4int) ((aProtone->GetPDGCharge())/eplus); 

   //   G4cout << "zincident:" << zIncident << G4endl;
    }
  else
    {
      if (massIncident == aAlpha->GetPDGMass())
	{
	  
	  zIncident  = (G4int) ((aAlpha->GetPDGCharge())/eplus); 
	  
	  //	  G4cout << "zincident:" << zIncident << G4endl;
	}
      else
	{ 
	  G4cout << "we can treat only Proton or Alpha incident particles " << G4endl;
	  massIncident =0.;
	}
    }

  G4double l2BindingEnergy = transitionManager->Shell(zTarget,2)->BindingEnergy();

  G4double massTarget = (massManager->GetAtomicMassAmu(zTarget))*amu_c2;
  G4double systemMass =((massIncident*massTarget)/(massIncident+massTarget))/electron_mass_c2;//the mass of the system (projectile, target) 
  const G4double zlshell= 4.15;
  G4double screenedzTarget = zTarget-zlshell;         // screenedzTarget is the screened nuclear charge of the target
  const G4double rydbergMeV= 13.6056923e-6;
  const G4double nl= 2.;        // nl is the quantum number of the L shell 

  G4double tetal2 = (l2BindingEnergy*nl*nl)/((screenedzTarget*screenedzTarget)*rydbergMeV); 

  const G4double bohrPow2Barn=(Bohr_radius*Bohr_radius)/barn ;        
  G4double sigma0 = 8.*pi*(zIncident*zIncident)*bohrPow2Barn*pow(screenedzTarget,-4.);     //sigma0 is the initial cross section of L shell at stable state

  G4double velocityl2 = CalculateVelocity(2,  zTarget, massIncident, energyIncident);

  const G4double l2AnalyticalApproximation= 1.25;

  G4double x2 = nl*l2AnalyticalApproximation/velocityl2;

  //---------------------------------------------------------------x of l2 sub shell------------------------------------------
 
G4double electrIonizationEnergyl2;
                                       
                                       
  if ( x2<0.035 && x2>= 0.)  
    {
      electrIonizationEnergyl2= 0.75*pi*(log(1./(x2*x2))-1.);  
    }
  else 
    { 
      if ( x2<3. && x2 >=0.035) 
	{ 
	  electrIonizationEnergyl2 =exp(-2.*x2)/(0.031+(0.213*pow(x2,0.5))+(0.005*x2)-(0.069*pow(x2,3./2.))+(0.324*x2*x2));
	}
     
      else  
	{
	  
	  if ( x2<=11.&& x2>=3.) {
	    electrIonizationEnergyl2 =2.*exp(-2.*x2)/pow(x2,1.6);
	  }
	  else {
	    electrIonizationEnergyl2=0.;
	  }
	  
	}
    }
  
  //----------------------------------------------------------------------------- h and g functions for l2  ----------------------------
  
  
  G4double hFunctionl2 =(electrIonizationEnergyl2*2.*nl)/(tetal2*pow(velocityl2,3)); //hFunction represents the correction for polarization effet
  
  G4double gFunctionl2 = (1.+(10.*velocityl2)+(45.*velocityl2*velocityl2)+(102.*pow(velocityl2,3.))+(331.*pow(velocityl2,4.))+(6.7*pow(velocityl2,5.))
			  +(58.*pow(velocityl2,6.))+(7.8*pow(velocityl2,7.))+ (0.888*pow(velocityl2,8.)) )/pow(1.+velocityl2,10.); //gFunction represents the correction for binding effet
  
  
  G4double sigmaPSS_l2 = 1.+(((2.*zIncident)/(screenedzTarget*tetal2))*(gFunctionl2-hFunctionl2)); 
  
  //----------------------------------------------------------------------------------------------------------------------------
  //const G4double cNaturalUnit= 1/fine_structure_const;  // it's the speed of light according to Atomic-Unit-System
  //--------------------------------------------------------------------------------------------------------------
  
  //G4double yl2Formula=0.15*(screenedzTarget/cNaturalUnit)*(screenedzTarget/cNaturalUnit)/(velocityl2/sigmaPSS_l2);
  
  //-----------------------------------------------Relativity effect correction L2 --------------------------------------------------------------
  
  //G4double relativityCorrectionl2 = pow((1.+(1.1*yl2Formula*yl2Formula)),0.5)+yl2Formula;// the relativistic correction parameter
  //  G4double reducedVelocityl2 = velocityl2*pow(relativityCorrectionl2,0.5);  // presents the reduced collision velocity parameter 
  
  G4double L2etaOverTheta2 = (energyIncident*electron_mass_c2)/(massIncident*rydbergMeV*screenedzTarget*screenedzTarget)/(sigmaPSS_l2*tetal2)/(sigmaPSS_l2*tetal2);
  
  
  G4double universalFunction_l2 ;
  
  
  universalFunction_l2 = FunctionFL2((sigmaPSS_l2*tetal2), L2etaOverTheta2);
  
  G4double sigmaPSSR_l2 = (sigma0/(sigmaPSS_l2*tetal2))*universalFunction_l2;
  G4double pssDeltal2 = (4./(systemMass*sigmaPSS_l2*tetal2))*(sigmaPSS_l2/velocityl2)*(sigmaPSS_l2/velocityl2);
  G4double energyLossl2 = pow(1-pssDeltal2,0.5);  
  
  G4double coulombDeflectionl2 = (4.*pi*zIncident/systemMass)*pow(tetal2*sigmaPSS_l2,-2.)*pow(velocityl2/sigmaPSS_l2,-3.)*(zTarget/screenedzTarget); //incorporates Coulomb deflection parameter 
  G4double cParameterl2 = 2.*coulombDeflectionl2/(energyLossl2*(energyLossl2+1.));
  G4double coulombDeflectionFunction_l2 = 11.*ExpIntFunction(12,cParameterl2);
  
  G4double crossSection_L2 = coulombDeflectionFunction_l2 * sigmaPSSR_l2 ;         
  
  
  if (crossSection_L2 >= 0) {
    return crossSection_L2;
  }
  else {return 0;}
  
  
}


G4double G4ecpssrLiCrossSection::CalculateL3CrossSection(G4int zTarget,G4double massIncident, G4double energyIncident)
 
 //this L-CrossSection calculation method is done according to W.Brandt and G.Lapicki, Phys.Rev.A23(1981)//					                

{


  G4NistManager* massManager = G4NistManager::Instance();   

  G4AtomicTransitionManager*  transitionManager =  G4AtomicTransitionManager::Instance();


  G4int zIncident = 0; 
  G4Proton* aProtone = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();

  if (massIncident == aProtone->GetPDGMass() )
  {

   zIncident = (G4int) ((aProtone->GetPDGCharge())/eplus); 

   //   G4cout << "zincident:" << zIncident << G4endl;
    }
  else
    {
      if (massIncident == aAlpha->GetPDGMass())
	{
	  
	  zIncident  =(G4int) ((aAlpha->GetPDGCharge())/eplus); 
	  
	  //	  G4cout << "zincident:" << zIncident << G4endl;
	}
      else
	{ 
	  G4cout << "we can treat only Proton or Alpha incident particles " << G4endl;
	  massIncident =0.;
	}
    }

  G4double l3BindingEnergy = transitionManager->Shell(zTarget,3)->BindingEnergy();
        
  G4double massTarget = (massManager->GetAtomicMassAmu(zTarget))*amu_c2;
 
  G4double systemMass =((massIncident*massTarget)/(massIncident+massTarget))/electron_mass_c2;//the mass of the system (projectile, target) 

  const G4double zlshell= 4.15;

  G4double screenedzTarget = zTarget-zlshell;         // screenedzTarget is the screened nuclear charge of the target

  const G4double rydbergMeV= 13.6056923e-6;

  const G4double nl= 2.;        // nl is the quantum number of the L shell 
     
  G4double tetal3 = (l3BindingEnergy*nl*nl)/((screenedzTarget*screenedzTarget)*rydbergMeV); 



  const G4double bohrPow2Barn=(Bohr_radius*Bohr_radius)/barn ;        
  
  G4double sigma0 = 8.*pi*(zIncident*zIncident)*bohrPow2Barn*pow(screenedzTarget,-4.);     //sigma0 is the initial cross section of L shell at stable state

  //---------------------------------------------------------------------------------------------------------------------

  G4double velocityl3 = CalculateVelocity(3, zTarget, massIncident, energyIncident);

  const G4double l3AnalyticalApproximation= 1.25;

  G4double x3 = nl*l3AnalyticalApproximation/velocityl3;


  //--------------------------------------------------------------------- x of l3 sub shell---------------------------------


 G4double electrIonizationEnergyl3;
                                     
                                       
  if ( x3<0.035 && x3>=0. )  
    {
      electrIonizationEnergyl3= 0.75*pi*(log(1./(x3*x3))-1.);  
    }
  else 
    { 
      if ( x3<3. && x3 >= 0.035) 
	{ 
	  electrIonizationEnergyl3 =exp(-2.*x3)/(0.031+(0.213*pow(x3,0.5))+(0.005*x3)-(0.069*pow(x3,3./2.))+(0.324*x3*x3));
	}
     
      else  
	{
	  if ( x3<=11.&& x3>=3.) {
	    electrIonizationEnergyl3 =2.*exp(-2.*x3)/pow(x3,1.6);
	  }
	  else {
            electrIonizationEnergyl3=0.;
	  }
	  
	}
    }

 //------------------------------------------------------------ h and g function  for l3  ---------------------------------------------


  G4double hFunctionl3 =(electrIonizationEnergyl3*2.*nl)/(tetal3*pow(velocityl3,3)); //hFunction represents the correction for polarization effet
    
 
  G4double gFunctionl3 = (1.+(10.*velocityl3)+(45.*velocityl3*velocityl3)+(102.*pow(velocityl3,3.))+(331.*pow(velocityl3,4.))+(6.7*pow(velocityl3,5.))
			  +(58.*pow(velocityl3,6.))+(7.8*pow(velocityl3,7.))+ (0.888*pow(velocityl3,8.)) )/pow(1.+velocityl3,10.); //gFunction represents the correction for binding effet
  //-----------------------------------------------------------------------------------------------------------------------------

  G4double sigmaPSS_l3 = 1.+(((2.*zIncident)/(screenedzTarget*tetal3))*(gFunctionl3-hFunctionl3)); 


  //----------------------------------------------------------------------------------------------------------------------------
  //const G4double cNaturalUnit= 1/fine_structure_const;  // it's the speed of light according to Atomic-Unit-System
  //--------------------------------------------------------------------------------------------------------------

  //G4double yl3Formula=0.15*(screenedzTarget/cNaturalUnit)*(screenedzTarget/cNaturalUnit)/(velocityl3/sigmaPSS_l3);

 //-----------------------------------------------Relativity effect correction L3 --------------------------------------------------------------

 //G4double relativityCorrectionl3 = pow((1.+(1.1*yl3Formula*yl3Formula)),0.5)+yl3Formula;// the relativistic correction parameter
 //G4double reducedVelocityl3 = velocityl3*pow(relativityCorrectionl3,0.5);  // presents the reduced collision velocity parameter 

  G4double universalFunction_l3 ;

  G4double L3etaOverTheta2 = (energyIncident*electron_mass_c2)/(massIncident*rydbergMeV*screenedzTarget*screenedzTarget)/(sigmaPSS_l3*tetal3)/(sigmaPSS_l3*tetal3);
  
  universalFunction_l3 = 2*FunctionFL2((sigmaPSS_l3*tetal3), L3etaOverTheta2);
  
  
  G4double sigmaPSSR_l3 = (sigma0/(sigmaPSS_l3*tetal3))*universalFunction_l3;
  
  G4double pssDeltal3 = (4./(systemMass*sigmaPSS_l3*tetal3))*(sigmaPSS_l3/velocityl3)*(sigmaPSS_l3/velocityl3);

  G4double energyLossl3 = pow(1-pssDeltal3,0.5); 

  G4double coulombDeflectionl3 = (4.*pi*zIncident/systemMass)*pow(tetal3*sigmaPSS_l3,-2.)*pow(velocityl3/sigmaPSS_l3,-3.)*(zTarget/screenedzTarget); //incorporates Coulomb deflection parameter 

  G4double cParameterl3 = 2.*coulombDeflectionl3/(energyLossl3*(energyLossl3+1.));
  
  G4double coulombDeflectionFunction_l3 = 11.*ExpIntFunction(12,cParameterl3);
  
  G4double crossSection_L3 = coulombDeflectionFunction_l3 * sigmaPSSR_l3;              
  
  if (crossSection_L3 >= 0) {
    return crossSection_L3;
  }
  else {return 0;}
  
}



G4double G4ecpssrLiCrossSection::CalculateVelocity(G4int subShell, G4int zTarget, G4double massIncident,  G4double energyIncident) 
			                     
{  

  G4AtomicTransitionManager*  transitionManager =  G4AtomicTransitionManager::Instance();

  G4double liBindingEnergy = transitionManager->Shell(zTarget,subShell)->BindingEnergy();

  
  G4Proton* aProtone = G4Proton::Proton();
  G4Alpha* aAlpha = G4Alpha::Alpha();
  
  if (!(massIncident == aProtone->GetPDGMass() || massIncident == aAlpha->GetPDGMass()))
    {
      G4cout << "we can treat only Proton or Alpha incident particles " << G4endl;
      return 0;
    }
  
  
  const G4double zlshell= 4.15;      
  
  G4double screenedzTarget = zTarget- zlshell;                                  
  
  const G4double rydbergMeV= 13.6e-6;
  
  const G4double nl= 2.;        // nl is the quantum number of the L shell 
  
  G4double tetali = (liBindingEnergy*nl*nl)/(screenedzTarget*screenedzTarget*rydbergMeV);     

   G4double velocity =(2.*nl/(tetali*screenedzTarget))*pow(((energyIncident*electron_mass_c2)/(massIncident*rydbergMeV)),0.5);
  
  return velocity;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double G4ecpssrLiCrossSection::FunctionFL1(G4double k, G4double theta) 							  
{

  G4double sigma = 0.;
  G4double valueT1 = 0;
  G4double valueT2 = 0;
  G4double valueE21 = 0;
  G4double valueE22 = 0;
  G4double valueE12 = 0;
  G4double valueE11 = 0;
  G4double xs11 = 0;   
  G4double xs12 = 0; 
  G4double xs21 = 0; 
  G4double xs22 = 0; 

  // PROTECTION TO ALLOW INTERPOLATION AT MINIMUM AND MAXIMUM EtaK/Theta2 values 
  // (in particular for FK computation at 95 for high velocity formula)
  
  if (
  theta==9.5e-2 ||
  theta==9.5e-1 ||
  theta==9.5e+00 ||
  theta==9.5e+01
  ) theta=theta-1e-12;

  if (
  theta==1.e-2 ||
  theta==1.e-1 ||
  theta==1.e+00 ||
  theta==1.e+01
  ) theta=theta+1e-12;

  // END PROTECTION
  
  {
    std::vector<double>::iterator t2 = std::upper_bound(dummyVec.begin(),dummyVec.end(), k);
    std::vector<double>::iterator t1 = t2-1;
 
    std::vector<double>::iterator e12 = std::upper_bound(aVecMap[(*t1)].begin(),aVecMap[(*t1)].end(), theta);
    std::vector<double>::iterator e11 = e12-1; 
	  
    std::vector<double>::iterator e22 = std::upper_bound(aVecMap[(*t2)].begin(),aVecMap[(*t2)].end(), theta);
    std::vector<double>::iterator e21 = e22-1;
	  	
    valueT1  =*t1;
    valueT2  =*t2;
    valueE21 =*e21;
    valueE22 =*e22;
    valueE12 =*e12;
    valueE11 =*e11;

    xs11 = FL1Data[valueT1][valueE11];
    xs12 = FL1Data[valueT1][valueE12];
    xs21 = FL1Data[valueT2][valueE21];
    xs22 = FL1Data[valueT2][valueE22];

}
     
  G4double xsProduct = xs11 * xs12 * xs21 * xs22;
  
  if (xs11==0 || xs12==0 ||xs21==0 ||xs22==0) return (0.);
     
  if (xsProduct != 0.)
  {
    sigma = QuadInterpolator(  valueE11, valueE12, 
    			       valueE21, valueE22, 
			       xs11, xs12, 
			       xs21, xs22, 
			       valueT1, valueT2, 
			       k, theta );
  }

  return sigma;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4double G4ecpssrLiCrossSection::FunctionFL2(G4double k, G4double theta) 							  
{

  G4double sigma = 0.;
  G4double valueT1 = 0;
  G4double valueT2 = 0;
  G4double valueE21 = 0;
  G4double valueE22 = 0;
  G4double valueE12 = 0;
  G4double valueE11 = 0;
  G4double xs11 = 0;   
  G4double xs12 = 0; 
  G4double xs21 = 0; 
  G4double xs22 = 0; 

  // PROTECTION TO ALLOW INTERPOLATION AT MINIMUM AND MAXIMUM EtaK/Theta2 values 
  // (in particular for FK computation at 95 for high velocity formula)
  
  if (
  theta==9.5e-2 ||
  theta==9.5e-1 ||
  theta==9.5e+00 ||
  theta==9.5e+01
  ) theta=theta-1e-12;

  if (
  theta==1.e-2 ||
  theta==1.e-1 ||
  theta==1.e+00 ||
  theta==1.e+01
  ) theta=theta+1e-12;

  // END PROTECTION
  
  {
    std::vector<double>::iterator t2 = std::upper_bound(dummyVec.begin(),dummyVec.end(), k);
    std::vector<double>::iterator t1 = t2-1;
 
    std::vector<double>::iterator e12 = std::upper_bound(aVecMap[(*t1)].begin(),aVecMap[(*t1)].end(), theta);
    std::vector<double>::iterator e11 = e12-1; 
	  
    std::vector<double>::iterator e22 = std::upper_bound(aVecMap[(*t2)].begin(),aVecMap[(*t2)].end(), theta);
    std::vector<double>::iterator e21 = e22-1;
	  	
    valueT1  =*t1;
    valueT2  =*t2;
    valueE21 =*e21;
    valueE22 =*e22;
    valueE12 =*e12;
    valueE11 =*e11;

    xs11 = FL2Data[valueT1][valueE11];
    xs12 = FL2Data[valueT1][valueE12];
    xs21 = FL2Data[valueT2][valueE21];
    xs22 = FL2Data[valueT2][valueE22];


}
     
  G4double xsProduct = xs11 * xs12 * xs21 * xs22;
  
  if (xs11==0 || xs12==0 ||xs21==0 ||xs22==0) return (0.);
     
  if (xsProduct != 0.)
  {
    sigma = QuadInterpolator(  valueE11, valueE12, 
    			       valueE21, valueE22, 
			       xs11, xs12, 
			       xs21, xs22, 
			       valueT1, valueT2, 
			       k, theta );
  }

  return sigma;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ecpssrLiCrossSection::LinLogInterpolate(G4double e1, 
						        G4double e2, 
						        G4double e, 
						        G4double xs1, 
						        G4double xs2)
{
  G4double d1 = std::log(xs1);
  G4double d2 = std::log(xs2);
  G4double value = std::exp(d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ecpssrLiCrossSection::LogLogInterpolate(G4double e1, 
						        G4double e2, 
						        G4double e, 
						        G4double xs1, 
						        G4double xs2)
{
  G4double a = (std::log10(xs2)-std::log10(xs1)) / (std::log10(e2)-std::log10(e1));
  G4double b = std::log10(xs2) - a*std::log10(e2);
  G4double sigma = a*std::log10(e) + b;
  G4double value = (std::pow(10.,sigma));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ecpssrLiCrossSection::QuadInterpolator(G4double e11, G4double e12, 
						       G4double e21, G4double e22, 
						       G4double xs11, G4double xs12, 
						       G4double xs21, G4double xs22, 
						       G4double t1, G4double t2, 
						       G4double t, G4double e)
{
// Log-Log
/*
  G4double interpolatedvalue1 = LogLogInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LogLogInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LogLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
*/

// Lin-Log
  G4double interpolatedvalue1 = LinLogInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LinLogInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LinLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
  return value;
}

