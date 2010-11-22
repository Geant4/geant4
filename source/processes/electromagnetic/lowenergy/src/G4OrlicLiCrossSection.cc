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
//$Id: G4OrlicLiCrossSection.cc,v 1.6 2010-11-22 18:32:00 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Haifa Ben Abdelouahed
//         
//
// History:
// -----------
//  23 Apr 2008   H. Ben Abdelouahed   1st implementation
//  28 Apr 2008   MGP        Major revision according to a design iteration
//  21 Apr 2009	  ALF Some correction for compatibility to G4VShellCrossSection
//		  and changed name to G4OrlicLiCrossSection 
//
// -------------------------------------------------------------------
// Class description:
// Low Energy Electromagnetic Physics, Cross section, proton ionisation, L shell
// Further documentation available from http://www.ge.infn.it/geant4/lowE
// -------------------------------------------------------------------


#include "globals.hh"
#include "G4OrlicLiCrossSection.hh"
#include "G4Proton.hh"


G4OrlicLiCrossSection::G4OrlicLiCrossSection()
{ 

  transitionManager =  G4AtomicTransitionManager::Instance();

}

G4OrlicLiCrossSection::~G4OrlicLiCrossSection()
{ 

}

//this L-CrossSection calculation method is done according to 
//I.ORLIC, C.H.SOW and S.M.TANG,International Journal of PIXE.Vol.4(1997) 217-230	
 

//*****************************************************************************************************************************************

G4double G4OrlicLiCrossSection::CalculateL1CrossSection(G4int zTarget, G4double energyIncident)
  			       
{

  if ( (energyIncident < 0.1*MeV) || energyIncident > 10*MeV )

    {return 0;}



  G4double  massIncident; 

  G4Proton* aProtone = G4Proton::Proton();
    
   massIncident = aProtone->GetPDGMass(); 
 
  G4double l1BindingEnergy = (transitionManager->Shell(zTarget,1)->BindingEnergy())/keV;
 
  G4double lamda =  massIncident/electron_mass_c2;

  G4double normalizedEnergy =  (energyIncident/keV)/(lamda*l1BindingEnergy);

  G4double x = std::log(normalizedEnergy);

  G4double a0 = 0.;
  G4double a1 = 0.;
  G4double a2 = 0.;
  G4double a3 = 0.;
  G4double a4 = 0.;
  G4double a5 = 0.;
  G4double a6 = 0.;  
  G4double a7 = 0.;  
  G4double a8 = 0.; 
  G4double a9 = 0.;

                                       
   if ( zTarget>=14 && zTarget<=40)  
    {

      return 0;
      /*
      // parameters used for calculating total L cross section
      a0=12.5081;
      a1=0.2177;
      a2=-0.3758;
      a3=0.0096;
      a4=0.0073;
      a5=0.0022;
      a6=0.;
      a7=0.;
      a8=0.;
      a9=0.; */
    }
  else 
    { 
       if ( zTarget>=41 &&  zTarget<=50 )
       { 	  
      a0=11.274881;
      a1=-0.187401;
      a2=-0.943341;
      a3=-1.47817;
      a4=-1.282343;
      a5=-0.386544; 
      a6=-0.037932;
      a7=0.;
      a8=0.;
      a9=0.;
	}
     
      else  
	{
	  if ( zTarget>=51 &&  zTarget<=60 ) 
	    { 
	      a0=11.242637;
	      a1=-0.162515;
	      a2=1.035774;
	      a3=3.970908;
	      a4=3.968233;
	      a5=1.655714;
	      a6=0.058885;
	      a7=-0.155743;
	      a8=-0.042228;
	      a9=-0.003371; 
	    } 
	
	  else  
	    {
	      if ( zTarget>=61 &&  zTarget<=70 ) 
		{ 
		  a0=6.476722;
		  a1=-25.804787;
		  a2=-54.061629;
		  a3=-56.684589;
		  a4=-33.223367;
		  a5=-11.034979; 
		  a6=-2.042851;
		  a7=-0.194075;
		  a8=-0.007252;
		  a9=0.;
		} 
	      else  
		{
		  if ( zTarget>=71 &&  zTarget<=80 ) 
		    { 
		      a0=12.776794;
		      a1=6.562907;
		      a2=10.158703;
		      a3=7.432592;
		      a4=2.332036;
		      a5=0.317946; 
		      a6=0.014479;
		      a7=0.;
		      a8=0.;
		      a9=0.;
		    } 
		  else  
		    {
		      if ( zTarget>=81 &&  zTarget<=92 ) 
			{ 
			  a0=28.243087;
			  a1=50.199585;
			  a2=58.281684;
			  a3=34.130538;
			  a4=10.268531;
			  a5=1.525302; 
			  a6=0.08835;
			  a7=0.;
			  a8=0.;
			  a9=0.;
			}
		      else
			{ 
			  G4cout << "ERROR: L1 Cross-Section exist only for ZTarget between 14 and 92!!! " << G4endl;
			
			}
		    }
		}
	    }
	}
      }
      

G4double analyticalFunction = a0 + (a1*x)+(a2*x*x)+(a3*std::pow(x,3))+(a4*std::pow(x,4))+(a5*std::pow(x,5))+(a6*std::pow(x,6))+
	(a7*std::pow(x,7))+(a8*std::pow(x,8))+(a9*std::pow(x,9)); 



  G4double L1crossSection =  std::exp(analyticalFunction)/(l1BindingEnergy*l1BindingEnergy);  


  if (L1crossSection >= 0) {
    return L1crossSection * barn;
  }
  else {return 0;}

}

//*****************************************************************************************************************************************


G4double G4OrlicLiCrossSection::CalculateL2CrossSection(G4int zTarget, G4double energyIncident)
  			       
{


  if ( (energyIncident < 0.1*MeV) || energyIncident > 10*MeV )

    {return 0;}


  G4double  massIncident; 

  G4Proton* aProtone = G4Proton::Proton();
    
   massIncident = aProtone->GetPDGMass(); 
 
 G4double L2crossSection;

 if (zTarget<41)
   { 
     L2crossSection =0.;
   }
 else
   {

 G4double l2BindingEnergy = (transitionManager->Shell(zTarget,2)->BindingEnergy())/keV;
 
  G4double lamda =  massIncident/electron_mass_c2;

  G4double normalizedEnergy =  (energyIncident/keV)/(lamda*l2BindingEnergy);

  G4double x = std::log(normalizedEnergy);

  G4double a0 = 0.;
  G4double a1 = 0.;
  G4double a2 = 0.;
  G4double a3 = 0.;
  G4double a4 = 0.;
  G4double a5 = 0.;
 
      if ( zTarget>=41 &&  zTarget<=50 ) 
	{ 
      a0=11.194798;
      a1=0.178807;
      a2=-0.449865;
      a3=-0.063528;
      a4=-0.015364;
      a5=0.; 
	}
     
      else  
	{
	  if ( zTarget>=51 &&  zTarget<=60 ) 
	    { 
	      a0=11.241409;
	      a1=0.149635;
	      a2=-0.633269;
	      a3=-0.17834;
	      a4=-0.034743;
	      a5=0.006474;

	    } 
	
	  else  
	    {
	      if ( zTarget>=61 &&  zTarget<=70 ) 
		{ 
		  a0=11.247424;
		  a1=0.203051;
		  a2=-0.219083;
		  a3=0.164514;
		  a4=0.058692;
		  a5=0.007866; 
		} 
	      else  
		{
		  if ( zTarget>=71 &&  zTarget<=80 ) 
		    { 
		      a0=11.229924;
		      a1=-0.087241;
		      a2=-0.753908;
		      a3=-0.181546;
		      a4=-0.030406;
		      a5=0.; 
		    } 
		  else  
		    {
		      if ( zTarget>=81 &&  zTarget<=92 ) 
			{ 
			  a0=11.586671;
			  a1=0.730838;
			  a2=-0.056713;
			  a3=0.053262;
			  a4=-0.003672;
			  a5=0.; 
			}
		      else
			{ 
			  G4cout << "ERROR: L2 Cross-Section exist only for ZTarget between 14 and 92!!! " << G4endl;
		       
			}
		    }
		}
	    }
	}
      

 G4double analyticalFunction = a0 + (a1*x)+(a2*x*x)+(a3*std::pow(x,3))+(a4*std::pow(x,4))+(a5*std::pow(x,5)); 


   L2crossSection =  std::exp(analyticalFunction)/(l2BindingEnergy*l2BindingEnergy); 
                                                                            
   } 

  if (L2crossSection >= 0) {
    return L2crossSection * barn;
  }
  else {return 0;}
   
}

//*****************************************************************************************************************************************


G4double G4OrlicLiCrossSection::CalculateL3CrossSection(G4int zTarget, G4double energyIncident)
  			       
{

  if ( (energyIncident < 0.1*MeV) || energyIncident > 10*MeV )

    {return 0;}



  G4double  massIncident; 

  G4Proton* aProtone = G4Proton::Proton();
    
  massIncident = aProtone->GetPDGMass(); 

 
 G4double L3crossSection;

 if (zTarget<41)
   { 
     L3crossSection =0.;
   }
 else
   {

 G4double l3BindingEnergy = (transitionManager->Shell(zTarget,3)->BindingEnergy())/keV;
 
 
  G4double lamda =  massIncident/electron_mass_c2;

  G4double normalizedEnergy =  (energyIncident/keV)/(lamda*l3BindingEnergy);

  G4double x = std::log(normalizedEnergy);


  G4double a0 = 0.;
  G4double a1 = 0.;
  G4double a2 = 0.;
  G4double a3 = 0.;
  G4double a4 = 0.;
 
      if ( zTarget>=41 &&  zTarget<=50 ) 
	{ 
      a0=11.91837;
      a1=0.03064;
      a2=-0.657644;
      a3=-0.14532;
      a4=-0.026059;      
	}
     
      else  
	{
	  if ( zTarget>=51 &&  zTarget<=60 ) 
	    { 
	      a0=11.909485;
	      a1=0.15918;
	      a2=-0.588004;
	      a3=-0.159466;
	      a4=-0.033184;	      
	    } 
	
	  else  
	    {
	      if ( zTarget>=61 &&  zTarget<=70 ) 
		{ 
		  a0=11.878472;
		  a1=-0.137007;
		  a2=-0.959475;
		  a3=-0.316505;
		  a4=-0.054154;	 
		} 
	      else  
		{
		  if ( zTarget>=71 &&  zTarget<=80 ) 
		    { 
		      a0=11.802538;
		      a1=-0.371796;
		      a2=-1.052238;
		      a3=-0.28766;
		      a4=-0.042608;	     
		    } 
		  else  
		    {
		      if ( zTarget>=81 &&  zTarget<=92 ) 
			{ 
			  a0=11.423712;
			  a1=-1.428823;
			  a2=-1.946979;
			  a3=-0.585198;
			  a4=-0.076467;			 
			}
		      else
			{ 
			  G4cout << "ERROR: L3 Cross-Section exist only for ZTarget between 14 and 92!!! " << G4endl;
		       
			}
		    }
		}
	    }
	}
      

 G4double analyticalFunction = a0 + (a1*x)+(a2*x*x)+(a3*std::pow(x,3))+(a4*std::pow(x,4)); 


   L3crossSection =  std::exp(analyticalFunction)/(l3BindingEnergy*l3BindingEnergy); 
                                                                            
   } 
  if (L3crossSection >= 0) {
    return L3crossSection * barn;
  }
  else {return 0;}


}
