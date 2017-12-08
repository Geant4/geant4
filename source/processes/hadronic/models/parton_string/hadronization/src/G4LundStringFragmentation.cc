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
//
// $Id: G4LundStringFragmentation.cc 106967 2017-10-31 08:41:49Z gcosmo $
// GEANT4 tag $Name:  $ 1.8
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Maxim Komogorov, 10-Jul-1998
// -----------------------------------------------------------------------------
#include "G4LundStringFragmentation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4FragmentingString.hh"
#include "G4DiQuarks.hh"
#include "G4Quarks.hh"

#include "G4Exp.hh"
#include "G4Pow.hh"

//#define debug_LUNDfragmentation

// Class G4LundStringFragmentation 
//*************************************************************************************

G4LundStringFragmentation::G4LundStringFragmentation()
{
  // ------ For estimation of a minimal string mass ---------------
  Mass_of_light_quark    =140.*MeV;
  Mass_of_heavy_quark    =500.*MeV;
  Mass_of_string_junction=720.*MeV;
  // ------ An estimated minimal string mass ----------------------
  MinimalStringMass  = 0.;              
  MinimalStringMass2 = 0.;              
  // ------ Minimal invariant mass used at a string fragmentation -
  WminLUND = 0.45*GeV; //0.23*GeV;                   // Uzhi 0.7 -> 0.23 3.8.10 //0.8 1.5
  // ------ Smooth parameter used at a string fragmentation for ---
  // ------ smearing sharp mass cut-off ---------------------------
  SmoothParam  = 0.2;                   

  SetStringTensionParameter(1.);                         
  SetDiquarkBreakProbability(0.05); 
  SetStrangenessSuppression(0.46); //(0.447);
  SetDiquarkSuppression(0.05);

  // For treating of small string decays
  for(G4int i=0; i<3; i++)
  { for(G4int j=0; j<3; j++)
    { for(G4int k=0; k<6; k++)
      {
        Meson[i][j][k]=0; MesonWeight[i][j][k]=0.;
      }
    }
  }
  //--------------------------
        Meson[0][0][0]=111;                       // dbar-d Pi0
  MesonWeight[0][0][0]=(1.-pspin_meson)*(1.-scalarMesonMix[0]);

        Meson[0][0][1]=221;                       // dbar-d Eta
  MesonWeight[0][0][1]=(1.-pspin_meson)*(scalarMesonMix[0]-scalarMesonMix[1]);

        Meson[0][0][2]=331;                       // dbar-d EtaPrime
  MesonWeight[0][0][2]=(1.-pspin_meson)*(scalarMesonMix[1]);

        Meson[0][0][3]=113;                       // dbar-d Rho0
  MesonWeight[0][0][3]=pspin_meson*(1.-vectorMesonMix[0]);

        Meson[0][0][4]=223;                       // dbar-d Omega
  MesonWeight[0][0][4]=pspin_meson*(vectorMesonMix[0]);
  //--------------------------

        Meson[0][1][0]=211;                       // dbar-u Pi+
  MesonWeight[0][1][0]=(1.-pspin_meson);

        Meson[0][1][1]=213;                       // dbar-u Rho+
  MesonWeight[0][1][1]=pspin_meson;
  //--------------------------

        Meson[0][2][0]=311;                      // dbar-s K0bar
  MesonWeight[0][2][0]=(1.-pspin_meson);

        Meson[0][2][1]=313;                       // dbar-s K*0bar
  MesonWeight[0][2][1]=pspin_meson;
  //--------------------------

        Meson[1][0][0]=211;                       // ubar-d Pi-
  MesonWeight[1][0][0]=(1.-pspin_meson);

        Meson[1][0][1]=213;                       // ubar-d Rho-
  MesonWeight[1][0][1]=pspin_meson;
  //--------------------------

        Meson[1][1][0]=111;                       // ubar-u Pi0
  MesonWeight[1][1][0]=(1.-pspin_meson)*(1.-scalarMesonMix[0]);

        Meson[1][1][1]=221;                       // ubar-u Eta
  MesonWeight[1][1][1]=(1.-pspin_meson)*(scalarMesonMix[0]-scalarMesonMix[1]);

        Meson[1][1][2]=331;                       // ubar-u EtaPrime
  MesonWeight[1][1][2]=(1.-pspin_meson)*(scalarMesonMix[1]);

        Meson[1][1][3]=113;                       // ubar-u Rho0
  MesonWeight[1][1][3]=pspin_meson*(1.-vectorMesonMix[0]);

        Meson[1][1][4]=223;                       // ubar-u Omega
  //MesonWeight[1][1][4]=pspin_meson*(scalarMesonMix[0]);
  MesonWeight[1][1][4]=pspin_meson*(vectorMesonMix[0]);  // Uzhi 2015 scalar -> vector
  //--------------------------

        Meson[1][2][0]=321;                      // ubar-s K-
  MesonWeight[1][2][0]=(1.-pspin_meson);

        Meson[1][2][1]=323;                      // ubar-s K*-bar -
  MesonWeight[1][2][1]=pspin_meson;
  //--------------------------

        Meson[2][0][0]=311;                       // sbar-d K0
  MesonWeight[2][0][0]=(1.-pspin_meson);

        Meson[2][0][1]=313;                       // sbar-d K*0
  MesonWeight[2][0][1]=pspin_meson;
  //--------------------------

        Meson[2][1][0]=321;                        // sbar-u K+
  MesonWeight[2][1][0]=(1.-pspin_meson);

        Meson[2][1][1]=323;                       // sbar-u K*+
  MesonWeight[2][1][1]=pspin_meson;
  //--------------------------

        Meson[2][2][0]=221;                       // sbar-s Eta
  MesonWeight[2][2][0]=(1.-pspin_meson)*(1.-scalarMesonMix[5]);

        Meson[2][2][1]=331;                       // sbar-s EtaPrime
  MesonWeight[2][2][1]=(1.-pspin_meson)*(1.-scalarMesonMix[5]);

        Meson[2][2][3]=333;                       // sbar-s EtaPrime
  MesonWeight[2][2][3]=pspin_meson*(vectorMesonMix[5]);
  //--------------------------

  for(G4int i=0; i<3; i++)
  { for(G4int j=0; j<3; j++)
    { for(G4int k=0; k<3; k++)
      { for(G4int l=0; l<4; l++)
        { Baryon[i][j][k][l]=0; BaryonWeight[i][j][k][l]=0.;}
      }
    }
  }

  G4double pspin_barion_in=pspin_barion;
  //pspin_barion=0.75;
  //---------------------------------------

        Baryon[0][0][0][0]=1114;         // Delta-
  BaryonWeight[0][0][0][0]=1.;
  //---------------------------------------

        Baryon[0][0][1][0]=2112;         // neutron
  BaryonWeight[0][0][1][0]=1.-pspin_barion;

        Baryon[0][0][1][1]=2114;         // Delta0
  BaryonWeight[0][0][1][1]=pspin_barion;
  //---------------------------------------

        Baryon[0][0][2][0]=3112;         // Sigma-
  BaryonWeight[0][0][2][0]=1.-pspin_barion;

        Baryon[0][0][2][1]=3114;         // Sigma*-
  BaryonWeight[0][0][2][1]=pspin_barion;
  //---------------------------------------

        Baryon[0][1][0][0]=2112;         // neutron
  BaryonWeight[0][1][0][0]=1.-pspin_barion;

        Baryon[0][1][0][1]=2114;         // Delta0
  BaryonWeight[0][1][0][1]=pspin_barion;
  //---------------------------------------

        Baryon[0][1][1][0]=2212;         // proton
  BaryonWeight[0][1][1][0]=1.-pspin_barion;

        Baryon[0][1][1][1]=2214;         // Delta+
  BaryonWeight[0][1][1][1]=pspin_barion;
  //---------------------------------------

        Baryon[0][1][2][0]=3122;         // Lambda
  BaryonWeight[0][1][2][0]=(1.-pspin_barion)*0.5;

        Baryon[0][1][2][1]=3212;         // Sigma0
  BaryonWeight[0][1][2][1]=(1.-pspin_barion)*0.5;

        Baryon[0][1][2][2]=3214;         // Sigma*0
  BaryonWeight[0][1][2][2]=pspin_barion;
  //---------------------------------------

        Baryon[0][2][0][0]=3112;         // Sigma-
  BaryonWeight[0][2][0][0]=1.-pspin_barion;

        Baryon[0][2][0][1]=3114;         // Sigma*-
  BaryonWeight[0][2][0][1]=pspin_barion;
  //---------------------------------------

        Baryon[0][2][1][0]=3122;         // Lambda
  BaryonWeight[0][2][1][0]=(1.-pspin_barion)*0.5;

        Baryon[0][2][1][1]=3212;         // Sigma0
  BaryonWeight[0][2][1][1]=(1.-pspin_barion)*0.5;

        Baryon[0][2][1][2]=3214;         // Sigma*0
  BaryonWeight[0][2][1][2]=pspin_barion;
  //---------------------------------------

        Baryon[0][2][2][0]=3312;         // Theta-
  BaryonWeight[0][2][2][0]=1.-pspin_barion;

        Baryon[0][2][2][1]=3314;         // Theta*-
  BaryonWeight[0][2][2][1]=pspin_barion;
  //---------------------------------------

        Baryon[1][0][0][0]=2112;         // neutron
  BaryonWeight[1][0][0][0]=1.-pspin_barion;

        Baryon[1][0][0][1]=2114;         // Delta0
  BaryonWeight[1][0][0][1]=pspin_barion;
  //---------------------------------------

        Baryon[1][0][1][0]=2212;         // proton
  BaryonWeight[1][0][1][0]=1.-pspin_barion;          

        Baryon[1][0][1][1]=2214;         // Delta+
  BaryonWeight[1][0][1][1]=pspin_barion;
  //---------------------------------------

        Baryon[1][0][2][0]=3122;         // Lambda
  BaryonWeight[1][0][2][0]=(1.-pspin_barion)*0.5;

        Baryon[1][0][2][1]=3212;         // Sigma0
  BaryonWeight[1][0][2][1]=(1.-pspin_barion)*0.5;

        Baryon[1][0][2][2]=3214;         // Sigma*0
  BaryonWeight[1][0][2][2]=pspin_barion;
  //---------------------------------------

        Baryon[1][1][0][0]=2212;         // proton
  BaryonWeight[1][1][0][0]=1.-pspin_barion;

        Baryon[1][1][0][1]=2214;         // Delta+
  BaryonWeight[1][1][0][1]=pspin_barion;
  //---------------------------------------

        Baryon[1][1][1][0]=2224;         // Delta++
  BaryonWeight[1][1][1][0]=1.;
  //---------------------------------------

        Baryon[1][1][2][0]=3222;         // Sigma+
  BaryonWeight[1][1][2][0]=1.-pspin_barion;

        Baryon[1][1][2][1]=3224;         // Sigma*+
  BaryonWeight[1][1][2][1]=pspin_barion;
  //---------------------------------------

        Baryon[1][2][0][0]=3122;         // Lambda
  BaryonWeight[1][2][0][0]=(1.-pspin_barion)*0.5;

        Baryon[1][2][0][1]=3212;         // Sigma0
  BaryonWeight[1][2][0][1]=(1.-pspin_barion)*0.5;

        Baryon[1][2][0][2]=3214;         // Sigma*0
  BaryonWeight[1][2][0][2]=pspin_barion;
  //---------------------------------------

        Baryon[1][2][1][0]=3222;         // Sigma+
  BaryonWeight[1][2][1][0]=1.-pspin_barion;

        Baryon[1][2][1][1]=3224;         // Sigma*+
  BaryonWeight[1][2][1][1]=pspin_barion;
  //---------------------------------------

        Baryon[1][2][2][0]=3322;         // Theta0
  BaryonWeight[1][2][2][0]=1.-pspin_barion;

        Baryon[1][2][2][1]=3324;         // Theta*0
  BaryonWeight[1][2][2][1]=pspin_barion;
  //---------------------------------------

        Baryon[2][0][0][0]=3112;         // Sigma-
  BaryonWeight[2][0][0][0]=1.-pspin_barion;

        Baryon[2][0][0][1]=3114;         // Sigma*-
  BaryonWeight[2][0][0][1]=pspin_barion;
  //---------------------------------------

        Baryon[2][0][1][0]=3122;         // Lambda
  BaryonWeight[2][0][1][0]=(1.-pspin_barion)*0.5;          

        Baryon[2][0][1][1]=3212;         // Sigma0
  BaryonWeight[2][0][1][1]=(1.-pspin_barion)*0.5; 

        Baryon[2][0][1][2]=3214;         // Sigma*0
  BaryonWeight[2][0][1][2]=pspin_barion;
  //---------------------------------------

        Baryon[2][0][2][0]=3312;         // Sigma-
  BaryonWeight[2][0][2][0]=1.-pspin_barion;

        Baryon[2][0][2][1]=3314;         // Sigma*-
  BaryonWeight[2][0][2][1]=pspin_barion;
  //---------------------------------------

        Baryon[2][1][0][0]=3122;         // Lambda
  BaryonWeight[2][1][0][0]=(1.-pspin_barion)*0.5;

        Baryon[2][1][0][1]=3212;         // Sigma0
  BaryonWeight[2][1][0][1]=(1.-pspin_barion)*0.5;

        Baryon[2][1][0][2]=3214;         // Sigma*0
  BaryonWeight[2][1][0][2]=pspin_barion;
  //---------------------------------------

        Baryon[2][1][1][0]=3222;         // Sigma+
  BaryonWeight[2][1][1][0]=1.-pspin_barion;

        Baryon[2][1][1][1]=3224;         // Sigma*+
  BaryonWeight[2][1][1][1]=pspin_barion;
  //---------------------------------------

        Baryon[2][1][2][0]=3322;         // Theta0
  BaryonWeight[2][1][2][0]=1.-pspin_barion;

        Baryon[2][1][2][1]=3324;         // Theta*0
  BaryonWeight[2][1][2][1]=pspin_barion;
  //---------------------------------------

        Baryon[2][2][0][0]=3312;         // Theta-
  BaryonWeight[2][2][0][0]=1.-pspin_barion;

        Baryon[2][2][0][1]=3314;         // Theta*-
  BaryonWeight[2][2][0][1]=pspin_barion;
  //---------------------------------------

        Baryon[2][2][1][0]=3322;         // Theta0
  BaryonWeight[2][2][1][0]=1.-pspin_barion;

        Baryon[2][2][1][1]=3324;         // Theta*0
  BaryonWeight[2][2][1][1]=pspin_barion;
  //---------------------------------------

        Baryon[2][2][2][0]=3334;         // Omega
  BaryonWeight[2][2][2][0]=1.;
  //---------------------------------------

  pspin_barion=pspin_barion_in;
  /*
  for(G4int i=0; i<3; i++)
  { for(G4int j=0; j<3; j++)
    { for(G4int k=0; k<3; k++)
      { for(G4int l=0; l<4; l++)
	{ G4cout<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<Baryon[i][j][k][l]<<G4endl;}
      }
    }
  }
  G4int Uzhi;
  G4cin>>Uzhi;
  */

  SetStrangenessSuppression(0.375);
  Prob_QQbar[0]=StrangeSuppress;         // Probability of ddbar production
  Prob_QQbar[1]=StrangeSuppress;         // Probability of uubar production
  Prob_QQbar[2]=1.0-2.*StrangeSuppress;  // Probability of ssbar production 
  SetStrangenessSuppression(0.46); //(0.447);

  //A.R. 25-Jul-2012 : Coverity fix.
  for ( G4int i=0 ; i<35 ; i++ ) { 
    FS_LeftHadron[i] = 0;
    FS_RightHadron[i] = 0;
    FS_Weight[i] = 0.0; 
  }
  NumberOf_FS = 0;

}

// --------------------------------------------------------------
G4LundStringFragmentation::~G4LundStringFragmentation()
{}


//--------------------------------------------------------------------------------------
void G4LundStringFragmentation::SetMinimalStringMass(const G4FragmentingString  * const string)  
{
  G4double EstimatedMass=0.;
  G4int Number_of_quarks=0;
  G4int Number_of_squarks=0;
        
  G4double StringM=string->Get4Momentum().mag();

  G4int Qleft =std::abs(string->GetLeftParton()->GetPDGEncoding());

  #ifdef debug_LUNDfragmentation
  //G4cout<<"MinStringMass// Input String mass "<<string->Get4Momentum().mag()<<" Qleft "<<Qleft<<G4endl;
  #endif

  if( Qleft > 1000)
  {
    Number_of_quarks+=2;
    G4int q1=Qleft/1000;
    if( q1 < 3) {EstimatedMass +=Mass_of_light_quark;}
    if( q1 > 2) {EstimatedMass +=Mass_of_heavy_quark; Number_of_squarks++;}

    G4int q2=(Qleft/100)%10;
    if( q2 < 3) {EstimatedMass +=Mass_of_light_quark;}
    if( q2 > 2) {EstimatedMass +=Mass_of_heavy_quark; Number_of_squarks++;}
    //		EstimatedMass +=Mass_of_string_junction;
    //if((q1 > 2)||(q2 > 2)) EstimatedMass -= 120*MeV;
  } 
  else 
  {
    Number_of_quarks++;
    if( Qleft < 3) {EstimatedMass +=Mass_of_light_quark;}
    if( Qleft > 2) {EstimatedMass +=Mass_of_heavy_quark; Number_of_squarks++;}
  }

  #ifdef debug_LUNDfragmentation
  //G4cout<<"Min mass with Qleft "<<Qleft<<" "<<EstimatedMass<<G4endl;
  #endif
  G4int Qright=std::abs(string->GetRightParton()->GetPDGEncoding());
  if( Qright > 1000)
  {
    Number_of_quarks+=2;
    G4int q1=Qright/1000;
    if( q1 < 3) {EstimatedMass +=Mass_of_light_quark;}
    if( q1 > 2) {EstimatedMass +=Mass_of_heavy_quark; Number_of_squarks++;}

    G4int q2=(Qright/100)%10;
    if( q2 < 3) {EstimatedMass +=Mass_of_light_quark;}
    if( q2 > 2) {EstimatedMass +=Mass_of_heavy_quark; Number_of_squarks++;}
		//EstimatedMass +=Mass_of_string_junction;
    //if((q1 > 2)||(q2 > 2)) EstimatedMass -= 120*MeV;
  }
  else 
  {
    Number_of_quarks++;
    if( Qright < 3) {EstimatedMass +=Mass_of_light_quark;}
    if( Qright > 2) {EstimatedMass +=Mass_of_heavy_quark; Number_of_squarks++;}
  }

  #ifdef debug_LUNDfragmentation
  //G4cout<<"Min mass with Qleft and Qright "<<Qright<<" "<<EstimatedMass<<G4endl;
  //G4cout<<"Number_of_quarks "<<Number_of_quarks<<" Number_of_squarks "<<Number_of_squarks<<G4endl;
  #endif

  if(Number_of_quarks==2){EstimatedMass += 70.*MeV;} //100.*MeV;}
  //if(Number_of_quarks==3){EstimatedMass += 20.*MeV;}
  if(Number_of_quarks==3)
  { 
    if(Number_of_squarks==0) {EstimatedMass += 740.*MeV;} // 700
    if(Number_of_squarks==1) {EstimatedMass += 740.*MeV;} // 740
    if(Number_of_squarks==2) {EstimatedMass += 400.*MeV;}        
    if(Number_of_squarks==3) {EstimatedMass += 382.*MeV;}
  }
  if(Number_of_quarks==4)
  {
    if((StringM > 1880.) && ( EstimatedMass < 2100))     {EstimatedMass = 2020.;}//1880.;}
    //if((StringM > 1880.) && ( EstimatedMass < 2100))     {EstimatedMass = 2051.;}
    else if((StringM > 2232.) && ( EstimatedMass < 2730)){EstimatedMass = 2570.;}
    else if((StringM > 5130.) && ( EstimatedMass < 3450)){EstimatedMass = 5130.;}
    else {
      // EstimatedMass -=2.*Mass_of_string_junction;
      if(EstimatedMass <= 1600.*MeV){EstimatedMass-=200.*MeV;}
      else                          {EstimatedMass+=100.*MeV;}
    }
  }

  #ifdef debug_LUNDfragmentation
  //G4cout<<"EstimatedMass -------------------- "<<EstimatedMass <<G4endl;
  #endif
  MinimalStringMass=EstimatedMass;
  SetMinimalStringMass2(EstimatedMass);
}

//--------------------------------------------------------------------------------------
void G4LundStringFragmentation::SetMinimalStringMass2(const G4double aValue)
{
  MinimalStringMass2=aValue * aValue;
}

//--------------------------------------------------------------------------------------
G4KineticTrackVector* G4LundStringFragmentation::FragmentString(const G4ExcitedString& theString)
{
  // Can no longer modify Parameters for Fragmentation.
  PastInitPhase=true;

  SetMassCut(160.*MeV); // For LightFragmentationTest it is required that no one pi-meson can be produced

  G4FragmentingString  aString(theString);
  SetMinimalStringMass(&aString);

  #ifdef debug_LUNDfragmentation
  G4cout<<G4endl<<"LUND StringFragm: String Mass " <<theString.Get4Momentum().mag()<<" Pz "
                                                   <<theString.Get4Momentum().pz()
                <<"------------------------------------"<<G4endl;
  G4cout<<"String ends Direct "<<theString.GetLeftParton()->GetPDGcode()<<" "
                               <<theString.GetRightParton()->GetPDGcode()<<" "
                               <<theString.GetDirection()<< G4endl;
  G4cout<<"Left  mom "<<theString.GetLeftParton()->Get4Momentum()<<G4endl;
  G4cout<<"Right mom "<<theString.GetRightParton()->Get4Momentum()<<G4endl;
  G4cout<<"Check for Fragmentation "<<G4endl;
  #endif

  G4KineticTrackVector * LeftVector(0);

  //if(!IsFragmentable(&aString)) // produce 1 hadron   True ===============
  if(!aString.FourQuarkString() && !IsFragmentable(&aString))
  {
    #ifdef debug_LUNDfragmentation
    G4cout<<"Non fragmentable - the string is converted to one hadron "<<G4endl;
    #endif

    SetMassCut(10000.*MeV);
    LeftVector=LightFragmentationTest(&theString);
    SetMassCut(160.*MeV);

    LeftVector->operator[](0)->SetFormationTime(theString.GetTimeOfCreation());
    LeftVector->operator[](0)->SetPosition(theString.GetPosition());

    if(LeftVector->size() > 1)
    {
      // 2 hadrons created from qq-qqbar are stored
      LeftVector->operator[](1)->SetFormationTime(theString.GetTimeOfCreation());
      LeftVector->operator[](1)->SetPosition(theString.GetPosition());
    }
    return LeftVector;
  }

  #ifdef debug_LUNDfragmentation
  G4cout<<"The string will be fragmented. "<<G4endl;
  #endif

  // The string can fragment. At least two particles can be produced.
  LeftVector =new G4KineticTrackVector;
  G4KineticTrackVector * RightVector=new G4KineticTrackVector;

  G4ExcitedString *theStringInCMS=CPExcited(theString);

  #ifdef debug_LUNDfragmentation
  G4cout<<"CMS Left  mom "<<theStringInCMS->GetLeftParton()->Get4Momentum()<<G4endl;
  G4cout<<"CMS Right mom "<<theStringInCMS->GetRightParton()->Get4Momentum()<<G4endl;
  #endif

  G4LorentzRotation toCms=theStringInCMS->TransformToAlignedCms();

  #ifdef debug_LUNDfragmentation
  G4cout<<"aligCMS Left  mom "<<theStringInCMS->GetLeftParton()->Get4Momentum()<<G4endl;
  G4cout<<"aligCMS Right mom "<<theStringInCMS->GetRightParton()->Get4Momentum()<<G4endl;
  G4cout<<G4endl<<"LUND StringFragm: String Mass " <<theStringInCMS->Get4Momentum().mag()<<" Pz Lab "
                                                   <<theStringInCMS->Get4Momentum().pz()
                <<"------------------------------------"<<G4endl;
  G4cout<<"String ends and Direction "<<theStringInCMS->GetLeftParton()->GetPDGcode()<<" "
                                      <<theStringInCMS->GetRightParton()->GetPDGcode()<<" "
                                      <<theStringInCMS->GetDirection()<< G4endl;
  #endif

  G4bool success = Loop_toFragmentString(theStringInCMS, LeftVector, RightVector);

  delete theStringInCMS;

  if ( ! success )
  {
    std::for_each(LeftVector->begin(), LeftVector->end(), DeleteKineticTrack());
    LeftVector->clear();
    std::for_each(RightVector->begin(), RightVector->end(), DeleteKineticTrack());
    delete RightVector;
    return LeftVector;
  }

  // Join Left- and RightVector into LeftVector in correct order.
  while(!RightVector->empty())  /* Loop checking, 07.08.2015, A.Ribon */
  {
    LeftVector->push_back(RightVector->back());
    RightVector->erase(RightVector->end()-1);
  }
  delete RightVector;

  CalculateHadronTimePosition(theString.Get4Momentum().mag(), LeftVector);

  G4LorentzRotation toObserverFrame(toCms.inverse());

  G4double TimeOftheStringCreation=theString.GetTimeOfCreation();
  G4ThreeVector PositionOftheStringCreation(theString.GetPosition());

  for(size_t C1 = 0; C1 < LeftVector->size(); C1++)
  {
    G4KineticTrack* Hadron = LeftVector->operator[](C1);
    G4LorentzVector Momentum = Hadron->Get4Momentum();
    //G4cout<<"Hadron "<<Hadron->GetDefinition()->GetParticleName()<<" "<<Momentum<<G4endl;
    Momentum = toObserverFrame*Momentum;
    Hadron->Set4Momentum(Momentum);

    G4LorentzVector Coordinate(Hadron->GetPosition(), Hadron->GetFormationTime());
    Momentum = toObserverFrame*Coordinate;
    Hadron->SetFormationTime(TimeOftheStringCreation + Momentum.e() - fermi/c_light);
    G4ThreeVector aPosition(Momentum.vect());
    Hadron->SetPosition(PositionOftheStringCreation+aPosition);
  }

  return LeftVector;
}

//----------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::IsFragmentable(const G4FragmentingString * const string)
{
  SetMinimalStringMass(string);
  //return sqr(MinimalStringMass + WminLUND) < string->Get4Momentum().mag2();
  //G4cout<<"MinM StrM "<<MinimalStringMass<<" "<< string->Get4Momentum().mag()<<G4endl;
  return MinimalStringMass < string->Get4Momentum().mag();
}

//----------------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::StopFragmenting(const G4FragmentingString * const string)
{
  SetMinimalStringMass(string);

  if (string->FourQuarkString())
  {
    return G4UniformRand() < G4Exp(-0.0005*(string->Mass() - MinimalStringMass));
  } else {
    // Uzhi 23 Jan. 2015 0.88 -> 0.66      
    return G4UniformRand() < G4Exp(-0.66e-6*(string->Mass()*string->Mass() -
				             MinimalStringMass*MinimalStringMass));
  }
}

//----------------------------------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::SplitLast(G4FragmentingString * string,
                                            G4KineticTrackVector * LeftVector,
                                            G4KineticTrackVector * RightVector)
{
  //... perform last cluster decay

  #ifdef debug_LUNDfragmentation
  G4cout<<"Split last-----------------------------------------"<<G4endl;
  #endif

  G4LorentzVector Str4Mom=string->Get4Momentum();
  G4ThreeVector ClusterVel=string->Get4Momentum().boostVector();
  G4double StringMass=string->Mass();

  G4ParticleDefinition * LeftHadron(0), * RightHadron(0);

  NumberOf_FS=0;
  for(G4int i=0; i<35; i++) {FS_Weight[i]=0.;}

  #ifdef debug_LUNDfragmentation
  G4cout<<"StrMass "<<StringMass<<" q "
        <<string->GetLeftParton()->GetParticleName()<<" "
        <<string->GetRightParton()->GetParticleName()<<G4endl;
  #endif

  string->SetLeftPartonStable(); // to query quark contents..

  if (string->FourQuarkString() )
  {
    // The string is qq-qqbar type. Diquarks are on the string ends
    //G4cout<<"The string is qq-qqbar type. Diquarks are on the string ends"<<G4endl;

    if(StringMass-MinimalStringMass < 0.)
    {
      if (! Diquark_AntiDiquark_belowThreshold_lastSplitting(string, LeftHadron, RightHadron) ) 
      {
        return false;
      }
    } else {
      Diquark_AntiDiquark_aboveThreshold_lastSplitting(string, LeftHadron, RightHadron);

      if(NumberOf_FS == 0) return false;

      G4int sampledState = SampleState();
      if(string->GetLeftParton()->GetPDGEncoding() < 0)
      {
	LeftHadron =FS_LeftHadron[sampledState];
	RightHadron=FS_RightHadron[sampledState];
      } else {
	LeftHadron =FS_RightHadron[sampledState];
	RightHadron=FS_LeftHadron[sampledState];
      }
      //G4cout<<"Selected "<<SampledState<<" "<<LeftHadron->GetParticleName()<<" "<<RightHadron->GetParticleName()<<G4endl;
    }
  } else {
    if (string->DecayIsQuark() && string->StableIsQuark() ) {  //... there are quarks on cluster ends

      #ifdef debug_LUNDfragmentation
      G4cout<<"Q Q string LastSplit"<<G4endl;
      #endif

      Quark_AntiQuark_lastSplitting(string, LeftHadron, RightHadron);
     } else {  //... there is a Diquark on one of the cluster ends

       #ifdef debug_LUNDfragmentation
       G4cout<<"DiQ Q string Last Split"<<G4endl;
       #endif

       Quark_Diquark_lastSplitting(string, LeftHadron, RightHadron);
     }
		
     if(NumberOf_FS == 0) return false;
     G4int sampledState = SampleState();
     LeftHadron =FS_LeftHadron[sampledState];
     RightHadron=FS_RightHadron[sampledState];

     #ifdef debug_LUNDfragmentation
     G4cout<<"Selected LeftHad RightHad "<<sampledState<<" "
           <<LeftHadron->GetParticleName()<<" "<<RightHadron->GetParticleName()<<G4endl;
     #endif

  }

  G4LorentzVector  LeftMom, RightMom;
  G4ThreeVector    Pos;

  Sample4Momentum(&LeftMom,  LeftHadron->GetPDGMass(), &RightMom, RightHadron->GetPDGMass(), StringMass);

  LeftMom.boost(ClusterVel);
  RightMom.boost(ClusterVel);

  LeftVector->push_back(new G4KineticTrack(LeftHadron, 0, Pos, LeftMom));
  RightVector->push_back(new G4KineticTrack(RightHadron, 0, Pos, RightMom));

  return true;
}

//----------------------------------------------------------------------------------------------------------
void G4LundStringFragmentation::Sample4Momentum(G4LorentzVector*     Mom, G4double     Mass, 
                                                G4LorentzVector* AntiMom, G4double AntiMass, 
                                                G4double InitialMass) 
{
  // ------ Sampling of momenta of 2 last produced hadrons --------------------
  G4ThreeVector Pt;
  G4double MassMt2, AntiMassMt2;
  G4double AvailablePz, AvailablePz2;
  //AR-Oct2017  G4double ProbIsotropy = sqr(sqr(938.0/InitialMass));

  #ifdef debug_LUNDfragmentation
  G4cout<<"Sampling of momenta of 2 last produced hadrons ----------------"<<G4endl;
  G4cout<<"Masses "<<InitialMass<<" "<<Mass<<" "<<AntiMass<<G4endl; 
  #endif

  //AR-Oct2017  if((Mass > 930. || AntiMass > 930.) && (G4UniformRand() < ProbIsotropy)) {  //If there is a baryon
  //AR-Oct2017    // ----------------- Isotropic decay ------------------------------------
  //AR-Oct2017    G4double r_val = sqr(InitialMass*InitialMass - Mass*Mass - AntiMass*AntiMass) - sqr(2.*Mass*AntiMass);
  //AR-Oct2017    G4double Pabs = (r_val > 0.)? std::sqrt(r_val)/(2.*InitialMass) : 0;
  //AR-Oct2017    G4cout<<"P for isotr decay "<<Pabs<<G4endl;
  //AR-Oct2017    //... sample unit vector
  //AR-Oct2017    G4double pz =1. - 2.*G4UniformRand();
  //AR-Oct2017    G4double st = std::sqrt(1. - pz * pz)*Pabs;
  //AR-Oct2017    G4double phi = 2.*pi*G4UniformRand();
  //AR-Oct2017    G4double px = st*std::cos(phi);
  //AR-Oct2017    G4double py = st*std::sin(phi);
  //AR-Oct2017    pz *= Pabs;
  //AR-Oct2017    Mom->setPx(px); Mom->setPy(py); Mom->setPz(pz);
  //AR-Oct2017    Mom->setE(std::sqrt(Pabs*Pabs + Mass*Mass));
  //AR-Oct2017    AntiMom->setPx(-px); AntiMom->setPy(-py); AntiMom->setPz(-pz);
  //AR-Oct2017    AntiMom->setE (std::sqrt(Pabs*Pabs + AntiMass*AntiMass));
  //AR-Oct2017    //G4int Uzhi; G4cin>>Uzhi;
  //AR-Oct2017  } else {

  const G4int maxNumberOfLoops = 1000;

  G4double SigmaQTw = SigmaQT;  //AR-Oct2017
  if ( Mass > 930. && AntiMass > 930.) SigmaQT *=(1.0-0.45*sqr((Mass+AntiMass)/InitialMass));  //AR-Oct2017
  // In the above condition, "&&" is used instead of "||" (as in development) to be conservative
  // i.e. to be applied (for the time being) only in annihilations.

  G4int loopCounter = 0;
  do
  {
    // GF 22-May-09, limit sampled pt to allowed range

    G4double termD = InitialMass*InitialMass -Mass*Mass - AntiMass*AntiMass;
    G4double termab = 4*sqr(Mass*AntiMass);
    G4double termN = 2*termD + 4*Mass*Mass + 4*AntiMass*AntiMass;
    G4double pt2max=(termD*termD - termab )/ termN ;
    //G4cout<<"Anis "<<pt2max<<" "<<(termD*termD-termab)/(4.*InitialMass*InitialMass)<<G4endl;

    Pt=SampleQuarkPt(std::sqrt(pt2max)); Pt.setZ(0); G4double Pt2=Pt.mag2();
    //G4cout<<"Sampl pt2 "<<Pt2<<G4endl;
    MassMt2    =     Mass *     Mass + Pt2;
    AntiMassMt2= AntiMass * AntiMass + Pt2;

    AvailablePz2= sqr(InitialMass*InitialMass - MassMt2 - AntiMassMt2) - 4.*MassMt2*AntiMassMt2;
    
  } while( (AvailablePz2 < 0.) &&  // GF will occur only for numerical precision problem with limit in sampled pt
           ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */

  if ( loopCounter >= maxNumberOfLoops ) {
    AvailablePz2 = 0.0;
  }

  AvailablePz2 /=(4.*InitialMass*InitialMass);
  AvailablePz = std::sqrt(AvailablePz2);

  G4double Px=Pt.getX();
  G4double Py=Pt.getY();

  Mom->setPx(Px); Mom->setPy(Py); Mom->setPz(AvailablePz);
  Mom->setE(std::sqrt(MassMt2+AvailablePz2));

  AntiMom->setPx(-Px); AntiMom->setPy(-Py); AntiMom->setPz(-AvailablePz);
  AntiMom->setE (std::sqrt(AntiMassMt2+AvailablePz2));

  if ( Mass > 930. && AntiMass > 930.) SigmaQT = SigmaQTw;  //AR-Oct2017
  // In the above condition, "&&" is used instead of "||" (as in development) to be conservative
  // i.e. to be applied (for the time being) only in annihilations.

  //AR-Oct2017  }
}

//-----------------------------------------------------------------------------
G4LorentzVector * G4LundStringFragmentation::SplitEandP(G4ParticleDefinition * pHadron,
                                                        G4FragmentingString *   string, 
                                                        G4FragmentingString * newString)
{ 
  G4LorentzVector String4Momentum=string->Get4Momentum();
  G4double StringMT2=string->MassT2();
  G4double StringMT =std::sqrt(StringMT2);

  G4double HadronMass = pHadron->GetPDGMass();

  SetMinimalStringMass(newString);

  #ifdef debug_LUNDfragmentation
  G4cout<<G4endl<<"Start LUND SplitEandP "<<G4endl;
  G4cout<<"String 4 mom, String M and Mt "<<String4Momentum<<" "<<String4Momentum.mag()<<" "<<std::sqrt(StringMT2)<<G4endl;
  G4cout<<"Hadron "<<pHadron->GetParticleName()<<G4endl;
  G4cout<<"HadM MinimalStringMassLeft StringM hM+sM "<<HadronMass<<" "<<MinimalStringMass<<" "
        <<String4Momentum.mag()<<" "<<HadronMass+MinimalStringMass<<G4endl;
  #endif

  if(HadronMass + MinimalStringMass > string->Mass()) 
  {
    #ifdef debug_LUNDfragmentation
    G4cout<<"Mass of the string is not sufficient to produce the hadron!"<<G4endl;
    #endif

    return 0;
  } // have to start all over!

  String4Momentum.setPz(0.);
  G4ThreeVector StringPt=String4Momentum.vect();

  // calculate and assign hadron transverse momentum component HadronPx and HadronPy
  G4ThreeVector HadronPt    , RemSysPt; 
  G4double      HadronMassT2, ResidualMassT2;

  //...  sample Pt of the hadron
  G4int attempt=0;
  do
  {
    attempt++; if(attempt > StringLoopInterrupt) return 0;

    HadronPt =SampleQuarkPt()  + string->DecayPt();	
    HadronPt.setZ(0);
    RemSysPt = StringPt - HadronPt;

    HadronMassT2 = sqr(HadronMass) + HadronPt.mag2();
    ResidualMassT2=sqr(MinimalStringMass) + RemSysPt.mag2();

  } while(std::sqrt(HadronMassT2) + std::sqrt(ResidualMassT2) > StringMT);  /* Loop checking, 07.08.2015, A.Ribon */

  //...  sample z to define hadron longitudinal momentum and energy
  //... but first check the available phase space

  G4double Pz2 = (sqr(StringMT2 - HadronMassT2 - ResidualMassT2) -
		  4*HadronMassT2 * ResidualMassT2)/4./StringMT2;

  if(Pz2 < 0 ) {return 0;}          // have to start all over!

  //... then compute allowed z region  z_min <= z <= z_max

  G4double Pz = std::sqrt(Pz2);
  G4double zMin = (std::sqrt(HadronMassT2+Pz2) - Pz)/std::sqrt(StringMT2);
  //G4double zMin = (std::sqrt(HadronMassT2+Pz2) - 0.)/std::sqrt(StringMT2);
  G4double zMax = (std::sqrt(HadronMassT2+Pz2) + Pz)/std::sqrt(StringMT2);

  if (zMin >= zMax) return 0;		// have to start all over!

  G4double z = GetLightConeZ(zMin, zMax, string->GetDecayParton()->GetPDGEncoding(),
                             pHadron, HadronPt.x(), HadronPt.y());

  //... now compute hadron longitudinal momentum and energy
  // longitudinal hadron momentum component HadronPz

  HadronPt.setZ(0.5* string->GetDecayDirection() *
		(z * string->LightConeDecay() - HadronMassT2/(z * string->LightConeDecay())));
  G4double HadronE  = 0.5* (z * string->LightConeDecay() +
			    HadronMassT2/(z * string->LightConeDecay()));

  G4LorentzVector * a4Momentum= new G4LorentzVector(HadronPt,HadronE);

  #ifdef debug_LUNDfragmentation
  G4cout<<"string->LightConeDecay() "<<string->LightConeDecay()<<G4endl;
  G4cout<<"HadronPt,HadronE "<<HadronPt<<" "<<HadronE<<G4endl;
  G4cout<<"String4Momentum "<<String4Momentum<<G4endl;
  //G4int Uzhi; G4cin>>Uzhi;
  G4cout<<"Out of LUND SplitEandP "<<G4endl;
  #endif

  return a4Momentum;
}

//-----------------------------------------------------------------------------------------
G4double G4LundStringFragmentation::
GetLightConeZ(G4double zmin, G4double zmax, G4int PDGEncodingOfDecayParton,
	      G4ParticleDefinition* pHadron, G4double Px, G4double Py)
{
  G4double Mass = pHadron->GetPDGMass();
  //G4int HadronEncoding=std::abs(pHadron->GetPDGEncoding());

  G4double Mt2 = Px*Px + Py*Py + Mass*Mass;

  G4double  alund;
  G4double zOfMaxyf(0.), maxYf(1.), z(0.), yf(1.);
  if(std::abs(PDGEncodingOfDecayParton) < 1000)
  { // ---------------- Quark fragmentation ----------------------
    alund=0.7/GeV/GeV;
    //    If blund get restored, you MUST adapt the calculation of zOfMaxyf.
    //    const G4double  blund = 1;

    zOfMaxyf=alund*Mt2/(alund*Mt2 + 1.);
    maxYf=(1-zOfMaxyf)/zOfMaxyf * G4Exp(-alund*Mt2/zOfMaxyf);

    const G4int maxNumberOfLoops = 1000;
    G4int loopCounter = 0;
    do
    {
      z = zmin + G4UniformRand()*(zmax-zmin);
      yf = (1-z)/z * G4Exp(-alund*Mt2/z);
      //yf = G4Pow::GetInstance()->powA(1.0-z,blund)/z*G4Exp(-alund*Mt2/z
    } while ( (G4UniformRand()*maxYf > yf) && ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */
    if ( loopCounter >= maxNumberOfLoops ) {
      z = 0.5*(zmin + zmax);  // Just a value between zmin and zmax, no physics considerations at all! 
    }
    return z;
  }

  if(std::abs(PDGEncodingOfDecayParton) > 1000)
  {
    /*
    if(HadronEncoding < 3000)
    {
      maxYf=(zmax-zmin);
      do
      {
        z = zmin + G4UniformRand()*(zmax-zmin);
	//yf=G4Exp(-sqr(z-Zc)/2/sqr(0.28));  // 0.42 0.632 0.28 a'la UrQMD
	yf =(z-zmin);
      } while (G4UniformRand()*maxYf > yf);
    } else { // Strange baryons
      z = zmin + G4UniformRand()*(zmax-zmin);
    }
    */
    z = zmin + G4UniformRand()*(zmax-zmin);
  }
  return z;
}

//------------------------------------------------------------------------
G4double G4LundStringFragmentation::lambda(G4double S, G4double m1_Sqr, G4double m2_Sqr)
{ 
  G4double lam = sqr(S - m1_Sqr - m2_Sqr) - 4.*m1_Sqr*m2_Sqr;
  return lam;
}


//------------------------------------------------------------------------
//------------------------------------------------------------------------
// Internal methods introduced to improve the code structure (AR Nov 2011)
//------------------------------------------------------------------------
//------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::Loop_toFragmentString(G4ExcitedString * & theStringInCMS, 
                                                        G4KineticTrackVector * & LeftVector, 
                                                        G4KineticTrackVector * & RightVector)
{
  #ifdef debug_LUNDfragmentation
  G4cout<<"Loop_toFrag "<<theStringInCMS->GetLeftParton()->GetPDGcode()<<" "
                        <<theStringInCMS->GetRightParton()->GetPDGcode()<<" "
                        <<theStringInCMS->GetDirection()<< G4endl;
  #endif

  G4bool final_success=false;
  G4bool inner_success=true;
  G4int attempt=0;
  while ( ! final_success && attempt++ < StringLoopInterrupt )  /* Loop checking, 07.08.2015, A.Ribon */
  { // If the string fragmentation do not be happend, repeat the fragmentation.

    G4FragmentingString *currentString=new G4FragmentingString(*theStringInCMS);
    //G4cout<<"Main loop start whilecounter "<<attempt<<G4endl;
    // Cleaning up the previously produced hadrons
    std::for_each(LeftVector->begin(), LeftVector->end(), DeleteKineticTrack());
    LeftVector->clear();
    std::for_each(RightVector->begin(), RightVector->end(), DeleteKineticTrack());
    RightVector->clear();

    // Main fragmentation loop until the string will not be able to fragment
    inner_success=true;  // set false on failure.
    const G4int maxNumberOfLoops = 1000;
    G4int loopCounter = -1;
    while ( (! StopFragmenting(currentString)) && ++loopCounter < maxNumberOfLoops )  /* Loop checking, 07.08.2015, A.Ribon */
    {       // Split current string into hadron + new string
      #ifdef debug_LUNDfragmentation
      G4cout<<"The string can fragment. "<<G4endl;;
      //G4cout<<"1 "<<currentString->GetDecayDirection()<<G4endl;
      #endif
      G4FragmentingString *newString=0;  // used as output from SplitUp.
      G4KineticTrack * Hadron=Splitup(currentString,newString);
      if ( Hadron != 0 )  // Store the hadron                               
      {
        #ifdef debug_LUNDfragmentation
        G4cout<<"Hadron prod at fragm. "<<Hadron->GetDefinition()->GetParticleName()<<G4endl;
        //G4cout<<"2 "<<currentString->GetDecayDirection()<<G4endl;
        #endif

	if ( currentString->GetDecayDirection() > 0 ) {
	  LeftVector->push_back(Hadron); 
        } else {
	  RightVector->push_back(Hadron);
        }
	delete currentString;
	currentString=newString;
      }
    };
    if ( loopCounter >= maxNumberOfLoops ) {
      inner_success=false;
    }

    // Split remaining string into 2 final hadrons.
    #ifdef debug_LUNDfragmentation
    G4cout<<"Split remaining string into 2 final hadrons."<<G4endl;
    #endif

    if ( inner_success && SplitLast(currentString, LeftVector, RightVector) )
    {
      final_success=true;
    }

    //final_success=true;
    delete currentString;
  }  // End of the loop where we try to fragment the string.
  return final_success;
}

//----------------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::
Diquark_AntiDiquark_belowThreshold_lastSplitting(G4FragmentingString * & string,
                                                 G4ParticleDefinition * & LeftHadron,
                                                 G4ParticleDefinition * & RightHadron)
{
  G4double StringMass   = string->Mass();
  G4int cClusterInterrupt = 0;
  do
  {
    //G4cout<<"cClusterInterrupt "<<cClusterInterrupt<<G4endl;
    if (cClusterInterrupt++ >= ClusterLoopInterrupt) return false;

    G4int LeftQuark1= string->GetLeftParton()->GetPDGEncoding()/1000;
    G4int LeftQuark2=(string->GetLeftParton()->GetPDGEncoding()/100)%10;

    G4int RightQuark1= string->GetRightParton()->GetPDGEncoding()/1000;
    G4int RightQuark2=(string->GetRightParton()->GetPDGEncoding()/100)%10;

    if(G4UniformRand()<0.5) {
      LeftHadron =hadronizer->Build(FindParticle( LeftQuark1), FindParticle(RightQuark1));
      RightHadron=hadronizer->Build(FindParticle( LeftQuark2), FindParticle(RightQuark2));
    } else {
      LeftHadron =hadronizer->Build(FindParticle( LeftQuark1), FindParticle(RightQuark2));
      RightHadron=hadronizer->Build(FindParticle( LeftQuark2), FindParticle(RightQuark1));
    }

    //... repeat procedure, if mass of cluster is too low to produce hadrons
    //... ClusterMassCut = 0.15*GeV model parameter
  }
  while ((StringMass <= LeftHadron->GetPDGMass() + RightHadron->GetPDGMass()));  /* Loop checking, 07.08.2015, A.Ribon */

  return true;
}

//----------------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::
Diquark_AntiDiquark_aboveThreshold_lastSplitting(G4FragmentingString * & string,
                                                 G4ParticleDefinition * & LeftHadron,
                                                 G4ParticleDefinition * & RightHadron)
{
  // StringMass-MinimalStringMass > 0. Creation of 2 baryons is possible ----

  G4double StringMass   = string->Mass();
  G4double StringMassSqr= sqr(StringMass); 
  G4ParticleDefinition * Di_Quark;
  G4ParticleDefinition * Anti_Di_Quark;

  if(string->GetLeftParton()->GetPDGEncoding() < 0) {
    Anti_Di_Quark   =string->GetLeftParton();
    Di_Quark=string->GetRightParton();
  } else {
    Anti_Di_Quark   =string->GetRightParton();
    Di_Quark=string->GetLeftParton();
  }

  G4int IDAnti_di_quark    =Anti_Di_Quark->GetPDGEncoding();
  G4int AbsIDAnti_di_quark =std::abs(IDAnti_di_quark);
  G4int IDdi_quark         =Di_Quark->GetPDGEncoding();
  G4int AbsIDdi_quark      =std::abs(IDdi_quark);

  G4int ADi_q1=AbsIDAnti_di_quark/1000;
  G4int ADi_q2=(AbsIDAnti_di_quark-ADi_q1*1000)/100;

  G4int Di_q1=AbsIDdi_quark/1000;
  G4int Di_q2=(AbsIDdi_quark-Di_q1*1000)/100;

  NumberOf_FS=0;
  for(G4int ProdQ=1; ProdQ < 4; ProdQ++)
  {
    G4int StateADiQ=0;
    const G4int maxNumberOfLoops = 1000;
    G4int loopCounter = 0;
    do
    {
      LeftHadron=G4ParticleTable::GetParticleTable()->FindParticle(
			-Baryon[ADi_q1-1][ADi_q2-1][ProdQ-1][StateADiQ]);
      G4double LeftHadronMass=LeftHadron->GetPDGMass();

      //G4cout<<"Anti Bar "<<LeftHadron->GetParticleName()<<G4endl;

      G4int StateDiQ=0;
      const G4int maxNumberOfInternalLoops = 1000;
      G4int internalLoopCounter = 0;
      do
      {
        RightHadron=G4ParticleTable::GetParticleTable()->FindParticle(
			+Baryon[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]);
	G4double RightHadronMass=RightHadron->GetPDGMass();

	if(StringMass > LeftHadronMass + RightHadronMass)
	{
          if ( NumberOf_FS > 34 ) {
            G4ExceptionDescription ed;
            ed << " NumberOf_FS exceeds its limit: NumberOf_FS=" << NumberOf_FS << G4endl;
            G4Exception( "G4LundStringFragmentation::Diquark_AntiDiquark_aboveThreshold_lastSplitting ",
                         "HAD_LUND_001", JustWarning, ed );
            NumberOf_FS = 34;
          }

	  G4double FS_Psqr=lambda(StringMassSqr,sqr(LeftHadronMass), sqr(RightHadronMass));
	  //FS_Psqr=1.;
	  FS_Weight[NumberOf_FS]=std::sqrt(FS_Psqr)*FS_Psqr*
				 BaryonWeight[ADi_q1-1][ADi_q2-1][ProdQ-1][StateADiQ]*
				 BaryonWeight[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]*
				 Prob_QQbar[ProdQ-1];

	  FS_LeftHadron[NumberOf_FS] = LeftHadron;
	  FS_RightHadron[NumberOf_FS]= RightHadron;

	  NumberOf_FS++;
	} 

	StateDiQ++;

      } while( (Baryon[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]!=0) && 
               ++internalLoopCounter < maxNumberOfInternalLoops );  /* Loop checking, 07.08.2015, A.Ribon */
      if ( internalLoopCounter >= maxNumberOfInternalLoops ) {
        return false;
      }
 
      StateADiQ++;
    } while( (Baryon[ADi_q1-1][ADi_q2-1][ProdQ-1][StateADiQ]!=0) &&
             ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */
    if ( loopCounter >= maxNumberOfLoops ) {
      return false;
    }
  }

  return true;
}

//----------------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::
Quark_Diquark_lastSplitting(G4FragmentingString * & string,
                            G4ParticleDefinition * & LeftHadron,
                            G4ParticleDefinition * & RightHadron)
{
  G4double StringMass   = string->Mass();
  G4double StringMassSqr= sqr(StringMass);

  G4ParticleDefinition * Di_Quark;
  G4ParticleDefinition * Quark;

  if(string->GetLeftParton()->GetParticleSubType()== "quark") {
    Quark   =string->GetLeftParton();
    Di_Quark=string->GetRightParton();
  } else {
    Quark   =string->GetRightParton();
    Di_Quark=string->GetLeftParton();
  }

  G4int IDquark        =Quark->GetPDGEncoding();
  G4int AbsIDquark     =std::abs(IDquark);
  G4int IDdi_quark   =Di_Quark->GetPDGEncoding();
  G4int AbsIDdi_quark=std::abs(IDdi_quark);
  G4int Di_q1=AbsIDdi_quark/1000;
  G4int Di_q2=(AbsIDdi_quark-Di_q1*1000)/100;

  G4int              SignDiQ= 1;
  if(IDdi_quark < 0) SignDiQ=-1;

  NumberOf_FS=0;
  for(G4int ProdQ=1; ProdQ < 4; ProdQ++)
  {
    G4int SignQ;
    if(IDquark > 0) {
      SignQ=-1;
      if(IDquark == 2)                   SignQ= 1;
      if((IDquark == 1) && (ProdQ == 3)) SignQ= 1; // K0
      if((IDquark == 3) && (ProdQ == 1)) SignQ=-1; // K0bar
    } else {
      SignQ= 1;
      if(IDquark == -2)                  SignQ=-1;
      if((IDquark ==-1) && (ProdQ == 3)) SignQ=-1; // K0bar
      if((IDquark ==-3) && (ProdQ == 1)) SignQ= 1; // K0
    }

    if(AbsIDquark == ProdQ)            SignQ= 1;

    //G4cout<<G4endl;
    //G4cout<<"Insert QQbar "<<ProdQ<<" Sign "<<SignQ<<G4endl;

    G4int StateQ=0;
    const G4int maxNumberOfLoops = 1000;
    G4int loopCounter = 0;
    do
    {
      LeftHadron=G4ParticleTable::GetParticleTable()->FindParticle(
                       SignQ*Meson[AbsIDquark-1][ProdQ-1][StateQ]);
      G4double LeftHadronMass=LeftHadron->GetPDGMass();

      G4int StateDiQ=0;
      const G4int maxNumberOfInternalLoops = 1000;
      G4int internalLoopCounter = 0;
      do
      {
	RightHadron=G4ParticleTable::GetParticleTable()->FindParticle(
                          SignDiQ*Baryon[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]);
	G4double RightHadronMass=RightHadron->GetPDGMass();

	if(StringMass > LeftHadronMass + RightHadronMass)
	{
          if ( NumberOf_FS > 34 ) {
            G4ExceptionDescription ed;
            ed << " NumberOf_FS exceeds its limit: NumberOf_FS=" << NumberOf_FS << G4endl;
            G4Exception( "G4LundStringFragmentation::Quark_Diquark_lastSplitting ",
                         "HAD_LUND_002", JustWarning, ed );
            NumberOf_FS = 34;
          }

	  G4double FS_Psqr=lambda(StringMassSqr,sqr(LeftHadronMass),sqr(RightHadronMass));
	  FS_Weight[NumberOf_FS]=std::sqrt(FS_Psqr)*MesonWeight[AbsIDquark-1][ProdQ-1][StateQ]*
			                            BaryonWeight[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]*
					            Prob_QQbar[ProdQ-1];

	  FS_LeftHadron[NumberOf_FS] = LeftHadron;
	  FS_RightHadron[NumberOf_FS]= RightHadron;

	  NumberOf_FS++;
        }

	StateDiQ++;

      } while( (Baryon[Di_q1-1][Di_q2-1][ProdQ-1][StateDiQ]!=0) &&
               ++internalLoopCounter < maxNumberOfInternalLoops );  /* Loop checking, 07.08.2015, A.Ribon */
      if ( internalLoopCounter >= maxNumberOfInternalLoops ) {
        return false;
      }

      StateQ++;
    } while( (Meson[AbsIDquark-1][ProdQ-1][StateQ]!=0) &&
             ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */
    if ( loopCounter >= maxNumberOfLoops ) {
      return false;
    }
  }

  return true;
}

//----------------------------------------------------------------------------------------
G4bool G4LundStringFragmentation::
Quark_AntiQuark_lastSplitting(G4FragmentingString * & string,
                              G4ParticleDefinition * & LeftHadron,
                              G4ParticleDefinition * & RightHadron)
{
  G4double StringMass   = string->Mass();
  G4double StringMassSqr= sqr(StringMass);

  G4ParticleDefinition * Quark;
  G4ParticleDefinition * Anti_Quark;

  if(string->GetLeftParton()->GetPDGEncoding()>0) {
    Quark     =string->GetLeftParton();
    Anti_Quark=string->GetRightParton();
  } else {
    Quark     =string->GetRightParton();
    Anti_Quark=string->GetLeftParton();
  }

  G4int IDquark        =Quark->GetPDGEncoding();
  G4int AbsIDquark     =std::abs(IDquark);
  G4int IDanti_quark   =Anti_Quark->GetPDGEncoding();
  G4int AbsIDanti_quark=std::abs(IDanti_quark);

  NumberOf_FS=0;
  for(G4int ProdQ=1; ProdQ < 4; ProdQ++)
  {
    G4int                              SignQ=-1;
    if(IDquark == 2)                   SignQ= 1;
    if((IDquark == 1) && (ProdQ == 3)) SignQ= 1; // K0
    if((IDquark == 3) && (ProdQ == 1)) SignQ=-1; // K0bar
    if(IDquark == ProdQ)               SignQ= 1;

    G4int                                   SignAQ= 1;
    if(IDanti_quark == -2)                  SignAQ=-1;
    if((IDanti_quark ==-1) && (ProdQ == 3)) SignAQ=-1; // K0bar
    if((IDanti_quark ==-3) && (ProdQ == 1)) SignAQ= 1; // K0
    if(AbsIDanti_quark == ProdQ)            SignAQ= 1;

    G4int StateQ=0;
    const G4int maxNumberOfLoops = 1000;
    G4int loopCounter = 0;
    do
    {
      LeftHadron=G4ParticleTable::GetParticleTable()->FindParticle(
                       SignQ*Meson[AbsIDquark-1][ProdQ-1][StateQ]);
      G4double LeftHadronMass=LeftHadron->GetPDGMass();

      G4int StateAQ=0;
      const G4int maxNumberOfInternalLoops = 1000;
      G4int internalLoopCounter = 0;
      do
      {
        RightHadron=G4ParticleTable::GetParticleTable()->FindParticle(
                          SignAQ*Meson[AbsIDanti_quark-1][ProdQ-1][StateAQ]);
	G4double RightHadronMass=RightHadron->GetPDGMass();

	if(StringMass > LeftHadronMass + RightHadronMass)
	{
          if ( NumberOf_FS > 34 ) {
            G4ExceptionDescription ed;
            ed << " NumberOf_FS exceeds its limit: NumberOf_FS=" << NumberOf_FS << G4endl;
            G4Exception( "G4LundStringFragmentation::Quark_AntiQuark_lastSplitting ",
                         "HAD_LUND_003", JustWarning, ed );
            NumberOf_FS = 34;
          }

          G4double FS_Psqr=lambda(StringMassSqr,sqr(LeftHadronMass),sqr(RightHadronMass));
	  //FS_Psqr=1.;
	  FS_Weight[NumberOf_FS]=std::sqrt(FS_Psqr)*MesonWeight[AbsIDquark-1][ProdQ-1][StateQ]*
					            MesonWeight[AbsIDanti_quark-1][ProdQ-1][StateAQ]*
					            Prob_QQbar[ProdQ-1];

	  if(string->GetLeftParton()->GetPDGEncoding()>0) {
	    FS_LeftHadron[NumberOf_FS] = RightHadron;
	    FS_RightHadron[NumberOf_FS]= LeftHadron;
	  } else {
	    FS_LeftHadron[NumberOf_FS] = LeftHadron;
	    FS_RightHadron[NumberOf_FS]= RightHadron;
          }
	  NumberOf_FS++;

        }

	StateAQ++;
      } while( (Meson[AbsIDanti_quark-1][ProdQ-1][StateAQ]!=0) &&
               ++internalLoopCounter < maxNumberOfInternalLoops );  /* Loop checking, 07.08.2015, A.Ribon */
      if ( internalLoopCounter >= maxNumberOfInternalLoops ) {
        return false;
      }

      StateQ++;
    } while( (Meson[AbsIDquark-1][ProdQ-1][StateQ]!=0) &&
             ++loopCounter < maxNumberOfLoops );  /* Loop checking, 07.08.2015, A.Ribon */
    if ( loopCounter >= maxNumberOfLoops ) {
      return false;
    }
  } 
  return true;
}

//----------------------------------------------------------------------------------------------------------
G4int G4LundStringFragmentation::SampleState(void) 
{
  if ( NumberOf_FS > 34 ) {
    G4ExceptionDescription ed;
    ed << " NumberOf_FS exceeds its limit: NumberOf_FS=" << NumberOf_FS << G4endl;
    G4Exception( "G4LundStringFragmentation::SampleState ", "HAD_LUND_004", JustWarning, ed );
    NumberOf_FS = 34;
  }
 
  G4double SumWeights=0.;

  for(G4int i=0; i<NumberOf_FS; i++) {SumWeights+=FS_Weight[i];}// G4cout<<i<<" "<<FS_Weight[i]<<G4endl;}

  G4double ksi=G4UniformRand();
  G4double Sum=0.;
  G4int indexPosition = 0;

  for(G4int i=0; i<NumberOf_FS; i++)
  {
    Sum+=(FS_Weight[i]/SumWeights);
    indexPosition=i;
    if(Sum >= ksi) break;
  }
  return indexPosition;
}

// Uzhi June 2014 Insert from G4ExcitedStringDecay.cc
//-----------------------------------------------------------------------------

G4ParticleDefinition * G4LundStringFragmentation::
DiQuarkSplitup( G4ParticleDefinition* decay, G4ParticleDefinition *&created )
{
  //... can Diquark break or not?
  if (G4UniformRand() < DiquarkBreakProb )
  {
    //... Diquark break

    G4int stableQuarkEncoding = decay->GetPDGEncoding()/1000;
    G4int decayQuarkEncoding = (decay->GetPDGEncoding()/100)%10;
    if (G4UniformRand() < 0.5)
    {
      G4int Swap = stableQuarkEncoding;
      stableQuarkEncoding = decayQuarkEncoding;
      decayQuarkEncoding = Swap;
    }

    G4int IsParticle=(decayQuarkEncoding>0) ? -1 : +1; // if we have a quark, we need antiquark)

    //G4cout<<"GetStrangeSuppress() "<<GetStrangeSuppress()<<G4endl;
    //G4int Uzhi; G4cin>>Uzhi;

    //G4double StrSup=GetStrangeSuppress(); 
    //StrangeSuppress=0.34;
    pDefPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
    //StrangeSuppress=StrSup;
    //... Build new Diquark
    G4int QuarkEncoding=QuarkPair.second->GetPDGEncoding();
    G4int i10  = std::max(std::abs(QuarkEncoding), std::abs(stableQuarkEncoding));
    G4int i20  = std::min(std::abs(QuarkEncoding), std::abs(stableQuarkEncoding));
    G4int spin = (i10 != i20 && G4UniformRand() <= 0.5)? 1 : 3;
    G4int NewDecayEncoding = -1*IsParticle*(i10 * 1000 + i20 * 100 + spin);
    created = FindParticle(NewDecayEncoding);
    G4ParticleDefinition * decayQuark=FindParticle(decayQuarkEncoding);
    G4ParticleDefinition * had=hadronizer->Build(QuarkPair.first, decayQuark);
    return had;
    // return hadronizer->Build(QuarkPair.first, decayQuark);

  } else {
    //... Diquark does not break

    G4int IsParticle=(decay->GetPDGEncoding()>0) ? +1 : -1;  // if we have a diquark, we need quark)

    G4double StrSup=GetStrangeSuppress();
    StrangeSuppress=0.43;   //0.42 0.38
    pDefPair QuarkPair = CreatePartonPair(IsParticle,false);  // no diquarks wanted
    StrangeSuppress=StrSup;

    created = QuarkPair.second;

    G4ParticleDefinition * had=hadronizer->Build(QuarkPair.first, decay);
    return had;
  }
}

// Uzhi June 2014 End of the inserting

