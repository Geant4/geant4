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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:          G4HIJING_Model.hh
//
// Version:        1.B
// Date:           10/09/2013
// Authors:        Khaled Abdel-Waged 
// Institute:      Umm Al-Qura University
// Country:        Saudi Arabia
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  
#include "G4HIJING_Model.hh"
#ifdef G4_USE_HIJING
#include "G4HIJING_Interface.hh"
//-------------------------------
#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4IonTable.hh"
#include "G4CollisionOutput.hh"
#include "G4V3DNucleus.hh"
#include "G4Track.hh"
#include "G4Nucleus.hh"
#include "G4LorentzRotation.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

//AND->
#include "G4Version.hh"
//AND<-
//----------------new_anti
#include "G4AntiHe3.hh"
#include "G4AntiDeuteron.hh"
#include "G4AntiTriton.hh"
#include "G4AntiAlpha.hh"
//---------------------------
#include <fstream>
#include <string>
#include "HistoManager.hh"  //newkhaled
#include "G4SystemOfUnits.hh"
///////////////////////////////////////////////////////////////////////////

//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4HIJING_Model::G4HIJING_Model(const G4String& nam)
  :G4VIntraNuclearTransportModel(nam), verbose(0)  
{


  if (verbose > 3) {
    G4cout << " >>> G4HIJING_Model default constructor" << G4endl;
  }

#ifdef G4ANALYSIS_USE
fHistoManager   = HistoManager::GetPointer();   //new_khaled
#endif

//
// Set the minimum and maximum range for the HIJING model

  SetMinEnergy(4.0*GeV);
//  SetMaxEnergy(2000.0*TeV);



//
  
// 
  WelcomeMessage();
//
  CurrentEvent=0;
  
//

InitialiseDataTables();


//
}
////////////////////////////////////////////////////////////////////////////////
//
// Destructor
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4HIJING_Model::~G4HIJING_Model (){}
////////////////////////////////////////////////////////////////////////////////

G4ReactionProductVector* G4HIJING_Model::Propagate(G4KineticTrackVector* , 
                              G4V3DNucleus* ) {
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
// ApplyYourself
//
// Member function to process an event, and get information about the products.

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  G4HadFinalState *G4HIJING_Model::ApplyYourself (
  const G4HadProjectile &theTrack, G4Nucleus &theTarget)
{
  G4cout<<"HERE I AM"<<G4endl;
//anti_new
//  ------------------define anti_light_nucleus
const G4ParticleDefinition* anti_deu =
  G4AntiDeuteron::AntiDeuteron();

const G4ParticleDefinition* anti_he3=
  G4AntiHe3::AntiHe3();

const G4ParticleDefinition* anti_tri=
  G4AntiTriton::AntiTriton();

const G4ParticleDefinition* anti_alp=
  G4AntiAlpha::AntiAlpha();

//---------------------------------------------------
//
// The secondaries will be returned in G4HadFinalState &theResult -
// initialise this.  The original track will always be discontinued and
// secondaries followed.
//
  theResult.Clear();
  theResult.SetStatusChange(stopAndKill);

  G4DynamicParticle* cascadeParticle=0; 
//
//
// Get relevant information about the projectile and target (A, Z, energy/nuc,
// momentum, etc).
//


  const G4ParticleDefinition *definitionP = theTrack.GetDefinition();
  const G4double AP        = definitionP->GetBaryonNumber();
  const G4double ZP        = definitionP->GetPDGCharge();
        G4int AT        = theTarget.GetN_asInt();
        G4int ZT        = theTarget.GetZ_asInt();
//  -----------------------------------------------
      G4int id=definitionP->GetPDGEncoding();  //get particle encoding

//      G4cout<<"particle id=========       "<<id<<G4endl;
// ------------------------------------------------
  G4int AP1 = G4lrint(AP);
  G4int ZP1 = G4lrint(ZP);
  G4int AT1 = AT;
  G4int ZT1 = ZT;


// ****************************************************************************
// The following is the parameters necessary to initiate HIJSET() and HIJING()
// ----------------------------------------------------------------------------
//           hiparnt_.ihpr2[3]=0;     //switch off(=0) /  on(=1) jet quenching
//           hiparnt_.ihpr2[2]=1;     //switch on triggered Jet production
// ---------------------------------------------------------------------------
//        hiparnt_.ihnt2[0]=AP1;  //Projectile
hiparnt_.ihnt2[1]=ZP1;
hiparnt_.ihnt2[2]=AT1;  //Target
hiparnt_.ihnt2[3]=ZT1;
hiparnt_.ihnt2[5]=0;    //Special Target
 
      
      if (AP1>1 ||definitionP==anti_deu ||definitionP==anti_he3 
            ||definitionP==anti_tri ||definitionP==anti_alp)
       
       {

        hiparnt_.ihnt2[0]=AP1; 
        hiparnt_.ihnt2[4]=0;    //Special Projectile

      }else if (id==2212) {                  //!proton 

      hiparnt_.ihnt2[0]=1; 
      hiparnt_.ihnt2[4]=2212;


      } else if(id==-2212){            //! anti-proton  

      hiparnt_.ihnt2[0]=1; 
      hiparnt_.ihnt2[4]=-2212;
        
      } else if(id==2112){              //! neutron      
        
      hiparnt_.ihnt2[0]=1; 
      hiparnt_.ihnt2[4]=2112;
        
      } else if(id==-2112){            //! anti-neutron 

     hiparnt_.ihnt2[0]=1; 
     hiparnt_.ihnt2[4]=-2112;
    

      } else if(id==211) {              //! pi+  
      hiparnt_.ihnt2[0]=1;          //needed by HIJING
      hiparnt_.ihnt2[4]=211;

    
      } else if(id==-211) {            //! pi-  

      hiparnt_.ihnt2[0]=1;         //needed by HIJING
      hiparnt_.ihnt2[4]=-211;
        
    } else if(id==321) {              //! K+   

     hiparnt_.ihnt2[0]=1;         //needed by HIJING
     hiparnt_.ihnt2[4]=321;
   
        
    } else if(id==-321) {              //! K-    

      hiparnt_.ihnt2[0]=1;          //needed by HIJING
      hiparnt_.ihnt2[4]=-321;
    
    }  else {

      G4cout << " Sorry, No definition for PROJECTLE for HIJING::"
      <<id<< "found" << G4endl;

      //AND->
#if G4VERSION_NUMBER>=950
      //New signature (9.5) for G4Exception
      //Using G4HadronicException
      throw G4HadronicException(__FILE__,__LINE__,
          "Sorry, no definition for PROJECTILE for HIJING");
#else
    G4Exception(" "); 
#endif
    //AND<-
      }  //end if id

//-------------------------------------------------------
//  -------------identify mass -------------------------
    
     G4int id_n=2112;
     G4int id_p=2212;

     hiparnt_.hint1[7]=std::max(ulmass_ (&id_n),ulmass_ (&id_p));

   
     hiparnt_.hint1[8]=hiparnt_.hint1[7];


     if (hiparnt_.ihnt2[4]!=0) 
     hiparnt_.hint1[7]=ulmass_ (&hiparnt_.ihnt2[4]);  
     //rest mass of the projectile HIJING

//----------------------------------------------------
//  identify Energy
//


G4double m= hiparnt_.hint1[7];   //mass in GeV

G4ThreeVector P3= theTrack.Get4Momentum().vect()/GeV;    
// momentum in GeV

G4double Pbeam=P3.z();      
//momentum in z-direction
       
G4double Ebeam=Eplab(m, Pbeam);  
//calculate Energy of beam

//G4cout<<"mass= "<<m<<"  P3= "<<P3<<endl;

//---------------------------Beam ---------------------------------------

//Lab frame: beam moves in negative z-direction

G4LorentzVector lab= G4LorentzVector(0.0,0.0,-1.0*Pbeam,Ebeam+m);

G4double TotalPbefore=-1.0*lab.z();     
//Calculate z-Momentum before collision
//        
G4double TotalEbefore = lab.e();  
//Calculate Energy before collision

 
//   --------------------------------------------------------
//                     Turn to CM frame: 
//   ---------------------------------------------------------

G4LorentzVector cms = G4LorentzVector(0.0,0.0,0.0,lab.mag());

// ----------------------Get relative speed between frames---------
// ----------------------------------------------------------------
G4LorentzVector Psum=(lab+cms);  //4-Momentum sum  
G4double beta_rel=Psum.beta();
    
//---------------------Transform to equal frame--------------------
//-----------------------------------------------------------------
      
Psum.boost(0.0,0.0,-1.0*beta_rel);

//-----------------Get equal speed velocity between frames--------
G4double betann=Psum.beta();
//G4double gama= Psum.gamma(); 

// ----------Colliding CM Energy per nucleon-nucleon for HIJING-
// ----------------------------------------------------

G4double Ecms=lab.mag();    //CM energy for HIJING 
efrm=Ecms;                 //units are in GeV for HIJING

///////////////////////// initialise/////////////////////

if (CurrentEvent==0)
{


G4cout << "\n initialise HIJING, wait-------"<<G4endl;      

G4cout << "\n"<<G4endl;


//hijset_ (&efrm,&AP1,&ZP1,&AT1,&ZT1);

hijset_ (&efrm);

G4cout << "\n end initialize "<<G4endl;

CurrentEvent=1;
}
////////////////////////////////////////////////////////
//------------------------------------------------------------
// identify impact parameter
   bmin=0.0;
//   bmax=0.5;

bmax=hiparnt_.hipr1[33]+hiparnt_.hipr1[34];
  
//----------------------------------------------

      do 
       {

       G4cout <<"HIJING_Model running-------------" <<G4endl;

      hijing_ (&bmin,&bmax);

       Nproduce=himain1_.natt;  //no of produced particles


    if (Nproduce<2)
      {



G4cout <<"===============Warning====================================="<<G4endl;
G4cout <<"-----------------------------------------------------------"<<G4endl;
G4cout <<"Number of produced particles is very low:  " <<himain1_.natt<<G4endl;
G4cout <<"------------------------------------------------------------"<<G4endl;
G4cout <<"============================================================"<<G4endl;
      }
      }
   while (Nproduce<2);
// =============================================================================
         
  G4double BB=hiparnt_.hint1[18];      //impact parameter HINT1(19) 
//     cout<<"HIJING=====impact parameter= "<<BB<<endl;

  for (G4int i=0; i<Nproduce; i++)
  {
    


  G4int pid=himain2_.katt[0][i];

// Particle is a final state secondary and not a nucleus.
// Determine what this secondary particle is, and if valid, load dynamic
// parameters.
//
//   G4cout<<"pid================"<<pid<<G4endl;

  G4ParticleDefinition* pd=
  G4ParticleTable::GetParticleTable()->FindParticle(pid);
///////////////////////////////////////////////////////////////
//  exclude beam nucleons as produced particles
//     cout<<" himain2_.katt[1][i]== "<<himain2_.katt[1][i]<<endl;
//    if(himain2_.katt[1][i]==0 || himain2_.katt[1][i]==10) continue;
//  -----------------------------------------------------------
//      --------------reject neutral particles by calling luchge <new>
//         G4int chg_HIJ=luchge_ (&pid);
//        if (chg_HIJ==0) continue;

      if (pd)
      {

//units are in MeV/c for G4

        G4double px        = himain2_.patt[0][i]*GeV;
        G4double py        = himain2_.patt[1][i]*GeV;
        G4double pz        = himain2_.patt[2][i]*GeV;
        G4double et        = himain2_.patt[3][i]*GeV;


//    ------------------------------Use  "Lorentz vector"----------
G4LorentzVector lorenzCM = G4LorentzVector(px,py,pz,et);
//    Move to the lab frame
lorenzCM.boost(0.0,0.0,-1.0*betann);
G4LorentzVector lorenzLab = G4LorentzVector(lorenzCM.px(),lorenzCM.py(),
                            -1.0*lorenzCM.pz(),lorenzCM.e()); 
//-------------------------------------------------------------------
cascadeParticle = new G4DynamicParticle(pd, lorenzLab);   

theResult.AddSecondary(cascadeParticle);
  
}  //if pd

}  //for

#ifdef G4ANALYSIS_USE      //khaled new
fHistoManager->StoreSecondaries(BB, theResult);
#endif
//} //if warning

//


//=======================================================================
if (verbose >= 3) {

//
    G4double TotalEafter = 0.0;
    G4ThreeVector TotalPafter;
    G4double charge    = 0.0;
    G4int baryon        = 0;
    G4int nSecondaries  = theResult.GetNumberOfSecondaries();

    for (G4int j=0; j<nSecondaries; j++) {
      TotalEafter += theResult.GetSecondary(j)->
        GetParticle()->GetTotalEnergy()/GeV;

      TotalPafter += theResult.GetSecondary(j)->
        GetParticle()->GetMomentum()/GeV;

      G4ParticleDefinition *pd = theResult.GetSecondary(j)->
        GetParticle()->GetDefinition();

      charge += pd->GetPDGCharge();
      baryon += pd->GetBaryonNumber();

    }  //for secondaries

    G4cout <<"----------------------------------------"
          <<"----------------------------------------"
          <<G4endl;
    G4cout <<"Total energy before collision  = " <<TotalEbefore    ///GeV
          <<" GeV" <<G4endl;
    G4cout <<"Total energy after collision    = " <<TotalEafter/nSecondaries
        //GeV
          <<" GeV" <<G4endl;

    G4cout <<"----------------------------------------"<<G4endl;

    G4cout <<"Total momentum before collision = " <<TotalPbefore        
    //GeV
<<" GeV/c" <<G4endl;
    G4cout <<"Total momentum after collision  = " <<TotalPafter.z()/nSecondaries
    //GeV
    <<" GeV/c" <<G4endl;
    G4cout <<"----------------------------------------"<<G4endl;

    if (verbose >= 4) {
      G4cout <<"Total charge before collision  = " <<(ZP+ZT)      //
            <<G4endl;
      G4cout <<"Total charge after collision    = " <<charge
            <<G4endl;

      G4cout <<"----------------------------------------"<<G4endl;

      G4cout <<"Total baryon number before collision = "<<AP+AT
            <<G4endl;
      G4cout <<"Total baryon number after collision  = "<<baryon
            <<G4endl;
      G4cout <<"----------------------------------------"<<G4endl;

    }  //if verbose4

    G4cout <<"----------------------------------------"
          <<"----------------------------------------"
          <<G4endl;

  }  //if verbose3


return &theResult;
}  //G4hadfinal


//---------------------------------------------------------------------

//---------------------------------------------------------------------
//
// WelcomeMessage
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4HIJING_Model::WelcomeMessage () const
{
G4cout <<G4endl;
G4cout <<" *****************************************************************"
<<G4endl;
G4cout <<" Interface to        G4HIJING_Model                      activated"
<<G4endl;
G4cout <<" Version number : 01.00.0B          File date : 10/09/2013" <<G4endl;
G4cout <<"  Interface written by    Khaled Abdel-Waged              "
<<G4endl;   
G4cout <<"                       Umm Al-Qura University             "
<<G4endl;
G4cout <<"                         SAUDI ARABIA                     "
<<G4endl;
G4cout <<G4endl;
G4cout <<" *****************************************************************"
<<G4endl;
G4cout << G4endl;
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4HIJING_Model::InitialiseDataTables ()
{
//
//
// The next line is to make sure the block data statements are
// executed.
//

g4hijingblockdata_ ();
   
}

G4double G4HIJING_Model::Eplab(G4double m, G4double P)
{
G4double Eb= std::sqrt(P*P+m*m);
return Eb;
}
#endif

