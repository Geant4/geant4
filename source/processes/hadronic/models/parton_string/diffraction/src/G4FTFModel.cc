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
// $Id$
// GEANT4 tag $Name:  $
//

// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4FTFModel ----------------
//             by Gunter Folger, May 1998.
//       class implementing the excitation in the FTF Parton String Model
// ------------------------------------------------------------

#include <utility> 

#include "G4FTFModel.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4FTFParameters.hh"
#include "G4FTFParticipants.hh"
#include "G4DiffractiveSplitableHadron.hh"
#include "G4InteractionContent.hh"
#include "G4LorentzRotation.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

// Class G4FTFModel 

G4FTFModel::G4FTFModel(const G4String& modelName):G4VPartonStringModel(modelName),
                         theExcitation(new G4DiffractiveExcitation()),
                         theElastic(new G4ElasticHNScattering()),
                         theAnnihilation(new G4FTFAnnihilation())
{
	G4VPartonStringModel::SetThisPointer(this);
        theParameters=0;
	NumberOfInvolvedNucleon=0;
        NumberOfInvolvedNucleonOfProjectile=0;
    SetEnergyMomentumCheckLevels(2*perCent, 150*MeV);
}

struct DeleteVSplitableHadron { void operator()(G4VSplitableHadron * aH){ delete aH;} };

G4FTFModel::~G4FTFModel()
{
// Because FTF model can be called for various particles
// theParameters must be erased at the end of each call.
// Thus the delete is also in G4FTFModel::GetStrings() method
   if( theParameters   != 0 ) delete theParameters; 
   if( theExcitation   != 0 ) delete theExcitation;
   if( theElastic      != 0 ) delete theElastic; 
   if( theAnnihilation != 0 ) delete theAnnihilation;

// Erasing of strings created at annihilation
   if(theAdditionalString.size() != 0)
   {
    std::for_each(theAdditionalString.begin(), theAdditionalString.end(), 
                  DeleteVSplitableHadron());
   }
   theAdditionalString.clear();

// Erasing of target involved nucleons
   if( NumberOfInvolvedNucleon != 0)
   {
    for(G4int i=0; i < NumberOfInvolvedNucleon; i++)
    {
     G4VSplitableHadron * aNucleon = TheInvolvedNucleon[i]->GetSplitableHadron();
     if(aNucleon) delete aNucleon;
    }
   }

// Erasing of projectile involved nucleons
/*
   if( NumberOfInvolvedNucleonOfProjectile != 0)
   {
    for(G4int i=0; i < NumberOfInvolvedNucleonOfProjectile; i++)
    {
     G4VSplitableHadron * aNucleon = TheInvolvedNucleonOfProjectile[i]->GetSplitableHadron();
     if(aNucleon) delete aNucleon;
    }
   }
*/
}

// ------------------------------------------------------------
void G4FTFModel::Init(const G4Nucleus & aNucleus, const G4DynamicParticle & aProjectile)
{
	theProjectile = aProjectile;  

        G4double PlabPerParticle(0.);  // Laboratory momentum Pz per particle/nucleon

/*
G4cout<<"FTF init Pro Name "<<theProjectile.GetDefinition()->GetParticleName()<<G4endl;
G4cout<<"FTF init Pro Mass "<<theProjectile.GetMass()<<" "<<theProjectile.GetMomentum()<<G4endl;
G4cout<<"FTF init Pro B Q  "<<theProjectile.GetDefinition()->GetBaryonNumber()<<" "<<(G4int) theProjectile.GetDefinition()->GetPDGCharge()<<G4endl; 
G4cout<<"FTF init A Z "<<aNucleus.GetA_asInt()<<" "<<aNucleus.GetZ_asInt()<<G4endl;
G4cout<<"             "<<aNucleus.GetN()<<" "<<aNucleus.GetZ()<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
*/

        theParticipants.SetProjectileNucleus(0);
        theParticipants.Init(aNucleus.GetA_asInt(),aNucleus.GetZ_asInt());

        if(std::abs(theProjectile.GetDefinition()->GetBaryonNumber()) <= 1) 
        { // Projectile is a hadron
         PlabPerParticle=theProjectile.GetMomentum().z();

//         S = sqr( theProjectile.GetMass() ) + sqr( ProtonMass ) +
//                 2*ProtonMass*theProjectile.GetTotalEnergy();
        }


        if(theProjectile.GetDefinition()->GetBaryonNumber() > 1) 
        { // Projectile is a nucleus
         theParticipants.InitProjectileNucleus(
                      theProjectile.GetDefinition()->GetBaryonNumber(),
              (G4int) theProjectile.GetDefinition()->GetPDGCharge()    );

         G4ThreeVector BoostVector=theProjectile.GetMomentum()/theProjectile.GetTotalEnergy();
         theParticipants.theProjectileNucleus->DoLorentzBoost(BoostVector);

         PlabPerParticle=theProjectile.GetMomentum().z()/
                         theProjectile.GetDefinition()->GetBaryonNumber();

//         S =         2.*sqr( ProtonMass ) + 2*ProtonMass*
//             theProjectile.GetTotalEnergy()/theProjectile.GetDefinition()->GetBaryonNumber();
        }

        if(theProjectile.GetDefinition()->GetBaryonNumber() < -1) 
        { // Projectile is an anti-nucleus
         theParticipants.InitProjectileNucleus(
             std::abs(        theProjectile.GetDefinition()->GetBaryonNumber()),
             std::abs((G4int) theProjectile.GetDefinition()->GetPDGCharge())    );

         G4ThreeVector BoostVector=theProjectile.GetMomentum()/theProjectile.GetTotalEnergy();

         theParticipants.theProjectileNucleus->StartLoop();
         G4Nucleon * aNucleon;
         while ( (aNucleon = theParticipants.theProjectileNucleus->GetNextNucleon()) )
         {
          if(aNucleon->GetDefinition()->GetPDGEncoding() == 2212) // Proton
          {aNucleon->SetParticleType(G4AntiProton::AntiProton());} 

          if(aNucleon->GetDefinition()->GetPDGEncoding() == 2112) // Neutron
          {aNucleon->SetParticleType(G4AntiNeutron::AntiNeutron());} 
         }   // end of while (theParticipant.theProjectileNucleus->GetNextNucleon())

         theParticipants.theProjectileNucleus->DoLorentzBoost(BoostVector);

         PlabPerParticle=         theProjectile.GetMomentum().z()/
                         std::abs(theProjectile.GetDefinition()->GetBaryonNumber());

//         S =        2.*sqr( ProtonMass ) + 2*ProtonMass*
//                      theProjectile.GetTotalEnergy()/
//             std::abs(theProjectile.GetDefinition()->GetBaryonNumber());
        }

// ------------------------------------------------------------------------
      if( theParameters != 0 ) delete theParameters;
      theParameters = new G4FTFParameters(theProjectile.GetDefinition(),
                                          aNucleus.GetA_asInt(),aNucleus.GetZ_asInt(),
                                          PlabPerParticle);
//G4cout<<" End Init "<<theProjectile.GetMomentum()<<G4endl;
// To turn on/off (1/0) elastic scattering close/open ...
//theParameters->SetProbabilityOfElasticScatt(0.); 
//G4cout<<" etProbabilityOfElasticScatt "<<theParameters->GetProbabilityOfElasticScatt()<<G4endl;
//G4cout<<" INIT ";
//G4int Uzhi; G4cin>>Uzhi;

   if(theAdditionalString.size() != 0)
   {
    std::for_each(theAdditionalString.begin(), theAdditionalString.end(), 
                  DeleteVSplitableHadron());
   }
   theAdditionalString.clear();
//G4cout<<" End Init theProjectile.GetMomentum()"<<theProjectile.GetMomentum()<<G4endl;
}

// ------------------------------------------------------------
G4ExcitedStringVector * G4FTFModel::GetStrings()
{ 
        G4ExcitedStringVector * theStrings(0);

	theParticipants.GetList(theProjectile,theParameters);
//        StoreInvolvedNucleon();
//G4cout<<" GetList theProjectile.GetMomentum()  GetBaryonNumber() "<<theProjectile.GetMomentum()<<" "<<theProjectile.GetDefinition()->GetBaryonNumber()<<G4endl;
        G4bool Success(true);

        if((std::abs(theProjectile.GetDefinition()->GetBaryonNumber()) <= 1) &&
                    (theProjectile.GetDefinition()->GetBaryonNumber() != -1)   )
        { // Standard variant of FTF for projectile hadron/nucleon
//G4cout<<"Standard variant of FTF for projectile hadron/nucleon"<<G4endl;
         ReggeonCascade(); 
//G4cout<<"Success after Reggeon "<<Success<<" PutOnMasShell"<<G4endl;
         Success=PutOnMassShell(); 
//G4cout<<"Success after PutOn "<<Success<<" GetResid"<<G4endl;
         GetResidualNucleus();
        } 
//G4cout<<"Success after GetN "<<Success<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
        if(theProjectile.GetDefinition()->GetBaryonNumber() > 1)
        { // Variant of FTF for projectile nuclei
//G4cout<<"Variant of FTF for projectile nuclei"<<G4endl;
         StoreInvolvedNucleon();
         ReggeonCascade(); 
         Success=PutOnMassShell(); 
         GetResidualNucleus();
        } 

//        G4bool LowE_Anti_Ion(false);
        if(theProjectile.GetDefinition()->GetBaryonNumber() <= -1) 
        { // Projectile is Anti-baryon or Anti-Nucleus
//G4cout<<"Projectile is Anti-baryon or Anti-Nucleus "<<G4endl;
//G4cout<<"Be4 Store"<<G4endl;
         StoreInvolvedNucleon();
         if(theProjectile.GetTotalMomentum()/
            std::abs(theProjectile.GetDefinition()->GetBaryonNumber()) > 5000.*MeV)
         {// High energy interaction
//G4cout<<"High energy interaction "<<G4endl;
//G4cout<<"Regeon "<<G4endl;
          ReggeonCascade(); 
//G4cout<<"Put on mass "<<G4endl;
          Success=PutOnMassShell(); 
//G4cout<<"Residual "<<G4endl;
          GetResidualNucleus();
         }
         else
         {
//G4cout<<"Low energy interaction "<<G4endl;
//          LowE_Anti_Ion=true;
          Success=true;
         }
        }
//G4cout<<"Before Excite Success "<<Success<<G4endl;
        Success=Success && ExciteParticipants();
//G4cout<<"Success ExciteParticipants()? "<<Success<<G4endl;
//        if(LowE_Anti_Ion) Success=Success && GetResidualNucleusAfterAnnihilation();

        if( Success )
        {       
	  theStrings = BuildStrings();
//G4cout<<"BuildStrings OK"<<G4endl;
          if( theParameters != 0 )
          {
           delete theParameters;
           theParameters=0;
          }
         }
/*
        if( Success )
        {
         if( ExciteParticipants() )
         {
//G4cout<<"Excite partic OK"<<G4endl;
	  theStrings = BuildStrings();
//G4cout<<"Build String OK"<<G4endl;
          if(LowE_Anti_Ion) GetResidualNucleusAfterAnnihilation();

          if( theParameters != 0 )
          {
           delete theParameters;
           theParameters=0;
          }
         } else                      // if( ExciteParticipants() )
         {     Success=false;}
        } else                       // if( Success )
        {      Success=false;}
*/
        if(!Success)   
        {
           // -------------- Erase the projectile ----------------
//G4cout<<"Erase Proj"<<G4endl;
	 std::vector<G4VSplitableHadron *> primaries;

	 theParticipants.StartLoop();    // restart a loop 
         while ( theParticipants.Next() ) 
	 {
	    const G4InteractionContent & interaction=theParticipants.GetInteraction();

                	 //  do not allow for duplicates ...
	    if ( primaries.end() == std::find(primaries.begin(), primaries.end(),
                                                   interaction.GetProjectile()) )
	    	primaries.push_back(interaction.GetProjectile());                
         }
         std::for_each(primaries.begin(), primaries.end(), DeleteVSplitableHadron());
         primaries.clear();
        }
// -------------- Cleaning of the memory --------------
        G4VSplitableHadron * aNucleon = 0;
// -------------- Erase the projectile nucleon --------
/*
        G4VSplitableHadron * aNucleon = 0;
        for(G4int i=0; i < NumberOfInvolvedNucleonOfProjectile; i++)
        {
           aNucleon = TheInvolvedNucleonOfProjectile[i]->GetSplitableHadron();
           if(aNucleon) delete aNucleon;
        } 

        NumberOfInvolvedNucleonOfProjectile=0;
*/  // Maybe it will be needed latter------------------

// -------------- Erase the target nucleons -----------
//G4cout<<"Erase Target Ninv "<<NumberOfInvolvedNucleon<<G4endl;
        for(G4int i=0; i < NumberOfInvolvedNucleon; i++)
        {
           aNucleon = TheInvolvedNucleon[i]->GetSplitableHadron();
           if(aNucleon) delete aNucleon;
        } 

        NumberOfInvolvedNucleon=0;
//G4cout<<"Go to fragmentation"<<G4endl;
        return theStrings;

}

//-------------------------------------------------------------------
void G4FTFModel::StoreInvolvedNucleon()                             
{ //--- To store nucleons involved in low energy interaction  -------
if(std::abs(theProjectile.GetDefinition()->GetBaryonNumber()) <= 1)
{ // the projectile is a hadron -----------
//G4cout<<"the projectile is a hadron"<<G4endl;
        NumberOfInvolvedNucleon=0;

        theParticipants.StartLoop();

	while (theParticipants.Next())
	{   
//G4cout<<"theParticipants.Next()"<<G4endl;
	   const G4InteractionContent & collision=theParticipants.GetInteraction();
//G4cout<<"collision=theParticipants.GetInteraction()"<<G4endl;
           G4Nucleon * TargetNucleon=collision.GetTargetNucleon();
//G4cout<<"TargetNucleon=collision.GetTargetNucleon()"<<G4endl;

           TheInvolvedNucleon[NumberOfInvolvedNucleon]=TargetNucleon;
           NumberOfInvolvedNucleon++;
//G4cout<<G4endl<<"Prim NumberOfInvolvedNucleon "<<NumberOfInvolvedNucleon<<G4endl;
	}      // end of while (theParticipants.Next())

        NumberOfInvolvedTargetNucleon=NumberOfInvolvedNucleon;
// ---------------- Calculation of creation time for each target nucleon -----------
//G4cout<<"theParticipants.StartLoop() "<<G4endl;
	theParticipants.StartLoop();    // restart a loop
//G4cout<<"theParticipants.Next();"<<G4endl;
        theParticipants.Next();
	G4VSplitableHadron * primary = theParticipants.GetInteraction().GetProjectile();
//G4cout<<"primary->Get4Momentum() "<<primary->Get4Momentum()<<G4endl;
//G4cout<<"primary->Get4Momentum().pz() "<<primary->Get4Momentum().pz()<<G4endl;
//G4cout<<"primary->Get4Momentum().e() "<<primary->Get4Momentum().e()<<G4endl;

        G4double betta_z=primary->Get4Momentum().pz()/primary->Get4Momentum().e();
//G4cout<<"betta_z "<<betta_z<<G4endl;
        primary->SetTimeOfCreation(0.);

        G4double ZcoordinateOfPreviousCollision(0.);
        G4double ZcoordinateOfCurrentInteraction(0.);
        G4double TimeOfPreviousCollision(0.);
        G4double TimeOfCurrentCollision(0);

        theParticipants.theNucleus->StartLoop();
        G4Nucleon * aNucleon;
        G4bool theFirstInvolvedNucleon(true);
	while ( (aNucleon = theParticipants.theNucleus->GetNextNucleon()) )
        {
//G4cout<<"aNucleon->AreYouHit() "<<aNucleon->AreYouHit()<<G4endl;
          if(aNucleon->AreYouHit())
          {
            if(theFirstInvolvedNucleon)
            {
              ZcoordinateOfPreviousCollision=aNucleon->GetPosition().z();
//G4cout<<"ZcoordinateOfPreviousCollision "<<ZcoordinateOfPreviousCollision/fermi<<G4endl;
              theFirstInvolvedNucleon=false;
            }

            ZcoordinateOfCurrentInteraction=aNucleon->GetPosition().z();
//G4cout<<"ZcoordinateOfCurrentInteraction "<<ZcoordinateOfCurrentInteraction/fermi<<G4endl;
//G4cout<<"TimeOfPreviousCollision "<<TimeOfPreviousCollision<<G4endl;

            // A.R. 18-Oct-2011 : Protection needed for nuclear capture of
            //                    anti-proton at rest.
	    if ( betta_z > 1.0e-10 ) {
              TimeOfCurrentCollision=TimeOfPreviousCollision+ 
              (ZcoordinateOfCurrentInteraction-ZcoordinateOfPreviousCollision)/betta_z;
            } else {
              TimeOfCurrentCollision=TimeOfPreviousCollision;
            } 

//G4cout<<"TimeOfCurrentCollision "<<TimeOfCurrentCollision<<G4endl;
// It is assumed that the nucleons are ordered on increasing z-coordinate ------------
            aNucleon->GetSplitableHadron()->SetTimeOfCreation(TimeOfCurrentCollision);

            ZcoordinateOfPreviousCollision=ZcoordinateOfCurrentInteraction;
            TimeOfPreviousCollision=TimeOfCurrentCollision;
          }  // end of if(aNucleon->AreYouHit())
	}   // end of while (theParticipant.theNucleus->GetNextNucleon())
//
        return;
} // end of if(std::abs(theProjectile.GetDefinition()->GetBaryonNumber()) <= 1)

// The projectile is a nucleus or an anti-nucleus
//G4cout<<"FTF The projectile is a nucleus or an anti-nucleus"<<G4endl;
        NumberOfInvolvedNucleonOfProjectile=0;

        G4V3DNucleus * ProjectileNucleus =theParticipants.GetProjectileNucleus();
	ProjectileNucleus->StartLoop();

        G4Nucleon *    ProjectileNucleon;
        while ( (ProjectileNucleon=ProjectileNucleus->GetNextNucleon()) )
        {
         if ( ProjectileNucleon->AreYouHit() )
         {  // Projectile nucleon was involved in the interaction.
           TheInvolvedNucleonOfProjectile[NumberOfInvolvedNucleonOfProjectile]=
                                    ProjectileNucleon;
           NumberOfInvolvedNucleonOfProjectile++;
         }
        } // End of while ( (ProjectileNucleon=ProjectileNucleus->GetNextNucleon()) )

//------------------
        NumberOfInvolvedNucleon=0;

        G4V3DNucleus * TargetNucleus =theParticipants.GetWoundedNucleus();
	TargetNucleus->StartLoop();

        G4Nucleon *    TargetNucleon;
        while ( (TargetNucleon=TargetNucleus->GetNextNucleon()) )
        {
         if ( TargetNucleon->AreYouHit() )
         {  // Target nucleon was involved in the interaction.
           TheInvolvedNucleon[NumberOfInvolvedNucleon]=TargetNucleon;
           NumberOfInvolvedNucleon++;
         }
        } // End of while ( (TargetNucleon=TargetNucleus->GetNextNucleon()) )
//G4cout<<"Store inv "<<NumberOfInvolvedNucleonOfProjectile<<" "<<NumberOfInvolvedNucleon<<G4endl;

        NumberOfInvolvedTargetNucleon=NumberOfInvolvedNucleon;
        return;
}                                                             // Uzhi 10 Feb. 2011

//-------------------------------------------------------------------
void G4FTFModel::ReggeonCascade()                             
{ //--- Implementation of the reggeon theory inspired model-------
//G4cout<<"In reggeon"<<G4endl;

        if(std::abs(theProjectile.GetDefinition()->GetBaryonNumber()) > 1) return;
//      For Nucleus-nucleus or Anti-nucleus - nucleus interactions 
//      the cascading will be implemented latter.

        NumberOfInvolvedNucleon=0;

        theParticipants.StartLoop();
	while (theParticipants.Next())
	{   
	   const G4InteractionContent & collision=theParticipants.GetInteraction();

           G4Nucleon * TargetNucleon=collision.GetTargetNucleon();

           TheInvolvedNucleon[NumberOfInvolvedNucleon]=TargetNucleon;
           NumberOfInvolvedNucleon++;
//G4cout<<"Prim NumberOfInvolvedNucleon "<<NumberOfInvolvedNucleon<<G4endl;
           G4double XofWoundedNucleon = TargetNucleon->GetPosition().x();
           G4double YofWoundedNucleon = TargetNucleon->GetPosition().y();

           theParticipants.theNucleus->StartLoop();
           G4Nucleon * Neighbour(0);

	   while ( (Neighbour = theParticipants.theNucleus->GetNextNucleon()) )
           {
            if(!Neighbour->AreYouHit())
            {
    	     G4double impact2= sqr(XofWoundedNucleon - Neighbour->GetPosition().x()) +
                               sqr(YofWoundedNucleon - Neighbour->GetPosition().y());

             if(G4UniformRand() < theParameters->GetCofNuclearDestruction()*
                std::exp(-impact2/theParameters->GetR2ofNuclearDestruction()))  
             { // The neighbour nucleon is involved in the reggeon cascade

              TheInvolvedNucleon[NumberOfInvolvedNucleon]=Neighbour;
              NumberOfInvolvedNucleon++;
//G4cout<<"Seco NumberOfInvolvedNucleon "<<NumberOfInvolvedNucleon<<G4endl;

              G4VSplitableHadron *targetSplitable; 
              targetSplitable = new G4DiffractiveSplitableHadron(*Neighbour); 

              Neighbour->Hit(targetSplitable);
              targetSplitable->SetStatus(2);     
             }
            }  // end of if(!Neighbour->AreYouHit())
	   }   // end of while (theParticipant.theNucleus->GetNextNucleon())
	}      // end of while (theParticipants.Next())

        NumberOfInvolvedTargetNucleon=NumberOfInvolvedNucleon;

// ---------------- Calculation of creation time for each target nucleon -----------
	theParticipants.StartLoop();    // restart a loop
        theParticipants.Next();
	G4VSplitableHadron * primary = theParticipants.GetInteraction().GetProjectile();
        G4double betta_z=primary->Get4Momentum().pz()/primary->Get4Momentum().e();
        primary->SetTimeOfCreation(0.);

        G4double ZcoordinateOfPreviousCollision(0.);
        G4double ZcoordinateOfCurrentInteraction(0.);
        G4double TimeOfPreviousCollision(0.);
        G4double TimeOfCurrentCollision(0);

        theParticipants.theNucleus->StartLoop();
        G4Nucleon * aNucleon;
        G4bool theFirstInvolvedNucleon(true);
	while ( (aNucleon = theParticipants.theNucleus->GetNextNucleon()) )
        {
          if(aNucleon->AreYouHit())
          {
            if(theFirstInvolvedNucleon)
            {
              ZcoordinateOfPreviousCollision=aNucleon->GetPosition().z();
              theFirstInvolvedNucleon=false;
            }

            ZcoordinateOfCurrentInteraction=aNucleon->GetPosition().z();
            TimeOfCurrentCollision=TimeOfPreviousCollision+ 
            (ZcoordinateOfCurrentInteraction-ZcoordinateOfPreviousCollision)/betta_z; 
// It is assumed that the nucleons are ordered on increasing z-coordinate ------------
            aNucleon->GetSplitableHadron()->SetTimeOfCreation(TimeOfCurrentCollision);

            ZcoordinateOfPreviousCollision=ZcoordinateOfCurrentInteraction;
            TimeOfPreviousCollision=TimeOfCurrentCollision;
          }  // end of if(aNucleon->AreYouHit())
	}   // end of while (theParticipant.theNucleus->GetNextNucleon())
//
// The algorithm can be improved, but it will be more complicated, and will require
// changes in G4DiffractiveExcitation.cc and G4ElasticHNScattering.cc
}                                                             // Uzhi 26 July 2009

// ------------------------------------------------------------
G4bool G4FTFModel::PutOnMassShell()
{
//G4cout<<"PutOnMassShell start "<<G4endl;
        if(std::abs(theProjectile.GetDefinition()->GetBaryonNumber()) > 1) // !!!!
        { // The projectile is a nucleus or an anti-nucleus
//G4cout<<"PutOnMassShell AA "<<G4endl;
         G4LorentzVector P_total(0.,0.,0.,0.);

         G4LorentzVector P_participants(0.,0.,0.,0.);
         G4int MultiplicityOfParticipants(0);

         G4V3DNucleus * ProjectileNucleus =theParticipants.GetProjectileNucleus();
	 ProjectileNucleus->StartLoop();

         G4Nucleon *    ProjectileNucleon;
         while ( (ProjectileNucleon=ProjectileNucleus->GetNextNucleon()) )
         {
          if ( ProjectileNucleon->AreYouHit() )
          {  // Projectile nucleon was involved in the interaction.
           P_total+=ProjectileNucleon->Get4Momentum();
           MultiplicityOfParticipants++;
           P_participants+=ProjectileNucleon->Get4Momentum();
          }
         } // End of while ( (ProjectileNucleon=ProjectileNucleus->GetNextNucleon()) )
//G4cout<<"MultParts mom Pr "<<MultiplicityOfParticipants<<" "<<P_participants<<G4endl;
//------------------
         G4int ResidualMassNumber(0);
         G4int ResidualCharge(0);
         ResidualExcitationEnergy=0.;
         G4LorentzVector PnuclearResidual(0.,0.,0.,0.);

         G4double ExcitationEnergyPerWoundedNucleon=
                  theParameters->GetExcitationEnergyPerWoundedNucleon();

         G4V3DNucleus * TargetNucleus =theParticipants.GetWoundedNucleus();
	 TargetNucleus->StartLoop();

         G4Nucleon *    TargetNucleon;
         while ( (TargetNucleon=TargetNucleus->GetNextNucleon()) )
         {
          P_total+=TargetNucleon->Get4Momentum();
          if ( TargetNucleon->AreYouHit() )
          {  // Target nucleon was involved in the interaction.
           MultiplicityOfParticipants++;
           P_participants+=TargetNucleon->Get4Momentum();
           ResidualExcitationEnergy+=ExcitationEnergyPerWoundedNucleon; 
          }
          else
          {
           ResidualMassNumber+=1;
           ResidualCharge+=(G4int) TargetNucleon->GetDefinition()->GetPDGCharge();
           PnuclearResidual+=TargetNucleon->Get4Momentum();
          }
         } // End of while ( (TargetNucleon=TargetNucleus->GetNextNucleon()) )

         G4double ResidualMass(0.);
         if(ResidualMassNumber == 0)
         {
          ResidualMass=0.;
          ResidualExcitationEnergy=0.;
         }
         else
         {
          ResidualMass=G4ParticleTable::GetParticleTable()->GetIonTable()->
                              GetIonMass(ResidualCharge ,ResidualMassNumber);
          if(ResidualMassNumber == 1) {ResidualExcitationEnergy=0.;}
         }
//      ResidualMass +=ResidualExcitationEnergy; // Will be given after checks

//G4cout<<"MultPars mom+Tr  "<<MultiplicityOfParticipants<<" "<<P_participants<<G4endl;
//G4cout<<"Res              "<<ResidualMassNumber<<" "<<PnuclearResidual<<G4endl;
//G4cout<<"P_total "<<P_total<<G4endl;

//G4cout<<"Rez A Z M E* "<<ResidualMassNumber<<" "<<ResidualCharge<<" "<<ResidualMass<<" "<<ResidualExcitationEnergy<<G4endl;
//--------------
         G4double SqrtS=P_total.mag();
         G4double     S=P_total.mag2();

//         if(theProjectile.GetDefinition()->GetBaryonNumber() > 1)
         { // For nucleus-nucleus interactions
          G4double MassOfParticipants=P_participants.mag();
          G4double MassOfParticipants2=sqr(MassOfParticipants);

//G4cout<<"ResidualMass + MassOfParticipants "<<ResidualMass + MassOfParticipants<<G4endl;
          if(SqrtS < ResidualMass + MassOfParticipants) {return false;}

          if(SqrtS < ResidualMass+ResidualExcitationEnergy + MassOfParticipants)
          {ResidualExcitationEnergy=0.;}

          ResidualMass +=ResidualExcitationEnergy; 
//G4cout<<"Rez A Z M E* "<<ResidualMassNumber<<" "<<ResidualCharge<<" "<<ResidualMass<<" "<<ResidualExcitationEnergy<<G4endl;

          G4double ResidualMass2=sqr(ResidualMass);

          G4LorentzRotation toCms(-1*P_total.boostVector());

          G4LorentzVector Pcms=toCms*P_participants;
//G4cout<<"Ppart in CMS "<<Ptmp<<G4endl;

          if ( Pcms.pz() <= 0. )                                
          {  // "String" moving backwards in  CMS, abort collision !!
           //G4cout << " abort ColliDeleteVSplitableHadronsion!! " << G4endl;
           return false; 
          }

          toCms.rotateZ(-1*Pcms.phi());              // Uzhi 5.12.09
          toCms.rotateY(-1*Pcms.theta());            // Uzhi 5.12.09

//G4cout<<"Mpa M res "<<ResidualMass<<" "<<MassOfParticipants<<" "<<SqrtS<<G4endl;
          G4LorentzRotation toLab(toCms.inverse());

          G4double DecayMomentum2= 
                      sqr(S)+sqr(MassOfParticipants2)+ sqr(ResidualMass2) -
                          2.*S*MassOfParticipants2 - 2.*S*ResidualMass2 
                              -2.*MassOfParticipants2*ResidualMass2;

          if(DecayMomentum2 < 0.) return false;

          DecayMomentum2/=(4.*S);
          G4double DecayMomentum = std::sqrt(DecayMomentum2);
//G4cout<<"DecayMomentum "<<DecayMomentum<<G4endl;

          G4LorentzVector New_P_participants(0.,0.,DecayMomentum,
                                             std::sqrt(DecayMomentum2+MassOfParticipants2));
          G4LorentzVector New_PnuclearResidual(0.,0.,-DecayMomentum,
                                             std::sqrt(DecayMomentum2+ResidualMass2));

//G4cout<<"MultPars mom     "<<MultiplicityOfParticipants<<" "<<New_P_participants<<G4endl;
//G4cout<<"Res              "<<ResidualMassNumber<<" "<<New_PnuclearResidual<<G4endl;
          New_P_participants.transform(toLab);
          New_PnuclearResidual.transform(toLab);

//G4cout<<"MultPars mom     "<<MultiplicityOfParticipants<<" "<<New_P_participants<<G4endl;
//G4cout<<"Res              "<<ResidualMassNumber<<" "<<New_PnuclearResidual<<G4endl;

          G4LorentzVector DeltaP_participants=(New_P_participants - P_participants)/
                                              ((G4double) MultiplicityOfParticipants);

//G4cout<<"DeltaP_participants "<<DeltaP_participants<<G4endl;
//-------------
          ProjectileNucleus->StartLoop();
          while ( (ProjectileNucleon=ProjectileNucleus->GetNextNucleon()) )
          {
           if ( ProjectileNucleon->AreYouHit() )
           {  // Projectile nucleon was involved in the interaction.
            G4LorentzVector Ptmp=ProjectileNucleon->Get4Momentum() + DeltaP_participants;
            ProjectileNucleon->GetSplitableHadron()->Set4Momentum(Ptmp);
            ProjectileNucleon->SetMomentum(Ptmp);
           }
          } // End of while ( (ProjectileNucleon=ProjectileNucleus->GetNextNucleon()) )

//-------------
          TargetNucleus->StartLoop();
          while ( (TargetNucleon=TargetNucleus->GetNextNucleon()) )
          {
           if ( TargetNucleon->AreYouHit() )
           {  // Target nucleon was involved in the interaction.
            G4LorentzVector Ptmp=TargetNucleon->Get4Momentum() + DeltaP_participants;
            TargetNucleon->GetSplitableHadron()->Set4Momentum(Ptmp);
           }
          } // End of while ( (TargetNucleon=TargetNucleus->GetNextNucleon()) )

          Residual4Momentum = New_PnuclearResidual;        
//          return true;
         } // End of if(theProjectile.GetDefinition()->GetBaryonNumber() > 1)

         return true;
        }
//---------------------------------------------------------------------
// -------- The projectile is hadron, or baryon, or anti-baryon -------
// -------------- Properties of the projectile ------------------------
	theParticipants.StartLoop();    // restart a loop
        theParticipants.Next();
	G4VSplitableHadron * primary = theParticipants.GetInteraction().GetProjectile();
	G4LorentzVector Pprojectile=primary->Get4Momentum();

// 13.06.2012 G4bool ProjectileIsAntiBaryon = primary->GetDefinition()->GetBaryonNumber() < 0;

//G4cout<<"PutOnMass Pprojectile "<<Pprojectile<<G4endl;
// To get original projectile particle

        if(Pprojectile.z() < 0.){return false;}

        G4double Mprojectile  = Pprojectile.mag();
        G4double M2projectile = Pprojectile.mag2();
//-------------------------------------------------------------
	G4LorentzVector Psum      = Pprojectile;

        G4double        SumMasses = Mprojectile + 20.*MeV; // 13.12.09
                                               // Separation energy for projectile
//        if(ProjectileIsAntiBaryon) {SumMasses = Mprojectile;}
//G4cout<<"SumMasses Pr "<<SumMasses<<G4endl;
//--------------- Target nucleus ------------------------------
        G4V3DNucleus *theNucleus = GetWoundedNucleus();
        G4int ResidualMassNumber=theNucleus->GetMassNumber();
        G4int ResidualCharge    =theNucleus->GetCharge();

        ResidualExcitationEnergy=0.;
	G4LorentzVector Ptarget(0.,0.,0.,0.);
	G4LorentzVector PnuclearResidual(0.,0.,0.,0.); // Uzhi 12.06.2012

        G4double ExcitationEnergyPerWoundedNucleon=
                  theParameters->GetExcitationEnergyPerWoundedNucleon();

//G4cout<<"ExcitationEnergyPerWoundedNucleon "<<ExcitationEnergyPerWoundedNucleon<<G4endl;

        theNucleus->StartLoop();

	while (G4Nucleon * aNucleon = theNucleus->GetNextNucleon())
        {
         Ptarget+=aNucleon->Get4Momentum();

         if(aNucleon->AreYouHit())
         {   // Involved nucleons
//G4cout<<"PutOn Tr "<<aNucleon->Get4Momentum()<<G4endl;
//        Psum += aNucleon->Get4Momentum();                       // Uzhi 20 Sept.
//          if(!ProjectileIsAntiBaryon)                           // Uzhi 13.06.2012
//          {
           SumMasses += std::sqrt(sqr(aNucleon->GetDefinition()->GetPDGMass()) //Uzhi 12.06.2012
                     +  aNucleon->Get4Momentum().perp2());                     //Uzhi 12.06.2012
           SumMasses += 20.*MeV;   // 13.12.09 Separation energy for a nucleon
           ResidualExcitationEnergy+=ExcitationEnergyPerWoundedNucleon;
/*                                                                             //Uzhi 13.06.2012
          } else 
          {
           SumMasses += aNucleon->Get4Momentum().mag();           // 4.12.2010
           G4LorentzVector tmp=aNucleon->Get4Momentum();
           tmp.setE(aNucleon->Get4Momentum().mag());   // It is need to save mass 6.12.2011
           aNucleon->SetMomentum(tmp);
          }
*/                                                                            //Uzhi 13.06.2012

          ResidualMassNumber--;
          ResidualCharge-=(G4int) aNucleon->GetDefinition()->GetPDGCharge();
         }
         else
         {   // Spectator nucleons
          PnuclearResidual += aNucleon->Get4Momentum();          // Uzhi 12.06.2012
         }  // end of if(!aNucleon->AreYouHit())
	}   // end of while (theNucleus->GetNextNucleon())

        Psum += Ptarget;   
        PnuclearResidual.setPz(0.); PnuclearResidual.setE(0.);   // Uzhi 12.06.2012
//G4cout<<"ResidualCharge ,ResidualMassNumber "<<ResidualCharge<<" "<<ResidualMassNumber<<" "<<Ptarget<<" "<<PnuclearResidual<<G4endl;

        G4double ResidualMass(0.);
        if(ResidualMassNumber == 0)
        {
         ResidualMass=0.;
         ResidualExcitationEnergy=0.;
        }
        else
        {
         ResidualMass=G4ParticleTable::GetParticleTable()->GetIonTable()->
                              GetIonMass(ResidualCharge ,ResidualMassNumber);
         if(ResidualMassNumber == 1) {ResidualExcitationEnergy=0.;}
        }
 
//      ResidualMass +=ResidualExcitationEnergy; // Will be given after checks
//G4cout<<"SumMasses And ResidualMass "<<SumMasses<<" "<<ResidualMass<<G4endl;
        SumMasses += std::sqrt(sqr(ResidualMass)+PnuclearResidual.perp2()); // Uzhi 16.05.2013
//G4cout<<"SumMasses + ResM "<<SumMasses<<G4endl;
//G4cout<<"Psum "<<Psum<<G4endl;
//-------------------------------------------------------------

        G4double SqrtS=Psum.mag();
        G4double     S=Psum.mag2();

//G4cout<<"SqrtS < SumMasses "<<SqrtS<<" "<<SumMasses<<G4endl;

        if(SqrtS < SumMasses)      {return false;} // It is impossible to simulate
                                                   // after putting nuclear nucleons
                                                   // on mass-shell

        SumMasses -= std::sqrt(sqr(ResidualMass)+PnuclearResidual.perp2()); // Uzhi 16.05.2013
        SumMasses += std::sqrt(sqr(ResidualMass+ResidualExcitationEnergy)
                    +PnuclearResidual.perp2());                             // Uzhi 16.05.2013
        if(SqrtS < SumMasses)                                              // Uzhi 12.06.2012
        {
         SumMasses -= std::sqrt(sqr(ResidualMass+ResidualExcitationEnergy)
                     +PnuclearResidual.perp2());                            // Uzhi 16.05.2013
         SumMasses += std::sqrt(sqr(ResidualMass)+PnuclearResidual.perp2());// Uzhi 16.05.2013
         ResidualExcitationEnergy=0.;
        }

        ResidualMass +=ResidualExcitationEnergy;
//      SumMasses    +=ResidualExcitationEnergy;                           // Uzhi 12.06.2012
//G4cout<<"ResidualMass SumMasses ResidualExcitationEnergy "<<ResidualMass<<" "<<SumMasses<<" "<<ResidualExcitationEnergy<<G4endl;
//-------------------------------------------------------------
// Sampling of nucleons what are transfered to delta-isobars --
        G4int MaxNumberOfDeltas = (int)((SqrtS - SumMasses)/(400.*MeV));
        G4int NumberOfDeltas(0);
//G4cout<<"MaxNumberOfDeltas "<<MaxNumberOfDeltas<<G4endl;
        if(theNucleus->GetMassNumber() != 1)
        {
//G4cout<<"NumberOfInvolvedNucleon "<<NumberOfInvolvedNucleon<<G4endl;
          G4double ProbDeltaIsobar(0.05);                                  // Uzhi 6.07.2012
	  for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
          {
            if((G4UniformRand() < ProbDeltaIsobar)&&(NumberOfDeltas < MaxNumberOfDeltas))
            {
              NumberOfDeltas++;
              G4VSplitableHadron * targetSplitable=TheInvolvedNucleon[i]->GetSplitableHadron();
              G4double MassDec=std::sqrt(sqr(targetSplitable->GetDefinition()->GetPDGMass())
                                           + targetSplitable->Get4Momentum().perp2());

              G4int PDGcode = targetSplitable->GetDefinition()->GetPDGEncoding();
              G4ParticleDefinition* Old_def = targetSplitable->GetDefinition();

              G4int newPDGcode = PDGcode/10; newPDGcode=newPDGcode*10+4; // Delta
              G4ParticleDefinition* ptr = 
                 G4ParticleTable::GetParticleTable()->FindParticle(newPDGcode);
              targetSplitable->SetDefinition(ptr);
              G4double MassInc=std::sqrt(sqr(targetSplitable->GetDefinition()->GetPDGMass())
                                           + targetSplitable->Get4Momentum().perp2());
              if(SqrtS < SumMasses + MassInc - MassDec)                    // Uzhi 12.06.2012
              { // Change cannot be acsepted!
               targetSplitable->SetDefinition(Old_def);
               ProbDeltaIsobar = 0.;
              } else
              { // Change is acsepted.
               SumMasses += (MassInc - MassDec);
              }
            } 
          }   // end of for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
        }   // end of if(theNucleus.GetMassNumber() != 1)

//-------------------------------------------------------------

        G4LorentzRotation toCms(-1*Psum.boostVector());
        G4LorentzVector Ptmp=toCms*Pprojectile;
        if ( Ptmp.pz() <= 0. )                                
        {  // "String" moving backwards in  CMS, abort collision !!
           //G4cout << " abort ColliDeleteVSplitableHadronsion!! " << G4endl;
         return false; 
        }

        G4LorentzRotation toLab(toCms.inverse());

        Ptmp=toCms*Ptarget;                      
        G4double YtargetNucleus=Ptmp.rapidity();
//G4cout<<"YtargetNucleus "<<YtargetNucleus<<G4endl;
//-------------------------------------------------------------
//------- Ascribing of the involved nucleons Pt and Xminus ----
        G4double Dcor        = theParameters->GetDofNuclearDestruction()/
                                               theNucleus->GetMassNumber();

        G4double AveragePt2  = theParameters->GetPt2ofNuclearDestruction();
        G4double maxPtSquare = theParameters->GetMaxPt2ofNuclearDestruction();
//G4cout<<"Dcor "<<theParameters->GetDofNuclearDestruction()<<" "<<Dcor<<" AveragePt2 "<<AveragePt2<<G4endl;
        G4double M2target(0.);
        G4double WminusTarget(0.);
        G4double WplusProjectile(0.);

        G4int    NumberOfTries(0);
        G4double ScaleFactor(1.);
        G4bool OuterSuccess(true);
        do    // while (!OuterSuccess)
        {
          OuterSuccess=true;

          do    // while (SqrtS < Mprojectile + std::sqrt(M2target))
          {     // while (DecayMomentum < 0.)

            NumberOfTries++;

            if(NumberOfTries == 100*(NumberOfTries/100))   // 100
            { // At large number of tries it would be better to reduce the values
              ScaleFactor/=2.;
              Dcor       *=ScaleFactor;
              AveragePt2 *=ScaleFactor;
            }

            G4ThreeVector PtSum(0.,0.,0.);
            G4double XminusSum(0.);
            G4double Xminus(0.);
            G4bool InerSuccess=true;

            do      // while(!InerSuccess);
            {
             InerSuccess=true;

             PtSum    =G4ThreeVector(0.,0.,0.);
             XminusSum=0.;

	     for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
             {
               G4Nucleon * aNucleon = TheInvolvedNucleon[i];

               G4ThreeVector tmpPt = GaussianPt(AveragePt2, maxPtSquare);
               PtSum += tmpPt;

               G4ThreeVector tmpX=GaussianPt(Dcor*Dcor, 1.);
               Xminus=tmpX.x();
               XminusSum+=Xminus;

               G4LorentzVector tmp(tmpPt.x(),tmpPt.y(),Xminus,     // 6 Dec.2010
                                   aNucleon->Get4Momentum().e());  // 6 Dec.2010

               aNucleon->SetMomentum(tmp);
             }   // end of for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )

//---------------------------------------------------------------------------
             G4double DeltaX = (PtSum.x()-PnuclearResidual.x())/NumberOfInvolvedNucleon;
             G4double DeltaY = (PtSum.y()-PnuclearResidual.y())/NumberOfInvolvedNucleon;
             G4double DeltaXminus(0.);

             if(ResidualMassNumber == 0)
             {
              DeltaXminus = (XminusSum-1.)/NumberOfInvolvedNucleon;
             }
             else
             {
              DeltaXminus = -1./theNucleus->GetMassNumber();
             }

             XminusSum=1.;
             M2target =0.;

	     for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
             {
               G4Nucleon * aNucleon = TheInvolvedNucleon[i];

               Xminus = aNucleon->Get4Momentum().pz() - DeltaXminus;
               XminusSum-=Xminus;               

               if(ResidualMassNumber == 0)               
               {
                if((Xminus <= 0.)   || (Xminus > 1.))    {InerSuccess=false; break;}
               } else
               {
                if((Xminus <= 0.)   || (Xminus > 1.) || 
                   (XminusSum <=0.) || (XminusSum > 1.)) {InerSuccess=false; break;}
               }                                          

               G4double Px=aNucleon->Get4Momentum().px() - DeltaX;
               G4double Py=aNucleon->Get4Momentum().py() - DeltaY;

//               if(!ProjectileIsAntiBaryon)                          // 4.12.2010
//               {
                M2target +=(aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass()*
                            aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass()  + 
                            Px*Px + Py*Py)/Xminus;
/*
               } else
               {
                M2target +=(aNucleon->Get4Momentum().e() *
                            aNucleon->Get4Momentum().e()  +      // 6.12.2010
                            Px*Px + Py*Py)/Xminus;
               }
*/
               G4LorentzVector tmp(Px,Py,Xminus,aNucleon->Get4Momentum().e()); // 6.12.2010
               aNucleon->SetMomentum(tmp);
             }   // end of for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )

             if(InerSuccess && (ResidualMassNumber != 0))
             {
              M2target +=(sqr(ResidualMass) + PnuclearResidual.perp2())/XminusSum; // Uzhi 16.05.2013
             }
//G4cout<<"InerSuccess "<<InerSuccess<<" "<<SqrtS<<" "<<Mprojectile + std::sqrt(M2target)<<" "<<std::sqrt(M2target)<<G4endl;
//G4int Uzhi;G4cin>>Uzhi;
            } while(!InerSuccess);
          } while (SqrtS < Mprojectile + std::sqrt(M2target));
//-------------------------------------------------------------
          G4double DecayMomentum2= S*S+M2projectile*M2projectile+M2target*M2target
                                    -2.*S*M2projectile - 2.*S*M2target 
                                         -2.*M2projectile*M2target;

          WminusTarget=(S-M2projectile+M2target+std::sqrt(DecayMomentum2))/2./SqrtS;
          WplusProjectile=SqrtS - M2target/WminusTarget;

          G4double Pzprojectile=WplusProjectile/2. - M2projectile/2./WplusProjectile;// 8.06.11
          G4double Eprojectile =WplusProjectile/2. + M2projectile/2./WplusProjectile;// 8.06.11
          G4double Yprojectile=0.5*std::log((Eprojectile+Pzprojectile)/
                                            (Eprojectile-Pzprojectile));            // 1.07.11

//G4cout<<"Yprojectile "<<Yprojectile<<G4endl;
//G4LorentzVector TestPprojectile=Pprojectile;
//TestPprojectile.setPz(Pzprojectile);  TestPprojectile.setE(Eprojectile);

//G4cout<<"DecayMomentum2 "<<DecayMomentum2<<G4endl;
//G4cout<<"WminusTarget WplusProjectile "<<WminusTarget<<" "<<WplusProjectile<<G4endl;
//G4int Uzhi;G4cin>>Uzhi;
//-------------------------------------------------------------
	  for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
          {
           G4Nucleon * aNucleon = TheInvolvedNucleon[i];
           G4LorentzVector tmp=aNucleon->Get4Momentum();

           G4double Mt2(0.);

//           if(!ProjectileIsAntiBaryon)                          // 4.12.2010
//           {
            Mt2 = sqr(tmp.x())+sqr(tmp.y())+
                  sqr(aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass());
/*
           } else
           {
            Mt2 = sqr(tmp.x())+sqr(tmp.y())+                   // 4.12.2010
                  sqr(aNucleon->Get4Momentum().e());    // sqr
           }
*/
           G4double Xminus=tmp.z();

           G4double Pz=-WminusTarget*Xminus/2. + Mt2/(2.*WminusTarget*Xminus);
           G4double E = WminusTarget*Xminus/2. + Mt2/(2.*WminusTarget*Xminus);
           G4double YtargetNucleon=0.5*std::log((E+Pz)/(E-Pz));   // 1.07.11 //Uzhi 20 Sept.

//G4cout<<"YtargetNucleon "<<YtargetNucleon<<G4endl;
//G4cout<<"YtargetNucleon-YtargetNucleus "<<YtargetNucleon-YtargetNucleus<<G4endl;
//G4cout<<"Yprojectile  YtargetNucleon "<<Yprojectile<<" "<<YtargetNucleon<<G4endl;
if((std::abs(YtargetNucleon-YtargetNucleus) > 2) || 
            (Yprojectile  < YtargetNucleon))        {OuterSuccess=false; break;} // 1.07.11

          }   // end of for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
//if(ProjectileIsAntiBaryon) {G4int Uzhi;G4cin>>Uzhi;}
        } while(!OuterSuccess);

//-------------------------------------------------------------
        G4double Pzprojectile=WplusProjectile/2. - M2projectile/2./WplusProjectile;
        G4double Eprojectile =WplusProjectile/2. + M2projectile/2./WplusProjectile;
        Pprojectile.setPz(Pzprojectile);  Pprojectile.setE(Eprojectile);
//G4cout<<"Proj after in CMS "<<Pprojectile<<G4endl;

        Pprojectile.transform(toLab);       // The work with the projectile
        primary->Set4Momentum(Pprojectile); // is finished at the moment.
//G4cout<<"Final proj mom "<<primary->Get4Momentum()<<G4endl;

//-------------------------------------------------------------
        G4ThreeVector Residual3Momentum(0.,0.,1.);

	for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
        {
           G4Nucleon * aNucleon = TheInvolvedNucleon[i];
           G4LorentzVector tmp=aNucleon->Get4Momentum();
//G4cout<<"trg "<<aNucleon->Get4Momentum()<<G4endl;
           Residual3Momentum-=tmp.vect();

           G4double Mt2(0.);

//           if(!ProjectileIsAntiBaryon)                          // 4.12.2010
//           {
            Mt2 = sqr(tmp.x())+sqr(tmp.y())+
                  sqr(aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass());
/*
           } else
           {
            Mt2 = sqr(tmp.x())+sqr(tmp.y())+                   // 4.12.2010
                  aNucleon->Get4Momentum().e()*aNucleon->Get4Momentum().e();
           }
*/
           G4double Xminus=tmp.z();

           G4double Pz=-WminusTarget*Xminus/2. + Mt2/(2.*WminusTarget*Xminus);
           G4double E = WminusTarget*Xminus/2. + Mt2/(2.*WminusTarget*Xminus);

           tmp.setPz(Pz); 
           tmp.setE(E);
//G4cout<<"Targ after in CMS "<<tmp<<G4endl;
           tmp.transform(toLab);

           aNucleon->SetMomentum(tmp);
//G4cout<<"Targ after in LAB "<<aNucleon->Get4Momentum()<<G4endl;
           G4VSplitableHadron * targetSplitable=aNucleon->GetSplitableHadron();
           targetSplitable->Set4Momentum(tmp);
           
        }   // end of for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )

//G4cout<<"ResidualMassNumber and Mom "<<ResidualMassNumber<<" "<<Residual3Momentum<<G4endl;
        G4double Mt2Residual=sqr(ResidualMass) +
                                 sqr(Residual3Momentum.x())+sqr(Residual3Momentum.y());
//==========================
//G4cout<<"WminusTarget Residual3Momentum.z() "<<WminusTarget<<" "<<Residual3Momentum.z()<<G4endl;
        G4double PzResidual=0.;
        G4double EResidual =0.;
        if(ResidualMassNumber != 0)
        {
         PzResidual=-WminusTarget*Residual3Momentum.z()/2. + 
                             Mt2Residual/(2.*WminusTarget*Residual3Momentum.z());
         EResidual = WminusTarget*Residual3Momentum.z()/2. + 
                             Mt2Residual/(2.*WminusTarget*Residual3Momentum.z());
        }
//==========================
        Residual4Momentum.setPx(Residual3Momentum.x());
        Residual4Momentum.setPy(Residual3Momentum.y());
        Residual4Momentum.setPz(PzResidual); 
        Residual4Momentum.setE(EResidual);
//G4cout<<"Residual4Momentum in CMS Y "<<Residual4Momentum.beta()<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
        Residual4Momentum.transform(toLab);
//G4cout<<"Residual4Momentum in Lab "<<Residual4Momentum<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
//-------------------------------------------------------------
 return true;
}

// ------------------------------------------------------------
G4bool G4FTFModel::ExciteParticipants()
{
//G4cout<<"G4FTFModel::ExciteParticipants() "<<G4endl;
        G4bool Successfull(true);  //(false); // 1.07.11

        theParticipants.StartLoop();
        G4int CurrentInteraction(0);   // Uzhi Feb26

        G4int MaxNumOfInelCollisions=G4int(theParameters->GetMaxNumberOfCollisions());
//G4cout<<"MaxNumOfInelCollisions "<<MaxNumOfInelCollisions<<G4endl;
        G4double NumberOfInel(0.);
//
        if(MaxNumOfInelCollisions > 0)  
        {   //  Plab > Pbound, Normal application of FTF is possible
         G4double ProbMaxNumber=theParameters->GetMaxNumberOfCollisions()-
                                                  MaxNumOfInelCollisions;
         if(G4UniformRand() < ProbMaxNumber) {MaxNumOfInelCollisions++;}
         NumberOfInel=MaxNumOfInelCollisions;
//G4cout<<"Plab > Pbound MaxNumOfInelCollisions "<<MaxNumOfInelCollisions<<G4endl;
        } else
        {   //  Plab < Pbound, Normal application of FTF is impossible, 
            //                 low energy corrections applied.
         if(theParticipants.theNucleus->GetMassNumber() > 1)
         {
          NumberOfInel = theParameters->GetProbOfInteraction();
          MaxNumOfInelCollisions = 1;
         } else
         { // Special case for hadron-nucleon interactions
          NumberOfInel = 1.;
          MaxNumOfInelCollisions = 1;
//G4cout<<"Plab < Pbound MaxNumOfInelCollisions "<<MaxNumOfInelCollisions<<G4endl;
         }
        }  // end of if(MaxNumOfInelCollisions > 0)
//
//G4cout<<"MaxNumOfInelCollisions MaxNumOfInelCollisions "<<MaxNumOfInelCollisions<<" "<<MaxNumOfInelCollisions<<G4endl;

	while (theParticipants.Next())
	{	   
           CurrentInteraction++;
	   const G4InteractionContent & collision=theParticipants.GetInteraction();

	   G4VSplitableHadron * projectile=collision.GetProjectile();
	   G4VSplitableHadron * target=collision.GetTarget();
//
//G4cout<<"Ppr tr "<<projectile<<" "<<target<<G4endl;
//G4cout<<"theInter Time "<<collision.GetInteractionTime()/fermi<<G4endl;
//G4cout<<"Int Status    "<<collision.GetStatus()<<" "<<CurrentInteraction<<G4endl;
//G4cout<<"Proj M "<<projectile->Get4Momentum()<<" "<<projectile->Get4Momentum().mag()<<G4endl;
//G4cout<<"Targ M "<<target->Get4Momentum()<<" "<<target->Get4Momentum().mag()<<G4endl;
//G4cout<<"ProbabilityOfElasticScatt "<<theParameters->GetProbabilityOfElasticScatt()<<G4endl;
//G4cout<<"ProbabilityOfAnnihilation "<<theParameters->GetProbabilityOfAnnihilation()<<G4endl;
//G4cout<<"projectile->GetStatus target->GetStatus "<<projectile->GetStatus()<<" "<<target->GetStatus()<<G4endl;
//if((projectile->GetStatus() == 1) && (target->GetStatus() ==1))
//
if(collision.GetStatus())                                          // Uzhi Feb26
{

//theParameters->SetProbabilityOfElasticScatt(1.);
//G4cout<<"before pro "<<projectile->Get4Momentum()<<" "<<projectile->Get4Momentum().mag()<<G4endl;
//G4cout<<"before tar "<<target->Get4Momentum()<<" "<<target->Get4Momentum().mag()<<G4endl;
//G4cout<<"Prob el "<<theParameters->GetProbabilityOfElasticScatt()<<G4endl;
//G4cout<<"Prob an "<<theParameters->GetProbabilityOfAnnihilation()<<G4endl;

//G4cout<<"Pr Tr "<<projectile->GetDefinition()->GetPDGEncoding()<<" "<<target->GetDefinition()->GetPDGEncoding()<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;

           if(G4UniformRand()< theParameters->GetProbabilityOfElasticScatt())
           { //   Elastic scattering -------------------------
//G4cout<<"Elastic FTF"<<G4endl;
            Successfull = theElastic->ElasticScattering(projectile, target, theParameters)
                       || Successfull;              // 9.06.2012
           }
           else if(G4UniformRand() > theParameters->GetProbabilityOfAnnihilation())
           { //   Inelastic scattering ---------------------- 
//G4cout<<"Inelastic FTF"<<G4endl;
//G4cout<<"NumberOfInel MaxNumOfInelCollisions "<<NumberOfInel<<" "<<MaxNumOfInelCollisions<<G4endl;
            if(G4UniformRand()< NumberOfInel/MaxNumOfInelCollisions)
            {
             if(theExcitation->ExciteParticipants(projectile, target, 
                                                 theParameters, theElastic))
             {
              NumberOfInel--;
//G4cout<<"Excitation FTF  Successfull "<<G4endl;
//G4cout<<"After  pro "<<projectile->Get4Momentum()<<" "<<projectile->Get4Momentum().mag()<<G4endl;
//G4cout<<"After  tar "<<target->Get4Momentum()<<" "<<target->Get4Momentum().mag()<<G4endl;
             } else
             {
//G4cout<<"Excitation FTF  Non Successfull -> Elastic scattering "<<Successfull<<G4endl;

              Successfull = theElastic->ElasticScattering(projectile, target, theParameters)
                         || Successfull;              // 9.06.2012

/*
              if(NumberOfInvolvedTargetNucleon > 1)
              {
               NumberOfInvolvedTargetNucleon--;
               target->SetStatus(0);  // 1->0 return nucleon to the target  VU 10.04.2012
              }
*/
             }
            } else // If NumOfInel
            {
//G4cout<<"Elastic at rejected inelastic scattering"<<G4endl;
             Successfull = theElastic->ElasticScattering(projectile, target, theParameters)
                        || Successfull;              // 9.06.2012
/*
             if(theElastic->ElasticScattering(projectile, target, theParameters))
             {
//              Successfull = Successfull || true;
             } else
             {
//            Successfull = Successfull || false;
//Successfull = Successfull && false;
Successfull = false; break;                         // 1.07.11

              if(NumberOfInvolvedTargetNucleon > 1)
              {
               NumberOfInvolvedTargetNucleon--;
               target->SetStatus(4); // 1->0 return nucleon to the target  VU 10.04.2012
              }
             }
*/
            }   // end if NumOfInel
           } 
           else  // Annihilation
           {
//G4cout<<"Annihilation"<<G4endl;
//G4cout<<"Before  pro "<<projectile->Get4Momentum()<<G4endl;
//G4cout<<"Before  tar "<<target->Get4Momentum()<<G4endl;
//G4cout<<"Mom pro "<<theProjectile.GetTotalMomentum()<<G4endl;
//if(theProjectile.GetTotalMomentum() < 2000.*MeV)           
{ 
            while (theParticipants.Next())
            {   
            G4InteractionContent & acollision=theParticipants.GetInteraction();//Uzhi Feb26

	     G4VSplitableHadron * NextProjectileNucleon=acollision.GetProjectile(); // Uzhi Feb26
	     G4VSplitableHadron * NextTargetNucleon    =acollision.GetTarget();
             if((projectile == NextProjectileNucleon) ||
                (target     == NextTargetNucleon        )) acollision.SetStatus(0);
//             if(target != NextTargetNucleon) NextTargetNucleon->SetStatus(0); // Uzhi Feb26
            }

            theParticipants.StartLoop(); 
            for(G4int I=0; I < CurrentInteraction; I++) theParticipants.Next();
            
//-----------------------------------------
// 1Nov2011            AjustTargetNucleonForAnnihilation(projectile, target);
//-----------------------------------------
//G4cout<<"After Ajsd pro "<<projectile->Get4Momentum()<<G4endl;
//G4cout<<"After Ajst tar "<<target->Get4Momentum()<<G4endl;
}
            G4VSplitableHadron *AdditionalString=0;
            if(theAnnihilation->Annihilate(projectile, target, AdditionalString, theParameters))
            {
             Successfull = Successfull || true;
//G4cout<<G4endl<<"*AdditionalString "<<AdditionalString<<G4endl;
//G4cout<<"After  pro "<<projectile->Get4Momentum()<<G4endl;
//G4cout<<"After  tar "<<target->Get4Momentum()<<G4endl;

             if(AdditionalString != 0) theAdditionalString.push_back(AdditionalString);

//             break;

            } else
            {
              //A.R. 25-Jul-2012 : commenting the next line to fix a Coverity
              //                   "logically dead code".
              //Successfull = Successfull || false;

//             target->SetStatus(2);
            }
           } 
//
} // End of if((projectile->GetStatus() == 1) && (target->GetStatus() ==1))

        }       // end of while (theParticipants.Next())

//Successfull=true;
//G4cout<<"G4FTFModel::ExciteParticipants() Successfull "<<Successfull<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
	return Successfull;
}

//-------------------------------------------------------------------
void G4FTFModel::AjustTargetNucleonForAnnihilation(G4VSplitableHadron *SelectedAntiBaryon,
                                                   G4VSplitableHadron *SelectedTargetNucleon)
{
        G4LorentzVector Pparticipants=SelectedAntiBaryon->Get4Momentum()+
                                      SelectedTargetNucleon->Get4Momentum();

        G4V3DNucleus *theNucleus = GetWoundedNucleus();
//G4cout<<"Init A mass "<<theNucleus->GetMass()<<" "<<theNucleus->GetMassNumber()<<" "<<theNucleus->GetCharge()<<G4endl;

        ResidualExcitationEnergy=0.;
        G4int ResidualCharge    =theNucleus->GetCharge();
        G4int ResidualMassNumber=theNucleus->GetMassNumber();

        G4ThreeVector   P3nuclearResidual(0.,0.,0.);
        G4LorentzVector PnuclearResidual(0.,0.,0.,0.);


        G4double ExcitationEnergyPerWoundedNucleon=
                 theParameters->GetExcitationEnergyPerWoundedNucleon();
//-------
        G4Nucleon * aNucleon;
        theNucleus->StartLoop();
        G4int NumberOfHoles(0);
//G4cout<<"Start loop"<<G4endl;
	while ((aNucleon = theNucleus->GetNextNucleon()))
        {
         G4int CurrentStatus=0;
         if(aNucleon->AreYouHit()) 
         {
          if(aNucleon->GetSplitableHadron() == SelectedTargetNucleon)
          {CurrentStatus=1;}
          else
          {
           if(aNucleon->GetSplitableHadron()->GetSoftCollisionCount() == 0) 
           {CurrentStatus=0;}
           else {CurrentStatus=1;}
          }
         }
//G4cout<<"CurrentStatus "<<CurrentStatus<<G4endl;
         if(CurrentStatus != 0)
         {   // Participating nucleons
//G4cout<<"              Partic "<<aNucleon->GetSplitableHadron()->GetStatus()<<" "<<aNucleon->GetSplitableHadron()->GetSoftCollisionCount()<<G4endl;
          NumberOfHoles++;
          ResidualExcitationEnergy+=ExcitationEnergyPerWoundedNucleon;
          ResidualCharge-=(G4int) aNucleon->GetDefinition()->GetPDGCharge();
          ResidualMassNumber--;
         }
         else
         {   // Spectator nucleon
          PnuclearResidual+=aNucleon->Get4Momentum();
         }
	}   // end of while (theNucleus->GetNextNucleon())

//G4cout<<"Res Z M "<<ResidualCharge<<" "<<ResidualMassNumber<<G4endl;
//-------------------------------
        G4LorentzVector Psum=Pparticipants + PnuclearResidual;  // 4-momentum in CMS

// Transform momenta to cms and then rotate parallel to z axis;
        G4LorentzRotation toCms(-1*Psum.boostVector());

        G4LorentzVector Ptmp=toCms*Psum;

        toCms.rotateZ(-1*Ptmp.phi());
        toCms.rotateY(-1*Ptmp.theta());

        G4LorentzRotation toLab(toCms.inverse());

//-------------------------------
        G4double SqrtS=Psum.mag();
        G4double S    =sqr(SqrtS);

        G4double ResidualMass(0.);
        if(ResidualMassNumber != 0) 
        {
         ResidualMass=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(
                                                   ResidualCharge,ResidualMassNumber);
        } else {return;}

//G4cout<<"Res Mass E* "<<ResidualMass<<" "<<ResidualExcitationEnergy<<G4endl;

        if(ResidualMass > SqrtS) {return;}
        else 
        {
         if(ResidualMass+ResidualExcitationEnergy > SqrtS)
           ResidualExcitationEnergy = SqrtS-ResidualMass;
        }

        ResidualMass+=ResidualExcitationEnergy;
        G4double ResidualMass2=sqr(ResidualMass);
//G4cout<<"New Res Mass E* "<<ResidualMass<<" "<<ResidualExcitationEnergy<<G4endl;

//-------
        G4double ParticipantMass=Pparticipants.mag();
        G4double ParticipantMass2=sqr(ParticipantMass);

        if(ResidualMass + ParticipantMass > SqrtS) ParticipantMass=SqrtS-ResidualMass;

//G4cout<<"Parts P "<<Pparticipants<<G4endl;
//G4cout<<"Res Nuc "<<PnuclearResidual<<G4endl;

        G4double DecayMomentum2= 
                      sqr(S)+sqr(ParticipantMass2)+ sqr(ResidualMass2) -
                          2.*S*ParticipantMass2 - 2.*S*ResidualMass2 
                              -2.*ParticipantMass2*ResidualMass2;

        if(DecayMomentum2 < 0.) return;

        DecayMomentum2/=(4.*S);
        G4double DecayMomentum = std::sqrt(DecayMomentum2);
//G4cout<<"DecayMomentum "<<DecayMomentum<<G4endl;

        G4LorentzVector New_Pparticipants(0.,0.,DecayMomentum,
                                std::sqrt(DecayMomentum2+ParticipantMass2));

        G4LorentzVector New_PnuclearResidual(0.,0.,-DecayMomentum,
                                std::sqrt(DecayMomentum2+ResidualMass2));

//G4cout<<"New part P "<<New_Pparticipants<<" "<<New_Pparticipants.mag()<<G4endl;
//G4cout<<"New resd P "<<New_PnuclearResidual<<" "<<New_PnuclearResidual.mag()<<G4endl;

        New_Pparticipants.transform(toLab);
        New_PnuclearResidual.transform(toLab);
//G4cout<<"New part P "<<New_Pparticipants<<" "<<New_Pparticipants.mag()<<G4endl;
//G4cout<<"New resd P "<<New_PnuclearResidual<<" "<<New_PnuclearResidual.mag()<<G4endl;

        G4LorentzVector DeltaP_participants=(Pparticipants - New_Pparticipants)/2.;
        G4LorentzVector DeltaP_nuclearResidual=(PnuclearResidual - New_PnuclearResidual)/
                                               (G4double) ResidualMassNumber;
//------------------

        Ptmp=SelectedAntiBaryon->Get4Momentum() - DeltaP_participants;
        SelectedAntiBaryon->Set4Momentum(Ptmp);

        Ptmp=SelectedTargetNucleon->Get4Momentum() - DeltaP_participants;
        SelectedTargetNucleon->Set4Momentum(Ptmp);
//-----------

        //A.R. 25-Jul-2012 : Fix to Covery "Division by zero"
        //G4double DeltaExcitationEnergy=ResidualExcitationEnergy/((G4double) NumberOfHoles);
        G4double DeltaExcitationEnergy = 0.0;
        if ( NumberOfHoles != 0 ) {
          DeltaExcitationEnergy = ResidualExcitationEnergy / ((G4double) NumberOfHoles);
        }
 
// Re-definition of the wounded nucleon momenta
        theNucleus->StartLoop();
	while ((aNucleon = theNucleus->GetNextNucleon()))
        {
         G4int CurrentStatus=0;
         if(aNucleon->AreYouHit()) 
         {
          if(aNucleon->GetSplitableHadron() == SelectedTargetNucleon)
          {CurrentStatus=1;}
          else
          {
           if(aNucleon->GetSplitableHadron()->GetSoftCollisionCount() == 0) 
           {CurrentStatus=0;}
           else {CurrentStatus=1;}
          }
         }
//G4cout<<"CurrentStatus "<<CurrentStatus<<G4endl;
         if(CurrentStatus != 0)
         {   // Participating nucleons
//G4cout<<"              Partic "<<aNucleon->GetSplitableHadron()->GetStatus()<<" "<<aNucleon->GetSplitableHadron()->GetSoftCollisionCount()<<G4endl;
          aNucleon->SetBindingEnergy(DeltaExcitationEnergy);
         }
         else
         { // Spectator nucleon of nuclear residual
          Ptmp=aNucleon->Get4Momentum() - DeltaP_nuclearResidual;
          aNucleon->SetMomentum(Ptmp);
         }
	}   // end of while (theNucleus->GetNextNucleon())

//-------------------------------
        return;
}

// ------------------------------------------------------------
G4ExcitedStringVector * G4FTFModel::BuildStrings()
{	
// Loop over all collisions; find all primaries, and all target ( targets may 
//  be duplicate in the List ( to unique G4VSplitableHadrons)

	G4ExcitedStringVector * strings;
	strings = new G4ExcitedStringVector();
	
	std::vector<G4VSplitableHadron *> primaries;
	
        G4ExcitedString * FirstString(0);     // If there will be a kink,
        G4ExcitedString * SecondString(0);    // two strings will be produced.

	theParticipants.StartLoop();    // restart a loop 
//
	while ( theParticipants.Next() ) 
	{
	    const G4InteractionContent & interaction=theParticipants.GetInteraction();
//         if(interaction.GetStatus() != 0)   // Uzhi Feb26
         {
                 //  do not allow for duplicates ...
	    if ( primaries.end() == std::find(primaries.begin(), primaries.end(),
                                                interaction.GetProjectile()) )
	    	primaries.push_back(interaction.GetProjectile());     
         } // Uzhi Feb26
	}

//G4cout<<G4endl<<"primaries.size() "<<primaries.size()<<G4endl;
	for (unsigned int ahadron=0; ahadron < primaries.size() ; ahadron++)
	{
            G4bool isProjectile(0);

            if(primaries[ahadron]->GetStatus() <= 1) {isProjectile=true; } // VU 10.04.2012
//            if(primaries[ahadron]->GetStatus() == 3) {isProjectile=false;}

            FirstString=0; SecondString=0;
            theExcitation->CreateStrings(primaries[ahadron], isProjectile,
                                         FirstString, SecondString,
                                         theParameters);

	    if(FirstString  != 0) strings->push_back(FirstString);
            if(SecondString != 0) strings->push_back(SecondString);
//G4cout<<"Quarks in the string in FTF"<<FirstString->GetRightParton()->GetPDGcode()<<" "<<FirstString->GetLeftParton()->GetPDGcode()<<G4endl;

//G4cout<<FirstString<<" "<<SecondString<<G4endl;
	}

//G4cout<<"Check "<<strings->operator[](0)->GetRightParton()->GetPDGcode()<<" "<<strings->operator[](0)->GetLeftParton()->GetPDGcode()<<G4endl;
//
// Looking for spectator nucleons of the projectile-----------
        G4V3DNucleus * ProjectileNucleus =theParticipants.GetProjectileNucleus();
        if(ProjectileNucleus)
        {
         ProjectileNucleus->StartLoop();

         G4Nucleon *    ProjectileNucleon;
         while ( (ProjectileNucleon=ProjectileNucleus->GetNextNucleon()) )
         {
          if ( !ProjectileNucleon->AreYouHit() )
          {  // Projectile nucleon was involved in the interaction.

           G4VSplitableHadron * ProjectileSplitable=0;
           ProjectileSplitable= new G4DiffractiveSplitableHadron(*ProjectileNucleon);
           ProjectileNucleon->Hit(0);

	   G4bool isProjectile=true;
           FirstString=0; SecondString=0;
           theExcitation->CreateStrings(ProjectileSplitable,
                                        isProjectile,
                                        FirstString, SecondString,
                                        theParameters);
           if(FirstString  != 0) strings->push_back(FirstString);
           if(SecondString != 0) strings->push_back(SecondString);

           delete  ProjectileSplitable;     
          }
         } // End of while ( (ProjectileNucleon=ProjectileNucleus->GetNextNucleon()) )
        }  // End of if(ProjectileNucleus)

//G4cout<<G4endl<<"theAdditionalString.size() "<<theAdditionalString.size()<<G4endl;
        if(theAdditionalString.size() != 0)
        {
	 for (unsigned int  ahadron=0; ahadron < theAdditionalString.size() ; ahadron++)
	 {
            G4bool isProjectile(0);

            if(theAdditionalString[ahadron]->GetStatus() <= 1) {isProjectile=true; } // VU 10.04.2012
//            if(theAdditionalString[ahadron]->GetStatus() == 3) {isProjectile=false;}

            FirstString=0; SecondString=0;
            theExcitation->CreateStrings(theAdditionalString[ahadron], isProjectile,
                                         FirstString, SecondString,
                                         theParameters);

	    if(FirstString  != 0) strings->push_back(FirstString);
            if(SecondString != 0) strings->push_back(SecondString);
//G4cout<<"Quarks in the string in FTF"<<FirstString->GetRightParton()->GetPDGcode()<<" "<<FirstString->GetLeftParton()->GetPDGcode()<<G4endl;
//G4cout<<FirstString<<" "<<SecondString<<G4endl;
	 }
        }
//G4cout<<"Check "<<strings->operator[](0)->GetRightParton()->GetPDGcode()<<" "<<strings->operator[](0)->GetLeftParton()->GetPDGcode()<<G4endl;
//G4cout<<"Check "<<strings->operator[](1)->GetRightParton()->GetPDGcode()<<" "<<strings->operator[](1)->GetLeftParton()->GetPDGcode()<<G4endl;
//
//G4cout<<G4endl<<"NumberOfInvolvedNucleon "<<NumberOfInvolvedNucleon<<G4endl;
	for (G4int ahadron=0; ahadron < NumberOfInvolvedNucleon ; ahadron++)
	{
//G4cout<<"Nucleon status & int# "<<ahadron<<" "<<TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetStatus()<<" "<<TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetSoftCollisionCount()<<G4endl;
            if(TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetStatus() ==4)
            { // A nucleon is returned back to the nucleus after annihilation act for example
//G4cout<<" Delete 0"<<G4endl;
             delete TheInvolvedNucleon[ahadron]->GetSplitableHadron();
             G4VSplitableHadron *aHit=0; 
             TheInvolvedNucleon[ahadron]->Hit(aHit);
            }
            else if((TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetStatus() <=1)  && // VU 10.04.2012
            (TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetSoftCollisionCount() ==0))
            { // A nucleon is returned back to the nucleus after rejected interactions
              // due to an annihilation before
//G4cout<<" Delete int# 0"<<G4endl;
             delete TheInvolvedNucleon[ahadron]->GetSplitableHadron();
             G4VSplitableHadron *aHit=0; 
             TheInvolvedNucleon[ahadron]->Hit(aHit);
            }
            else if((TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetStatus() <=1)  && // VU 10.04.2012
            (TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetSoftCollisionCount() !=0))
            { // Nucleon which participate in the interactions, 
//G4cout<<"Taken 1 !=0"<<G4endl;
	     G4bool isProjectile=false;
             FirstString=0; SecondString=0;
             theExcitation->CreateStrings(
                            TheInvolvedNucleon[ahadron]->GetSplitableHadron(),
                                          isProjectile,
                                          FirstString, SecondString,
                                          theParameters);
	     if(FirstString  != 0) strings->push_back(FirstString);
             if(SecondString != 0) strings->push_back(SecondString);
            }
            else if(TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetStatus() ==2)
            { // Nucleon which was involved in the Reggeon cascading
//G4cout<<"Taken st 2"<<G4endl;
	     G4bool isProjectile=false;
             FirstString=0; SecondString=0;
             theExcitation->CreateStrings(
                            TheInvolvedNucleon[ahadron]->GetSplitableHadron(),
                                          isProjectile,
                                          FirstString, SecondString,
                                          theParameters);
	     if(FirstString  != 0) strings->push_back(FirstString);
             if(SecondString != 0) strings->push_back(SecondString);
//G4cout<<FirstString<<" "<<SecondString<<G4endl;
            }
            else if(TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetStatus() ==3)
            { // Nucleon which has participated in annihilation and disappiered
//G4cout<<"Status 3 "<<G4endl;
              TheInvolvedNucleon[ahadron]->SetBindingEnergy(theParameters->GetExcitationEnergyPerWoundedNucleon());
            }
            else {}

	}
/*
G4cout<<"Check "<<strings->operator[](0)->GetRightParton()->GetPDGcode()<<" "<<strings->operator[](0)->GetLeftParton()->GetPDGcode()<<G4endl;
G4cout<<"Check "<<strings->operator[](1)->GetRightParton()->GetPDGcode()<<" "<<strings->operator[](1)->GetLeftParton()->GetPDGcode()<<G4endl;
//G4cout<<"Check "<<strings->operator[](2)->GetRightParton()->GetPDGcode()<<" "<<strings->operator[](2)->GetLeftParton()->GetPDGcode()<<G4endl;

G4cout<<"*** "<<strings->operator[](0)->GetRightParton()<<" "<<strings->operator[](0)->GetLeftParton()<<G4endl;
G4cout<<"*** "<<strings->operator[](1)->GetRightParton()<<" "<<strings->operator[](1)->GetLeftParton()<<G4endl;
//G4cout<<"*** "<<strings->operator[](2)->GetRightParton()<<" "<<strings->operator[](2)->GetLeftParton()<<G4endl;
*/
	std::for_each(primaries.begin(), primaries.end(), DeleteVSplitableHadron());
	primaries.clear();
/*
G4cout<<"*** "<<strings->operator[](0)->GetRightParton()<<" "<<strings->operator[](0)->GetLeftParton()<<G4endl;
G4cout<<"*** "<<strings->operator[](1)->GetRightParton()<<" "<<strings->operator[](1)->GetLeftParton()<<G4endl;
G4cout<<"*** "<<strings->operator[](2)->GetRightParton()<<" "<<strings->operator[](2)->GetLeftParton()<<G4endl;

G4cout<<"Check "<<strings->operator[](0)->GetRightParton()->GetPDGcode()<<" "<<strings->operator[](0)->GetLeftParton()->GetPDGcode()<<G4endl;
G4cout<<"Check "<<strings->operator[](1)->GetRightParton()->GetPDGcode()<<" "<<strings->operator[](1)->GetLeftParton()->GetPDGcode()<<G4endl;
G4cout<<"Check "<<strings->operator[](2)->GetRightParton()->GetPDGcode()<<" "<<strings->operator[](2)->GetLeftParton()->GetPDGcode()<<G4endl;
*/

/*
for (unsigned int ahadron=0; ahadron < strings->size() ; ahadron++)
{
G4cout<<ahadron<<" "<<strings->operator[](ahadron)->GetRightParton()->GetPDGcode()<<" "<<strings->operator[](ahadron)->GetLeftParton()->GetPDGcode()<<G4endl;
}
G4cout<<"------------------------"<<G4endl;
*/
//G4int Uzhi; G4cin >> Uzhi;
	return strings;
}
// ------------------------------------------------------------
void G4FTFModel::GetResidualNucleus()
{ // This method is needed for the correct application of G4PrecompoundModelInterface
        G4double DeltaExcitationE=ResidualExcitationEnergy/
                                  (G4double) NumberOfInvolvedNucleon;
        G4LorentzVector DeltaPResidualNucleus = Residual4Momentum/
                                  (G4double) NumberOfInvolvedNucleon;

	for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
        {
         G4Nucleon * aNucleon = TheInvolvedNucleon[i];
//         G4LorentzVector tmp=aNucleon->Get4Momentum()-DeltaPResidualNucleus;
         G4LorentzVector tmp=-DeltaPResidualNucleus;
         aNucleon->SetMomentum(tmp);
         aNucleon->SetBindingEnergy(DeltaExcitationE);
        }   // end of for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )

}

// ------------------------------------------------------------
G4ThreeVector G4FTFModel::GaussianPt(G4double AveragePt2, G4double maxPtSquare) const
{            //  @@ this method is used in FTFModel as well. Should go somewhere common!
	
	G4double Pt2(0.);
        if(AveragePt2 <= 0.) {Pt2=0.;}
        else
        {
         Pt2 = -AveragePt2 * std::log(1. + G4UniformRand() * 
                (std::exp(-maxPtSquare/AveragePt2)-1.)); 
	}
	G4double Pt=std::sqrt(Pt2);
	G4double phi=G4UniformRand() * twopi;
	
	return G4ThreeVector (Pt*std::cos(phi), Pt*std::sin(phi), 0.);    
}

void G4FTFModel::ModelDescription(std::ostream& desc) const
{
	desc << "please add description here" << G4endl;
}
