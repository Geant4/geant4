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
// $Id: G4FTFModel.cc,v 1.38 2010/12/07 10:42:40 vuzhinsk Exp $
// GEANT4 tag $Name:  $
//

// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4FTFModel ----------------
//             by Gunter Folger, May 1998.
//       class implementing the excitation in the FTF Parton String Model
// ------------------------------------------------------------

#include "G4FTFModel.hh"
#include "G4FTFParameters.hh"
#include "G4FTFParticipants.hh"
#include "G4DiffractiveSplitableHadron.hh"
#include "G4InteractionContent.hh"
#include "G4LorentzRotation.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"
#include <utility> 
#include "G4IonTable.hh"

// Class G4FTFModel 

G4FTFModel::G4FTFModel():theExcitation(new G4DiffractiveExcitation()),
                         theElastic(new G4ElasticHNScattering()),
                         theAnnihilation(new G4FTFAnnihilation())
{
	G4VPartonStringModel::SetThisPointer(this);
        theParameters=0;
	NumberOfInvolvedNucleon=0;
        NumberOfInvolvedNucleonOfProjectile=0;
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
   if( NumberOfInvolvedNucleonOfProjectile != 0)
   {
    for(G4int i=0; i < NumberOfInvolvedNucleonOfProjectile; i++)
    {
     G4VSplitableHadron * aNucleon = TheInvolvedNucleonOfProjectile[i]->GetSplitableHadron();
     if(aNucleon) delete aNucleon;
    }
   }
}

const G4FTFModel & G4FTFModel::operator=(const G4FTFModel &)
{
	throw G4HadronicException(__FILE__, __LINE__, "G4FTFModel::operator= is not meant to be accessed ");
	return *this;
}

int G4FTFModel::operator==(const G4FTFModel &right) const
{
	return this==&right;
}

int G4FTFModel::operator!=(const G4FTFModel &right) const
{
	return this!=&right;
}

// ------------------------------------------------------------
void G4FTFModel::Init(const G4Nucleus & aNucleus, const G4DynamicParticle & aProjectile)
{
	theProjectile = aProjectile;  

        G4double ProtonMass=G4Proton::Proton()->GetPDGMass();
        G4double S(0.);  // CMS energy squered
/*
G4cout<<"FTF init Pro Name "<<theProjectile.GetDefinition()->GetParticleName()<<G4endl;
G4cout<<"FTF init Pro Mass "<<theProjectile.GetMass()<<" "<<theProjectile.GetMomentum()<<G4endl;
G4cout<<"FTF init Pro B Q  "<<theProjectile.GetDefinition()->GetBaryonNumber()<<" "<<(G4int) theProjectile.GetDefinition()->GetPDGCharge()<<G4endl; 
G4cout<<"FTF init A Z "<<aNucleus.GetA_asInt()<<" "<<aNucleus.GetZ_asInt()<<G4endl;
G4cout<<"             "<<aNucleus.GetN()<<" "<<aNucleus.GetZ()<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
*/

        theParticipants.Init(aNucleus.GetA_asInt(),aNucleus.GetZ_asInt());

        if(std::abs(theProjectile.GetDefinition()->GetBaryonNumber()) <= 1) 
        { // Projectile is a hadron
         S = sqr( theProjectile.GetMass() ) + sqr( ProtonMass ) +
                 2*ProtonMass*theProjectile.GetTotalEnergy();
        }


        if(theProjectile.GetDefinition()->GetBaryonNumber() > 1) 
        { // Projectile is a nucleus
         theParticipants.InitProjectileNucleus(
                      theProjectile.GetDefinition()->GetBaryonNumber(),
              (G4int) theProjectile.GetDefinition()->GetPDGCharge()    );
         S =         2.*sqr( ProtonMass ) + 2*ProtonMass*
             theProjectile.GetTotalEnergy()/theProjectile.GetDefinition()->GetBaryonNumber();

         G4ThreeVector BoostVector=theProjectile.GetMomentum()/theProjectile.GetTotalEnergy();
         theParticipants.theProjectileNucleus->DoLorentzBoost(BoostVector);
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

         S =        2.*sqr( ProtonMass ) + 2*ProtonMass*
                      theProjectile.GetTotalEnergy()/
             std::abs(theProjectile.GetDefinition()->GetBaryonNumber());
        }

// ------------------------------------------------------------------------
      if( theParameters != 0 ) delete theParameters;
      theParameters = new G4FTFParameters(theProjectile.GetDefinition(),
                                          aNucleus.GetA_asInt(),aNucleus.GetZ_asInt(),
                                          S);

// To turn on/off (1/0) elastic scattering close/open ...
//theParameters->SetProbabilityOfElasticScatt(0.); 
//G4cout<<theParameters->GetProbabilityOfElasticScatt()<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;

   if(theAdditionalString.size() != 0)
   {
    std::for_each(theAdditionalString.begin(), theAdditionalString.end(), 
                  DeleteVSplitableHadron());
   }
   theAdditionalString.clear();
}

// ------------------------------------------------------------
G4ExcitedStringVector * G4FTFModel::GetStrings()
{ 
        G4ExcitedStringVector * theStrings(0);

	theParticipants.GetList(theProjectile,theParameters);
//        StoreInvolvedNucleon();

        G4bool Success(true);

        if((std::abs(theProjectile.GetDefinition()->GetBaryonNumber()) <= 1) &&
                    (theProjectile.GetDefinition()->GetBaryonNumber() != -1)   )
        { // Standard variant of FTF for projectile hadron/nucleon
//G4cout<<"Standard variant of FTF for projectile hadron/nucleon"<<G4endl;
         ReggeonCascade(); 
         Success=PutOnMassShell(); 
         GetResidualNucleus();
        } 

        if(theProjectile.GetDefinition()->GetBaryonNumber() > 1)
        { // Variant of FTF for projectile nuclei
//G4cout<<"Variant of FTF for projectile nuclei"<<G4endl;
         StoreInvolvedNucleon();
         ReggeonCascade(); 
         Success=PutOnMassShell(); 
         GetResidualNucleus();
        } 

        G4bool LowE_Anti_Ion(false);
        if(theProjectile.GetDefinition()->GetBaryonNumber() <= -1) 
        { // Projectile is Anti-baryon or Anti-Nucleus
//G4cout<<"Projectile is Anti-baryon or Anti-Nucleus "<<G4endl;
         StoreInvolvedNucleon();
         if(theProjectile.GetTotalMomentum()/
            std::abs(theProjectile.GetDefinition()->GetBaryonNumber()) > 5000.*MeV)
         {// High energy interaction
//G4cout<<"High energy interaction "<<G4endl;
          ReggeonCascade(); 
          Success=PutOnMassShell(); 
          GetResidualNucleus();
         }
         else
         {
//G4cout<<"Low energy interaction "<<G4endl;
          LowE_Anti_Ion=true;
          Success=true;
         }
        }

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
        NumberOfInvolvedNucleon=0;

        theParticipants.StartLoop();

	while (theParticipants.Next())
	{   
	   const G4InteractionContent & collision=theParticipants.GetInteraction();
           G4Nucleon * TargetNucleon=collision.GetTargetNucleon();

           TheInvolvedNucleon[NumberOfInvolvedNucleon]=TargetNucleon;
           NumberOfInvolvedNucleon++;
//G4cout<<G4endl<<"Prim NumberOfInvolvedNucleon "<<NumberOfInvolvedNucleon<<G4endl;
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

          G4LorentzVector Ptmp=toCms*P_participants;
//G4cout<<"Ppart in CMS "<<Ptmp<<G4endl;

          if ( Ptmp.pz() <= 0. )                                
          {  // "String" moving backwards in  CMS, abort collision !!
           //G4cout << " abort ColliDeleteVSplitableHadronsion!! " << G4endl;
           return false; 
          }

          toCms.rotateZ(-1*Ptmp.phi());              // Uzhi 5.12.09
          toCms.rotateY(-1*Ptmp.theta());            // Uzhi 5.12.09

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

        G4bool ProjectileIsAntiBaryon = primary->GetDefinition()->GetBaryonNumber() < 0;

//G4cout<<"PutOnMass Pprojectile "<<Pprojectile<<G4endl;
// To get original projectile particle

        if(Pprojectile.z() < 0.){return false;}

        G4double Mprojectile  = Pprojectile.mag();
        G4double M2projectile = Pprojectile.mag2();
//-------------------------------------------------------------
	G4LorentzVector Psum      = Pprojectile;

        G4double        SumMasses = Mprojectile + 20.*MeV; // 13.12.09
                                               // Separation energy for projectile
        if(ProjectileIsAntiBaryon) {SumMasses = Mprojectile;}
//G4cout<<"SumMasses Pr "<<SumMasses<<G4endl;
//--------------- Target nucleus ------------------------------
        G4V3DNucleus *theNucleus = GetWoundedNucleus();
        G4Nucleon * aNucleon;
        G4int ResidualMassNumber=theNucleus->GetMassNumber();
        G4int ResidualCharge    =theNucleus->GetCharge();

        ResidualExcitationEnergy=0.;
	G4LorentzVector PnuclearResidual(0.,0.,0.,0.);

        G4double ExcitationEnergyPerWoundedNucleon=
                  theParameters->GetExcitationEnergyPerWoundedNucleon();

//G4cout<<"ExcitationEnergyPerWoundedNucleon "<<ExcitationEnergyPerWoundedNucleon<<G4endl;

        theNucleus->StartLoop();

	while ((aNucleon = theNucleus->GetNextNucleon()))
        {
         if(aNucleon->AreYouHit())
         {   // Involved nucleons
//G4cout<<"PutOn Tr "<<aNucleon->Get4Momentum()<<G4endl;
          Psum += aNucleon->Get4Momentum();
          if(!ProjectileIsAntiBaryon)
          {
           SumMasses += aNucleon->GetDefinition()->GetPDGMass();  
           SumMasses += 20.*MeV;   // 13.12.09 Separation energy for a nucleon
           ResidualExcitationEnergy+=ExcitationEnergyPerWoundedNucleon;
          } else 
          {
           SumMasses += aNucleon->Get4Momentum().mag();           // 4.12.2010
           G4LorentzVector tmp=aNucleon->Get4Momentum();
           tmp.setE(aNucleon->Get4Momentum().mag());   // It is need to save mass 6.12.2011
           aNucleon->SetMomentum(tmp);
          }

//G4cout<<"SumMasses Tr "<<SumMasses<<G4endl;
//if(SumMasses+ResidualExcitationEnergy > Psum.mag())
//{
// SetStatus ???
// if(!ProjectileIsAntiBaryon)
// {
//  SumMasses -= aNucleon->GetDefinition()->GetPDGMass();
//  SumMasses -= 20.*MeV;
//  ResidualExcitationEnergy-=ExcitationEnergyPerWoundedNucleon;
// } else
// {
//  SumMasses -= aNucleon->Get4Momentum().mag();
//  tmp ???
// }
//}
          ResidualMassNumber--;
          ResidualCharge-=(G4int) aNucleon->GetDefinition()->GetPDGCharge();
         }
         else
         {   // Spectator nucleons
          PnuclearResidual += aNucleon->Get4Momentum();
         }  // end of if(!aNucleon->AreYouHit())
	}   // end of while (theNucleus->GetNextNucleon())

        Psum += PnuclearResidual;
//G4cout<<"ResidualCharge ,ResidualMassNumber "<<ResidualCharge<<" "<<ResidualMassNumber<<G4endl;
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
        SumMasses += ResidualMass;
//G4cout<<"SumMasses + ResM "<<SumMasses<<G4endl;
//G4cout<<"Psum "<<Psum<<G4endl;
//-------------------------------------------------------------

        G4double SqrtS=Psum.mag();
        G4double     S=Psum.mag2();

//G4cout<<"SqrtS < SumMasses "<<SqrtS<<" "<<SumMasses<<G4endl;

        if(SqrtS < SumMasses)      {return false;} // It is impossible to simulate
                                                   // after putting nuclear nucleons
                                                   // on mass-shell

        if(SqrtS < SumMasses+ResidualExcitationEnergy) {ResidualExcitationEnergy=0.;}

        ResidualMass +=ResidualExcitationEnergy;
        SumMasses    +=ResidualExcitationEnergy;
//G4cout<<"ResidualMass "<<ResidualMass<<" "<<SumMasses<<G4endl;
//-------------------------------------------------------------
// Sampling of nucleons what are transfered to delta-isobars --
        G4int MaxNumberOfDeltas = (int)((SqrtS - SumMasses)/(400.*MeV));
        G4int NumberOfDeltas(0);

        if(theNucleus->GetMassNumber() != 1)
        {
          G4double ProbDeltaIsobar(0.);  // 1. *** Can be set if it is needed
	  for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
          {
            if((G4UniformRand() < ProbDeltaIsobar)&&(NumberOfDeltas < MaxNumberOfDeltas))
            {
              NumberOfDeltas++;
              G4VSplitableHadron * targetSplitable=TheInvolvedNucleon[i]->GetSplitableHadron();
              SumMasses-=targetSplitable->GetDefinition()->GetPDGMass();

              G4int PDGcode = targetSplitable->GetDefinition()->GetPDGEncoding();
              G4int newPDGcode = PDGcode/10; newPDGcode=newPDGcode*10+4; // Delta
              G4ParticleDefinition* ptr = 
                 G4ParticleTable::GetParticleTable()->FindParticle(newPDGcode);
              targetSplitable->SetDefinition(ptr);
              SumMasses+=targetSplitable->GetDefinition()->GetPDGMass();
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

//        toCms.rotateZ(-1*Ptmp.phi());              // Uzhi 5.12.09
//        toCms.rotateY(-1*Ptmp.theta());            // Uzhi 5.12.09
	
        G4LorentzRotation toLab(toCms.inverse());

//-------------------------------------------------------------
//------- Ascribing of the involved nucleons Pt and Xminus ----
        G4double Dcor        = theParameters->GetDofNuclearDestruction()/
                                               theNucleus->GetMassNumber();

        G4double AveragePt2  = theParameters->GetPt2ofNuclearDestruction();
        G4double maxPtSquare = theParameters->GetMaxPt2ofNuclearDestruction();
//G4cout<<"Dcor "<<Dcor<<" AveragePt2 "<<AveragePt2<<G4endl;
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
//G4cout<<"NumberOfTries "<<NumberOfTries<<G4endl;
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

//             G4LorentzVector tmp(tmpPt.x(),tmpPt.y(),Xminus,0.); // 6 Dec.2010
               G4LorentzVector tmp(tmpPt.x(),tmpPt.y(),Xminus,     // 6 Dec.2010
                                   aNucleon->Get4Momentum().e());// 6 Dec.2010

//G4cout<<"Inv i mom "<<i<<" "<<tmp<<G4endl;
               aNucleon->SetMomentum(tmp);
             }   // end of for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )

//---------------------------------------------------------------------------
             G4double DeltaX(0.);
             G4double DeltaY(0.);
             G4double DeltaXminus(0.);

//G4cout<<"ResidualMassNumber "<<ResidualMassNumber<<" "<<PtSum<<G4endl;
             if(ResidualMassNumber == 0)
             {
              DeltaX      = PtSum.x()/NumberOfInvolvedNucleon;
              DeltaY      = PtSum.y()/NumberOfInvolvedNucleon;
              DeltaXminus = (XminusSum-1.)/NumberOfInvolvedNucleon;
             }
             else
             {
              DeltaXminus = -1./theNucleus->GetMassNumber();
             }
//G4cout<<"Dx y xmin "<<DeltaX<<" "<<DeltaY<<" "<<DeltaXminus<<G4endl;
             XminusSum=1.;
             M2target =0.;

	     for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
             {
               G4Nucleon * aNucleon = TheInvolvedNucleon[i];

               Xminus = aNucleon->Get4Momentum().pz() - DeltaXminus;
               XminusSum-=Xminus;               
//G4cout<<" i X-sum "<<i<<" "<<Xminus<<" "<<XminusSum<<G4endl;
               if(ResidualMassNumber == 0)                // Uzhi 5.07.10
               {
                if((Xminus <= 0.)   || (Xminus > 1.))    {InerSuccess=false; break;}
               } else
               {
                if((Xminus <= 0.)   || (Xminus > 1.) || 
                   (XminusSum <=0.) || (XminusSum > 1.)) {InerSuccess=false; break;}
               }                                          // Uzhi 5.07.10

               G4double Px=aNucleon->Get4Momentum().px() - DeltaX;
               G4double Py=aNucleon->Get4Momentum().py() - DeltaY;

               if(!ProjectileIsAntiBaryon)                          // 4.12.2010
               {
                M2target +=(aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass()*
                            aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass()  + 
                            Px*Px + Py*Py)/Xminus;
               } else
               {
                M2target +=(aNucleon->Get4Momentum().e() *
                            aNucleon->Get4Momentum().e()  +      // 6.12.2010
                            Px*Px + Py*Py)/Xminus;
               }

               G4LorentzVector tmp(Px,Py,Xminus,aNucleon->Get4Momentum().e()); // 6.12.2010
               aNucleon->SetMomentum(tmp);
             }   // end of for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )

//G4cout<<"Rescale O.K."<<G4endl;

             if(InerSuccess && (ResidualMassNumber != 0))
             {
              M2target +=(ResidualMass*ResidualMass + PtSum.mag2())/XminusSum;
             }
//G4cout<<"InerSuccess "<<InerSuccess<<G4endl;
//G4int Uzhi;G4cin>>Uzhi;
            } while(!InerSuccess);
          } while (SqrtS < Mprojectile + std::sqrt(M2target));
//-------------------------------------------------------------
          G4double DecayMomentum2= S*S+M2projectile*M2projectile+M2target*M2target
                                    -2.*S*M2projectile - 2.*S*M2target 
                                         -2.*M2projectile*M2target;

          WminusTarget=(S-M2projectile+M2target+std::sqrt(DecayMomentum2))/2./SqrtS;
          WplusProjectile=SqrtS - M2target/WminusTarget;
//G4cout<<"DecayMomentum2 "<<DecayMomentum2<<G4endl;
//-------------------------------------------------------------
	  for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
          {
           G4Nucleon * aNucleon = TheInvolvedNucleon[i];
           G4LorentzVector tmp=aNucleon->Get4Momentum();
//G4cout<<"Invol Nucl "<<tmp<<G4endl;
           G4double Mt2(0.);

           if(!ProjectileIsAntiBaryon)                          // 4.12.2010
           {
            Mt2 = sqr(tmp.x())+sqr(tmp.y())+
                  aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass()*
                  aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass();
           } else
           {
            Mt2 = sqr(tmp.x())+sqr(tmp.y())+                   // 4.12.2010
                  aNucleon->Get4Momentum().e();
           }
           G4double Xminus=tmp.z();

           G4double Pz=-WminusTarget*Xminus/2. + Mt2/(2.*WminusTarget*Xminus);
           G4double E = WminusTarget*Xminus/2. + Mt2/(2.*WminusTarget*Xminus);

           if( E+Pz > WplusProjectile ){OuterSuccess=false; break;}
          }   // end of for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
//G4int Uzhi;G4cin>>Uzhi;
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

           if(!ProjectileIsAntiBaryon)                          // 4.12.2010
           {
            Mt2 = sqr(tmp.x())+sqr(tmp.y())+
                  aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass()*
                  aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass();
           } else
           {
            Mt2 = sqr(tmp.x())+sqr(tmp.y())+                   // 4.12.2010
                  aNucleon->Get4Momentum().e()*aNucleon->Get4Momentum().e();
           }

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
//G4cout<<"Residual4Momentum in CMS "<<Residual4Momentum<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
        Residual4Momentum.transform(toLab);
//G4cout<<"Residual4Momentum in Lab "<<Residual4Momentum<<G4endl;
//-------------------------------------------------------------
 return true;
}

// ------------------------------------------------------------
G4bool G4FTFModel::ExciteParticipants()
{
//G4cout<<"G4FTFModel::ExciteParticipants() "<<G4endl;
        G4bool Successfull(false);

        theParticipants.StartLoop();

        G4int MaxNumOfInelCollisions=G4int(theParameters->GetMaxNumberOfCollisions());

        G4double NumberOfInel(0.);
//
        if(MaxNumOfInelCollisions > 0)  
        {   //  Plab > Pbound, Normal application of FTF is possible
         G4double ProbMaxNumber=theParameters->GetMaxNumberOfCollisions()-
                                                  MaxNumOfInelCollisions;
         if(G4UniformRand() < ProbMaxNumber) {MaxNumOfInelCollisions++;}
         NumberOfInel=MaxNumOfInelCollisions;
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
         }
        }  // end of if(MaxNumOfInelCollisions > 0)
//
//G4cout<<"MaxNumOfInelCollisions MaxNumOfInelCollisions "<<MaxNumOfInelCollisions<<" "<<MaxNumOfInelCollisions<<G4endl;

	while (theParticipants.Next())
	{	   
	   const G4InteractionContent & collision=theParticipants.GetInteraction();

	   G4VSplitableHadron * projectile=collision.GetProjectile();
	   G4VSplitableHadron * target=collision.GetTarget();
//
//G4cout<<"Ppr tr "<<projectile<<" "<<target<<G4endl;
//G4cout<<"theInter Time "<<collision.GetInteractionTime()/fermi<<G4endl;
//G4cout<<"Proj M "<<projectile->Get4Momentum()<<G4endl;
//G4cout<<"Targ M "<<target->Get4Momentum()<<G4endl;
//G4cout<<"ProbabilityOfElasticScatt "<<theParameters->GetProbabilityOfElasticScatt()<<G4endl;
//G4cout<<"ProbabilityOfAnnihilation "<<theParameters->GetProbabilityOfAnnihilation()<<G4endl;

if((projectile->GetStatus() == 1) && (target->GetStatus() ==1))
{

//theParameters->SetProbabilityOfElasticScatt(1.);
//G4cout<<"before pro "<<projectile->Get4Momentum()<<" "<<projectile->Get4Momentum().mag()<<G4endl;
//G4cout<<"before tar "<<target->Get4Momentum()<<" "<<target->Get4Momentum().mag()<<G4endl;
           if(G4UniformRand()< theParameters->GetProbabilityOfElasticScatt())
           { //   Elastic scattering -------------------------
//G4cout<<"Elastic FTF"<<G4endl;
            if(theElastic->ElasticScattering(projectile, target, theParameters))
            {
//G4cout<<"Elastic FTF  Successfull "<<target->GetStatus()<<G4endl;
//G4cout<<"After  pro "<<projectile->Get4Momentum()<<" "<<projectile->Get4Momentum().mag()<<G4endl;
//G4cout<<"After  tar "<<target->Get4Momentum()<<" "<<target->Get4Momentum().mag()<<G4endl;
            Successfull = Successfull || true;
            } else
            {
//G4cout<<"Elastic FTF  Not Successfull "<<target->GetStatus()<<G4endl;
             Successfull = Successfull || false;
             if(NumberOfInvolvedTargetNucleon > 1)
             {
              NumberOfInvolvedTargetNucleon--;
              target->SetStatus(0); // 1->0 return nucleon to the target VU 18.02.11
             }
            }
           }
           else if(G4UniformRand() > theParameters->GetProbabilityOfAnnihilation())
           { //   Inelastic scattering ---------------------- 
//G4cout<<"Inelastic FTF"<<G4endl;
//G4cout<<"MaxNumOfInelCollisions MaxNumOfInelCollisions "<<MaxNumOfInelCollisions<<" "<<MaxNumOfInelCollisions<<G4endl;
            if(G4UniformRand()< NumberOfInel/MaxNumOfInelCollisions)
            {
             if(theExcitation->ExciteParticipants(projectile, target, 
                                                 theParameters, theElastic))
             {
              Successfull = Successfull || true; 
              NumberOfInel--;
             } else
             {
              Successfull = Successfull || false;
              if(NumberOfInvolvedTargetNucleon > 1)
              {
               NumberOfInvolvedTargetNucleon--;
               target->SetStatus(0);  // 1->0 return nucleon to the target VU 18.02.11
              }
             }
            } else // If NumOfInel
            {
             if(theElastic->ElasticScattering(projectile, target, theParameters))
             {
              Successfull = Successfull || true;
             } else
             {
              Successfull = Successfull || false;
              if(NumberOfInvolvedTargetNucleon > 1)
              {
               NumberOfInvolvedTargetNucleon--;
               target->SetStatus(0); // 1->0 return nucleon to the target VU 18.02.11
              }
             }
            }   // end if NumOfInel
           } 
           else  // Annihilation
           {
//G4cout<<"Annihilation"<<G4endl;
//G4cout<<"After  pro "<<projectile->Get4Momentum()<<G4endl;
//G4cout<<"After  tar "<<target->Get4Momentum()<<G4endl;
//G4cout<<"Mom pro "<<theProjectile.GetTotalMomentum()<<G4endl;
if(theProjectile.GetTotalMomentum() < 2000.*MeV)
{ 
            while (theParticipants.Next())
            {   
             const G4InteractionContent & collision=theParticipants.GetInteraction();
	     G4VSplitableHadron * NextTargetNucleon=collision.GetTarget();
             NextTargetNucleon->SetStatus(0);
            }
//-----------------------------------------
            AjustTargetNucleonForAnnihilation(target);
//-----------------------------------------
//G4cout<<"After  pro "<<projectile->Get4Momentum()<<G4endl;
//G4cout<<"After  tar "<<target->Get4Momentum()<<G4endl;
}
            G4VSplitableHadron *AdditionalString=0;
            if(theAnnihilation->Annihilate(projectile, target, AdditionalString, theParameters))
            {
             Successfull = Successfull || true;
//G4cout<<G4endl<<"*AdditionalString "<<AdditionalString<<G4endl;
//G4cout<<"After  pro "<<projectile->Get4Momentum()<<G4endl;
//G4cout<<"After  tar "<<target->Get4Momentum()<<G4endl;

             if(AdditionalString != 0) theAdditionalString.push_back(AdditionalString);

             break;

            } else
            {
             Successfull = Successfull || false;
//             target->SetStatus(2);
            }
           } 
//
} // End of if((projectile->GetStatus() == 1) && (target->GetStatus() ==1))
        }       // end of while (theParticipants.Next())

//Successfull=true;
//G4cout<<"Successfull "<<Successfull<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
	return Successfull;
}

//-------------------------------------------------------------------
void G4FTFModel::AjustTargetNucleonForAnnihilation(G4VSplitableHadron *SelectedTargetNucleon)
{
        G4V3DNucleus *theNucleus = GetWoundedNucleus();
//G4cout<<"Init A mass "<<theNucleus->GetMass()<<" "<<theNucleus->GetMassNumber()<<" "<<theNucleus->GetCharge()<<G4endl;

        G4double SqrtS=theNucleus->GetMass();

        ResidualExcitationEnergy=0.;
        G4int ResidualCharge    =theNucleus->GetCharge();
        G4int ResidualMassNumber=theNucleus->GetMassNumber();

        G4ThreeVector P3nuclearResidual(0.,0.,0.);
        G4LorentzVector Pparticipant(0.,0.,0.,0.);

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
         if(aNucleon->AreYouHit()) CurrentStatus=aNucleon->GetSplitableHadron()->GetStatus();
         if(CurrentStatus != 0)
         {   // Participating nucleons
//G4cout<<"Partic "<<aNucleon->GetSplitableHadron()->GetStatus()<<G4endl;
          NumberOfHoles++;
          ResidualExcitationEnergy+=ExcitationEnergyPerWoundedNucleon;
          ResidualCharge-=(G4int) aNucleon->GetDefinition()->GetPDGCharge();
          ResidualMassNumber--;
          P3nuclearResidual-=aNucleon->Get4Momentum().vect();
          if(aNucleon->GetSplitableHadron() != SelectedTargetNucleon)
             Pparticipant+=aNucleon->Get4Momentum();
         }
	}   // end of while (theNucleus->GetNextNucleon())

//G4cout<<"Res Z M "<<ResidualCharge<<" "<<ResidualMassNumber<<G4endl;

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
//G4cout<<"New Res Mass E* "<<ResidualMass<<" "<<ResidualExcitationEnergy<<G4endl;

//-------
	G4LorentzVector P_SelectedTargetNucleon=SelectedTargetNucleon->Get4Momentum();
        G4LorentzVector PnuclearResidual(P3nuclearResidual,
                                         std::sqrt(P3nuclearResidual.mag2()+sqr(ResidualMass)));

//G4cout<<"Sel N P "<<P_SelectedTargetNucleon<<G4endl;
//G4cout<<"Res Nuc "<<PnuclearResidual<<G4endl;

        G4double NewNucleonEnergy=SqrtS-(PnuclearResidual+Pparticipant).mag();

        P_SelectedTargetNucleon.setE(NewNucleonEnergy);
        SelectedTargetNucleon->Set4Momentum(P_SelectedTargetNucleon);

//G4cout<<"NewSelP "<<P_SelectedTargetNucleon<<G4endl;

        G4double DeltaExcitationEnergy=ResidualExcitationEnergy/((G4double) NumberOfHoles);

// Re-definition of the wounded nucleon momenta
        theNucleus->StartLoop();
	while ((aNucleon = theNucleus->GetNextNucleon()))
        {
         if(aNucleon->AreYouHit())
         {   // Participating nucleons
          aNucleon->SetBindingEnergy(DeltaExcitationEnergy);
if(aNucleon->GetSplitableHadron() == SelectedTargetNucleon)  aNucleon->SetMomentum(P_SelectedTargetNucleon);
         }
        }
//
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
                 //  do not allow for duplicates ...

	    if ( primaries.end() == std::find(primaries.begin(), primaries.end(),
                                                interaction.GetProjectile()) )
	    	primaries.push_back(interaction.GetProjectile());     
	}

	unsigned int ahadron;
//G4cout<<G4endl<<"primaries.size() "<<primaries.size()<<G4endl;
	for ( ahadron=0; ahadron < primaries.size() ; ahadron++)
	{
            G4bool isProjectile(0);

            if(primaries[ahadron]->GetStatus() == 1) {isProjectile=true; }
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
	 for ( ahadron=0; ahadron < theAdditionalString.size() ; ahadron++)
	 {
            G4bool isProjectile(0);

            if(theAdditionalString[ahadron]->GetStatus() == 1) {isProjectile=true; }
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
//G4cout<<"Nucleon status "<<ahadron<<" "<<TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetStatus()<<G4endl;
            if(TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetStatus() ==0)
            { // A nucleon is returned back to the nucleus after annihilation act for example
             delete TheInvolvedNucleon[ahadron]->GetSplitableHadron();
             G4VSplitableHadron *aHit=0; 
             TheInvolvedNucleon[ahadron]->Hit(aHit);
            }
            else if((TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetStatus() ==1)||
                    (TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetStatus() ==2)  )
            { // Nucleon which participate in the interactions, 
              // or nucleon which is involved in the Reggeon cascading
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
            { // Nucleon which has participated in annihilation
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
