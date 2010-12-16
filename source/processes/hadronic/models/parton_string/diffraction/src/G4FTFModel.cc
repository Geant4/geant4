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
// $Id: G4FTFModel.cc,v 1.37 2010-11-15 10:02:38 vuzhinsk Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
                         theElastic(new G4ElasticHNScattering())
{
	G4VPartonStringModel::SetThisPointer(this);
        theParameters=0;
	NumberOfInvolvedNucleon=0;
}


G4FTFModel::~G4FTFModel()
{
// Because FTF model can be called for various particles
// theParameters must be erased at the end of each call.
// Thus the delete is also in G4FTFModel::GetStrings() method
   if( theParameters != 0 ) delete theParameters; 
   if( theExcitation != 0 ) delete theExcitation;
   if( theElastic    != 0 ) delete theElastic; 

   if( NumberOfInvolvedNucleon != 0)
   {
    for(G4int i=0; i < NumberOfInvolvedNucleon; i++)
    {
     G4VSplitableHadron * aNucleon = TheInvolvedNucleon[i]->GetSplitableHadron();
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
//G4cout<<"FTF init Pro "<<theProjectile.GetMass()<<" "<<theProjectile.GetMomentum()<<G4endl;
//G4cout<<"FTF init A Z "<<aNucleus.GetA_asInt()<<" "<<aNucleus.GetZ_asInt()<<G4endl;
//G4cout<<"             "<<aNucleus.GetN()<<" "<<aNucleus.GetZ()<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;

	theParticipants.Init(aNucleus.GetA_asInt(),aNucleus.GetZ_asInt()); 
//G4cout<<"End nucl init"<<G4endl;
// ----------- N-mass number Z-charge -------------------------

// --- cms energy
        G4double s = sqr( theProjectile.GetMass() ) +
                     sqr( G4Proton::Proton()->GetPDGMass() ) +
                     2*theProjectile.GetTotalEnergy()*G4Proton::Proton()->GetPDGMass();

      if( theParameters != 0 ) delete theParameters;
      theParameters = new G4FTFParameters(theProjectile.GetDefinition(),
                                          aNucleus.GetA_asInt(),aNucleus.GetZ_asInt(),
                                          s);
//      theParameters = new G4FTFParameters(theProjectile.GetDefinition(),
//                                          aNucleus.GetN(),aNucleus.GetZ(),
//                                          s);

//theParameters->SetProbabilityOfElasticScatt(0.); 
//G4cout<<theParameters->GetProbabilityOfElasticScatt()<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
// To turn on/off (1/0) elastic scattering

}

// ------------------------------------------------------------
struct DeleteVSplitableHadron { void operator()(G4VSplitableHadron * aH){ delete aH;} };


// ------------------------------------------------------------
G4ExcitedStringVector * G4FTFModel::GetStrings()
{ 
        G4ExcitedStringVector * theStrings(0);
//G4cout<<"GetString"<<G4endl;
	theParticipants.GetList(theProjectile,theParameters);
//G4cout<<"Reggeon"<<G4endl;
        ReggeonCascade(); 

        G4bool Success(true);
        if( PutOnMassShell() )
        {
//G4cout<<"PutOn mass Shell OK"<<G4endl;
         if( ExciteParticipants() )
         {
//G4cout<<"Excite partic OK"<<G4endl;
	  theStrings = BuildStrings();
//G4cout<<"Build String OK"<<G4endl;
          GetResidualNucleus();

          if( theParameters != 0 )
          {
           delete theParameters;
           theParameters=0;
          }
         } else                      // if( ExciteParticipants() )
         {     Success=false;}
        } else                       // if( PutOnMassShell() )
        {      Success=false;}

        if(!Success)   
        {
           // -------------- Erase the projectile ----------------
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
// -------------- Erase the target nucleons -----------
        G4VSplitableHadron * aNucleon = 0;
        for(G4int i=0; i < NumberOfInvolvedNucleon; i++)
        {
           aNucleon = TheInvolvedNucleon[i]->GetSplitableHadron();
           if(aNucleon) delete aNucleon;
        } 

        NumberOfInvolvedNucleon=0;

        return theStrings;

}
//-------------------------------------------------------------------
void G4FTFModel::ReggeonCascade()                             
{ //--- Implementation of reggeon theory inspired model-------
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
// -------------- Properties of the projectile ----------------
	theParticipants.StartLoop();    // restart a loop
        theParticipants.Next();
	G4VSplitableHadron * primary = theParticipants.GetInteraction().GetProjectile();
	G4LorentzVector Pprojectile=primary->Get4Momentum();

//G4cout<<"Pprojectile "<<Pprojectile<<G4endl;
// To get original projectile particle

        if(Pprojectile.z() < 0.){return false;}

        G4double Mprojectile  = Pprojectile.mag();
        G4double M2projectile = Pprojectile.mag2();
//-------------------------------------------------------------
	G4LorentzVector Psum      = Pprojectile;
        G4double        SumMasses = Mprojectile + 20.*MeV; // 13.12.09
                                               // Separation energy for projectile
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

        theNucleus->StartLoop();

	while ((aNucleon = theNucleus->GetNextNucleon()))
        {
         if(aNucleon->AreYouHit())
         {   // Involved nucleons
          Psum += aNucleon->Get4Momentum();
          SumMasses += aNucleon->GetDefinition()->GetPDGMass();
          SumMasses += 20.*MeV;   // 13.12.09 Separation energy for a nucleon
//G4cout<<"SumMasses Tr "<<SumMasses<<G4endl;
          ResidualMassNumber--;
          ResidualCharge-=(G4int) aNucleon->GetDefinition()->GetPDGCharge();
          ResidualExcitationEnergy+=ExcitationEnergyPerWoundedNucleon;
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

               G4LorentzVector tmp(tmpPt.x(),tmpPt.y(),Xminus,0.);
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

               M2target +=(aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass()*
                           aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass()  + 
                           Px*Px + Py*Py)/Xminus;

               G4LorentzVector tmp(Px,Py,Xminus,0.);
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

           G4double Mt2 = sqr(tmp.x())+sqr(tmp.y())+
                          aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass()*
                          aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass();
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

        Pprojectile.transform(toLab);       // The work with the projectile
        primary->Set4Momentum(Pprojectile); // is finished at the moment.
//G4cout<<"Final proj mom "<<primary->Get4Momentum()<<G4endl;

//-------------------------------------------------------------
        G4ThreeVector Residual3Momentum(0.,0.,1.);

	for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
        {
           G4Nucleon * aNucleon = TheInvolvedNucleon[i];
           G4LorentzVector tmp=aNucleon->Get4Momentum();
           Residual3Momentum-=tmp.vect();

           G4double Mt2 = sqr(tmp.x())+sqr(tmp.y())+
                          aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass()*
                          aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass();
           G4double Xminus=tmp.z();

           G4double Pz=-WminusTarget*Xminus/2. + Mt2/(2.*WminusTarget*Xminus);
           G4double E = WminusTarget*Xminus/2. + Mt2/(2.*WminusTarget*Xminus);

           tmp.setPz(Pz); 
           tmp.setE(E);

           tmp.transform(toLab);

           aNucleon->SetMomentum(tmp);

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
//G4cout<<"Residual4Momentum "<<Residual4Momentum<<G4endl;
        Residual4Momentum.transform(toLab);
//-------------------------------------------------------------
 return true;
}

// ------------------------------------------------------------
G4bool G4FTFModel::ExciteParticipants()
{
    G4bool Successfull(false);
//    do {                           // } while (Successfull == false) // Closed 15.12.09
        Successfull=false;
        theParticipants.StartLoop();

G4int MaxNumOfInelCollisions=G4int(theParameters->GetMaxNumberOfCollisions());
G4double NumberOfInel(0.);
//
if(MaxNumOfInelCollisions > 0)  
{   //  Plab > Pbound, Normal application of FTF is possible
 G4double ProbMaxNumber=theParameters->GetMaxNumberOfCollisions()-MaxNumOfInelCollisions;
 if(G4UniformRand() < ProbMaxNumber) {MaxNumOfInelCollisions++;}
 NumberOfInel=MaxNumOfInelCollisions;
} else
{   //  Plab < Pbound, Normal application of FTF is impossible, low energy corrections
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
	while (theParticipants.Next())
	{	   
	   const G4InteractionContent & collision=theParticipants.GetInteraction();

	   G4VSplitableHadron * projectile=collision.GetProjectile();
	   G4VSplitableHadron * target=collision.GetTarget();
//G4cout<<"ProbabilityOfElasticScatt "<<theParameters->GetProbabilityOfElasticScatt()<<G4endl;
           if(G4UniformRand()< theParameters->GetProbabilityOfElasticScatt())
           { //   Elastic scattering -------------------------
//G4cout<<"Elastic FTF"<<G4endl;
            if(theElastic->ElasticScattering(projectile, target, theParameters))
            {
            Successfull = Successfull || true;
            } else
            {
             Successfull = Successfull || false;
             target->SetStatus(2);
            }
           }
           else
           { //   Inelastic scattering ---------------------- 
/*
            if(theExcitation->ExciteParticipants(projectile, target, 
                                                 theParameters, theElastic))
            {
             Successfull = Successfull || true; 
            } else
            {
             Successfull = Successfull || false;
             target->SetStatus(2);
            }
*/
//G4cout<<"InElastic FTF"<<G4endl;
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
              target->SetStatus(2);
             }
            } else // If NumOfInel
            {
             if(theElastic->ElasticScattering(projectile, target, theParameters))
             {
              Successfull = Successfull || true;
             } else
             {
              Successfull = Successfull || false;
              target->SetStatus(2);
             }
            }   // end if NumOfInel
           }
        }       // end of while (theParticipants.Next())
//       } while (Successfull == false);                        // Closed 15.12.09
	return Successfull;
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
	for ( ahadron=0; ahadron < primaries.size() ; ahadron++)
	{
            G4bool isProjectile(0);
            if(primaries[ahadron]->GetStatus() == 1) {isProjectile=true; }
            if(primaries[ahadron]->GetStatus() == 3) {isProjectile=false;}

            FirstString=0; SecondString=0;
            theExcitation->CreateStrings(primaries[ahadron], isProjectile,
                                         FirstString, SecondString,
                                         theParameters);

	    if(FirstString  != 0) strings->push_back(FirstString);
            if(SecondString != 0) strings->push_back(SecondString);
	}
//
	for (G4int ahadron=0; ahadron < NumberOfInvolvedNucleon ; ahadron++)
	{
            if(TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetStatus() !=0) //== 2)
            {
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
	}

	std::for_each(primaries.begin(), primaries.end(), DeleteVSplitableHadron());
	primaries.clear();
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
