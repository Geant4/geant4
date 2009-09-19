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
// $Id: G4FTFModel.cc,v 1.25 2009-09-19 19:02:48 vuzhinsk Exp $
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

// Class G4FTFModel 

G4FTFModel::G4FTFModel():theExcitation(new G4DiffractiveExcitation()),
                         theElastic(new G4ElasticHNScattering())
{
	G4VPartonStringModel::SetThisPointer(this);
        theParameters=0;
}


G4FTFModel::~G4FTFModel()
{
   if( theParameters != 0 ) delete theParameters; 
// Because FTF model can be called for various particles
// theParameters must be erased at the end of each call.
// Thus the delete is also in G4FTFModel::GetStrings() method
   if( theExcitation != 0 ) delete theExcitation;
   if( theElastic    != 0 ) delete theElastic; 
}

const G4FTFModel & G4FTFModel::operator=(const G4FTFModel &)
{
	throw G4HadronicException(__FILE__, __LINE__, "G4FTFModel::operator= is not meant to be accessed ");
	return *this;
}
#include "G4ParticleTable.hh"

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
	theParticipants.Init(aNucleus.GetN(),aNucleus.GetZ()); 
// ----------- N-mass number Z-charge -------------------------

// --- cms energy

        G4double s = sqr( theProjectile.GetMass() ) +
                     sqr( G4Proton::Proton()->GetPDGMass() ) +
                     2*theProjectile.GetTotalEnergy()*G4Proton::Proton()->GetPDGMass();
/*
G4cout << " primary Total E (GeV): " << theProjectile.GetTotalEnergy()/GeV << G4endl;
G4cout << " primary Mass    (GeV): " << theProjectile.GetMass() /GeV << G4endl;
G4cout << " primary 3Mom           " << theProjectile.GetMomentum() << G4endl;
G4cout << " primary space position " << theProjectile.GetPositionInNucleus() << G4endl;
G4cout << "cms std::sqrt(s) (GeV) = " << std::sqrt(s) / GeV << G4endl;
*/

      if( theParameters != 0 ) delete theParameters;
      theParameters = new G4FTFParameters(theProjectile.GetDefinition(),
                                          aNucleus.GetN(),aNucleus.GetZ(),
                                          s);
//theParameters->SetProbabilityOfElasticScatt(0.); 
// To turn on/off (1/0) elastic scattering

}

// ------------------------------------------------------------
G4ExcitedStringVector * G4FTFModel::GetStrings()
{ 
//G4cout<<"GetList"<<G4endl;
	theParticipants.GetList(theProjectile,theParameters);
//G4cout<<"Regge"<<G4endl;
        ReggeonCascade();                                     // Uzhi 26 July 09 
//G4cout<<"On mass shell"<<G4endl;
        if (! PutOnMassShell()    ) return NULL;              // Uzhi 26 July 09
//G4cout<<"Excite"<<G4endl;
	if (! ExciteParticipants()) return NULL;
//G4cout<<"Strings"<<G4endl;
	G4ExcitedStringVector * theStrings = BuildStrings();
//G4cout<<"Out FTF N strings "<<theStrings->size()<<G4endl;
        if( theParameters != 0 )
        {
          delete theParameters;
          theParameters=0;
        }
	return theStrings;
}

// ------------------------------------------------------------
struct DeleteVSplitableHadron { void operator()(G4VSplitableHadron * aH){ delete aH;} };

// ------------------------------------------------------------
void G4FTFModel::ReggeonCascade()                             // Uzhi 26 July 2009
{ //--- Implementation of reggeon theory inspired model-------
        NumberOfInvolvedNucleon=0;

//G4int PrimInt(0);
        theParticipants.StartLoop();
	while (theParticipants.Next())
	{
//PrimInt++; 	   
	   const G4InteractionContent & collision=theParticipants.GetInteraction();
           G4Nucleon * TargetNucleon=collision.GetTargetNucleon();
//G4cout<<"Prim Nnucl "<<TargetNucleon->Get4Momentum()<<G4endl;
           TheInvolvedNucleon[NumberOfInvolvedNucleon]=TargetNucleon;
           NumberOfInvolvedNucleon++;

           G4double XofWoundedNucleon = TargetNucleon->GetPosition().x();
           G4double YofWoundedNucleon = TargetNucleon->GetPosition().y();

           theParticipants.theNucleus->StartLoop();
           G4Nucleon * Neighbour;
//G4int NucleonNum(0);
	   while ( (Neighbour = theParticipants.theNucleus->GetNextNucleon()) )
           {
            if(!Neighbour->AreYouHit())
            {
    	     G4double impact2= sqr(XofWoundedNucleon - Neighbour->GetPosition().x()) +
                               sqr(YofWoundedNucleon - Neighbour->GetPosition().y());

             if(G4UniformRand() < theParameters->GetCofNuclearDestruction()*
                std::exp(-impact2/theParameters->GetR2ofNuclearDestruction()))  
             { // The neighbour nucleon is involved in the reggeon cascade
//G4cout<<" involved "<<NucleonNum<<" "<<Neighbour->Get4Momentum()<<G4endl;
              TheInvolvedNucleon[NumberOfInvolvedNucleon]=Neighbour;
              NumberOfInvolvedNucleon++;
//G4cout<<" PrimInt"<<" "<<NumberOfInvolvedNucleon<<G4endl;

              G4VSplitableHadron *targetSplitable; 
              targetSplitable = new G4DiffractiveSplitableHadron(*Neighbour); 

              Neighbour->Hit(targetSplitable);
              targetSplitable->SetStatus(2);     
             }
            }  // end of if(!Neighbour->AreYouHit())
//NucleonNum++;
	   }   // end of while (theParticipant.theNucleus->GetNextNucleon())
//G4cout<<"Prim Int N Ninvolv "<<PrimInt<<" "<<NumberOfInvolvedNucleon<<G4endl;
	}      // end of while (theParticipants.Next())
//G4cout<<"At end "<<PrimInt<<" "<<NumberOfInvolvedNucleon<<G4endl;

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
//G4cout<<"Put on Mass-shell"<<G4endl;
	theParticipants.StartLoop();    // restart a loop
        theParticipants.Next();
	G4VSplitableHadron * primary = theParticipants.GetInteraction().GetProjectile();
	G4LorentzVector Pprojectile=primary->Get4Momentum();
//G4cout<<"Proj "<<Pprojectile<<G4endl;
        if(Pprojectile.z() < 0.){return false;}

        G4double Mprojectile  = Pprojectile.mag();
        G4double M2projectile = Pprojectile.mag2();
//-------------------------------------------------------------
	G4LorentzVector Psum      = Pprojectile;
        G4double        SumMasses = Mprojectile;

//--------------- Target nucleus ------------------------------
G4int Ninvolv(0);
	for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
        {
          G4Nucleon * aNucleon = TheInvolvedNucleon[i];
//G4cout<<"          "<<Ninvolv<<" "<<aNucleon->Get4Momentum()<<G4endl;
Ninvolv++;
          Psum += aNucleon->Get4Momentum();
          SumMasses += aNucleon->GetDefinition()->GetPDGMass();                
	}   // end of for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
//G4cout<<"Ninvolv "<<Ninvolv<<G4endl;
//-------------------------------------------------------------

        G4double SqrtS=Psum.mag();
        G4double     S=Psum.mag2();
//G4cout<<" SqrtS SumMasses "<<SqrtS <<" "<< SumMasses<<G4endl;
        if(SqrtS < SumMasses)      {return false;} // It is impossible to simulate
                                                   // after putting nuclear nucleons
                                                   // on mass-shell
//-------------------------------------------------------------
// Sampling of nucleons what are transfered to delta-isobars --
        G4int MaxNumberOfDeltas = (int)((SqrtS - SumMasses)/(400.*MeV));
        G4int NumberOfDeltas(0);
SumMasses=Mprojectile;
        if(GetWoundedNucleus()->GetMassNumber() != 1)
        {
          G4double ProbDeltaIsobar(1.);
	  for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
          {
            if((G4UniformRand() < ProbDeltaIsobar)&&(NumberOfDeltas < MaxNumberOfDeltas))
            {
              NumberOfDeltas++;
              G4VSplitableHadron * targetSplitable=TheInvolvedNucleon[i]->GetSplitableHadron();
              G4int PDGcode = targetSplitable->GetDefinition()->GetPDGEncoding();
              G4int newPDGcode = PDGcode/10; newPDGcode=newPDGcode*10+4; // Delta
              G4ParticleDefinition* ptr = 
                 G4ParticleTable::GetParticleTable()->FindParticle(newPDGcode);
              targetSplitable->SetDefinition(ptr);
SumMasses+=targetSplitable->GetDefinition()->GetPDGMass();
//G4cout<<" i "<<i<<" Delta "<<targetSplitable->GetDefinition()->GetPDGMass()<<G4endl;
            } else {
SumMasses+=TheInvolvedNucleon[i]->GetSplitableHadron()->GetDefinition()->GetPDGMass();
//G4cout<<" i "<<i<<" Nuclo "<<TheInvolvedNucleon[i]->GetSplitableHadron()->GetDefinition()->GetPDGMass()<<G4endl;
}

          }   // end of for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
        }   // end of if(GetWoundedNucleus().GetN() != 1)
//G4cout<<"MaxNumberOfDeltas NumberOfDeltas "<<MaxNumberOfDeltas<<" "<<NumberOfDeltas<<G4endl;
//G4cout<<" SqrtS SumMasses "<<SqrtS <<" "<< SumMasses<<G4endl;
//-------------------------------------------------------------

        G4LorentzRotation toCms(-1*Psum.boostVector());
        G4LorentzVector Ptmp=toCms*Pprojectile;

        if ( Ptmp.pz() <= 0. )                                
        {  // "String" moving backwards in  CMS, abort collision !!
           //G4cout << " abort ColliDeleteVSplitableHadronsion!! " << G4endl;
         return false; 
        }

        toCms.rotateZ(-1*Ptmp.phi());
        toCms.rotateY(-1*Ptmp.theta());
	
        G4LorentzRotation toLab(toCms.inverse());

//        Pprojectile.transform(toCms);
//G4cout<<"Proj in CMS "<<Pprojectile<<G4endl;

//G4cout<<"Main work"<<G4endl;
//-------------------------------------------------------------
//------- Ascribing of the involved nucleons Pt and Xminus ----
        G4double Dcor        = theParameters->GetDofNuclearDestruction()/
                                                  NumberOfInvolvedNucleon;
        G4double AveragePt2  = theParameters->GetPt2ofNuclearDestruction();
        G4double maxPtSquare = theParameters->GetMaxPt2ofNuclearDestruction();
//G4cout<<" D Pt2 "<<Dcor<<" "<<AveragePt2<<G4endl;
        G4double M2target;
        G4int NumberOfTries=0;
        G4double ScaleFactor=1.;
        do    // while (SqrtS < Mprojectile + std::sqrt(M2target))
        {     // while (DecayMomentum < 0.)

          NumberOfTries++;
//G4cout<<"NumberOfTries "<<NumberOfTries<<G4endl;
          if(NumberOfTries == 100*(NumberOfTries/100))
          { // At large number of tries it would be better to reduce the values
            ScaleFactor/=2.;
            Dcor       *=ScaleFactor;
            AveragePt2 *=ScaleFactor;
//G4cout<<"NumberOfTries "<<NumberOfTries<<" "<<Dcor<<" "<<AveragePt2<<G4endl;
//G4int Uzhi; G4cin>>Uzhi;
          }
//G4cout<<"Start Decay "<<G4endl; G4int Uzhi; G4cin>>Uzhi;
          G4ThreeVector PtSum;
          G4double XminusSum(0.);
          G4double Xminus(0.);

          do      // while((Xminus <= 0.) || (Xminus >= 1.));
          {
//G4cout<<"Sample Pt and X"<<G4endl; // G4int Uzhi1; G4cin>>Uzhi1;
             PtSum    =G4ThreeVector(0.,0.,0.);
             XminusSum=0.;
             M2target =0.;

	     for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
             {
               G4Nucleon * aNucleon = TheInvolvedNucleon[i];

               G4ThreeVector tmpPt;
               tmpPt = GaussianPt(AveragePt2, maxPtSquare);
               PtSum += tmpPt;

               G4ThreeVector tmpX;
               tmpX =GaussianPt(Dcor*Dcor, 1.);
               Xminus=tmpX.x();
               XminusSum+=Xminus;

               G4LorentzVector tmp(tmpPt.x(),tmpPt.y(),Xminus,0.);
               aNucleon->SetMomentum(tmp);
//G4cout<<"Nucl mom "<<aNucleon->GetMomentum()<<G4endl;
             }   // end of for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )

//---------------------------------------------------------------------------
             G4double DeltaX      = PtSum.x()/NumberOfInvolvedNucleon;
             G4double DeltaY      = PtSum.y()/NumberOfInvolvedNucleon;
             G4double DeltaXminus = (XminusSum-1.)/NumberOfInvolvedNucleon;
//G4cout<<"Deltas "<<DeltaX<<" "<<DeltaY<<" "<<DeltaXminus<<G4endl;
             M2target =0.;

	     for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )
             {
               G4Nucleon * aNucleon = TheInvolvedNucleon[i];

               Xminus = aNucleon->Get4Momentum().pz() - DeltaXminus;
               if((Xminus <= 0.) || (Xminus > 1.)) {//G4cout<<"-X break "<<Xminus<<G4endl; 
                                                    break;}
 
               G4double Px=aNucleon->Get4Momentum().px() - DeltaX;
               G4double Py=aNucleon->Get4Momentum().py() - DeltaY;

               M2target +=(aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass()*
                           aNucleon->GetSplitableHadron()->GetDefinition()->GetPDGMass()  + 
                           Px*Px + Py*Py)/Xminus;

               G4LorentzVector tmp(Px,Py,Xminus,0.);
               aNucleon->SetMomentum(tmp);
             }   // end of for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )

          } while((Xminus <= 0.) || (Xminus >= 1.));

//-------------------------------------------------------------
//G4cout<<"SqrtS > Mprojectile + std::sqrt(M2target) "<<SqrtS<<" "<<Mprojectile<<" "<< std::sqrt(M2target)<<" "<<Mprojectile + std::sqrt(M2target)<<G4endl;
//G4int Uzhi3; G4cin>>Uzhi3;
        } while (SqrtS < Mprojectile + std::sqrt(M2target));

//-------------------------------------------------------------
        G4double DecayMomentum2= S*S+M2projectile*M2projectile+M2target*M2target
                                -2.*S*M2projectile - 2.*S*M2target 
                                       -2.*M2projectile*M2target;

        G4double WminusTarget=(S-M2projectile+M2target+std::sqrt(DecayMomentum2))/2./SqrtS;
        G4double WplusProjectile=SqrtS - M2target/WminusTarget;
//G4cout<<"DecM W+ W- "<<DecayMomentum2<<" "<<WplusProjectile<<" "<<WminusTarget<<G4endl;
//-------------------------------------------------------------
        G4double Pzprojectile=WplusProjectile/2. - M2projectile/2./WplusProjectile;
        G4double Eprojectile =WplusProjectile/2. + M2projectile/2./WplusProjectile;
        Pprojectile.setPz(Pzprojectile);  Pprojectile.setE(Eprojectile);

        Pprojectile.transform(toLab);       // The work with the projectile
        primary->Set4Momentum(Pprojectile); // is finished at the moment.
//G4cout<<"Proj After "<<Pprojectile<<G4endl;
//-------------------------------------------------------------
//Ninvolv=0;
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

           tmp.setPz(Pz); 
           tmp.setE(E);
           tmp.transform(toLab);
//G4cout<<"invol  "<<Ninvolv<<" "<<tmp<<G4endl;
//Ninvolv++;
           aNucleon->SetMomentum(tmp);
           G4VSplitableHadron * targetSplitable=aNucleon->GetSplitableHadron();
           targetSplitable->Set4Momentum(tmp);
//G4cout<<"nucleon M "<<aNucleon->Get4Momentum()<<G4endl;
        }   // end of for(G4int i=0; i < NumberOfInvolvedNucleon; i++ )

//G4cout<<"Return Nucleus"<<G4endl;
//-------------------------------------------------------------
//-------------------------------------------------------------
//-------------------------------------------------------------
//G4int Uzhi2; G4cin>>Uzhi2;

 return true;
}

// ------------------------------------------------------------
G4bool G4FTFModel::ExciteParticipants()
{
/*    // Uzhi 29.03.08                     For elastic Scatt.
G4cout<<"  In ExciteParticipants() "<<theParticipants.theInteractions.size()<<G4endl;
G4cout<<" test Params Tot "<<theParameters->GetTotalCrossSection()<<G4endl;
G4cout<<" test Params Ela "<<theParameters->GetElasticCrossSection()<<G4endl;
	
G4int counter=0;
*/   // Uzhi 29.03.08


//G4int InterNumber=0; // Vova

        G4bool Successfull=false;
        theParticipants.StartLoop();                          //Uzhi 26 July 09

//	while (theParticipants.Next()&& (InterNumber < 3)) // Vova
	while (theParticipants.Next())
	{	   
	   const G4InteractionContent & collision=theParticipants.GetInteraction();
/*
counter++;
G4cout<<" Int num "<<counter<<G4endl;
*/
	   G4VSplitableHadron * projectile=collision.GetProjectile();
	   G4VSplitableHadron * target=collision.GetTarget();
//         G4Nucleon * TargetNucleon=collision.GetTargetNucleon(); // Uzhi 16.07.09
// Uzhi 16.07.09 ----------------------------
           if(G4UniformRand()< theParameters->GetProbabilityOfElasticScatt())
           { //   Elastic scattering -------------------------
            if(theElastic->ElasticScattering(projectile, target, theParameters))
            {
//G4cout<<"Elastic"<<G4endl;
             Successfull = Successfull || true;
            } else
            {
//G4cout<<"Elastic Not succes"<<G4endl;
             Successfull = Successfull || false;
             target->SetStatus(2);
/*
             if(target->GetStatus() == 0)                         // Uzhi 17.07.09
             {
              G4VSplitableHadron * aHit=0;                        // Uzhi 16.07.09
              TargetNucleon->Hit(aHit);                           // Uzhi 16.07.09
             };
*/
            };
           }
           else
           { //   Inelastic scattering ---------------------- 
            if(theExcitation->ExciteParticipants(projectile, target, 
                                                 theParameters, theElastic))
            {
             Successfull = Successfull || true; 
//G4cout<<"InElastic"<<G4endl;
            } else
            {
//G4cout<<"InElastic Non succes"<<G4endl;
             Successfull = Successfull || false;
             target->SetStatus(2);
/*
             if(target->GetStatus() == 0)                         // Uzhi 16.06.09
             {
              G4VSplitableHadron * aHit=0;                        // Uzhi 16.07.09
              TargetNucleon->Hit(aHit);                           // Uzhi 16.07.09
             };
*/
            };
           }
        }       // end of the loop Uzhi 9.07.09
// Uzhi 16.07.09 ----------------------------

        if(!Successfull)
	{
//G4cout<<"Process not successfull"<<G4endl;
//           give up, clean up
	  std::vector<G4VSplitableHadron *> primaries;
//	  std::vector<G4VSplitableHadron *> targets;        // Uzhi 31.07.09
	  theParticipants.StartLoop();    // restart a loop 
	  while ( theParticipants.Next() ) 
	  {
	    const G4InteractionContent & interaction=theParticipants.GetInteraction();
                	 //  do not allow for duplicates ...
	    if ( primaries.end() == std::find(primaries.begin(), primaries.end(),
                                                   interaction.GetProjectile()) )
	    	primaries.push_back(interaction.GetProjectile());
/*  // Uzhi 31.07.09
	    if ( targets.end()   == std::find(targets.begin(), targets.end(),
                                                      interaction.GetTarget()) ) 
	    	targets.push_back(interaction.GetTarget());
*/  // Uzhi 31.07.09
	  }
	  std::for_each(primaries.begin(), primaries.end(), DeleteVSplitableHadron());
	  primaries.clear();
/*  // Uzhi 31.07.09	
          std::for_each(targets.begin(), targets.end(), DeleteVSplitableHadron());
	  targets.clear();
*/  // Uzhi 31.07.09
//          theParticipants.theNucleus->StartLoop();

//G4cout<<"NumberOfInvolvedNucleon "<<NumberOfInvolvedNucleon<<G4endl;
          G4VSplitableHadron * aNucleon = 0;
          for(G4int i=0; i < NumberOfInvolvedNucleon; i++)
          {
           aNucleon = TheInvolvedNucleon[i]->GetSplitableHadron();
           if(aNucleon)
           { 
             if(aNucleon->GetStatus() != 0) delete aNucleon;
//           if(aNucleon->GetStatus() == 2)  DeleteVSplitableHadron()(aNucleon);
           }
          } 

          NumberOfInvolvedNucleon=0;
//G4cout<<"Out of Excit"<<G4endl; G4int Uzhi; G4cin>>Uzhi;
	  return false;
	}  // End of if(!Successfull)

	return true;
}
// ------------------------------------------------------------
G4ExcitedStringVector * G4FTFModel::BuildStrings()
{	
// Loop over all collisions; find all primaries, and all target ( targets may 
//  be duplicate in the List ( to unique G4VSplitableHadrons)

//G4cout<<"In build string"<<G4endl;

	G4ExcitedStringVector * strings;
	strings = new G4ExcitedStringVector();
	
	std::vector<G4VSplitableHadron *> primaries;
	std::vector<G4VSplitableHadron *> targets;
	std::vector<G4Nucleon          *> TargetNucleons;     // Uzhi 16.07.09
	
        G4ExcitedString * FirstString(0);     // If there will be a kink,
        G4ExcitedString * SecondString(0);    // two strings will be prodused.

	theParticipants.StartLoop();    // restart a loop 
//G4int InterCount(0); // Uzhi
	while ( theParticipants.Next() ) 
	{
	    const G4InteractionContent & interaction=theParticipants.GetInteraction();
                 //  do not allow for duplicates ...

	    if ( primaries.end() == std::find(primaries.begin(), primaries.end(),
                                                interaction.GetProjectile()) )
	    	primaries.push_back(interaction.GetProjectile());
		
	    if ( targets.end()   == std::find(targets.begin(), targets.end(),
                                                interaction.GetTarget()) ) 
	    	targets.push_back(interaction.GetTarget());

	    if ( TargetNucleons.end()   == std::find(TargetNucleons.begin(),     
                                                     TargetNucleons.end(),       
                                                interaction.GetTargetNucleon()) )
	    	TargetNucleons.push_back(interaction.GetTargetNucleon());        
//InterCount++;
	}
	    
	
//G4cout << "BuildStrings prim/targ/woundN " << primaries.size() << " , " <<targets.size() <<", "<<TargetNucleons.size()<< G4endl;

	unsigned int ahadron;
// Only for hA-interactions Uzhi -------------------------------------
//G4int StringN(0);
//G4cout<<"Proj strings -----------------------"<<G4endl;
	for ( ahadron=0; ahadron < primaries.size() ; ahadron++)
	{
            G4bool isProjectile(0);
            if(primaries[ahadron]->GetStatus() == 1) {isProjectile=true; }
            if(primaries[ahadron]->GetStatus() == 3) {isProjectile=false;}

            FirstString=0; SecondString=0;
            theExcitation->CreateStrings(primaries[ahadron], isProjectile,
                                         FirstString, SecondString,
                                         theParameters);
//G4cout<<"1str 2str "<<FirstString<<" "<<SecondString<<G4endl;
	    if(FirstString  != 0) strings->push_back(FirstString);
            if(SecondString != 0) strings->push_back(SecondString);

//StringN++; G4cout<<"Proj string "<<strings->size()<<G4endl;
	}
//G4cout<<"Targ strings ------------------------------"<<G4endl;
	for ( ahadron=0; ahadron < targets.size() ; ahadron++)
	{
//G4cout<<"targets[ahadron]->GetStatus() "<<targets[ahadron]->GetStatus()<<G4endl;
            if(targets[ahadron]->GetStatus() == 1)   // Uzhi 17.07.09
            {
	     G4bool isProjectile=false;
             FirstString=0; SecondString=0;
             theExcitation->CreateStrings(targets[ahadron], isProjectile,
                                          FirstString, SecondString,
                                          theParameters);
	     if(FirstString  != 0) strings->push_back(FirstString);
             if(SecondString != 0) strings->push_back(SecondString);

//StringN++; G4cout<<"Targ string"<<StringN<<G4endl;
            } else
            {
             if(targets[ahadron]->GetStatus() == 0)// Uzhi 17.07.09 Nucleon was rejected
             {
              G4VSplitableHadron * aHit=0;          // Uzhi 16.07.09
              TargetNucleons[ahadron]->Hit(aHit);   // Uzhi 16.07.09
             }
            }
	}

//G4cout<<"Proj + Targ string "<<strings->size()<<G4endl;
//G4cout<<"Involv strings NumberOfInvolvedNucleon "<<NumberOfInvolvedNucleon<<G4endl;
	for (G4int ahadron=0; ahadron < NumberOfInvolvedNucleon ; ahadron++)
	{
/*
G4cout<<"ahadron "<<ahadron<<" Status "<<
TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetStatus()<<
TheInvolvedNucleon[ahadron]->GetSplitableHadron()->Get4Momentum()<<G4endl;
*/
            if(TheInvolvedNucleon[ahadron]->GetSplitableHadron()->GetStatus() == 2)
            {
//StringN++; G4cout<<"Invol string "<<StringN<<G4endl;
	     G4bool isProjectile=false;
             FirstString=0; SecondString=0;
             theExcitation->CreateStrings(
                            TheInvolvedNucleon[ahadron]->GetSplitableHadron(),
                                          isProjectile,
                                          FirstString, SecondString,
                                          theParameters);
	     if(FirstString  != 0) strings->push_back(FirstString);
             if(SecondString != 0) strings->push_back(SecondString);

//	     strings->push_back(theExcitation->String(
//                      TheInvolvedNucleon[ahadron]->GetSplitableHadron(), 
//                                                           isProjectile));
            } 
//G4cout<<"Proj + Targ+Involved string "<<strings->size()<<G4endl;
/*
else
            {
G4cout<<"Else ????????????"<<G4endl; 
             if(targets[ahadron]->GetStatus() == 0)// Uzhi 17.07.09 Nucleon was rejected
             {
              G4VSplitableHadron * aHit=0;          // Uzhi 16.07.09
              TargetNucleons[ahadron]->Hit(aHit);   // Uzhi 16.07.09
             }
            }
*/
	}

//G4int Uzhi; G4cin>>Uzhi;

	std::for_each(primaries.begin(), primaries.end(), DeleteVSplitableHadron());
	primaries.clear();

	std::for_each(targets.begin(), targets.end(), DeleteVSplitableHadron());
	targets.clear();

	return strings;
}
// ------------------------------------------------------------
G4ThreeVector G4FTFModel::GaussianPt(G4double AveragePt2, G4double maxPtSquare) const
{            //  @@ this method is used in FTFModel as well. Should go somewhere common!
	
	G4double Pt2;
        Pt2 = -AveragePt2 * std::log(1. + G4UniformRand() * 
                (std::exp(-maxPtSquare/AveragePt2)-1.)); 
	
	G4double Pt=std::sqrt(Pt2);
	G4double phi=G4UniformRand() * twopi;
	
	return G4ThreeVector (Pt*std::cos(phi), Pt*std::sin(phi), 0.);    
}
