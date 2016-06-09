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
//
//// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      ---------------- G4VPartonStringModel ----------------
//             by Gunter Folger, May 1998.
//      abstract class for all Parton String Models
// ------------------------------------------------------------
// debug switch
//#define debug_PartonStringModel


#include "G4VPartonStringModel.hh"
#include "G4ios.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"


G4VPartonStringModel::G4VPartonStringModel(const G4String& modelName)
    : G4VHighEnergyGenerator(modelName),
      stringFragmentationModel(0),
      theThis(0)
{
//  Make shure Shotrylived partyicles are constructed.
	G4ShortLivedConstructor ShortLived;
	ShortLived.ConstructParticle();
}

G4VPartonStringModel::~G4VPartonStringModel()
{
}

G4KineticTrackVector * G4VPartonStringModel::Scatter(const G4Nucleus &theNucleus, 
                                                const G4DynamicParticle &aPrimary)
{  
  G4ExcitedStringVector * strings = NULL;

  G4DynamicParticle thePrimary=aPrimary;
  
  G4LorentzRotation toZ;
  G4LorentzVector Ptmp=thePrimary.Get4Momentum();
  toZ.rotateZ(-1*Ptmp.phi());
  toZ.rotateY(-1*Ptmp.theta());
  thePrimary.Set4Momentum(toZ*Ptmp);
  G4LorentzRotation toLab(toZ.inverse());

  G4int attempts = 0, maxAttempts=20;
  while ( strings  == NULL )
  {
  	if (attempts++ > maxAttempts ) 
  	{
		throw G4HadronicException(__FILE__, __LINE__, "G4VPartonStringModel::Scatter(): fails to generate strings");
  	}

	theThis->Init(theNucleus,thePrimary);

  	strings = GetStrings();
  }
  
  G4double stringEnergy(0);
  G4LorentzVector SumStringMom(0.,0.,0.,0.);

  for ( unsigned int astring=0; astring < strings->size(); astring++)
  {
//    rotate string to lab frame, models have it aligned to z
    stringEnergy += (*strings)[astring]->GetLeftParton()->Get4Momentum().t();
    stringEnergy += (*strings)[astring]->GetRightParton()->Get4Momentum().t();
    (*strings)[astring]->LorentzRotate(toLab);
    SumStringMom+=(*strings)[astring]->Get4Momentum();
  }

  G4double InvMass=SumStringMom.mag();   

//#define debug_PartonStringModel
#ifdef debug_PartonStringModel

     G4V3DNucleus * fancynucleus=theThis->GetWoundedNucleus();
  
       // loop over wounded nucleus
     G4int hits(0), charged_hits(0);
     G4ThreeVector hitNucleonMomentum(0.,0.,0.);
     G4Nucleon * theCurrentNucleon = fancynucleus->StartLoop() ? fancynucleus->GetNextNucleon() : NULL;
     while( theCurrentNucleon )
     {
       if(theCurrentNucleon->AreYouHit()) 
       {
         hits++;
         hitNucleonMomentum += theCurrentNucleon->Get4Momentum().vect();
         if ( theCurrentNucleon->GetDefinition() == G4Proton::Proton() )  ++charged_hits;
       }
       theCurrentNucleon = fancynucleus->GetNextNucleon();
     }
     
     G4int initialZ=fancynucleus->GetCharge();
     G4int initialA=fancynucleus->GetMassNumber();
     G4double initial_mass=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(initialZ,initialA);
     G4double final_mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(initialZ-charged_hits, initialA-hits);
     G4cout << "G4VPSM: strE, hit nucleons, Primary, SumStringE, nucleus intial, nucleus final, excitation estimate "
            << stringEnergy << " "    << hits << ", " << Ptmp.e() << ", "<< SumStringMom.e() << ", "
	        << initial_mass<< ", " << final_mass<< ", "  << Ptmp.e() + initial_mass - final_mass - stringEnergy  << G4endl;
     G4cout << "momentum balance " <<  thePrimary.GetMomentum() + hitNucleonMomentum - SumStringMom.vect()<< G4endl;
#endif

//  Fragment strings

  G4KineticTrackVector * theResult = 0;
  G4double SumMass(0.); 
  attempts = 0; 
  maxAttempts=100;
  do 
  {	    
   attempts++;   
   if(theResult != 0)
   {
    std::for_each(theResult->begin(), theResult->end(), DeleteKineticTrack());
    delete theResult;
   }
   theResult = stringFragmentationModel->FragmentStrings(strings);
   if(attempts > maxAttempts ) break;

    //G4cout<<"G4endl<<"G4VPartonStringModel:: Final Result, Size "<<theResult->size()<<G4endl;

   SumMass=0.;
    //G4LorentzVector SumP(0.,0.,0.,0.);
   for ( unsigned int i=0; i < theResult->size(); i++)
   {
    SumMass+=(*theResult)[i]->GetDefinition()->GetPDGMass();
     //SumP+=(*theResult)[i]->Get4Momentum();
     //G4cout<<i<<" : "<<(*theResult)[i]->GetDefinition()->GetParticleName();
     //G4cout<<"p= " << (*theResult)[i]->Get4Momentum()<<" m= "<<(*theResult)[i]->Get4Momentum().mag()<<G4endl;
   }

     //G4cout<<"SumP "<<SumP<<G4endl;
  } while(SumMass > InvMass);

  std::for_each(strings->begin(), strings->end(), DeleteString() );
  delete strings;

  return theResult;
}

void G4VPartonStringModel::ModelDescription(std::ostream& outFile) const
{
	outFile << GetModelName() << " has no description yet.\n";
}
