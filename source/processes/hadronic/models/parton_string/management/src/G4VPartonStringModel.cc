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
// $Id: G4VPartonStringModel.cc,v 1.7 2009-12-06 11:29:33 vuzhinsk Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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


G4VPartonStringModel::G4VPartonStringModel()
{
//  Make shure Shotrylived partyicles are constructed.
	G4ShortLivedConstructor ShortLived;
	ShortLived.ConstructParticle();
}

G4VPartonStringModel::G4VPartonStringModel(const G4VPartonStringModel &) : G4VHighEnergyGenerator()
{
}


G4VPartonStringModel::~G4VPartonStringModel()
{
}


const G4VPartonStringModel & G4VPartonStringModel::operator=(const G4VPartonStringModel &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4VPartonStringModel::operator= meant to not be accessable");
  return *this;
}


int G4VPartonStringModel::operator==(const G4VPartonStringModel &) const
{
 return 0;
}

int G4VPartonStringModel::operator!=(const G4VPartonStringModel &) const
{
  return 1;
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
  
  G4KineticTrackVector * theResult = 0;
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

#ifdef debug_PartonStringModel
  G4V3DNucleus * fancynucleus=theThis->GetWoundedNucleus();
  
       // loop over wounded nucleus
     G4int hits(0);
     G4Nucleon * theCurrentNucleon = fancynucleus->StartLoop() ? fancynucleus->GetNextNucleon() : NULL;
     while(theCurrentNucleon != NULL)
     {
       if(theCurrentNucleon->AreYouHit()) 
       {
         hits++;
       }
       theCurrentNucleon = fancynucleus->GetNextNucleon();
     }
     
     G4cout << " strE, nucleons, inE " 
            << stringEnergy << " "    
	    << hits << " "
	    << Ptmp.e() << " " 
	    << stringEnergy - 939.*hits - Ptmp.e()<< G4endl;
#endif
  G4double SumMass(0.); 
  attempts = 0; 
  maxAttempts=100;
  do 
  {	    
   attempts++;   
   if(theResult != 0)
   {
    std::for_each(theResult->begin(), theResult->end(), DeleteKineticTrack());
    theResult->clear();
   }

   theResult = stringFragmentationModel->FragmentStrings(strings);

   if(attempts > maxAttempts ) break;

   SumMass=0.;
   for ( unsigned int i=0; i < theResult->size(); i++)
   {
    SumMass+=(*theResult)[i]->GetDefinition()->GetPDGMass();
   }
  } while(SumMass > InvMass);

  std::for_each(strings->begin(), strings->end(), DeleteString() );
  delete strings;

  return theResult;
}

