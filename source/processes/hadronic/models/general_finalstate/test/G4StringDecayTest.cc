//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4ExcitedStringVector.hh"
#include "G4Parton.hh"
#include "G4ExcitedString.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4BosonConstructor.hh"

main()
{
  G4ShortLivedConstructor aC;
  G4MesonConstructor bC;
  G4BaryonConstructor cC;
  G4LeptonConstructor dC;
  G4BosonConstructor eC;
  aC.ConstructParticle();
  bC.ConstructParticle();
  cC.ConstructParticle();
  dC.ConstructParticle();
  eC.ConstructParticle();
 
  G4QGSMFragmentation theFragmentation;
  G4ExcitedStringDecay theModel(&theFragmentation);
  G4ExcitedStringVector theStrings;
  G4ThreeVector aPosition( 0., 0., 0. );
  G4LorentzVector aMom(0, 0, 200*GeV, 200*GeV);
  G4LorentzVector bMom(0, 0, -200*GeV, 200*GeV);
  G4Parton p1(1);
  p1.SetPosition(aPosition);
  p1.Set4Momentum(aMom);
  G4Parton p2(-1);
  p2.SetPosition(aPosition);
  p2.Set4Momentum(bMom);
  G4ExcitedString aString(&p1, &p2);
  theStrings.insert(&aString);
  G4int counter = 0;
  for(;;)
  {
    counter++;
    G4KineticTrackVector * result = theModel.FragmentStrings(&theStrings);
    G4cout << "Event number "<<counter<<" "<<result->length()<<G4endl;
    for(G4int i=0; i<result->length(); i++)
    {
      delete result->at(i);
    }
    delete result;
  }
}
