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
#include "G4QGSModel.hh"
#include "G4Proton.hh"
#include "G4DynamicParticle.hh"

main()
{
  G4QGSModel theModel;
  G4Nucleus theNuc(63., 29.);
  G4ParticleMomentum theDirection( 0., 0., 1. );
  G4double incomingEnergy = 200*GeV;
  G4DynamicParticle aParticle(G4Proton::ProtonDefinition(), 
                              theDirection, 
			      incomingEnergy );
  G4int counter = 0;
  for(;;)
  {
    counter++;
    theModel.Init(theNuc, aParticle);
    G4ExcitedStringVector * result = theModel.GetStrings();
    G4cout << "Event number "<<counter<<G4endl;
    for(G4int i=0; i<result->length(); i++)
    {
      delete result->at(i);
    }
    delete result;
  }
}
