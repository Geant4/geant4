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
//
// $Id: G4EvaporationFactory.cc,v 1.1 2002/06/06 18:00:31 larazb Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4EvaporationFactory.hh"


#include "G4NeutronEvaporationChannel.hh"
#include "G4ProtonEvaporationChannel.hh"
#include "G4DeuteronEvaporationChannel.hh"
#include "G4TritonEvaporationChannel.hh"
#include "G4He3EvaporationChannel.hh"
#include "G4AlphaEvaporationChannel.hh"

#include "G4CompetitiveFission.hh"
#include "G4PhotonEvaporation.hh"


const G4EvaporationFactory & 
G4EvaporationFactory::operator=(const G4EvaporationFactory & val)
{
  G4Exception("G4EvaporationFactory::operator= meant to not be accessable.");
  return *this;
}

G4bool 
G4EvaporationFactory::operator==(const G4EvaporationFactory & val) const
{
  G4Exception("G4EvaporationFactory::operator== meant to not be accessable.");
  return false;
}

G4bool 
G4EvaporationFactory::operator!=(const G4EvaporationFactory & val) const
{
  G4Exception("G4EvaporationFactory::operator!= meant to not be accessable.");
  return true;
}



G4std::vector<G4VEvaporationChannel*> * 
G4EvaporationFactory::CreateChannel()
{
  G4std::vector<G4VEvaporationChannel*> * theChannel = 
    new G4std::vector<G4VEvaporationChannel*>;
  theChannel->reserve(8);

  theChannel->push_back( new G4NeutronEvaporationChannel() );  // n
  theChannel->push_back( new G4ProtonEvaporationChannel() );   // p
  theChannel->push_back( new G4DeuteronEvaporationChannel() ); // Deuteron
  theChannel->push_back( new G4TritonEvaporationChannel() );   // Triton
  theChannel->push_back( new G4He3EvaporationChannel() );      // He3
  theChannel->push_back( new G4AlphaEvaporationChannel() );    // Alpha

  theChannel->push_back( new G4CompetitiveFission() );         // Fission Channel
  theChannel->push_back( new G4PhotonEvaporation() );          // Photon Channel

  return theChannel;

}



