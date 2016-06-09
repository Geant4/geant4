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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4EvaporationFactory.cc,v 1.2 2003/11/03 17:53:02 hpw Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
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
G4EvaporationFactory::operator=(const G4EvaporationFactory & )
{
  throw G4HadronicException(__FILE__, __LINE__, "G4EvaporationFactory::operator= meant to not be accessable.");
  return *this;
}

G4bool 
G4EvaporationFactory::operator==(const G4EvaporationFactory & ) const
{
  throw G4HadronicException(__FILE__, __LINE__, "G4EvaporationFactory::operator== meant to not be accessable.");
  return false;
}

G4bool 
G4EvaporationFactory::operator!=(const G4EvaporationFactory & ) const
{
  throw G4HadronicException(__FILE__, __LINE__, "G4EvaporationFactory::operator!= meant to not be accessable.");
  return true;
}



std::vector<G4VEvaporationChannel*> * 
G4EvaporationFactory::CreateChannel()
{
  std::vector<G4VEvaporationChannel*> * theChannel = 
    new std::vector<G4VEvaporationChannel*>;
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



