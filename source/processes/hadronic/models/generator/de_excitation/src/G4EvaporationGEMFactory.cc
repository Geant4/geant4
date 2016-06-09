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
// $Id: G4EvaporationGEMFactory.cc,v 1.4 2003/06/16 17:06:22 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4EvaporationGEMFactory.hh"

#include "G4NeutronGEMChannel.hh"
#include "G4ProtonGEMChannel.hh"
#include "G4DeuteronGEMChannel.hh"
#include "G4TritonGEMChannel.hh"
#include "G4He3GEMChannel.hh"
#include "G4AlphaGEMChannel.hh"
#include "G4He6GEMChannel.hh"
#include "G4He8GEMChannel.hh"
#include "G4Li6GEMChannel.hh"
#include "G4Li7GEMChannel.hh"
#include "G4Li8GEMChannel.hh"
#include "G4Li9GEMChannel.hh"
#include "G4Be7GEMChannel.hh"
#include "G4Be9GEMChannel.hh"
#include "G4Be10GEMChannel.hh"
#include "G4Be11GEMChannel.hh"
#include "G4Be12GEMChannel.hh"
#include "G4B8GEMChannel.hh"
#include "G4B10GEMChannel.hh"
#include "G4B11GEMChannel.hh"
#include "G4B12GEMChannel.hh"
#include "G4B13GEMChannel.hh"
#include "G4C10GEMChannel.hh"
#include "G4C11GEMChannel.hh"
#include "G4C12GEMChannel.hh"
#include "G4C13GEMChannel.hh"
#include "G4C14GEMChannel.hh"
#include "G4C15GEMChannel.hh"
#include "G4C16GEMChannel.hh"
#include "G4N12GEMChannel.hh"
#include "G4N13GEMChannel.hh"
#include "G4N14GEMChannel.hh"
#include "G4N15GEMChannel.hh"
#include "G4N16GEMChannel.hh"
#include "G4N17GEMChannel.hh"
#include "G4O14GEMChannel.hh"
#include "G4O15GEMChannel.hh"
#include "G4O16GEMChannel.hh"
#include "G4O17GEMChannel.hh"
#include "G4O18GEMChannel.hh"
#include "G4O19GEMChannel.hh"
#include "G4O20GEMChannel.hh"
#include "G4F17GEMChannel.hh"
#include "G4F18GEMChannel.hh"
#include "G4F19GEMChannel.hh"
#include "G4F20GEMChannel.hh"
#include "G4F21GEMChannel.hh"
#include "G4Ne18GEMChannel.hh"
#include "G4Ne19GEMChannel.hh"
#include "G4Ne20GEMChannel.hh"
#include "G4Ne21GEMChannel.hh"
#include "G4Ne22GEMChannel.hh"
#include "G4Ne23GEMChannel.hh"
#include "G4Ne24GEMChannel.hh"
#include "G4Na21GEMChannel.hh"
#include "G4Na22GEMChannel.hh"
#include "G4Na23GEMChannel.hh"
#include "G4Na24GEMChannel.hh"
#include "G4Na25GEMChannel.hh"
#include "G4Mg22GEMChannel.hh"
#include "G4Mg23GEMChannel.hh"
#include "G4Mg24GEMChannel.hh"
#include "G4Mg25GEMChannel.hh"
#include "G4Mg26GEMChannel.hh"
#include "G4Mg27GEMChannel.hh"
#include "G4Mg28GEMChannel.hh"

#include "G4CompetitiveFission.hh"
#include "G4PhotonEvaporation.hh"


const G4EvaporationGEMFactory & 
G4EvaporationGEMFactory::operator=(const G4EvaporationGEMFactory & )
{
  G4Exception("G4EvaporationGEMFactory::operator= meant to not be accessable.");
  return *this;
}

G4bool 
G4EvaporationGEMFactory::operator==(const G4EvaporationGEMFactory & ) const
{
  G4Exception("G4EvaporationGEMFactory::operator== meant to not be accessable.");
  return false;
}

G4bool 
G4EvaporationGEMFactory::operator!=(const G4EvaporationGEMFactory & ) const
{
  G4Exception("G4EvaporationGEMFactory::operator!= meant to not be accessable.");
  return true;
}



std::vector<G4VEvaporationChannel*> * 
G4EvaporationGEMFactory::CreateChannel()
{
  std::vector<G4VEvaporationChannel*> * theChannel = 
    new std::vector<G4VEvaporationChannel*>;
  theChannel->reserve(68);

  theChannel->push_back( new G4NeutronGEMChannel() );  // n
  theChannel->push_back( new G4ProtonGEMChannel() );   // p
  theChannel->push_back( new G4DeuteronGEMChannel() ); // Deuteron
  theChannel->push_back( new G4TritonGEMChannel() );   // Triton
  theChannel->push_back( new G4He3GEMChannel() );      // He3
  theChannel->push_back( new G4AlphaGEMChannel() );    // Alpha
  theChannel->push_back( new G4He6GEMChannel() );      // He6
  theChannel->push_back( new G4He8GEMChannel() );      // He8
  theChannel->push_back( new G4Li6GEMChannel() );      // Li6
  theChannel->push_back( new G4Li7GEMChannel() );      // Li7
  theChannel->push_back( new G4Li8GEMChannel() );      // Li8
  theChannel->push_back( new G4Li9GEMChannel() );      // Li9
  theChannel->push_back( new G4Be7GEMChannel() );      // Be7    
  theChannel->push_back( new G4Be9GEMChannel() );      // Be9   
  theChannel->push_back( new G4Be10GEMChannel() );     // Be10   
  theChannel->push_back( new G4Be11GEMChannel() );     // Be11   
  theChannel->push_back( new G4Be12GEMChannel() );     // Be12
  theChannel->push_back( new G4B8GEMChannel() );       // B8
  theChannel->push_back( new G4B10GEMChannel() );      // B10
  theChannel->push_back( new G4B11GEMChannel() );      // B11
  theChannel->push_back( new G4B12GEMChannel() );      // B12
  theChannel->push_back( new G4B13GEMChannel() );      // B13
  theChannel->push_back( new G4C10GEMChannel() );      // C10
  theChannel->push_back( new G4C11GEMChannel() );      // C11
  theChannel->push_back( new G4C12GEMChannel() );      // C12
  theChannel->push_back( new G4C13GEMChannel() );      // C13
  theChannel->push_back( new G4C14GEMChannel() );      // C14
  theChannel->push_back( new G4C15GEMChannel() );      // C15
  theChannel->push_back( new G4C16GEMChannel() );      // C16
  theChannel->push_back( new G4N12GEMChannel() );      // N12
  theChannel->push_back( new G4N13GEMChannel() );      // N13
  theChannel->push_back( new G4N14GEMChannel() );      // N14
  theChannel->push_back( new G4N15GEMChannel() );      // N15
  theChannel->push_back( new G4N16GEMChannel() );      // N16
  theChannel->push_back( new G4N17GEMChannel() );      // N17
  theChannel->push_back( new G4O14GEMChannel() );      // O14
  theChannel->push_back( new G4O15GEMChannel() );      // O15
  theChannel->push_back( new G4O16GEMChannel() );      // O16
  theChannel->push_back( new G4O17GEMChannel() );      // O17
  theChannel->push_back( new G4O18GEMChannel() );      // O18
  theChannel->push_back( new G4O19GEMChannel() );      // O19
  theChannel->push_back( new G4O20GEMChannel() );      // O20
  theChannel->push_back( new G4F17GEMChannel() );      // F17
  theChannel->push_back( new G4F18GEMChannel() );      // F18
  theChannel->push_back( new G4F19GEMChannel() );      // F19
  theChannel->push_back( new G4F20GEMChannel() );      // F20
  theChannel->push_back( new G4F21GEMChannel() );      // F21
  theChannel->push_back( new G4Ne18GEMChannel() );     // Ne18
  theChannel->push_back( new G4Ne19GEMChannel() );     // Ne19
  theChannel->push_back( new G4Ne20GEMChannel() );     // Ne20
  theChannel->push_back( new G4Ne21GEMChannel() );     // Ne21
  theChannel->push_back( new G4Ne22GEMChannel() );     // Ne22
  theChannel->push_back( new G4Ne23GEMChannel() );     // Ne23
  theChannel->push_back( new G4Ne24GEMChannel() );     // Ne24
  theChannel->push_back( new G4Na21GEMChannel() );     // Na21
  theChannel->push_back( new G4Na22GEMChannel() );     // Na22
  theChannel->push_back( new G4Na23GEMChannel() );     // Na23
  theChannel->push_back( new G4Na24GEMChannel() );     // Na24
  theChannel->push_back( new G4Na25GEMChannel() );     // Na25
  theChannel->push_back( new G4Mg22GEMChannel() );     // Mg22
  theChannel->push_back( new G4Mg23GEMChannel() );     // Mg23
  theChannel->push_back( new G4Mg24GEMChannel() );     // Mg24
  theChannel->push_back( new G4Mg25GEMChannel() );     // Mg25
  theChannel->push_back( new G4Mg26GEMChannel() );     // Mg26
  theChannel->push_back( new G4Mg27GEMChannel() );     // Mg27
  theChannel->push_back( new G4Mg28GEMChannel() );     // Mg28
  
  theChannel->push_back( new G4CompetitiveFission() ); // Fission Channel
  theChannel->push_back( new G4PhotonEvaporation() );  // Photon Channel

  return theChannel;

}



