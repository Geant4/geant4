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
// $Id: G4EvaporationGEMFactory.cc,v 1.5 2005-06-28 11:09:37 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

const G4EvaporationGEMFactory & 
G4EvaporationGEMFactory::operator=(const G4EvaporationGEMFactory & )
{
  throw G4HadronicException(__FILE__, __LINE__, "G4EvaporationGEMFactory::operator= meant to not be accessable.");
  return *this;
}

G4bool 
G4EvaporationGEMFactory::operator==(const G4EvaporationGEMFactory & ) const
{
  throw G4HadronicException(__FILE__, __LINE__, "G4EvaporationGEMFactory::operator== meant to not be accessable.");
  return false;
}

G4bool 
G4EvaporationGEMFactory::operator!=(const G4EvaporationGEMFactory & ) const
{
  throw G4HadronicException(__FILE__, __LINE__, "G4EvaporationGEMFactory::operator!= meant to not be accessable.");
  return true;
}

typedef GROUP68(G4NeutronGEMChannel, G4ProtonGEMChannel, G4DeuteronGEMChannel, G4TritonGEMChannel, 
                G4He3GEMChannel, G4AlphaGEMChannel, G4He6GEMChannel, G4He8GEMChannel, G4Li6GEMChannel,
                G4Li7GEMChannel, G4Li8GEMChannel, G4Li9GEMChannel, G4Be7GEMChannel, G4Be9GEMChannel,
                G4Be10GEMChannel, G4Be11GEMChannel, G4Be12GEMChannel, G4B8GEMChannel, G4B10GEMChannel, 
                G4B11GEMChannel, G4B12GEMChannel, G4B13GEMChannel, G4C10GEMChannel, G4C11GEMChannel,
                G4C12GEMChannel, G4C13GEMChannel, G4C14GEMChannel, G4C15GEMChannel, G4C16GEMChannel,
                G4N12GEMChannel, G4N13GEMChannel, G4N14GEMChannel, G4N15GEMChannel, G4N16GEMChannel,
                G4N17GEMChannel, G4O14GEMChannel, G4O15GEMChannel, G4O16GEMChannel, G4O17GEMChannel,
                G4O18GEMChannel, G4O19GEMChannel, G4O20GEMChannel, G4F17GEMChannel, G4F18GEMChannel,
                G4F19GEMChannel, G4F20GEMChannel, G4F21GEMChannel, G4Ne18GEMChannel, G4Ne19GEMChannel,
                G4Ne20GEMChannel, G4Ne21GEMChannel, G4Ne22GEMChannel, G4Ne23GEMChannel, G4Ne24GEMChannel,
                G4Na21GEMChannel, G4Na22GEMChannel, G4Na23GEMChannel, G4Na24GEMChannel, G4Na25GEMChannel,
                G4Mg22GEMChannel, G4Mg23GEMChannel, G4Mg24GEMChannel, G4Mg25GEMChannel, G4Mg26GEMChannel,
                G4Mg27GEMChannel, G4Mg28GEMChannel, G4CompetitiveFission, G4PhotonEvaporation) the_channels;
                 
std::vector<G4VEvaporationChannel*> * G4EvaporationGEMFactory::CreateChannel()
{
  std::vector<G4VEvaporationChannel*> * theChannel = 
    new std::vector<G4VEvaporationChannel*>;
  theChannel->reserve(68);
  gem_push_one_new<the_channels>(theChannel);

  return theChannel;
}
