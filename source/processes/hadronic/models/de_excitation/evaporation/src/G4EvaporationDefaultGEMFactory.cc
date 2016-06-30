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
// $Id: G4EvaporationDefaultGEMFactory.cc 96634 2016-04-27 09:31:49Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by J. M. Quesada (July 2009) on base of V. Lara code
// V.Ivanchenko cleanup
//
// new hybrid Default-GEM evaoration model:
//      - default evaporation for n,p,d,t and alpha particles
//      - GEM evaporation for light nuclei evaporation (2<Z<13,4<A<29) 
//
//
// Modifications:
//
// 23 January 2012 V.Ivanchenko added pointer of G4VPhotonEvaporation 

#include "G4EvaporationDefaultGEMFactory.hh"

#include "G4NeutronEvaporationChannel.hh"
#include "G4ProtonEvaporationChannel.hh"
#include "G4DeuteronEvaporationChannel.hh"
#include "G4TritonEvaporationChannel.hh"
#include "G4He3EvaporationChannel.hh"
#include "G4AlphaEvaporationChannel.hh"
//JMQ 220709 new GEM channels
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

G4EvaporationDefaultGEMFactory::G4EvaporationDefaultGEMFactory(G4VEvaporationChannel* ptr)
  : G4VEvaporationFactory(ptr)
{}

G4EvaporationDefaultGEMFactory::~G4EvaporationDefaultGEMFactory()
{}

std::vector<G4VEvaporationChannel*>* G4EvaporationDefaultGEMFactory::GetChannel()
{
  std::vector<G4VEvaporationChannel*> * theChannel = 
    new std::vector<G4VEvaporationChannel*>;
  theChannel->reserve(68);

  theChannel->push_back( thePhotonEvaporation );        // Photon Channel
  theChannel->push_back( new G4CompetitiveFission() );  // Fission Channel

  // JMQ 220709 standard particle evaporation channels (Z<3,A<5)
  theChannel->push_back( new G4NeutronEvaporationChannel() );  // n
  theChannel->push_back( new G4ProtonEvaporationChannel() );   // p
  theChannel->push_back( new G4DeuteronEvaporationChannel() ); // Deuteron
  theChannel->push_back( new G4TritonEvaporationChannel() );   // Triton
  theChannel->push_back( new G4He3EvaporationChannel() );      // He3
  theChannel->push_back( new G4AlphaEvaporationChannel() );    // Alpha

  //JMQ 220709 new GEM particle avaporation channels channels (2<Z<13,4<A<29)
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

  size_t nn = theChannel->size();
  for(size_t i=1; i<nn; ++i) { 
    (*theChannel)[i]->SetPhotonEvaporation(thePhotonEvaporation);
  }
  return theChannel;
}



