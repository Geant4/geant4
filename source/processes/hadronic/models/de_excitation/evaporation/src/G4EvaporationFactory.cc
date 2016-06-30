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
// $Id: G4EvaporationFactory.cc 96634 2016-04-27 09:31:49Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modifications:
//
// 23 January 2012 V.Ivanchenko added pointer of G4VPhotonEvaporation 

#include "G4EvaporationFactory.hh"

#include "G4NeutronEvaporationChannel.hh"
#include "G4ProtonEvaporationChannel.hh"
#include "G4DeuteronEvaporationChannel.hh"
#include "G4TritonEvaporationChannel.hh"
#include "G4He3EvaporationChannel.hh"
#include "G4AlphaEvaporationChannel.hh"

#include "G4CompetitiveFission.hh"

G4EvaporationFactory::G4EvaporationFactory(G4VEvaporationChannel* ptr)
  : G4VEvaporationFactory(ptr)
{}

G4EvaporationFactory::~G4EvaporationFactory()
{}

std::vector<G4VEvaporationChannel*>* G4EvaporationFactory::GetChannel()
{
  std::vector<G4VEvaporationChannel*> * theChannel = 
    new std::vector<G4VEvaporationChannel*>;
  theChannel->reserve(8);

  theChannel->push_back( thePhotonEvaporation );          // Photon Channel
  theChannel->push_back( new G4CompetitiveFission() );    // Fission Channel

  theChannel->push_back( new G4NeutronEvaporationChannel() );  // n
  theChannel->push_back( new G4ProtonEvaporationChannel() );   // p
  theChannel->push_back( new G4DeuteronEvaporationChannel() ); // Deuteron
  theChannel->push_back( new G4TritonEvaporationChannel() );   // Triton
  theChannel->push_back( new G4He3EvaporationChannel() );      // He3
  theChannel->push_back( new G4AlphaEvaporationChannel() );    // Alpha

  size_t nn = theChannel->size();
  for(size_t i=1; i<nn; ++i) { 
    (*theChannel)[i]->SetPhotonEvaporation(thePhotonEvaporation);
  }

  return theChannel;

}



