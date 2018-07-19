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
// $Id: G4EvaporationGEMFactoryVI.cc 92422 2015-09-01 08:38:35Z vnivanch $
//
//
// GEM de-excitation model
// by V. Ivanchenko (July 2016)
//

#include "G4EvaporationGEMFactoryVI.hh"
#include "G4GEMChannelVI.hh"
#include "G4CompetitiveFission.hh"

G4EvaporationGEMFactoryVI::G4EvaporationGEMFactoryVI(G4VEvaporationChannel* ptr)
  : G4VEvaporationFactory(ptr)
{} 
  
G4EvaporationGEMFactoryVI::~G4EvaporationGEMFactoryVI() 
{}
                 
std::vector<G4VEvaporationChannel*>* G4EvaporationGEMFactoryVI::GetChannel()
{
  std::vector<G4VEvaporationChannel*> * theChannel = 
    new std::vector<G4VEvaporationChannel*>;
  theChannel->reserve(81);

  theChannel->push_back( thePhotonEvaporation );  // Photon Channel
  theChannel->push_back( new G4CompetitiveFission() ); // Fission Channel

  theChannel->push_back( new G4GEMChannelVI( 1, 0) );// n
  theChannel->push_back( new G4GEMChannelVI( 1, 1) );// p
  theChannel->push_back( new G4GEMChannelVI( 2, 1) );// Deuteron
  theChannel->push_back( new G4GEMChannelVI( 3, 1) );// Triton
  theChannel->push_back( new G4GEMChannelVI( 3, 2) );// He3
  theChannel->push_back( new G4GEMChannelVI( 4, 2) );// Alpha
  theChannel->push_back( new G4GEMChannelVI( 5, 2) );// He5
  theChannel->push_back( new G4GEMChannelVI( 5, 3) );// Li5
  theChannel->push_back( new G4GEMChannelVI( 6, 3) );// Li6
  theChannel->push_back( new G4GEMChannelVI( 7, 3) );// Li7
  theChannel->push_back( new G4GEMChannelVI( 8, 3) );// Li8
  theChannel->push_back( new G4GEMChannelVI( 9, 3) );// Li9
  theChannel->push_back( new G4GEMChannelVI( 7, 4) );// Be7
  theChannel->push_back( new G4GEMChannelVI( 8, 4) );// Be8
  theChannel->push_back( new G4GEMChannelVI( 9, 4) );// Be9
  theChannel->push_back( new G4GEMChannelVI(10, 4) );// Be10
  theChannel->push_back( new G4GEMChannelVI(11, 4) );// Be11
  theChannel->push_back( new G4GEMChannelVI( 8, 5) );// B8
  theChannel->push_back( new G4GEMChannelVI( 9, 5) );// B9
  theChannel->push_back( new G4GEMChannelVI(10, 5) );// B10
  theChannel->push_back( new G4GEMChannelVI(11, 5) );// B11
  theChannel->push_back( new G4GEMChannelVI(12, 5) );// B12
  theChannel->push_back( new G4GEMChannelVI(13, 5) );// B13
  theChannel->push_back( new G4GEMChannelVI(10, 6) );// C10
  theChannel->push_back( new G4GEMChannelVI(11, 6) );// C11
  theChannel->push_back( new G4GEMChannelVI(12, 6) );// C12
  theChannel->push_back( new G4GEMChannelVI(13, 6) );// C13
  theChannel->push_back( new G4GEMChannelVI(14, 6) );// C14
  theChannel->push_back( new G4GEMChannelVI(15, 6) );// C15
  theChannel->push_back( new G4GEMChannelVI(16, 6) );// C16
  theChannel->push_back( new G4GEMChannelVI(13, 7) );// N13
  theChannel->push_back( new G4GEMChannelVI(14, 7) );// N14
  theChannel->push_back( new G4GEMChannelVI(15, 7) );// N15
  theChannel->push_back( new G4GEMChannelVI(16, 7) );// N16
  theChannel->push_back( new G4GEMChannelVI(17, 7) );// N17
  theChannel->push_back( new G4GEMChannelVI(18, 7) );// N17
  theChannel->push_back( new G4GEMChannelVI(15, 8) );// O15
  theChannel->push_back( new G4GEMChannelVI(16, 8) );// O16
  theChannel->push_back( new G4GEMChannelVI(17, 8) );// O17
  theChannel->push_back( new G4GEMChannelVI(18, 8) );// O18
  theChannel->push_back( new G4GEMChannelVI(19, 8) );// O19
  theChannel->push_back( new G4GEMChannelVI(20, 8) );// O20
  theChannel->push_back( new G4GEMChannelVI(21, 8) );// O21
  theChannel->push_back( new G4GEMChannelVI(22, 8) );// O22
  theChannel->push_back( new G4GEMChannelVI(17, 9) );// F17
  theChannel->push_back( new G4GEMChannelVI(18, 9) );// F18
  theChannel->push_back( new G4GEMChannelVI(19, 9) );// F19
  theChannel->push_back( new G4GEMChannelVI(20, 9) );// F20
  theChannel->push_back( new G4GEMChannelVI(21, 9) );// F21
  theChannel->push_back( new G4GEMChannelVI(22, 9) );// F22
  theChannel->push_back( new G4GEMChannelVI(23, 9) );// F23
  theChannel->push_back( new G4GEMChannelVI(24, 9) );// F24
  theChannel->push_back( new G4GEMChannelVI(25, 9) );// F25
  theChannel->push_back( new G4GEMChannelVI(26, 9) );// F26
  theChannel->push_back( new G4GEMChannelVI(27, 9) );// F27
  theChannel->push_back( new G4GEMChannelVI(18,10) );// Ne18
  theChannel->push_back( new G4GEMChannelVI(19,10) );// Ne19
  theChannel->push_back( new G4GEMChannelVI(20,10) );// Ne20
  theChannel->push_back( new G4GEMChannelVI(21,10) );// Ne21
  theChannel->push_back( new G4GEMChannelVI(22,10) );// Ne22
  theChannel->push_back( new G4GEMChannelVI(23,10) );// Ne23
  theChannel->push_back( new G4GEMChannelVI(24,10) );// Ne24
  theChannel->push_back( new G4GEMChannelVI(25,10) );// Ne25
  theChannel->push_back( new G4GEMChannelVI(26,10) );// Ne26
  theChannel->push_back( new G4GEMChannelVI(27,10) );// Ne27
  theChannel->push_back( new G4GEMChannelVI(28,10) );// Ne28
  theChannel->push_back( new G4GEMChannelVI(21,11) );// Na21
  theChannel->push_back( new G4GEMChannelVI(22,11) );// Na22
  theChannel->push_back( new G4GEMChannelVI(23,11) );// Na23
  theChannel->push_back( new G4GEMChannelVI(24,11) );// Na24
  theChannel->push_back( new G4GEMChannelVI(25,11) );// Na25
  theChannel->push_back( new G4GEMChannelVI(26,11) );// Na26
  theChannel->push_back( new G4GEMChannelVI(27,11) );// Na27
  theChannel->push_back( new G4GEMChannelVI(28,11) );// Na28
  theChannel->push_back( new G4GEMChannelVI(22,12) );// Mg22
  theChannel->push_back( new G4GEMChannelVI(23,12) );// Mg23
  theChannel->push_back( new G4GEMChannelVI(24,12) );// Mg24
  theChannel->push_back( new G4GEMChannelVI(25,12) );// Mg25
  theChannel->push_back( new G4GEMChannelVI(26,12) );// Mg26
  theChannel->push_back( new G4GEMChannelVI(27,12) );// Mg27
  theChannel->push_back( new G4GEMChannelVI(28,12) );// Mg28

  size_t nn = theChannel->size();
  for(size_t i=2; i<nn; ++i) { 
    (*theChannel)[i]->SetPhotonEvaporation(thePhotonEvaporation);
  }

  return theChannel;
}
