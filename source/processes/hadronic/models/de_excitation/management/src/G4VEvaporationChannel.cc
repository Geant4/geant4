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
// $Id: G4VEvaporationChannel.cc 86986 2014-11-21 13:00:05Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modified:
// 24.04.2010 (V.Ivanchenko) moved constructor and destructor to source; added two
//                          new virtual methods EmittedFragment(s) to allow more optimal
//                          work with G4Fragment objects; removed unnecesary exceptions
// 28.10.2010 V.Ivanchenko defined members in constructor and cleaned up

#include "G4VEvaporationChannel.hh"

G4VEvaporationChannel::G4VEvaporationChannel(const G4String & aName,
					     G4EvaporationChannelType timeType) 
  :sampleDecayTime(timeType),OPTxs(3),useSICB(false),Name(aName) 
{}

G4VEvaporationChannel::~G4VEvaporationChannel() 
{}

void G4VEvaporationChannel::Initialise()
{}

G4double G4VEvaporationChannel::GetLifeTime(G4Fragment*)
{
  return 0.0;
}

G4Fragment* G4VEvaporationChannel::EmittedFragment(G4Fragment*)
{
  return 0;
}

G4FragmentVector* G4VEvaporationChannel::BreakUpFragment(G4Fragment*)
{
  return 0;
}

G4bool G4VEvaporationChannel::BreakUpChain(G4FragmentVector*, G4Fragment*)
{
  return false;
}

void G4VEvaporationChannel::Dump() const
{}
