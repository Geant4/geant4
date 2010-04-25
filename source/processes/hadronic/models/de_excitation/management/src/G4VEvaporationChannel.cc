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
// $Id: G4VEvaporationChannel.cc,v 1.6 2010-04-25 18:43:08 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modified:
// 24.04.2010 (V.Ivanchenko) moved constructor and destructor to source; added two
//                          new virtual methods EmittedFragment(s) to allow more optimal
//                          work with G4Fragment objects; removed unnecesary exceptions

#include "G4VEvaporationChannel.hh"
#include "G4HadronicException.hh"

G4VEvaporationChannel::G4VEvaporationChannel(const G4String & aName) 
  : Name(aName) 
{}

G4VEvaporationChannel::~G4VEvaporationChannel() 
{}

//G4VEvaporationChannel::G4VEvaporationChannel(const G4VEvaporationChannel &)
//{
// throw G4HadronicException(__FILE__, __LINE__, "G4VEvaporationChannel::copy_constructor meant to not be accessable");
//}
//const G4VEvaporationChannel & G4VEvaporationChannel::operator=(const G4VEvaporationChannel &)
//{
//  throw G4HadronicException(__FILE__, __LINE__, "G4VEvaporationChannel::operator= meant to not be accessable");
//  return *this;
//}

G4bool G4VEvaporationChannel::operator==(const G4VEvaporationChannel &right) const
{
  return (this == (G4VEvaporationChannel *) &right);
}

G4bool G4VEvaporationChannel::operator!=(const G4VEvaporationChannel &right) const
{
  return (this != (G4VEvaporationChannel *) &right);
}

G4Fragment* G4VEvaporationChannel::EmittedFragment(G4Fragment*)
{
  return 0;
}

G4FragmentVector* G4VEvaporationChannel::BreakUpFragment(G4Fragment*)
{
  return 0;
}
