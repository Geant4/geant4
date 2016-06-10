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
// $Id: G4VEvaporation.cc 85443 2014-10-29 14:35:57Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998) writen from G4Evaporation.cc (May 1998)
//
// Modifications:
//
// 23 January 2012 V.Ivanchenko added pointer of G4VPhotonEvaporation 

#include "G4VEvaporation.hh"
#include "G4VEvaporationChannel.hh"

G4VEvaporation::G4VEvaporation()
  :thePhotonEvaporation(0),OPTxs(3),useSICB(false)
   ,theChannels(0),theChannelFactory(0)
{}

G4VEvaporation::~G4VEvaporation() 
{}

void G4VEvaporation::Initialise()
{}

void G4VEvaporation::SetPhotonEvaporation(G4VEvaporationChannel* ptr)
{
  if(thePhotonEvaporation != ptr) {
    delete thePhotonEvaporation;
    thePhotonEvaporation = ptr;
  }
}

void G4VEvaporation::BreakFragment(G4FragmentVector*, G4Fragment*)
{}




