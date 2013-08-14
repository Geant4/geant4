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
// $Id: $
// GEANT4 tag: $Name: not supported by cvs2svn $
//

#include "G4CascadeMuMinusPChannel.hh"
#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;

namespace {
  // mu- p : Outgoing particle types of a given multiplicity
  static const G4int mump2bfs[1][2] =
  {{neu,mnu}};

  static const G4int mump3bfs[1][3] =
  {{neu,mnu,gam}};

  static const G4int mump4bfs[1][4] =
  {{neu,mnu,gam,gam}};

  static const G4int mump5bfs[1][5] =
  {{neu,mnu,gam,gam,gam}};

  static const G4int mump6bfs[1][6] =
  {{neu,mnu,gam,gam,gam,gam}};

  static const G4int mump7bfs[1][7] =
  {{neu,mnu,gam,gam,gam,gam,gam}};

  static const G4int mump8bfs[1][8] =
  {{neu,mnu,gam,gam,gam,gam,gam,gam}};

  static const G4int mump9bfs[1][9] =
  {{neu,mnu,gam,gam,gam,gam,gam,gam,gam}};
}

namespace {
  // Total mu- p cross section as a function of kinetic energy
  static const G4double mumpTotXSec[30] =
   { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
     0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
     0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
/*
   {100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0,
    100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0,
    100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0};
*/
  static const G4double mumpCrossSections[8][30] = {
  //
  // multiplicity 2 (1 channel)
  //
  // nu_mu n
   { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
     0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
     0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01},
/*
   {100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0,
    100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0,
    100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0},
*/ 
  //
  // multiplicity 3 (1 channel)
  //
  // nu_mu n gam
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
  //
  // multiplicity 4 (1 channel)
  //
  // nu_mu n gam gam
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}, 
  //
  // multiplicity 5 (1 channel)
  //
  // nu_mu n gam gam gam
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
  //
  // multiplicity 6 (1 channel)
  //
  // nu_mu n gam gam gam gam
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
  //
  // multiplicity 7 (1 channel)
  //
  // nu_mu n gam gam gam gam gam
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
  //
  // multiplicity 8 (1 channel)
  //
  // nu_mu n gam gam gam gam gam gam
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0},
  //
  // multiplicity 9 (1 channel)
  //
  // nu_mu n gam gam gam gam gam gam gam 
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0}};
}


const G4CascadeMuMinusPChannelData::data_t
G4CascadeMuMinusPChannelData::data(mump2bfs, mump3bfs, mump4bfs, mump5bfs,
                               mump6bfs, mump7bfs, mump8bfs, mump9bfs, 
                               mumpCrossSections, mumpTotXSec, mum*pro,
                               "MuMinusProton");

