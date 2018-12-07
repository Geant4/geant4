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

#include "G4CascadePiZeroPChannel.hh"
#include "G4CascadePiZeroNChannel.hh"
#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;

namespace {
  // pi0 p : Outgoing particle types of a given multiplicity
  static const G4int pizP2bfs[5][2] =
  {{pi0,pro}, {pip,neu}, {kpl,lam}, {kpl,s0}, {k0,sp}};

  static const G4int pizP3bfs[13][3] =
  {{pip,pro,pim}, {pi0,pro,pi0}, {pi0,neu,pip}, {pi0,kpl,lam}, 
   {pip,k0,lam},  {pi0,kpl,s0},  {pip,k0,s0},   {pi0,k0,sp},
   {pim,kpl,sp},  {pip,kpl,sm},  {kpl,pro,kmi}, {k0,pro,k0b},
   {k0b,neu,kpl}};

  static const G4int pizP4bfs[22][4] =
  {{pi0,pro,pip,pim}, {pi0,pro,pi0,pi0}, {pip,neu,pip,pim},
   {pi0,neu,pip,pi0}, {pip,pim,kpl,lam}, {pip,pim,kpl,s0}, 
   {pi0,pi0,kpl,lam}, {pi0,pi0,kpl,s0},  {pi0,pi0,k0,sp}, 
   {pi0,pip,k0,lam},  {pi0,pip,k0,s0},   {pip,pim,k0,sp}, 
   {pi0,pim,kpl,sp},  {pi0,pip,kpl,sm},  {pip,pip,k0,sm},
   {pi0,pro,kpl,kmi}, {pi0,pro,k0,k0b},  {pip,pro,kmi,k0},
   {pim,pro,kpl,k0b}, {pip,neu,kpl,kmi}, {pip,neu,k0,k0b},
   {pi0,neu,kpl,k0b}}; 

  static const G4int pizP5bfs[31][5] =
  {{pip,pro,pip,pim,pim}, {pi0,pro,pip,pim,pi0}, {pi0,pro,pi0,pi0,pi0},
   {pi0,neu,pip,pip,pim}, {pi0,neu,pip,pi0,pi0}, {pi0,pip,pim,kpl,lam},
   {pi0,pi0,pi0,kpl,lam}, {pip,pip,pim,k0,lam},  {pi0,pip,pi0,k0,lam},
   {pi0,pim,pip,kpl,s0},  {pi0,pi0,pi0,kpl,s0},  {pip,pip,pim,k0,s0},
   {pi0,pi0,pip,k0,s0},   {pi0,pim,pip,k0,sp},   {pi0,pi0,pi0,k0,sp},
   {pip,pim,pim,kpl,sp},  {pi0,pi0,pim,kpl,sp},  {pip,pip,pim,kpl,sm}, 
   {pi0,pi0,pip,kpl,sm},  {pi0,pip,pip,k0,sm},   {pi0,pro,pi0,kpl,kmi},
   {pi0,pro,pi0,k0,k0b},  {pi0,pro,pim,kpl,k0b}, {pi0,pro,pip,kmi,k0},
   {pip,pro,pim,kpl,kmi}, {pip,pro,pim,k0,k0b},  {pi0,neu,pip,kpl,kmi},
   {pi0,neu,pip,k0,k0b},  {pi0,neu,pi0,kpl,k0b}, {pip,neu,pip,kmi,k0},
   {pip,neu,pim,kpl,k0b}};

  static const G4int pizP6bfs[39][6] =
  {{pi0,pro,pip,pip,pim,pim}, {pi0,pro,pip,pim,pi0,pi0},
   {pi0,pro,pi0,pi0,pi0,pi0}, {pip,neu,pip,pip,pim,pim},
   {pi0,neu,pip,pip,pim,pi0}, {pip,neu,pi0,pi0,pi0,pi0},
   {pip,pip,pim,pim,kpl,lam}, {pi0,pip,pim,pi0,kpl,lam},
   {pi0,pi0,pi0,pi0,kpl,lam}, {pi0,pip,pip,pim,k0,lam},
   {pi0,pi0,pi0,pip,k0,lam},  {pip,pip,pim,pim,kpl,s0},
   {pi0,pip,pim,pi0,kpl,s0},  {pi0,pi0,pi0,pi0,kpl,s0},
   {pi0,pip,pip,pim,k0,s0},   {pi0,pi0,pi0,pip,k0,s0},
   {pi0,pip,pim,pim,kpl,sp},  {pi0,pi0,pi0,pim,kpl,sp},
   {pip,pip,pim,pim,k0,sp},   {pi0,pim,pip,pi0,k0,sp},
   {pi0,pip,pip,pim,kpl,sm},  {pi0,pi0,pi0,pip,kpl,sm},
   {pip,pip,pip,pim,k0,sm},   {pi0,pip,pip,pi0,k0,sm},
   {pi0,pro,pim,pi0,kpl,k0b}, {pro,pip,pim,pim,kpl,k0b},
   {pi0,pro,pip,pim,kpl,kmi}, {pi0,pro,pi0,pi0,kpl,kmi},
   {pi0,pro,pip,pim,k0,k0b},  {pi0,pro,pi0,pi0,k0,k0b},
   {pro,pip,pip,pim,kmi,k0},  {pi0,pro,pip,pi0,kmi,k0},
   {pip,neu,pi0,pim,kpl,k0b}, {pi0,neu,pi0,pi0,kpl,k0b},
   {pip,neu,pip,pim,kpl,kmi}, {pip,neu,pi0,pi0,kpl,kmi},
   {pip,neu,pip,pim,k0,k0b},  {pip,neu,pi0,pi0,k0,k0b},
   {pip,neu,pip,pi0,kmi,k0}};

  static const G4int pizP7bfs[46][7] =
  {{pip,pro,pip,pip,pim,pim,pim}, {pi0,pro,pip,pip,pim,pim,pi0},
   {pi0,pro,pip,pim,pi0,pi0,pi0}, {pi0,pro,pi0,pi0,pi0,pi0,pi0},
   {pi0,neu,pip,pip,pip,pim,pim}, {pi0,neu,pip,pip,pim,pi0,pi0},
   {pi0,neu,pip,pi0,pi0,pi0,pi0}, {pi0,pip,pip,pim,pim,kpl,lam}, 
   {pi0,pim,pip,pi0,pi0,kpl,lam}, {pip,pip,pip,pim,pim,k0,lam},
   {pi0,pip,pip,pim,pi0,k0,lam},  {pi0,pip,pi0,pi0,pi0,k0,lam},
   {pi0,pip,pip,pim,pim,kpl,s0},  {pi0,pip,pim,pi0,pi0,kpl,s0},
   {pip,pip,pip,pim,pim,k0,s0},   {pi0,pip,pip,pim,pi0,k0,s0},
   {pi0,pip,pi0,pi0,pi0,k0,s0},   {pip,pip,pim,pim,pim,kpl,sp},  
   {pi0,pip,pim,pim,pi0,kpl,sp},  {pi0,pim,pi0,pi0,pi0,kpl,sp},
   {pi0,pip,pip,pim,pim,k0,sp},   {pi0,pip,pim,pi0,pi0,k0,sp},
   {pip,pip,pip,pim,pim,kpl,sm},  {pi0,pip,pip,pim,pi0,kpl,sm},
   {pi0,pip,pi0,pi0,pi0,kpl,sm},  {pi0,pip,pip,pip,pim,k0,sm},
   {pi0,pip,pip,pi0,pi0,k0,sm},   {pi0,pro,pim,pi0,pi0,kpl,k0b},
   {pi0,pro,pip,pim,pim,kpl,k0b}, {pi0,pro,pi0,pi0,pi0,kpl,kmi},
   {pi0,pro,pip,pim,pi0,kpl,kmi}, {pro,pip,pip,pim,pim,kpl,kmi},
   {pi0,pro,pi0,pi0,pi0,k0,k0b},  {pi0,pro,pip,pim,pi0,k0,k0b},
   {pip,pro,pip,pim,pim,k0,k0b},  {pi0,pro,pip,pi0,pi0,kmi,k0},
   {pi0,pro,pip,pip,pim,kmi,k0},  {pi0,neu,pi0,pi0,pi0,kpl,k0b},
   {pi0,neu,pip,pim,pi0,kpl,k0b}, {pip,neu,pip,pim,pim,kpl,k0b},
   {pi0,neu,pip,pi0,pi0,kpl,kmi}, {pi0,neu,pip,pip,pim,kpl,kmi},
   {pip,neu,pi0,pi0,pi0,k0,k0b},  {pip,neu,pip,pim,pi0,k0,k0b},
   {pip,neu,pip,pi0,pi0,kmi,k0},  {pip,neu,pip,pip,pim,kmi,k0}};

  static const G4int pizP8bfs[51][8] =
  {{pi0,pro,pip,pip,pip,pim,pim,pim}, {pi0,pro,pip,pip,pim,pim,pi0,pi0},
   {pi0,pro,pip,pim,pi0,pi0,pi0,pi0}, {pi0,pro,pi0,pi0,pi0,pi0,pi0,pi0},
   {pip,neu,pip,pip,pip,pim,pim,pim}, {pi0,neu,pip,pip,pip,pim,pim,pi0},
   {pi0,neu,pip,pip,pim,pi0,pi0,pi0}, {pi0,neu,pip,pi0,pi0,pi0,pi0,pi0},
   {pip,pip,pip,pim,pim,pim,kpl,lam}, {pi0,pip,pip,pim,pim,pi0,kpl,lam},
   {pi0,pip,pim,pi0,pi0,pi0,kpl,lam}, {pi0,pip,pip,pip,pim,pim,k0,lam},
   {pi0,pip,pip,pim,pi0,pi0,k0,lam},  {pi0,pip,pi0,pi0,pi0,pi0,k0,lam},
   {pip,pip,pip,pim,pim,pim,kpl,s0},  {pi0,pip,pip,pim,pim,pi0,kpl,s0},
   {pi0,pip,pim,pi0,pi0,pi0,kpl,s0},  {pi0,pip,pip,pip,pim,pim,k0,s0},
   {pi0,pip,pip,pim,pi0,pi0,k0,s0},   {pi0,pip,pip,pim,pim,pim,kpl,sp},
   {pi0,pip,pim,pim,pi0,pi0,kpl,sp},  {pi0,pim,pi0,pi0,pi0,pi0,kpl,sp},
   {pip,pip,pip,pim,pim,pim,k0,sp},   {pi0,pip,pip,pim,pim,pi0,k0,sp},
   {pi0,pip,pim,pi0,pi0,pi0,k0,sp},   {pi0,pip,pip,pip,pim,pim,kpl,sm},
   {pi0,pip,pip,pim,pi0,pi0,kpl,sm},  {pi0,pip,pi0,pi0,pi0,pi0,kpl,sm},
   {pi0,pip,pip,pip,pim,pi0,k0,sm},   {pip,pip,pip,pip,pim,pim,k0,sm},
   {pi0,pip,pip,pi0,pi0,pi0,k0,sm},   {pip,pro,pip,pim,pim,pim,kpl,k0b}, 
   {pi0,pro,pip,pim,pim,pi0,kpl,k0b}, {pi0,pro,pim,pi0,pi0,pi0,kpl,k0b},
   {pi0,pro,pip,pip,pim,pim,kpl,kmi}, {pi0,pro,pip,pim,pi0,pi0,kpl,kmi},
   {pi0,pro,pip,pip,pim,pim,k0,k0b},  {pi0,pro,pip,pim,pi0,pi0,k0,k0b},
   {pip,pro,pip,pip,pim,pim,kmi,k0},  {pi0,pro,pip,pip,pim,pi0,kmi,k0},
   {pi0,pro,pip,pi0,pi0,pi0,kmi,k0},  {pi0,neu,pip,pip,pim,pim,kpl,k0b},
   {pi0,neu,pip,pim,pi0,pi0,kpl,k0b}, {pip,neu,pip,pip,pim,pim,kpl,kmi},
   {pi0,neu,pip,pip,pim,pi0,kpl,kmi}, {pi0,neu,pip,pi0,pi0,pi0,kpl,kmi},
   {pip,neu,pip,pip,pim,pim,k0,k0b},  {pi0,neu,pip,pip,pim,pi0,k0,k0b},
   {pi0,neu,pip,pi0,pi0,pi0,k0,k0b},  {pi0,neu,pip,pip,pip,pim,kmi,k0},
   {pi0,neu,pip,pip,pi0,pi0,kmi,k0}};

  static const G4int pizP9bfs[58][9] =
  {{pip,pro,pip,pip,pip,pim,pim,pim,pim}, {pi0,pro,pip,pip,pip,pim,pim,pim,pi0},
   {pi0,pro,pip,pip,pim,pim,pi0,pi0,pi0}, {pi0,pro,pip,pim,pi0,pi0,pi0,pi0,pi0},
   {pi0,pro,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {pi0,neu,pip,pip,pip,pip,pim,pim,pim},
   {pi0,neu,pip,pip,pip,pim,pim,pi0,pi0}, {pi0,neu,pip,pip,pim,pi0,pi0,pi0,pi0},
   {pi0,neu,pip,pi0,pi0,pi0,pi0,pi0,pi0}, {pi0,pip,pip,pip,pim,pim,pim,kpl,lam},
   {pi0,pip,pip,pim,pim,pi0,pi0,kpl,lam}, {pi0,pip,pim,pi0,pi0,pi0,pi0,kpl,lam},
   {pip,pip,pip,pip,pim,pim,pim,k0,lam},  {pi0,pip,pip,pip,pim,pim,pi0,k0,lam},
   {pi0,pip,pip,pim,pi0,pi0,pi0,k0,lam},  {pi0,pip,pi0,pi0,pi0,pi0,pi0,k0,lam},
   {pi0,pip,pip,pip,pim,pim,pim,kpl,s0},  {pi0,pip,pip,pim,pim,pi0,pi0,kpl,s0},
   {pi0,pip,pim,pi0,pi0,pi0,pi0,kpl,s0},  {pip,pip,pip,pip,pim,pim,pim,k0,s0},
   {pi0,pip,pip,pip,pim,pim,pi0,k0,s0},   {pi0,pip,pip,pim,pi0,pi0,pi0,k0,s0},
   {pip,pip,pip,pim,pim,pim,pim,kpl,sp},  {pi0,pip,pip,pim,pim,pim,pi0,kpl,sp},
   {pi0,pip,pim,pim,pi0,pi0,pi0,kpl,sp},  {pi0,pip,pip,pip,pim,pim,pim,k0,sp},
   {pi0,pip,pip,pim,pim,pi0,pi0,k0,sp},   {pi0,pip,pim,pi0,pi0,pi0,pi0,k0,sp},
   {pip,pip,pip,pip,pim,pim,pim,kpl,sm},  {pi0,pip,pip,pip,pim,pim,pi0,kpl,sm},
   {pi0,pip,pip,pim,pi0,pi0,pi0,kpl,sm},  {pi0,pip,pip,pip,pip,pim,pim,k0,sm},
   {pi0,pip,pip,pip,pim,pi0,pi0,k0,sm},   {pi0,pip,pip,pi0,pi0,pi0,pi0,k0,sm},
   {pi0,pro,pip,pip,pim,pim,pim,kpl,k0b}, {pi0,pro,pip,pim,pim,pi0,pi0,kpl,k0b},
   {pi0,pro,pim,pi0,pi0,pi0,pi0,kpl,k0b}, {pip,pro,pip,pip,pim,pim,pim,kpl,kmi},
   {pi0,pro,pip,pip,pim,pim,pi0,kpl,kmi}, {pi0,pro,pip,pim,pi0,pi0,pi0,kpl,kmi},
   {pip,pro,pip,pip,pim,pim,pim,k0,k0b},  {pi0,pro,pip,pip,pim,pim,pi0,k0,k0b},
   {pi0,pro,pip,pim,pi0,pi0,pi0,k0,k0b},  {pi0,pro,pip,pip,pip,pim,pim,kmi,k0},
   {pi0,pro,pip,pip,pim,pi0,pi0,kmi,k0},  {pi0,pro,pip,pi0,pi0,pi0,pi0,kmi,k0},
   {pip,neu,pip,pip,pim,pim,pim,kpl,k0b}, {pi0,neu,pip,pip,pim,pim,pi0,kpl,k0b},
   {pi0,neu,pip,pim,pi0,pi0,pi0,kpl,k0b}, {pi0,neu,pip,pip,pip,pim,pim,kpl,kmi},
   {pi0,neu,pip,pip,pim,pi0,pi0,kpl,kmi}, {pi0,neu,pip,pi0,pi0,pi0,pi0,kpl,kmi},
   {pi0,neu,pip,pip,pip,pim,pim,k0,k0b},  {pi0,neu,pip,pip,pim,pi0,pi0,k0,k0b},
   {pi0,neu,pip,pi0,pi0,pi0,pi0,k0,k0b},  {pip,neu,pip,pip,pip,pim,pim,kmi,k0},
   {pi0,neu,pip,pip,pip,pim,pi0,kmi,k0},  {pi0,neu,pip,pip,pi0,pi0,pi0,kmi,k0}};
}

namespace {
  // pi0 n : Outgoing particle types of a given multiplicity
  static const G4int pizN2bfs[5][2] =
  {{pi0,neu}, {pim,pro}, {k0,lam}, {k0,s0}, {kpl,sm}};

  static const G4int pizN3bfs[13][3] =
  {{pip,neu,pim}, {pi0,neu,pi0}, {pi0,pro,pim}, {pi0,k0,lam},
   {pim,kpl,lam}, {pi0,k0,s0},   {pim,kpl,s0},  {pi0,kpl,sm},  
   {pip,k0,sm},   {pim,k0,sp},   {k0b,neu,k0},  {kpl,neu,kmi},
   {k0,pro,kmi}};

  static const G4int pizN4bfs[22][4] =
  {{pi0,neu,pip,pim}, {pi0,neu,pi0,pi0}, {pip,pro,pim,pim},
   {pi0,pro,pim,pi0}, {pip,pim,k0,lam},  {pip,pim,k0,s0},
   {pi0,pi0,k0,lam},  {pi0,pi0,k0,s0},   {pi0,pi0,kpl,sm},
   {pi0,pim,kpl,lam}, {pi0,pim,kpl,s0},  {pip,pim,kpl,sm},
   {pi0,pip,k0,sm},   {pi0,pim,k0,sp},   {pim,pim,kpl,sp},
   {pi0,neu,k0,k0b},  {pi0,neu,kpl,kmi}, {pim,neu,kpl,k0b},
   {pip,neu,k0,kmi},  {pim,pro,k0,k0b},  {pim,pro,kpl,kmi},
   {pi0,pro,k0,kmi}};

  static const G4int pizN5bfs[31][5] =
  {{pip,neu,pip,pim,pim}, {pi0,neu,pip,pim,pi0}, {pi0,neu,pi0,pi0,pi0},
   {pi0,pro,pip,pim,pim}, {pi0,pro,pim,pi0,pi0}, {pi0,pip,pim,k0,lam},
   {pi0,pi0,pi0,k0,lam},  {pip,pim,pim,kpl,lam}, {pi0,pi0,pim,kpl,lam},
   {pi0,pip,pim,k0,s0},   {pi0,pi0,pi0,k0,s0},   {pip,pim,pim,kpl,s0},
   {pi0,pi0,pim,kpl,s0},  {pi0,pip,pim,kpl,sm},  {pi0,pi0,pi0,kpl,sm},
   {pip,pip,pim,k0,sm},   {pi0,pip,pi0,k0,sm},   {pip,pim,pim,k0,sp},
   {pi0,pi0,pim,k0,sp},   {pi0,pim,pim,kpl,sp},  {pi0,neu,pi0,k0,k0b},
   {pi0,neu,pi0,kpl,kmi}, {pi0,neu,pip,kmi,k0},  {pi0,neu,pim,kpl,k0b},
   {pip,neu,pim,k0,k0b},  {pip,neu,pim,kpl,kmi}, {pi0,pro,pim,k0,k0b},
   {pi0,pro,pim,kpl,kmi}, {pi0,pro,pi0,kmi,k0},  {pim,pro,pim,kpl,k0b},
   {pip,pro,pim,kmi,k0}};

  static const G4int pizN6bfs[39][6] =
  {{pi0,neu,pip,pip,pim,pim}, {pi0,neu,pip,pim,pi0,pi0},
   {pi0,neu,pi0,pi0,pi0,pi0}, {pip,pro,pip,pim,pim,pim},
   {pi0,pro,pip,pim,pim,pi0}, {pi0,pro,pim,pi0,pi0,pi0},
   {pip,pip,pim,pim,k0,lam},  {pi0,pi0,pip,pim,k0,lam},
   {pi0,pi0,pi0,pi0,k0,lam},  {pi0,pip,pim,pim,kpl,lam},
   {pi0,pim,pi0,pi0,kpl,lam}, {pip,pip,pim,pim,k0,s0},
   {pi0,pip,pim,pi0,k0,s0},   {pi0,pi0,pi0,pi0,k0,s0},
   {pi0,pip,pim,pim,kpl,s0},  {pi0,pi0,pi0,pim,kpl,s0},
   {pi0,pip,pip,pim,k0,sm},   {pi0,pi0,pi0,pip,k0,sm},
   {pip,pip,pim,pim,kpl,sm},  {pi0,pip,pim,pi0,kpl,sm},
   {pi0,pip,pim,pim,k0,sp},   {pi0,pi0,pi0,pim,k0,sp},
   {pip,pim,pim,pim,kpl,sp},  {pi0,pim,pim,pi0,kpl,sp},
   {pi0,neu,pip,pi0,kmi,k0},  {pip,neu,pip,pim,kmi,k0},
   {pi0,neu,pip,pim,k0,k0b},  {pi0,neu,pi0,pi0,k0,k0b},
   {pi0,neu,pip,pim,kpl,kmi}, {pi0,neu,pi0,pi0,kpl,kmi},
   {pip,neu,pim,pim,kpl,k0b}, {pi0,neu,pim,pi0,kpl,k0b},
   {pi0,pro,pip,pim,kmi,k0},  {pi0,pro,pi0,pi0,kmi,k0},
   {pip,pro,pim,pim,k0,k0b},  {pi0,pro,pim,pi0,k0,k0b},
   {pip,pro,pim,pim,kpl,kmi}, {pi0,pro,pim,pi0,kpl,kmi},
   {pi0,pro,pim,pim,kpl,k0b}};

  static const G4int pizN7bfs[46][7] =
  {{pip,neu,pip,pip,pim,pim,pim}, {pi0,neu,pip,pip,pim,pim,pi0},
   {pi0,neu,pip,pim,pi0,pi0,pi0}, {pi0,neu,pi0,pi0,pi0,pi0,pi0},
   {pi0,pro,pip,pip,pim,pim,pim}, {pi0,pro,pip,pim,pim,pi0,pi0},
   {pi0,pro,pim,pi0,pi0,pi0,pi0}, {pi0,pip,pip,pim,pim,k0,lam},
   {pi0,pi0,pi0,pip,pim,k0,lam},  {pip,pip,pim,pim,pim,kpl,lam},
   {pi0,pi0,pip,pim,pim,kpl,lam}, {pi0,pi0,pi0,pi0,pim,kpl,lam},
   {pi0,pip,pip,pim,pim,k0,s0},   {pi0,pi0,pi0,pip,pim,k0,s0},
   {pip,pip,pim,pim,pim,kpl,s0},  {pi0,pi0,pip,pim,pim,kpl,s0},
   {pi0,pi0,pi0,pi0,pim,kpl,s0},  {pip,pip,pip,pim,pim,k0,sm},
   {pi0,pi0,pip,pip,pim,k0,sm},   {pi0,pi0,pi0,pi0,pip,k0,sm},
   {pi0,pip,pip,pim,pim,kpl,sm},  {pi0,pi0,pi0,pip,pim,kpl,sm},
   {pip,pip,pim,pim,pim,k0,sp},   {pi0,pi0,pip,pim,pim,k0,sp},
   {pi0,pi0,pi0,pi0,pim,k0,sp},   {pi0,pip,pim,pim,pim,kpl,sp},
   {pi0,pi0,pi0,pim,pim,kpl,sp},  {pi0,neu,pip,pi0,pi0,kmi,k0},
   {pi0,neu,pip,pip,pim,kmi,k0},  {pi0,neu,pi0,pi0,pi0,k0,k0b},
   {pi0,neu,pip,pim,pi0,k0,k0b},  {pip,neu,pip,pim,pim,k0,k0b},
   {pi0,neu,pi0,pi0,pi0,kpl,kmi}, {pi0,neu,pip,pim,pi0,kpl,kmi},
   {pip,neu,pip,pim,pim,kpl,kmi}, {pi0,neu,pim,pi0,pi0,kpl,k0b},
   {pi0,neu,pip,pim,pim,kpl,k0b}, {pi0,pro,pi0,pi0,pi0,kmi,k0},
   {pi0,pro,pip,pim,pi0,kmi,k0},  {pip,pro,pip,pim,pim,kmi,k0},
   {pi0,pro,pim,pi0,pi0,k0,k0b},  {pi0,pro,pip,pim,pim,k0,k0b},
   {pi0,pro,pim,pi0,pi0,kpl,kmi}, {pi0,pro,pip,pim,pim,kpl,kmi},
   {pi0,pro,pim,pim,pi0,kpl,k0b}, {pip,pro,pim,pim,pim,kpl,k0b}};

  static const G4int pizN8bfs[51][8] =
  {{pi0,neu,pip,pip,pip,pim,pim,pim}, {pi0,neu,pip,pip,pim,pim,pi0,pi0},
   {pi0,neu,pip,pim,pi0,pi0,pi0,pi0}, {pi0,neu,pi0,pi0,pi0,pi0,pi0,pi0},
   {pip,pro,pip,pip,pim,pim,pim,pim}, {pi0,pro,pip,pip,pim,pim,pim,pi0},
   {pi0,pro,pip,pim,pim,pi0,pi0,pi0}, {pi0,pro,pim,pi0,pi0,pi0,pi0,pi0},
   {pip,pip,pip,pim,pim,pim,k0,lam},  {pi0,pi0,pip,pip,pim,pim,k0,lam},
   {pi0,pi0,pi0,pi0,pip,pim,k0,lam},  {pi0,pip,pip,pim,pim,pim,kpl,lam},
   {pi0,pi0,pi0,pip,pim,pim,kpl,lam}, {pi0,pi0,pi0,pi0,pi0,pim,kpl,lam},
   {pip,pip,pip,pim,pim,pim,k0,s0},   {pi0,pi0,pip,pip,pim,pim,k0,s0},
   {pi0,pi0,pi0,pi0,pip,pim,k0,s0},   {pi0,pip,pip,pim,pim,pim,kpl,s0},
   {pi0,pi0,pi0,pip,pim,pim,kpl,s0},  {pi0,pip,pip,pip,pim,pim,k0,sm},
   {pi0,pi0,pi0,pip,pip,pim,k0,sm},   {pi0,pi0,pi0,pi0,pi0,pip,k0,sm},
   {pip,pip,pip,pim,pim,pim,kpl,sm},  {pi0,pi0,pip,pip,pim,pim,kpl,sm},
   {pi0,pi0,pi0,pi0,pip,pim,kpl,sm},  {pi0,pip,pip,pim,pim,pim,k0,sp},
   {pi0,pi0,pi0,pip,pim,pim,k0,sp},   {pi0,pi0,pi0,pi0,pi0,pim,k0,sp},
   {pi0,pi0,pip,pim,pim,pim,kpl,sp},  {pip,pip,pim,pim,pim,pim,kpl,sp},
   {pi0,pi0,pi0,pi0,pim,pim,kpl,sp},  {pip,neu,pip,pip,pim,pim,kmi,k0},
   {pi0,pi0,neu,pip,pip,pim,kmi,k0},  {pi0,neu,pip,pi0,pi0,pi0,kmi,k0},
   {pi0,neu,pip,pip,pim,pim,k0,k0b},  {pi0,neu,pip,pim,pi0,pi0,k0,k0b},
   {pi0,neu,pip,pip,pim,pim,kpl,kmi}, {pi0,neu,pip,pim,pi0,pi0,kpl,kmi},
   {pip,neu,pip,pim,pim,pim,kpl,k0b}, {pi0,neu,pip,pim,pim,pi0,kpl,k0b},
   {pi0,neu,pim,pi0,pi0,pi0,kpl,k0b}, {pi0,pro,pip,pip,pim,pim,kmi,k0},
   {pi0,pro,pip,pim,pi0,pi0,kmi,k0},  {pip,pro,pip,pim,pim,pim,k0,k0b},
   {pi0,pro,pip,pim,pim,pi0,k0,k0b},  {pi0,pro,pim,pi0,pi0,pi0,k0,k0b},
   {pip,pro,pip,pim,pim,pim,kpl,kmi}, {pi0,pro,pip,pim,pim,pi0,kpl,kmi},
   {pi0,pro,pim,pi0,pi0,pi0,kpl,kmi}, {pi0,pro,pip,pim,pim,pim,kpl,k0b},
   {pi0,pro,pim,pim,pi0,pi0,kpl,k0b}};

  static const G4int pizN9bfs[58][9] =
  {{pip,neu,pip,pip,pip,pim,pim,pim,pim}, {pi0,neu,pip,pip,pip,pim,pim,pim,pi0},
   {pi0,neu,pip,pip,pim,pim,pi0,pi0,pi0}, {pi0,neu,pip,pim,pi0,pi0,pi0,pi0,pi0},
   {pi0,neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {pi0,pro,pip,pip,pip,pim,pim,pim,pim},
   {pi0,pro,pip,pip,pim,pim,pim,pi0,pi0}, {pi0,pro,pip,pim,pim,pi0,pi0,pi0,pi0},
   {pi0,pro,pim,pi0,pi0,pi0,pi0,pi0,pi0}, {pi0,pip,pip,pip,pim,pim,pim,k0,lam},
   {pi0,pip,pip,pim,pim,pi0,pi0,k0,lam},  {pi0,pi0,pi0,pi0,pi0,pip,pim,k0,lam},
   {pip,pip,pip,pim,pim,pim,pim,kpl,lam}, {pi0,pi0,pip,pip,pim,pim,pim,kpl,lam},
   {pi0,pi0,pi0,pi0,pip,pim,pim,kpl,lam}, {pi0,pi0,pi0,pi0,pi0,pi0,pim,kpl,lam},
   {pi0,pip,pip,pip,pim,pim,pim,k0,s0},   {pi0,pi0,pi0,pip,pip,pim,pim,k0,s0},
   {pi0,pi0,pi0,pi0,pi0,pip,pim,k0,s0},   {pip,pip,pip,pim,pim,pim,pim,kpl,s0},
   {pi0,pi0,pip,pip,pim,pim,pim,kpl,s0},  {pi0,pi0,pi0,pi0,pip,pim,pim,kpl,s0},
   {pip,pip,pip,pip,pim,pim,pim,k0,sm},   {pi0,pi0,pip,pip,pip,pim,pim,k0,sm},
   {pi0,pi0,pi0,pi0,pip,pip,pim,k0,sm},   {pi0,pip,pip,pip,pim,pim,pim,kpl,sm},
   {pi0,pi0,pi0,pip,pip,pim,pim,kpl,sm},  {pi0,pi0,pi0,pi0,pi0,pip,pim,kpl,sm},
   {pip,pip,pip,pim,pim,pim,pim,k0,sp},   {pi0,pi0,pip,pip,pim,pim,pim,k0,sp},
   {pi0,pi0,pi0,pi0,pip,pim,pim,k0,sp},   {pi0,pip,pip,pim,pim,pim,pim,kpl,sp},
   {pi0,pi0,pi0,pip,pim,pim,pim,kpl,sp},  {pi0,pi0,pi0,pi0,pi0,pim,pim,kpl,sp},
   {pi0,neu,pip,pip,pip,pim,pim,kmi,k0},  {pi0,neu,pip,pip,pim,pi0,pi0,kmi,k0},
   {pi0,neu,pip,pi0,pi0,pi0,pi0,kmi,k0},  {pip,neu,pip,pip,pim,pim,pim,k0,k0b},
   {pi0,neu,pip,pip,pim,pim,pi0,k0,k0b},  {pi0,neu,pip,pim,pi0,pi0,pi0,k0,k0b},
   {pip,neu,pip,pip,pim,pim,pim,kpl,kmi}, {pi0,neu,pip,pip,pim,pim,pi0,kpl,kmi},
   {pi0,neu,pip,pim,pi0,pi0,pi0,kpl,kmi}, {pi0,neu,pip,pip,pim,pim,pim,kpl,k0b},
   {pi0,neu,pip,pim,pim,pi0,pi0,kpl,k0b}, {pi0,neu,pim,pi0,pi0,pi0,pi0,kpl,k0b},
   {pip,pro,pip,pip,pim,pim,pim,kmi,k0},  {pi0,pro,pip,pip,pim,pim,pi0,kmi,k0},
   {pi0,pro,pip,pim,pi0,pi0,pi0,kmi,k0},  {pi0,pro,pip,pip,pim,pim,pim,k0,k0b},
   {pi0,pro,pip,pim,pim,pi0,pi0,k0,k0b},  {pi0,pro,pim,pi0,pi0,pi0,pi0,k0,k0b},
   {pi0,pro,pip,pip,pim,pim,pim,kpl,kmi}, {pi0,pro,pip,pim,pim,pi0,pi0,kpl,kmi},
   {pi0,pro,pim,pi0,pi0,pi0,pi0,kpl,kmi}, {pro,pip,pip,pim,pim,pim,pim,kpl,k0b},
   {pi0,pro,pip,pim,pim,pim,pi0,kpl,k0b}, {pi0,pro,pim,pim,pi0,pi0,pi0,kpl,k0b}};
}

namespace {
  // Total pi0 p (pi0 n) cross section as a function of kinetic energy
  static const G4double pizPtotXSec[30] =
   { 6.43,   7.18,   7.54,   8.01,   8.52,   9.13,  10.22,  14.37,  20.96,  34.73,
    61.07,  98.23,  61.75,  31.85,  26.338, 30.49,  35.048, 39.756, 34.577, 30.562,
    29.595, 28.915, 28.459, 26.782, 25.637, 24.937, 24.391, 24.137, 23.895, 23.805};

  // pi0 p (pi0 n) cross sections as functions of kinetic energy and multiplicity
  static const G4double pizPCrossSections[265][30] = {
  //
  // multiplicity 2 (5 channels)
  //
  // p pi0 (n pi0) 
   { 1.73,  2.28,  2.44,  2.71,  3.02,   3.43,  3.92,  5.37,  9.96, 17.73,
    29.07, 50.23, 33.68, 16.69, 12.60,  11.89, 13.49, 15.24, 12.39,  9.59,
     7.92,  6.97,  6.13,  5.37,  4.82,   4.5,   4.3,   4.0,   3.8,   3.6},

  // n pi+ (p pi-)
   { 4.7,  4.9,  5.1,  5.3,  5.5,  5.7,  6.3,  9.0, 11.0,  17.0,
    32.0, 48.0, 28.0, 14.5, 11.04, 8.99, 4.79, 5.02, 2.08,  0.94,
     0.5,  0.25, 0.15, 0.09, 0.06, 0.04, 0.03, 0.02, 0.014, 0.01},

  // L K+ (L K0)
   { 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.65,  0.29,  0.17,
     0.15, 0.08, 0.05, 0.032, 0.021, 0.015, 0.011, 0.008, 0.006, 0.005},

  // S0 K+ (S0 K0)
   { 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,    0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.28,   0.19,  0.12,
     0.11, 0.06, 0.04, 0.026, 0.018, 0.013, 0.009, 0.006,  0.004, 0.003},

  // S+ K0 (S- K+)
   { 0.0,  0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.20, 0.25, 0.09,
     0.04, 0.015, 0.007, 0.003, 0.001, 0.0,  0.0,  0.0,  0.0,  0.0},
  //
  // multiplicity 3 (13 channels)
  //
  // p pi+ pi- (n pi+ pi-)
   {0.0,  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0, 0.03, 0.19, 0.73, 3.4,  7.01, 8.35, 8.9,  5.69,
    4.01, 2.7, 2.0,  1.30, 0.9,  0.68, 0.48, 0.34, 0.27, 0.19},

  // p pi0 pi0 (n pi0 pi0)
   { 0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.025, 0.32, 1.1,  2.2,   3.0,   2.3,   1.4,   0.7,
     0.39, 0.21, 0.11,  0.06, 0.03, 0.015, 0.008, 0.004, 0.002, 0.001},

  // n pi+ pi0 (p pi- pi0)
   { 0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,
     0.0,  0.0,  0.015, 0.15, 0.82, 3.4,  5.4,  5.5, 5.1,  5.0,
     3.3,  2.24, 1.65,  1.1,  0.8,  0.59, 0.42, 0.3, 0.22, 0.15},

  // L K+ pi0 (L K0 pi0)
   { 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.043, 0.143,
     0.179, 0.143, 0.104, 0.074, 0.051, 0.035, 0.026, 0.018, 0.013, 0.009},

  // L K0 pi+ (L K+ pi-)
   { 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.038, 0.11,
     0.135, 0.11, 0.08, 0.055, 0.039, 0.028, 0.02, 0.014, 0.01,  0.007},

  // S0 K+ pi0 (S0 K0 pi0)
   { 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.03,  0.08,
     0.07, 0.05, 0.03, 0.019, 0.012, 0.008, 0.005, 0.003, 0.002, 0.001},

  // S0 K0 pi+ (S0 K+ pi-)
   { 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.04,
     0.08, 0.05, 0.03, 0.019, 0.012, 0.008, 0.005, 0.003, 0.002, 0.001},

  // S+ K0 pi0 (S- K+ pi0)
  {  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0, 0.0, 0.0,
     0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0, 0.0, 0.04,
     0.07, 0.04, 0.02, 0.009, 0.004, 0.002, 0.001, 0.0, 0.0, 0.0},

  // S+ K+ pi- (S- K0 pi+)
   {0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.09,
    0.11, 0.07, 0.05, 0.04, 0.03, 0.02, 0.015, 0.01, 0.007, 0.005},

  // S- K+ pi+ (S+ K0 pi-)
  { 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0, 0.0, 0.0,
    0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0, 0.0, 0.04,
    0.06, 0.04, 0.02, 0.009, 0.004, 0.002, 0.001, 0.0, 0.0, 0.0},

  // p K+ K- (n K0 K0b)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.03,
     0.10, 0.15, 0.16, 0.10, 0.05, 0.025, 0.014, 0.006, 0.003, 0.002},

  // p K0 K0b (n K+ K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.03,
     0.12, 0.18, 0.18, 0.12, 0.07, 0.035, 0.02, 0.01, 0.005, 0.003},

  // n K+ K0b (p K0 K-)
   { 0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.03,
     0.07, 0.07, 0.055, 0.04, 0.03, 0.022, 0.018, 0.013, 0.01, 0.008},
  //
  // multiplicity 4 (22 channels)
  //
  // p pi+ pi- pi0 (n pi+ pi- pi0)
   {0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0, 0.0, 0.04, 0.14, 0.44, 1.11, 1.88,
    2.05, 2.07, 1.75, 1.5, 1.3, 1.1,  0.95, 0.80, 0.72, 0.66},

  // p pi0 pi0 pi0 (n pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.04, 0.5,  0.95, 0.9,  0.8,  0.71,
     0.6, 0.5, 0.43, 0.35, 0.3,  0.25, 0.21, 0.17, 0.14, 0.115},

  // n pi+ pi+ pi- (p pi+ pi- pi-)
   { 0.0, 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.004, 0.035, 0.13, 0.39, 0.75, 1.4,
     1.8, 1.85, 1.65, 1.45, 1.3,   1.2,   1.05, 0.95, 0.83, 0.75},

  // n pi+ pi0 pi0 (p pi- pi0 pi0)
   { 0.0, 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.004, 0.035, 0.13, 0.39, 0.75, 1.4,
     1.8, 1.85, 1.65, 1.45, 1.3,   1.2,   1.05, 0.95, 0.83, 0.75},

  // L K+ pi+ pi- (L K0 pi+ pi-)
   { 0.0,  0.0, 0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.012,
     0.06, 0.1, 0.09, 0.08, 0.065, 0.055, 0.047, 0.039, 0.035, 0.03},

  // S0 K+ pi+ pi- (S0 K0 pi+ pi-)
   {0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
    0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
    0.02, 0.03, 0.03, 0.025, 0.021, 0.018, 0.015, 0.013, 0.011, 0.009},

  // L K+ pi0 pi0 (L K0 pi0 pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.012,
     0.06, 0.1, 0.09, 0.08, 0.065, 0.055, 0.047, 0.039, 0.035, 0.03},

  // S0 K+ pi0 pi0 (S0 K0 pi0 pi0)
  {  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0,   0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.009,
     0.035, 0.05, 0.05, 0.043, 0.036, 0.03, 0.025, 0.02, 0.017, 0.014},

  // S+ K0 pi0 pi0 (S- K+ pi0 pi0)
   { 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.007,
     0.03, 0.045, 0.04, 0.035, 0.03, 0.025, 0.022, 0.018, 0.015, 0.012},

  // L K0 pi+ pi0 (L K+ pi- pi0)
   { 0.0,  0.0, 0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,  0.012,
     0.06, 0.1, 0.09, 0.075, 0.06, 0.05, 0.042, 0.035, 0.03, 0.025},

  // S0 K0 pi+ pi0 (S0 K+ pi- pi0)
   { 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.007,
     0.03, 0.045, 0.04, 0.035, 0.03, 0.025, 0.022, 0.018, 0.015, 0.012},

  // S+ K0 pi+ pi- (S- K+ pi+ pi-)
   {0.0,  0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
    0.0,  0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
    0.02, 0.045, 0.05, 0.045, 0.037, 0.03, 0.024, 0.019, 0.015, 0.012},

  // S+ K+ pi- pi0 (S- K0 pi+ pi0)
   {0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
    0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
    0.02, 0.05, 0.075, 0.06, 0.04, 0.03, 0.024, 0.019, 0.015, 0.012},

  // S- K+ pi+ pi0 (S+ K0 pi- pi0)
   { 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.007,
     0.03, 0.045, 0.04, 0.035, 0.03, 0.025, 0.022, 0.018, 0.015, 0.012},

  // S- K0 pi+ pi+  (S+ K+ pi- pi-)
   { 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.007,
     0.03, 0.045, 0.04, 0.035, 0.03, 0.025, 0.022, 0.018, 0.015, 0.012},

  // p pi0 K+ K- (n pi0 K0 K0b)
   { 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.02, 0.05, 0.07, 0.07, 0.065, 0.06, 0.055, 0.049, 0.045, 0.04},

  // p pi0 K0 K0b (n pi0 K+ K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.02, 0.05, 0.07, 0.07, 0.065, 0.06, 0.055, 0.049, 0.045, 0.04},

  // p pi+ K0 K- (n pi- K+ K0b)
   { 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.02, 0.05, 0.07, 0.07, 0.065, 0.06, 0.055, 0.049, 0.045, 0.04},

  // p pi- K+ K0b (n pi+ K0 K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.02, 0.05, 0.07, 0.07, 0.065, 0.06, 0.055, 0.049, 0.045, 0.04},

  // n pi+ K+ K- (p pi- K0 K0b)
   { 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.02, 0.05, 0.07, 0.07, 0.065, 0.06, 0.055, 0.049, 0.045, 0.04},

  // n pi+ K0 K0b (p pi- K+ K-)
   { 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.047, 0.09, 0.09, 0.083, 0.075, 0.068, 0.057, 0.05, 0.045},

  // n pi0 K+ K0b (p pi0 K0 K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.02, 0.05, 0.07, 0.07, 0.065, 0.06, 0.055, 0.049, 0.045, 0.04},
  //
  // multiplicity 5 (31 channels)
  //
  // p pi+ pi+ pi- pi- (n pi+ pi+ pi- pi-)
   { 0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0,   0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.003, 0.03, 0.12, 0.34,
     0.67, 0.93, 1.15, 1.2, 1.1, 0.94, 0.76,  0.59, 0.47, 0.37},

  // p pi+ pi- pi0 pi0 (n pi+ pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.0,   0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.002, 0.02, 0.10, 0.52,
     1.4, 1.9, 2.1, 1.85, 1.6, 1.35, 1.15,  0.91, 0.74, 0.62},

  // p pi0 pi0 pi0 pi0 (n pi0 pi0 pi0 pi0)
   { 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.003, 0.014, 0.073,
     0.198, 0.269, 0.292, 0.262, 0.227, 0.189, 0.162, 0.128, 0.102, 0.083},

  // n pi+ pi+ pi- pi0 (p pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.0,   0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.002, 0.02, 0.10, 0.61,
     1.5, 1.9, 2.1, 1.85, 1.6, 1.35, 1.15,  0.91, 0.74, 0.62},

  // n pi+ pi0 pi0 pi0 (p pi- pi0 pi0 pi0)
   { 0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  0.0,   0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  0.001, 0.01, 0.05, 0.26,
     0.7, 0.98, 1.05, 0.93, 0.8, 0.67, 0.56,  0.45, 0.37, 0.3},

  // L K+ pi+ pi- pi0 (L K0 pi+ pi- pi0)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
     0.01, 0.04, 0.07, 0.08, 0.08, 0.072, 0.065, 0.06, 0.055, 0.05},

  // L K+ pi0 pi0 pi0 (L K0 pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.006, 0.019, 0.025, 0.022, 0.02, 0.018, 0.015, 0.013, 0.011},

  //  L K0 pi+ pi+ pi- (L K+ pi+ pi- pi-)
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.01, 0.03, 0.04, 0.04, 0.035, 0.03, 0.024, 0.02, 0.017},

  // L K0 pi+ pi0 pi0 (L K+ pi- pi0 pi0)
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.01, 0.03, 0.04, 0.04, 0.035, 0.03, 0.024, 0.02, 0.017},

  // S0 K+ pi+ pi- pi0 (S0 K0 pi+ pi- pi0)
   { 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.01, 0.03, 0.04, 0.04, 0.035, 0.03, 0.024, 0.02, 0.017},

  // S0 K+ pi0 pi0 pi0 (S0 K0 pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.004, 0.006, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002},

  // S0 K0 pi+ pi+ pi- (S0 K+ pi+ pi- pi-)
   { 0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.005, 0.015, 0.02, 0.02, 0.017, 0.015, 0.012, 0.01, 0.008},

  // S0 K0 pi+ pi0 pi0 (S0 K+ pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.005, 0.015, 0.02, 0.02, 0.017, 0.015, 0.012, 0.01, 0.008},

  // S+ K0 pi+ pi- pi0 (S- K+ pi+ pi- pi0)
   { 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.006, 0.02, 0.027, 0.025, 0.023, 0.021, 0.019, 0.017, 0.015},

  // S+ K0 pi0 pi0 pi0 (S- K+ pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.003, 0.005, 0.005, 0.004, 0.003, 0.002, 0.001, 0.0},

  // S+ K+ pi+ pi- pi- (S- K0 pi+ pi+ pi-)
   { 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.003, 0.01, 0.013, 0.012, 0.011, 0.01, 0.009, 0.008, 0.007},

  // S+ K+ pi- pi0 pi0 (S- K0 pi+ pi0 pi0)
   { 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.003, 0.01, 0.013, 0.012, 0.011, 0.01, 0.009, 0.008, 0.007},

  // S- K+ pi+ pi+ pi- (S+ K0 pi+ pi- pi-)
   { 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.004, 0.01, 0.012, 0.011, 0.01, 0.009, 0.008, 0.007},

  // S- K+ pi+ pi0 pi0 (S+ K0 pi- pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.004, 0.01, 0.012, 0.011, 0.01, 0.009, 0.008, 0.007},

  // S- K0 pi+ pi+ pi0 (S+ K+ pi- pi- pi0)
   { 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0, 0.004, 0.01, 0.012, 0.011, 0.01, 0.009, 0.008, 0.007},

  // p K+ K- pi0 pi0 (n K0 K0b pi0 pi0)
   { 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.004, 0.02, 0.025, 0.025, 0.022, 0.02, 0.017, 0.015, 0.012},

  // p K0 K0b pi0 pi0 (n K+ K- pi0 pi0)
   { 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.004, 0.02, 0.025, 0.025, 0.022, 0.02, 0.017, 0.015, 0.012},
 
  // p K+ K0b pi- pi0 (n K0 K- pi+ pi0)
   { 0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.008, 0.04, 0.05, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025},

  // p K0 K- pi+ pi0 (n K+ K0b pi- pi0)
   { 0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.008, 0.04, 0.05, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025},

  // p K+ K- pi+ pi- (n K0 K0b pi+ pi-)
   { 0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.008, 0.04, 0.05, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025},

  // p K0 K0b pi+ pi- (n K+ K- pi+ pi-)
   { 0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
     0.0, 0.008, 0.04, 0.05, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025},

  // n K+ K- pi+ pi0 (p K0 K0b pi- pi0)
   { 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.008, 0.04, 0.065, 0.065, 0.06, 0.055, 0.05, 0.045, 0.04},

  // n K0 K0b pi+ pi0 (p K+ K- pi- pi0)
   { 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.008, 0.04, 0.065, 0.065, 0.06, 0.055, 0.05, 0.045, 0.04},

  // n K+ K0b pi0 pi0 (p K0 K- pi0 pi0)
   { 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.004, 0.02, 0.032, 0.032, 0.03, 0.027, 0.025, 0.022, 0.02},

  // n K0 K- pi+ pi+ (p K+ K0b pi- pi-)
   { 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.004, 0.02, 0.032, 0.032, 0.03, 0.027, 0.025, 0.022, 0.02},

  // n K+ K0b pi+ pi- (p K0 K- pi+ pi-)
   { 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0, 0.008, 0.04, 0.065, 0.065, 0.06, 0.055, 0.05, 0.045, 0.04},
  //
  // multiplicity 6 (39 channels)
  //
  // p pi+ pi+ pi- pi- pi0 (n pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0,   0.0, 0.0, 0.0,  0.0,  0.0,  0.0,   0.0,  0.0,
     0.0, 0.0,   0.0, 0.0, 0.0,  0.0,  0.0,  0.004, 0.02, 0.095,
     0.2, 0.321, 0.4, 0.5, 0.55, 0.57, 0.58, 0.59,  0.6,  0.6},

  // p pi+ pi- pi0 pi0 pi0 (n pi+ pi- pi0 pi0 pi0)
   { 0.0,   0.0,   0.0,  0.0, 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,  0.0, 0.0,  0.0,  0.0,  0.003, 0.013, 0.067,
     0.167, 0.333, 0.43, 0.5, 0.57, 0.62, 0.63, 0.63,  0.63,  0.63},

  // p pi0 pi0 pi0 pi0 pi0 (n pi0 pi0 pi0 pi0 pi0)
   { 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.001, 0.004,
     0.009, 0.017, 0.024, 0.03, 0.036, 0.04, 0.04, 0.04, 0.04,  0.04},

  // n pi+ pi+ pi+ pi- pi- (p pi+ pi+ pi- pi- pi-)
   { 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.001, 0.007, 0.033,
     0.088, 0.167, 0.22, 0.27, 0.33, 0.37, 0.38, 0.39,  0.4,   0.4},

  // n pi+ pi+ pi- pi0 pi0 (p pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0,   0.0, 0.0, 0.0,  0.0,  0.0,  0.0,   0.0,  0.0,
     0.0, 0.0,   0.0, 0.0, 0.0,  0.0,  0.0,  0.004, 0.02, 0.095,
     0.2, 0.321, 0.4, 0.5, 0.55, 0.57, 0.58, 0.59,  0.6,  0.6},

  // n pi+ pi0 pi0 pi0 pi0 (p pi- pi0 pi0 pi0 pi0)
   { 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.001, 0.004, 0.017,
     0.045, 0.083, 0.12, 0.145, 0.165, 0.185, 0.192, 0.195, 0.2,   0.2},

  // L K+ pi+ pi+ pi- pi- (L K0 pi+ pi+ pi- pi-)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.004, 0.012, 0.017, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011},

  // L K+ pi+ pi- pi0 (L K0 pi+ pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.008, 0.025, 0.033, 0.033, 0.03, 0.027, 0.024, 0.021, 0.019},

  // L K+ pi0 pi0 pi0 pi0 (L K0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003, 0.002, 0.002, 0.002},

  // L K0 pi+ pi+ pi- pi0 (L K+ pi+ pi- pi- pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.008, 0.025, 0.033, 0.033, 0.03, 0.027, 0.024, 0.021, 0.019},

  // L K0 pi+ pi0 pi0 pi0 (L K+ pi- pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.002, 0.007, 0.011, 0.011, 0.01, 0.009, 0.008, 0.007, 0.006},

  // S0 K+ pi+ pi+ pi- pi- (S0 K0 pi+ pi+ pi- pi-)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.002, 0.006, 0.009, 0.009, 0.009, 0.008, 0.007, 0.006, 0.006},

  // S0 K+ pi+ pi- pi0 pi0 (S0 K0 pi+ pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.004, 0.012, 0.017, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011},

  // S0 K+ pi0 pi0 pi0 pi0 (S0 K0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001},

  // S0 K0 pi+ pi+ pi- pi0 (S0 K+ pi+ pi- pi- pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.004, 0.012, 0.017, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011},

  // S0 K0 pi+ pi0 pi0 pi0 (S0 K+ pi- pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.003, 0.005, 0.006, 0.006, 0.006, 0.005, 0.004, 0.003},

  // S+ K+ pi- pi- pi+ pi0 (S- K0 pi+ pi+ pi- pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.004, 0.012, 0.016, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011},

  // S+ K+ pi- pi0 pi0 pi0 (S- K0 pi+ pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.004, 0.006, 0.006, 0.006, 0.005, 0.004, 0.004, 0.003},

  // S+ K0 pi+ pi+ pi- pi- (S- K+ pi+ pi+ pi- pi-)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.002, 0.006, 0.008, 0.008, 0.008, 0.007, 0.006, 0.005, 0.004},

  // S+ K0 pi+ pi- pi0 pi0 (S- K+ pi+ pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.004, 0.012, 0.016, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011},

  // S+ K0 pi0 pi0 pi0 pi0 (S- K+ pi0 pi0 pi0 pi0) (negligible) 

  // S- K+ pi+ pi+ pi- pi0 (S+ K0 pi+ pi- pi- pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.004, 0.012, 0.016, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011},

  // S- K+ pi+ pi0 pi0 pi0 (S+ K0 pi- pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.003, 0.005, 0.005, 0.005, 0.005, 0.004, 0.004, 0.004},

  // S- K0 pi+ pi+ pi+ pi- (S+ K+ pi+ pi- pi- pi-)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.003, 0.005, 0.005, 0.005, 0.005, 0.004, 0.004, 0.004},

  // S- K0 pi+ pi+ pi0 pi0 (S+ K+ pi- pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.002, 0.006, 0.008, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003},

  // p pi- pi0 pi0 K+ K0b (n pi+ pi0 pi0 K- K0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  // p pi+ pi- pi- K+ K0b (n pi+ pi+ pi- K- K0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  // p pi+ pi- pi0 K+ K- (n pi+ pi- pi0 K0 K0b)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.005, 0.026, 0.046, 0.054, 0.054, 0.05, 0.044, 0.041, 0.036},

  // p pi0 pi0 pi0 K+ K- (n pi0 pi0 pi0 K0 K0b)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.005, 0.008, 0.008, 0.008, 0.007, 0.006, 0.005, 0.004},

  // p pi+ pi- pi0 K0 K0b (n pi+ pi- pi0 K+ K-)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.005, 0.026, 0.046, 0.054, 0.054, 0.05, 0.044, 0.041, 0.036},

  // p pi0 pi0 pi0 K0 K0b (n pi0 pi0 pi0 K+ K-)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.005, 0.008, 0.008, 0.008, 0.007, 0.006, 0.005, 0.004},

  // p pi+ pi+ pi- K- K0 (n pi+ pi- pi- K+ K0b)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  // p pi+ pi0 pi0 K- K0 (n pi- pi0 pi0 K+ K0b)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  // n pi+ pi- pi0 K+ K0b (p pi+ pi- pi0 K- K0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.005, 0.026, 0.046, 0.054, 0.054, 0.05, 0.044, 0.041, 0.036},

  // n pi0 pi0 pi0 K+ K0b (p pi0 pi0 pi0 K- K0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.005, 0.008, 0.008, 0.008, 0.007, 0.006, 0.005, 0.004},

  // n pi+ pi+ pi- K+ K- (p pi+ pi- pi- K0 K0b)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  // n pi+ pi0 pi0 K+ K- (p pi- pi0 pi0 K0 K0b)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  // n pi+ pi+ pi- K0 K0b (p pi+ pi- pi- K+ K-)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  // n pi+ pi0 pi0 K0 K0b (p pi- pi0 pi0 K+ K-)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  // n pi+ pi+ pi0 K- K0 (p pi- pi- pi0 K+ K0b)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},
  //
  // multiplicity 7 (46 channels)
  //
  // p pi+ pi+ pi+ pi- pi- pi- (n pi+ pi+ pi+ pi- pi- pi-)
   { 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.001, 0.005,
     0.017, 0.04, 0.085, 0.145, 0.21, 0.265, 0.319, 0.365, 0.38,  0.38},

  // p pi+ pi+ pi- pi- pi0 pi0 (n pi+ pi+ pi- pi- pi0 pi0)
   { 0.0,   0.0,  0.0,   0.0,  0.0, 0.0,  0.0,  0.0,   0.0,   0.0,
     0.0,   0.0,  0.0,   0.0,  0.0, 0.0,  0.0,  0.0,   0.002, 0.013,
     0.045, 0.12, 0.203, 0.31, 0.4, 0.48, 0.53, 0.572, 0.593, 0.604},

  // p pi+ pi- pi0 pi0 pi0 pi0 (n pi+ pi- pi0 pi0 pi0 pi0)
   { 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.001, 0.005,
     0.017, 0.04, 0.085, 0.145, 0.21, 0.265, 0.319, 0.365, 0.38,  0.38},

  // p pi0 pi0 pi0 pi0 pi0 pi0 (n pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.001, 0.003, 0.006, 0.009, 0.012, 0.014, 0.016, 0.018, 0.019, 0.02},

  // n pi+ pi+ pi+ pi- pi- pi0 (p pi+ pi+ pi- pi- pi- pi0)
   { 0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.002, 0.012,
     0.045, 0.12, 0.213, 0.34, 0.45, 0.55, 0.641, 0.739, 0.784, 0.816},

  // n pi+ pi+ pi- pi0 pi0 pi0 (p pi+ pi- pi- pi0 pi0 pi0)
   { 0.0,  0.0,   0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,   0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.001, 0.006,
     0.02, 0.048, 0.1, 0.172, 0.25, 0.315, 0.377, 0.438, 0.451, 0.46},

  // n pi+ pi0 pi0 pi0 pi0 pi0 (p pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.001,
     0.004, 0.011, 0.025, 0.04, 0.057, 0.07, 0.08, 0.091, 0.095, 0.098},

  // L K+ pi+ pi+ pi- pi- pi0 (L K0 pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.002, 0.006, 0.012, 0.018, 0.021, 0.021, 0.021, 0.021, 0.021},

  // L K+ pi+ pi- pi0 pi0 pi0 (L K0 pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
     0.0, 0.001, 0.004, 0.008, 0.012, 0.014,  0.014,  0.014, 0.014, 0.014},

  // L K+ pi0 pi0 pi0 pi0 pi0 (L K0 pi0 pi0 pi0 pi0 pi0) (negligible)

  // L K0 pi+ pi+ pi+ pi- pi- (L K+ pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.002, 0.004, 0.006, 0.007, 0.007, 0.007, 0.007, 0.007},

  // L K0 pi+ pi+ pi- pi0 pi0 (L K+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.002, 0.006, 0.012, 0.018, 0.021, 0.021, 0.021, 0.021, 0.021},

  // L K0 pi+ pi0 pi0 pi0 pi0 (L K+ pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004, 0.004},

  // S0 K+ pi+ pi+ pi- pi- pi0 (S0 K0 pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
     0.0, 0.001, 0.003, 0.006, 0.009, 0.01,  0.01,  0.01, 0.01, 0.01},

  // S0 K+ pi+ pi- pi0 pi0 pi0 (S0 K0 pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.004, 0.006, 0.007,  0.007,  0.007, 0.007, 0.007},

  // S0 K+ pi0 pi0 pi0 pi0 pi0 (S0 K0 pi0 pi0 pi0 pi0 pi0) (negligible)

  // S0 K0 pi+ pi+ pi+ pi- pi- (S0 K+ pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004, 0.004},

  // S0 K0 pi+ pi+ pi- pi0 pi0 (S0 K+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
     0.0, 0.001, 0.003, 0.006, 0.009, 0.01,  0.01,  0.01, 0.01, 0.01},

  // S0 K0 pi+ pi0 pi0 pi0 pi0 (S0 K+ pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.002,  0.002,  0.002, 0.002, 0.002},

  // S+ K+ pi+ pi+ pi- pi- pi- (S- K0 pi+ pi+ pi+ pi- pi-)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004, 0.004},

  // S+ K+ pi+ pi- pi- pi0 pi0 (S- K0 pi+ pi+ pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.003, 0.007, 0.011, 0.012, 0.012, 0.012, 0.012, 0.012},

  // S+ K+ pi- pi0 pi0 pi0 pi0 (S- K0 pi+ pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002},

  // S+ K0 pi+ pi+ pi- pi- pi0 (S- K+ pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.003, 0.007, 0.011, 0.012, 0.012, 0.012, 0.012, 0.012},

  // S+ K0 pi+ pi- pi0 pi0 pi0 (S- K+ pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.002, 0.004, 0.007, 0.008, 0.008, 0.008, 0.008, 0.008},

  // S+ K0 pi0 pi0 pi0 pi0 pi0 (S- K+ pi0 pi0 pi0 pi0 pi0) (negligible) 

  // S- K+ pi+ pi+ pi+ pi- pi- (S+ K0 pi+ pi+ pi- pi- pi-) 
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004, 0.004},

  // S- K+ pi+ pi+ pi- pi0 pi0 (S+ K0 pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
     0.0, 0.001, 0.003, 0.007, 0.011, 0.012,  0.012,  0.012, 0.012, 0.012},

  // S- K+ pi+ pi0 pi0 pi0 pi0 (S+ K0 pi- pi0 pi0 pi0 pi0) 
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002},

  // S- K0 pi+ pi+ pi+ pi- pi0 (S+ K+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.002, 0.004, 0.007, 0.008, 0.008, 0.008, 0.008, 0.008},

  // S- K0 pi+ pi+ pi0 pi0 pi0 (S+ K+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004, 0.004},

  // p K+ K0b pi- pi0 pi0 pi0 (n K- K0 pi+ pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.001, 0.004, 0.008, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01},
 
  // p K+ K0b pi+ pi- pi- pi0 (n K- K0 pi+ pi+ pi- pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.003, 0.013, 0.023, 0.029, 0.03, 0.03, 0.03, 0.03, 0.03},

  // p K+ K- pi0 pi0 pi0 pi0 (n K0 K0b pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002},

  // p K+ K- pi+ pi- pi0 pi0 (n K0 K0b pi+ pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.003, 0.013, 0.023, 0.029, 0.029, 0.028, 0.027, 0.027, 0.027},

  // p K+ K- pi+ pi+ pi- pi- (n K0 K0b pi+ pi+ pi- pi-)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.002, 0.006, 0.012, 0.015, 0.015, 0.014, 0.014, 0.014, 0.014},

  // p K0 K0b pi0 pi0 pi0 pi0 (n K+ K- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002},

  // p K0 K0b pi+ pi- pi0 pi0 (n K+ K- pi+ pi- pi0 pi0) 
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.003, 0.013, 0.023, 0.029, 0.029, 0.027, 0.025, 0.025, 0.025},

  // p K0 K0b pi+ pi+ pi- pi- (n K+ K- pi+ pi+ pi- pi-) 
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.002, 0.006, 0.012, 0.015, 0.015, 0.014, 0.014, 0.014, 0.014},

  // p K- K0 pi+ pi0 pi0 pi0 (n K+ K0b pi- pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.004, 0.008, 0.01, 0.01, 0.009, 0.008, 0.008, 0.008},

  // p K- K0 pi+ pi+ pi- pi0 (n K+ K0b pi+ pi- pi- pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.003, 0.013, 0.023, 0.029, 0.029, 0.027, 0.025, 0.025, 0.025},

  // n K+ K0b pi0 pi0 pi0 pi0 (p K- K0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002},

  // n K+ K0b pi+ pi- pi0 pi0 (p K- K0 pi+ pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.003, 0.013, 0.023, 0.029, 0.029, 0.027, 0.025, 0.025, 0.025},

  // n K+ K0b pi+ pi+ pi- pi- (p K- K0 pi+ pi+ pi- pi-)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.002, 0.006, 0.012, 0.015, 0.015, 0.014, 0.013, 0.013, 0.013},

  // n K+ K- pi+ pi0 pi0 pi0 (p K0 K0b pi- pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.004, 0.008, 0.01, 0.01, 0.009, 0.008, 0.007, 0.007},

  // n K+ K- pi+ pi+ pi- pi0 (p K0 K0b pi+ pi- pi- pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.003, 0.013, 0.023, 0.029, 0.029, 0.027, 0.025, 0.025, 0.025},

  // n K0 K0b pi+ pi0 pi0 pi0 (p K+ K- pi- pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.004, 0.008, 0.01, 0.01, 0.009, 0.008, 0.007, 0.007},

  // n K0 K0b pi+ pi+ pi- pi0 (p K+ K- pi+ pi- pi- pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.003, 0.013, 0.023, 0.029, 0.029, 0.027, 0.025, 0.025, 0.025},

  // n K- K0 pi+ pi+ pi0 pi0 (p K+ K0b pi- pi- pi0 pi0) 
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.002, 0.006, 0.012, 0.015, 0.015, 0.014, 0.013, 0.013, 0.013},

  // n K- K0 pi+ pi+ pi+ pi- (p K+ K0b pi+ pi- pi- pi-)
   { 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.004, 0.008, 0.01, 0.01, 0.009, 0.008, 0.007, 0.007},
  //
  // multiplicity 8 (51 channels)
  //
  // p pi+ pi+ pi+ pi- pi- pi- pi0 (n pi+ pi+ pi+ pi- pi- pi- pi0)
   { 0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0, 0.0,  0.0,   0.0,
     0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0, 0.0,  0.0,   0.002,
     0.013, 0.05, 0.113, 0.17, 0.24, 0.32, 0.4, 0.51, 0.605, 0.69},

  // p pi+ pi+ pi- pi- pi0 pi0 pi0 (n pi+ pi+ pi- pi- pi0 pi0 pi0)
   { 0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
     0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.003,
     0.02, 0.065, 0.14, 0.22, 0.32, 0.45, 0.58, 0.778, 0.902, 1.04},

  // p pi+ pi- pi0 pi0 pi0 pi0 pi0 (n pi+ pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.001,
     0.005, 0.015, 0.034, 0.055, 0.08, 0.11, 0.139, 0.176, 0.196, 0.213},

  // p pi0 pi0 pi0 pi0 pi0 pi0 pi0 (n pi0 pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.006, 0.006},

  // n pi+ pi+ pi+ pi+ pi- pi- pi- (p pi+ pi+ pi+ pi- pi- pi- pi-)
   { 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,
     0.004, 0.018, 0.039, 0.06, 0.081, 0.11, 0.13, 0.152, 0.172, 0.189},

  // n pi+ pi+ pi+ pi- pi- pi0 pi0 (p pi+ pi+ pi- pi- pi- pi0 pi0)
   { 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,  0.0, 0.0,   0.0,   0.0,
     0.0,  0.0,   0.0,  0.0,   0.0,  0.0,  0.0, 0.0,   0.0,   0.003,
     0.02, 0.065, 0.14, 0.235, 0.34, 0.47, 0.6, 0.778, 0.902, 1.04},

  // n pi+ pi+ pi- pi0 pi0 pi0 pi0 (p pi+ pi- pi- pi0 pi0 pi0 pi0)
   { 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.001,
     0.009, 0.032, 0.07, 0.12, 0.18, 0.25, 0.31, 0.389, 0.451, 0.518},

  // n pi+ pi0 pi0 pi0 pi0 pi0 pi0 (p pi- pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.001, 0.003, 0.006, 0.01, 0.014, 0.019, 0.022, 0.028, 0.031, 0.034},

  // L K+ pi+ pi+ pi+ pi- pi- pi- (L K0 pi+ pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.005, 0.005, 0.005},

  // L K+ pi+ pi+ pi- pi- pi0 pi0 (L K0 pi+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.002, 0.005, 0.011, 0.017, 0.021, 0.024, 0.028, 0.032, 0.036},

  // L K+ pi+ pi- pi0 pi0 pi0 pi0 (L K0 pi+ pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.001, 0.003, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011},

  // L K+ pi0 pi0 pi0 pi0 pi0 pi0 (L K0 pi0 pi0 pi0 pi0 pi0 pi0) (negligible)

  // L K0 pi+ pi+ pi+ pi- pi- pi0 (L K+ pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.004, 0.007, 0.011, 0.014, 0.016, 0.019, 0.021, 0.024},

  // L K0 pi+ pi+ pi- pi0 pi0 pi0 (L K+ pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.004, 0.007, 0.011, 0.014, 0.016, 0.019, 0.021, 0.024},

  // L K0 pi+ pi0 pi0 pi0 pi0 pi0 (L K+ pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003, 0.003},

  // S0 K+ pi+ pi+ pi+ pi- pi- pi- (S0 K0 pi+ pi+ pi+ pi- pi- pi-) 
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003, 0.003},

  // S0 K+ pi+ pi+ pi- pi- pi0 pi0 (S0 K0 pi+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.001, 0.003, 0.006, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019},

  // S0 K+ pi+ pi- pi0 pi0 pi0 pi0 (S0 K0 pi+ pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.005, 0.005, 0.005},

  // S0 K+ pi0 pi0 pi0 pi0 pi0 pi0 (S0 K0 pi0 pi0 pi0 pi0 pi0 pi0) (negligible)

  // S0 K0 pi+ pi+ pi+ pi- pi- pi0 (S0 K+ pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.002, 0.004, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011},

  // S0 K0 pi+ pi+ pi- pi0 pi0 pi0 (S0 K+ pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.002, 0.004, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011},

  // S0 K0 pi+ pi0 pi0 pi0 pi0 pi0 (S0 K+ pi- pi0 pi0 pi0 pi0 pi0) (negligible)

  // S+ K+ pi+ pi+ pi- pi- pi- pi0 (S- K0 pi+ pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.005, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018},

  // S+ K+ pi+ pi- pi- pi0 pi0 pi0 (S- K0 pi+ pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.005, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018},

  // S+ K+ pi- pi0 pi0 pi0 pi0 pi0 (S- K0 pi+ pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002},

  // S+ K0 pi+ pi+ pi+ pi- pi- pi- (S- K+ pi+ pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.006, 0.006},

  // S+ K0 pi+ pi+ pi- pi- pi0 pi0 (S- K+ pi+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.003, 0.008, 0.012, 0.015, 0.018, 0.021, 0.024, 0.027},

  // S+ K0 pi+ pi- pi0 pi0 pi0 pi0 (S- K+ pi+ pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009},

  // S+ K0 pi0 pi0 pi0 pi0 pi0 pi0 (S- K+ pi0 pi0 pi0 pi0 pi0 pi0) (negligible)

  // S- K+ pi+ pi+ pi+ pi- pi- pi0 (S+ K0 pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.005, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018},

  // S- K+ pi+ pi+ pi- pi0 pi0 pi0 (S+ K0 pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.005, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018},

  // S- K+ pi+ pi0 pi0 pi0 pi0 pi0 (S+ K0 pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002},

  // S- K0 pi+ pi+ pi+ pi- pi0 pi0 (S+ K+ pi+ pi- pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.005, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018},

  // S- K0 pi+ pi+ pi+ pi+ pi- pi- (S+ K+ pi+ pi+ pi- pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004},

  // S- K0 pi+ pi+ pi0 pi0 pi0 pi0 (S+ K+ pi- pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004},

  // p K+ K0b pi+ pi+ pi- pi- pi- (n K- K0 pi+ pi+ pi+ pi- pi-)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.012, 0.013},

  // p K+ K0b pi+ pi- pi- pi0 pi0 (n K- K0 pi+ pi+ pi- pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.003, 0.01, 0.016, 0.021, 0.025, 0.029, 0.033, 0.037},

  // p K+ K0b pi- pi0 pi0 pi0 pi0 (n K- K0 pi+ pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.005, 0.005},
 
  // p K+ K- pi+ pi+ pi- pi- pi0 (n K0 K0b pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.003, 0.01, 0.016, 0.021, 0.025, 0.029, 0.033, 0.037},

  // p K+ K- pi+ pi- pi0 pi0 pi0 (n K0 K0b pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.002, 0.006, 0.011, 0.015, 0.018, 0.02, 0.022, 0.024},

  // p K+ K- pi0 pi0 pi0 pi0 pi0 (n K0 K0b pi0 pi0 pi0 pi0 pi0) (negligible)

  // p K0 K0b pi+ pi+ pi- pi- pi0 (n K+ K- pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.003, 0.009, 0.016, 0.021, 0.025, 0.029, 0.033, 0.037},

  // p K0 K0b pi+ pi- pi0 pi0 pi0 (n K+ K- pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.002, 0.006, 0.011, 0.015, 0.018, 0.02, 0.022, 0.024},

  // p K0 K0b pi0 pi0 pi0 pi0 pi0 (n K+ K- pi0 pi0 pi0 pi0 pi0) (negligible) 

  // p K- K0 pi+ pi+ pi+ pi- pi- (n K+ K0b pi+ pi+ pi- pi- pi-) 
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.012, 0.013},

  // p K- K0 pi+ pi+ pi- pi0 pi0 (n K+ K0b pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.003, 0.009, 0.016, 0.021, 0.025, 0.029, 0.033, 0.037},

  // p K- K0 pi+ pi0 pi0 pi0 pi0 (n K+ K0b pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.005, 0.005, 0.005},

  // n K+ K0b pi+ pi+ pi- pi- pi0 (p K- K0 pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.002, 0.006, 0.011, 0.015, 0.018, 0.02, 0.022, 0.024},

  // n K+ K0b pi+ pi- pi0 pi0 pi0 (p K- K0 pi+ pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.002, 0.006, 0.011, 0.015, 0.018, 0.02, 0.022, 0.024},

  // n K+ K0b pi0 pi0 pi0 pi0 pi0 (p K- K0 pi0 pi0 pi0 pi0 pi0) (negligible)

  // n K+ K- pi+ pi+ pi+ pi- pi- (p K0 K0b pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.01, 0.011, 0.012},

  // n K+ K- pi+ pi+ pi- pi0 pi0 (p K0 K0b pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.003, 0.01, 0.016, 0.021, 0.025, 0.029, 0.033, 0.037},

  // n K+ K- pi+ pi0 pi0 pi0 pi0 (p K0 K0b pi- pi0 pi0 pi0 pi0) 
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.005, 0.005, 0.005},

  // n K0 K0b pi+ pi+ pi+ pi- pi- (p K+ K- pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.01, 0.011, 0.012},

  // n K0 K0b pi+ pi+ pi- pi0 pi0 (p K+ K- pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.003, 0.01, 0.017, 0.021, 0.025, 0.029, 0.033, 0.037},

  // n K0 K0b pi+ pi0 pi0 pi0 pi0 (p K+ K- pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.003, 0.004, 0.005, 0.006, 0.006, 0.006},

  // n K- K0 pi+ pi+ pi+ pi- pi0 (p K+ K0b pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.002, 0.007, 0.011, 0.015, 0.018, 0.02, 0.022, 0.024},

  // n K- K0 pi+ pi+ pi0 pi0 pi0 (p K+ K0b pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.01, 0.011, 0.012},
  //
  // multiplicity 9 (58 channels)
  //
  // p pi+ pi+ pi+ pi+ pi- pi- pi- pi- (n pi+ pi+ pi+ pi+ pi- pi- pi- pi-)
   { 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0, 0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0, 0.0,   0.0,   0.0,
     0.001, 0.005, 0.018, 0.04, 0.08, 0.13, 0.2, 0.274, 0.332, 0.385},

  // p pi+ pi+ pi+ pi- pi- pi- pi0 pi0 (n pi+ pi+ pi+ pi- pi- pi- pi0 pi0)
   { 0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.002, 0.01, 0.036, 0.08, 0.16, 0.26, 0.39, 0.55, 0.66, 0.77},

  // p pi+ pi+ pi- pi- pi0 pi0 pi0 pi0 (n pi+ pi+ pi- pi- pi0 pi0 pi0 pi0)
   { 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.001, 0.007, 0.024, 0.06, 0.12, 0.195, 0.292, 0.413, 0.499, 0.575},

  // p pi+ pi- pi0 pi0 pi0 pi0 pi0 pi0 (n pi+ pi- pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.001, 0.004, 0.008, 0.016, 0.026, 0.04, 0.055, 0.066, 0.077},

  // p pi0 pi0 pi0 pi0 pi0 pi0 pi0 pi0 (n pi0 pi0 pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004},

  // n pi+ pi+ pi+ pi+ pi- pi- pi- pi0 (p pi+ pi+ pi+ pi- pi- pi- pi- pi0)
   { 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.001, 0.005, 0.018, 0.04, 0.075, 0.13, 0.192, 0.274, 0.332, 0.387},

  // n pi+ pi+ pi+ pi- pi- pi0 pi0 pi0 (p pi+ pi+ pi- pi- pi- pi0 pi0 pi0)
   { 0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.002, 0.01, 0.036, 0.08, 0.15, 0.26, 0.385, 0.547, 0.665, 0.77},

  // n pi+ pi+ pi- pi0 pi0 pi0 pi0 pi0 (p pi+ pi- pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.001, 0.003, 0.011, 0.024, 0.05, 0.087, 0.129, 0.182, 0.222, 0.26},

  // n pi+ pi0 pi0 pi0 pi0 pi0 pi0 pi0 (p pi- pi0 pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.006, 0.007, 0.008, 0.009},

  // L K+ pi+ pi+ pi+ pi- pi- pi- pi0 (L K0 pi+ pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.004, 0.008, 0.012, 0.015, 0.018, 0.021, 0.025},

  // L K+ pi+ pi+ pi- pi- pi0 pi0 pi0 (L K0 pi+ pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.006, 0.011, 0.017, 0.022, 0.027, 0.032, 0.037},

  // L K+ pi+ pi- pi0 pi0 pi0 pi0 pi0 (L K0 pi+ pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.005, 0.005},

  // L K+ pi0 pi0 pi0 pi0 pi0 pi0 pi0 (L K0 pi0 pi0 pi0 pi0 pi0 pi0 pi0) (negligible) 

  // L K0 pi+ pi+ pi+ pi+ pi- pi- pi- (L K+ pi+ pi+ pi+ pi- pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.005, 0.005},

  // L K0 pi+ pi+ pi+ pi- pi- pi0 pi0 (L K+ pi+ pi+ pi- pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.006, 0.011, 0.017, 0.022, 0.027, 0.032, 0.037},

  // L K0 pi+ pi+ pi- pi0 pi0 pi0 pi0 (L K+ pi+ pi- pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.012, 0.014, 0.016, 0.018},

  // L K0 pi+ pi0 pi0 pi0 pi0 pi0 pi0 (L K+ pi- pi0 pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002},

  // S0 K+ pi+ pi+ pi+ pi- pi- pi- pi0 (S0 K0 pi+ pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014},

  // S0 K+ pi+ pi+ pi- pi- pi0 pi0 pi0 (S0 K0 pi+ pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.011, 0.014, 0.016, 0.019},

  // S0 K+ pi+ pi- pi0 pi0 pi0 pi0 pi0 (S0 K0 pi+ pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004},

  // S0 K+ pi0 pi0 pi0 pi0 pi0 pi0 pi0 (S0 K0 pi0 pi0 pi0 pi0 pi0 pi0 pi0) (negligible) 

  // S0 K0 pi+ pi+ pi+ pi+ pi- pi- pi- (S0 K+ pi+ pi+ pi+ pi- pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004},

  // S0 K0 pi+ pi+ pi+ pi- pi- pi0 pi0 (S0 K+ pi+ pi+ pi- pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.011, 0.014, 0.016, 0.019},

  // S0 K0 pi+ pi+ pi- pi0 pi0 pi0 pi0 (S0 K+ pi+ pi- pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.006, 0.007, 0.008, 0.009},

  // S0 K0 pi+ pi0 pi0 pi0 pi0 pi0 pi0 (S0 K+ pi- pi0 pi0 pi0 pi0 pi0 pi0) (negligible)

  // S+ K+ pi+ pi+ pi+ pi- pi- pi- pi- (S- K0 pi+ pi+ pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.005, 0.005},

  // S+ K+ pi+ pi+ pi- pi- pi- pi0 pi0 (S- K0 pi+ pi+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.002, 0.005, 0.008, 0.012, 0.016, 0.02, 0.024, 0.028},

  // S+ K+ pi+ pi- pi- pi0 pi0 pi0 pi0 (S- K0 pi+ pi+ pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014},

  // S+ K+ pi- pi0 pi0 pi0 pi0 pi0 pi0 (S- K0 pi+ pi0 pi0 pi0 pi0 pi0 pi0) (negligible)

  // S+ K0 pi+ pi+ pi+ pi- pi- pi- pi0 (S- K+ pi+ pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.003, 0.005, 0.008, 0.011, 0.014, 0.016, 0.019},

  // S+ K0 pi+ pi+ pi- pi- pi0 pi0 pi0 (S- K+ pi+ pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.002, 0.005, 0.008, 0.012, 0.016, 0.02, 0.024, 0.028},

  // S+ K0 pi+ pi- pi0 pi0 pi0 pi0 pi0 (S- K+ pi+ pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.005, 0.005},

  // S+ K0 pi0 pi0 pi0 pi0 pi0 pi0 pi0 (S- K+ pi0 pi0 pi0 pi0 pi0 pi0 pi0) (negligible)

  // S- K+ pi+ pi+ pi+ pi+ pi- pi- pi- (S+ K0 pi+ pi+ pi+ pi- pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  // S- K+ pi+ pi+ pi+ pi- pi- pi0 pi0 (S+ K0 pi+ pi+ pi- pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.002, 0.006, 0.01, 0.014, 0.016, 0.02, 0.024, 0.028},

  // S- K+ pi+ pi+ pi- pi0 pi0 pi0 pi0 (S+ K0 pi+ pi- pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013, 0.015},

  // S- K+ pi+ pi0 pi0 pi0 pi0 pi0 pi0 (S+ K0 pi- pi0 pi0 pi0 pi0 pi0 pi0) (negligible) 

  // S- K0 pi+ pi+ pi+ pi+ pi- pi- pi0 (S+ K+ pi+ pi+ pi- pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013, 0.015},

  // S- K0 pi+ pi+ pi+ pi- pi0 pi0 pi0 (S+ K+ pi+ pi- pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.011, 0.014, 0.016, 0.019},

  // S- K0 pi+ pi+ pi0 pi0 pi0 pi0 pi0 (S+ K+ pi- pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  // p K+ K0b pi+ pi+ pi- pi- pi- pi0 (n K- K0 pi+ pi+ pi+ pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.011, 0.014, 0.018, 0.022, 0.026},

  // p K+ K0b pi+ pi- pi- pi0 pi0 pi0 (n K- K0 pi+ pi+ pi- pi0 pi0 pi0) 
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.007, 0.011, 0.014, 0.018, 0.022, 0.026},

  // p K+ K0b pi- pi0 pi0 pi0 pi0 pi0 (n K- K0 pi+ pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003},
 
  // p K+ K- pi+ pi+ pi+ pi- pi- pi- (n K0 K0b pi+ pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008},

  // p K+ K- pi+ pi+ pi- pi- pi0 pi0 (n K0 K0b pi+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.003, 0.01, 0.016, 0.021, 0.027, 0.033, 0.039},

  // p K+ K- pi+ pi- pi0 pi0 pi0 pi0 (n K0 K0b pi+ pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013},

  // p K+ K- pi0 pi0 pi0 pi0 pi0 pi0 (n K0 K0b pi0 pi0 pi0 pi0 pi0 pi0) (negligible)

  // p K0 K0b pi+ pi+ pi+ pi- pi- pi- (n K+ K- pi+ pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008},

  // p K0 K0b pi+ pi+ pi- pi- pi0 pi0 (n K+ K- pi+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.003, 0.009, 0.016, 0.021, 0.027, 0.033, 0.039},

  // p K0 K0b pi+ pi- pi0 pi0 pi0 pi0 (n K+ K- pi+ pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013},

  // p K0 K0b pi0 pi0 pi0 pi0 pi0 pi0 (n K+ K- pi0 pi0 pi0 pi0 pi0 pi0) (negligible)

  // p K- K0 pi+ pi+ pi+ pi- pi- pi0 (n K+ k0b pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.006, 0.011, 0.014, 0.018, 0.022, 0.026},

  // p K- K0 pi+ pi+ pi- pi0 pi0 pi0 (n K+ K0b pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.006, 0.011, 0.014, 0.018, 0.022, 0.026},

  // p K- K0 pi+ pi0 pi0 pi0 pi0 pi0 (n K+ K0b pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003},

  // n K+ K0b pi+ pi+ pi+ pi- pi- pi- (p K- K0 pi+ pi+ pi+ pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008},

  // n K+ K0b pi+ pi+ pi- pi- pi0 pi0 (p K- K0 pi+ pi+ pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.003, 0.009, 0.016, 0.021, 0.027, 0.033, 0.039},

  // n K+ K0b pi+ pi- pi0 pi0 pi0 pi0 (p K- K0 pi+ pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013},

  // n K+ K0b pi0 pi0 pi0 pi0 pi0 pi0 (p K- K0 pi0 pi0 pi0 pi0 pi0 pi0) (negligible) 

  // n K+ K- pi+ pi+ pi+ pi- pi- pi0 (p K0 K0b pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.006, 0.011, 0.014, 0.018, 0.022, 0.026},

  // n K+ K- pi+ pi+ pi- pi0 pi0 pi0 (p K0 K0b pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.006, 0.011, 0.014, 0.018, 0.022, 0.026},

  // n K+ K- pi+ pi0 pi0 pi0 pi0 pi0 (p K0 K0b pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003},

  // n K0 K0b pi+ pi+ pi+ pi- pi- pi0 (p K+ K- pi+ pi+ pi- pi- pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.006, 0.011, 0.014, 0.018, 0.022, 0.026},

  // n K0 K0b pi+ pi+ pi- pi0 pi0 pi0 (p K+ K- pi+ pi- pi- pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.006, 0.011, 0.014, 0.018, 0.022, 0.026},

  // n K0 K0b pi+ pi0 pi0 pi0 pi0 pi0 (p K+ K- pi- pi0 pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003},

  // n K- K0 pi+ pi+ pi+ pi+ pi- pi- (p K+ K0b pi+ pi+ pi- pi- pi- pi-)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007},

  // n K- K0 pi+ pi+ pi+ pi- pi0 pi0 (p K+ K0b pi+ pi- pi- pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.006, 0.011, 0.014, 0.018, 0.022, 0.026},

  // n K- K0 pi+ pi+ pi0 pi0 pi0 pi0 (p K+ K0b pi- pi- pi0 pi0 pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007}};

}

// Initialize both |T Tz> = |1/2 1/2> channels, using pizP cross-section table

const G4CascadePiZeroPChannelData::data_t
G4CascadePiZeroPChannelData::data(pizP2bfs, pizP3bfs, pizP4bfs,
				  pizP5bfs, pizP6bfs, pizP7bfs,
				  pizP8bfs, pizP9bfs, pizPCrossSections,
				  pizPtotXSec, pi0*pro, "PiZeroP");

const G4CascadePiZeroNChannelData::data_t
G4CascadePiZeroNChannelData::data(pizN2bfs, pizN3bfs, pizN4bfs,
				  pizN5bfs, pizN6bfs, pizN7bfs,
				  pizN8bfs, pizN9bfs, pizPCrossSections,
				  pizPtotXSec, pi0*neu, "PiZeroN");
