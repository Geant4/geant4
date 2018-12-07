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

#include "G4CascadePiMinusPChannel.hh"
#include "G4CascadePiPlusNChannel.hh"
#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;

namespace {
  // pi- p : Outgoing particle types of a given multiplicity
  static const G4int pimP2bfs[5][2] =
  {{pim,pro}, {pi0,neu}, {k0,lam}, {k0,s0}, {kpl,sm}};

  static const G4int pimP3bfs[13][3] =
  {{pim,pro,pi0}, {pim,neu,pip}, {pi0,neu,pi0}, {pi0,k0,lam}, 
   {pim,kpl,lam}, {pip,k0,sm},   {pi0,kpl,sm},  {pim,k0,sp},
   {pim,kpl,s0},  {pi0,k0,s0},   {kmi,pro,k0},  {kmi,neu,kpl},
   {k0b,neu,k0}};

  static const G4int pimP4bfs[22][4] =
  {{pim,pro,pip,pim}, {pim,pro,pi0,pi0}, {pi0,neu,pip,pim},
   {pi0,neu,pi0,pi0}, {pim,pip,k0,lam},  {pi0,pi0,k0,lam},
   {pim,pi0,kpl,lam}, {pim,pip,k0,s0},   {pi0,pi0,k0,s0},
   {pim,pi0,kpl,s0},  {pim,pim,kpl,sp},  {pim,pi0,k0,sp},
   {pim,pip,kpl,sm},  {pi0,pi0,kpl,sm},  {pip,pi0,k0,sm},
   {pim,pro,kpl,kmi}, {pim,pro,k0,k0b},  {pi0,pro,kmi,k0},
   {pip,neu,kmi,k0},  {pi0,neu,k0,k0b},  {pi0,neu,kpl,kmi},
   {pim,neu,kpl,k0b}};

  static const G4int pimP5bfs[31][5] =
  {{pim,pro,pip,pim,pi0}, {pim,pro,pi0,pi0,pi0}, {pim,neu,pip,pip,pim},
   {pi0,neu,pip,pim,pi0}, {pi0,neu,pi0,pi0,pi0}, {pim,pip,pi0,k0,lam},
   {pim,pi0,pi0,kpl,lam}, {pim,pip,pim,kpl,lam}, {pi0,pi0,pi0,k0,lam},
   {pim,pip,pim,kpl,s0},  {pim,pi0,pi0,kpl,s0},  {pim,pip,pi0,k0,s0},
   {pi0,pi0,pi0,k0,s0},   {pim,pip,pim,k0,sp},   {pim,pi0,pi0,k0,sp},
   {pim,pim,pi0,kpl,sp},  {pim,pip,pip,k0,sm},   {pip,pi0,pi0,k0,sm},
   {pim,pip,pi0,kpl,sm},  {pi0,pi0,pi0,kpl,sm},  {pim,pro,pi0,kpl,kmi},
   {pim,pro,pi0,k0,k0b},  {pim,pro,pip,kmi,k0},  {pi0,pro,pi0,kmi,k0},
   {pim,pro,pim,kpl,k0b}, {pim,neu,pip,kpl,kmi}, {pim,neu,pip,k0,k0b},
   {pi0,neu,pip,kmi,k0},  {pi0,neu,pim,kpl,k0b}, {pi0,neu,pi0,k0,k0b},
   {pi0,neu,pi0,kpl,kmi}};

  static const G4int pimP6bfs[39][6] =
  {{pim,pro,pip,pip,pim,pim}, {pim,pro,pip,pim,pi0,pi0},
   {pim,pro,pi0,pi0,pi0,pi0}, {pi0,neu,pip,pip,pim,pim},
   {pi0,neu,pip,pim,pi0,pi0}, {pi0,neu,pi0,pi0,pi0,pi0},
   {pim,pip,pim,pi0,kpl,lam}, {pi0,pim,pi0,pi0,kpl,lam},
   {pim,pip,pi0,pi0,k0,lam},  {pim,pip,pip,pim,k0,lam}, 
   {pi0,pi0,pi0,pi0,k0,lam},  {pim,pip,pim,pi0,kpl,s0},
   {pim,pi0,pi0,pi0,kpl,s0},  {pim,pip,pi0,pi0,k0,s0},
   {pim,pip,pip,pim,k0,s0},   {pi0,pi0,pi0,pi0,k0,s0},
   {pim,pim,pi0,pi0,kpl,sp},  {pim,pip,pim,pim,kpl,sp},
   {pim,pi0,pi0,pi0,k0,sp},   {pim,pip,pim,pi0,k0,sp},
   {pim,pip,pi0,pi0,kpl,sm},  {pim,pip,pip,pim,kpl,sm},
   {pim,pip,pip,pi0,k0,sm},   {pip,pi0,pi0,pi0,k0,sm},

   {pim,pro,pip,pi0,kmi,k0},  {pim,pro,pim,pi0,kpl,k0b},
   {pim,pro,pi0,pi0,kpl,kmi}, {pim,pro,pip,pim,kpl,kmi},
   {pim,pro,pi0,pi0,k0,k0b},  {pim,pro,pip,pim,k0,k0b},
   {pi0,pro,pi0,pi0,kmi,k0},  {pi0,neu,pip,pim,kpl,kmi},
   {pi0,neu,pip,pim,k0,k0b},  {pi0,neu,pim,pi0,kpl,k0b},
   {pim,neu,pip,pim,kpl,k0b}, {pi0,neu,pip,pi0,kmi,k0},
   {pim,neu,pip,pip,kmi,k0},  {pi0,neu,pi0,pi0,kpl,kmi},
   {pi0,neu,pi0,pi0,k0,k0b}};

  static const G4int pimP7bfs[46][7] =
  {{pim,pro,pip,pip,pim,pim,pi0}, {pim,pro,pip,pim,pi0,pi0,pi0},
   {pim,pro,pi0,pi0,pi0,pi0,pi0}, {pim,neu,pip,pip,pip,pim,pim},
   {pi0,neu,pip,pip,pim,pim,pi0}, {pi0,neu,pip,pim,pi0,pi0,pi0},
   {pi0,neu,pi0,pi0,pi0,pi0,pi0}, {pim,pip,pim,pi0,pi0,kpl,lam},
   {pim,pip,pip,pim,pim,kpl,lam}, {pim,pi0,pi0,pi0,pi0,kpl,lam},
   {pim,pip,pip,pim,pi0,k0,lam},  {pim,pip,pi0,pi0,pi0,k0,lam}, 
   {pim,pip,pim,pi0,pi0,kpl,s0},  {pim,pip,pip,pim,pim,kpl,s0},
   {pim,pi0,pi0,pi0,pi0,kpl,s0},  {pim,pip,pip,pim,pi0,k0,s0},
   {pim,pip,pi0,pi0,pi0,k0,s0},   {pim,pip,pim,pi0,pi0,k0,sp},
   {pim,pip,pim,pim,pi0,kpl,sp},  {pim,pim,pi0,pi0,pi0,kpl,sp},
   {pim,pip,pip,pim,pim,k0,sp},   {pim,pi0,pi0,pi0,pi0,k0,sp},
   {pim,pip,pip,pim,pi0,kpl,sm},  {pim,pip,pi0,pi0,pi0,kpl,sm},
   {pim,pip,pip,pi0,pi0,k0,sm},   {pim,pip,pip,pip,pim,k0,sm}, 
   {pip,pi0,pi0,pi0,pi0,k0,sm},   {pim,pro,pip,pim,pi0,kpl,kmi},
   {pim,pro,pi0,pi0,pi0,kpl,kmi}, {pim,pro,pim,pi0,pi0,kpl,k0b},
   {pim,pro,pip,pim,pim,kpl,k0b}, {pim,pro,pip,pim,pi0,k0,k0b},
   {pim,pro,pi0,pi0,pi0,k0,k0b},  {pim,pro,pip,pi0,pi0,kmi,k0},
   {pim,pro,pip,pip,pim,kmi,k0},  {pi0,pro,pi0,pi0,pi0,kmi,k0},
   {pi0,neu,pip,pim,pim,kpl,k0b}, {pi0,neu,pim,pi0,pi0,kpl,k0b},
   {pi0,neu,pi0,pi0,pi0,kpl,kmi}, {pi0,neu,pip,pim,pi0,kpl,kmi},
   {pim,neu,pip,pip,pim,kpl,kmi}, {pi0,neu,pip,pim,pi0,k0,k0b},
   {pim,neu,pip,pip,pim,k0,k0b},  {pi0,neu,pi0,pi0,pi0,k0,k0b},
   {pi0,neu,pip,pip,pim,kmi,k0},  {pi0,neu,pip,pi0,pi0,kmi,k0}};

  static const G4int pimP8bfs[51][8] =
  {{pim,pro,pip,pip,pip,pim,pim,pim}, {pim,pro,pip,pip,pim,pim,pi0,pi0},
   {pim,pro,pip,pim,pi0,pi0,pi0,pi0}, {pim,pro,pi0,pi0,pi0,pi0,pi0,pi0},
   {pi0,neu,pi0,pi0,pi0,pi0,pi0,pi0}, {pi0,neu,pip,pim,pi0,pi0,pi0,pi0},
   {pi0,neu,pip,pip,pim,pim,pi0,pi0}, {pi0,neu,pip,pip,pip,pim,pim,pim},
   {pim,pip,pip,pim,pim,pi0,kpl,lam}, {pim,pip,pim,pi0,pi0,pi0,kpl,lam},
   {pim,pi0,pi0,pi0,pi0,pi0,kpl,lam}, {pim,pip,pip,pim,pi0,pi0,k0,lam},
   {pim,pip,pi0,pi0,pi0,pi0,k0,lam},  {pim,pip,pip,pip,pim,pim,k0,lam},
   {pim,pip,pip,pim,pim,pi0,kpl,s0},  {pim,pip,pim,pi0,pi0,pi0,kpl,s0},
   {pim,pip,pip,pim,pi0,pi0,k0,s0},   {pim,pip,pi0,pi0,pi0,pi0,k0,s0}, 
   {pim,pip,pip,pip,pim,pim,k0,s0},   {pim,pip,pim,pim,pi0,pi0,kpl,sp},
   {pim,pim,pi0,pi0,pi0,pi0,kpl,sp},  {pim,pip,pip,pim,pim,pim,kpl,sp},
   {pim,pim,pip,pi0,pi0,pi0,k0,sp},   {pim,pip,pip,pim,pim,pi0,k0,sp},
   {pim,pi0,pi0,pi0,pi0,pi0,k0,sp},   {pim,pip,pip,pim,pi0,pi0,kpl,sm},
   {pim,pip,pi0,pi0,pi0,pi0,kpl,sm},  {pim,pip,pip,pip,pim,pim,kpl,sm},
   {pim,pip,pip,pi0,pi0,pi0,k0,sm},   {pim,pip,pip,pip,pim,pi0,k0,sm},
   {pip,pi0,pi0,pi0,pi0,pi0,k0,sm},   {pim,pro,pip,pim,pi0,pi0,kpl,kmi},
   {pim,pro,pip,pip,pim,pim,kpl,kmi}, {pim,pro,pi0,pi0,pi0,pi0,kpl,kmi},
   {pim,pro,pim,pi0,pi0,pi0,kpl,k0b}, {pim,pro,pip,pim,pim,pi0,kpl,k0b},
   {pim,pro,pi0,pi0,pi0,pi0,k0,k0b},  {pim,pro,pip,pim,pi0,pi0,k0,k0b},
   {pim,pro,pip,pip,pim,pim,k0,k0b},  {pim,pro,pip,pip,pim,pi0,kmi,k0},
   {pim,pro,pip,pi0,pi0,pi0,kmi,k0},  {pi0,neu,pip,pim,pim,pi0,kpl,k0b},
   {pim,neu,pip,pip,pim,pim,kpl,k0b}, {pi0,neu,pim,pi0,pi0,pi0,kpl,k0b},
   {pi0,neu,pip,pip,pim,pim,kpl,kmi}, {pi0,neu,pip,pim,pi0,pi0,kpl,kmi},
   {pi0,neu,pip,pip,pim,pim,k0,k0b},  {pi0,neu,pip,pim,pi0,pi0,k0,k0b},
   {pi0,neu,pip,pip,pim,pi0,kmi,k0},  {pim,neu,pip,pip,pip,pim,kmi,k0},
   {pi0,neu,pip,pi0,pi0,pi0,kmi,k0}};

  static const G4int pimP9bfs[58][9] =
  {{pim,pro,pip,pip,pip,pim,pim,pim,pi0}, {pim,pro,pip,pip,pim,pim,pi0,pi0,pi0},
   {pim,pro,pip,pim,pi0,pi0,pi0,pi0,pi0}, {pim,pro,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
   {pi0,neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {pi0,neu,pip,pim,pi0,pi0,pi0,pi0,pi0},
   {pi0,neu,pip,pip,pim,pim,pi0,pi0,pi0}, {pi0,neu,pip,pip,pip,pim,pim,pim,pi0},
   {pim,neu,pip,pip,pip,pip,pim,pim,pim}, {pim,pip,pip,pim,pim,pi0,pi0,kpl,lam},
   {pim,pip,pim,pi0,pi0,pi0,pi0,kpl,lam}, {pim,pip,pip,pip,pim,pim,pim,kpl,lam},
   {pim,pi0,pi0,pi0,pi0,pi0,pi0,kpl,lam}, {pim,pip,pip,pim,pi0,pi0,pi0,k0,lam},
   {pim,pip,pip,pip,pim,pim,pi0,k0,lam},  {pim,pip,pi0,pi0,pi0,pi0,pi0,k0,lam},
   {pim,pip,pip,pim,pim,pi0,pi0,kpl,s0},  {pim,pip,pim,pi0,pi0,pi0,pi0,kpl,s0},
   {pim,pip,pip,pip,pim,pim,pim,kpl,s0},  {pim,pip,pip,pim,pi0,pi0,pi0,k0,s0},
   {pim,pip,pip,pip,pim,pim,pi0,k0,s0},   {pim,pip,pi0,pi0,pi0,pi0,pi0,k0,s0},
   {pim,pip,pim,pim,pi0,pi0,pi0,kpl,sp},  {pim,pip,pip,pim,pim,pim,pi0,kpl,sp},
   {pim,pim,pi0,pi0,pi0,pi0,pi0,kpl,sp},  {pim,pip,pip,pim,pim,pi0,pi0,k0,sp},
   {pim,pip,pim,pi0,pi0,pi0,pi0,k0,sp},   {pim,pip,pip,pip,pim,pim,pim,k0,sp},
   {pim,pip,pip,pim,pi0,pi0,pi0,kpl,sm},  {pim,pip,pip,pip,pim,pim,pi0,kpl,sm},
   {pim,pip,pi0,pi0,pi0,pi0,pi0,kpl,sm},  {pim,pip,pip,pip,pim,pi0,pi0,k0,sm},
   {pim,pip,pip,pi0,pi0,pi0,pi0,k0,sm},   {pim,pip,pip,pip,pip,pim,pim,k0,sm},
   {pim,pro,pip,pim,pim,pi0,pi0,kpl,k0b}, {pim,pro,pip,pip,pim,pim,pim,kpl,k0b},
   {pim,pro,pim,pi0,pi0,pi0,pi0,kpl,k0b}, {pim,pro,pip,pim,pi0,pi0,pi0,kpl,kmi},
   {pim,pro,pip,pip,pim,pim,pi0,kpl,kmi}, {pim,pro,pi0,pi0,pi0,pi0,pi0,kpl,kmi},
   {pim,pro,pip,pim,pi0,pi0,pi0,k0,k0b},  {pim,pro,pip,pip,pim,pim,pi0,k0,k0b},
   {pim,pro,pi0,pi0,pi0,pi0,pi0,k0,k0b},  {pim,pro,pip,pip,pim,pi0,pi0,kmi,k0},
   {pim,pro,pip,pi0,pi0,pi0,pi0,kmi,k0},  {pim,pro,pip,pip,pip,pim,pim,kmi,k0},
   {pi0,neu,pip,pim,pim,pi0,pi0,kpl,k0b}, {pi0,neu,pip,pip,pim,pim,pim,kpl,k0b},
   {pi0,neu,pim,pi0,pi0,pi0,pi0,kpl,k0b}, {pi0,neu,pip,pip,pim,pim,pi0,kpl,kmi},
   {pi0,neu,pip,pim,pi0,pi0,pi0,kpl,kmi}, {pim,neu,pip,pip,pip,pim,pim,kpl,kmi},
   {pi0,neu,pip,pip,pim,pim,pi0,k0,k0b},  {pi0,neu,pip,pim,pi0,pi0,pi0,k0,k0b},
   {pim,neu,pip,pip,pip,pim,pim,k0,k0b},  {pi0,neu,pip,pip,pim,pi0,pi0,kmi,k0},
   {pi0,neu,pip,pip,pip,pim,pim,kmi,k0},  {pi0,neu,pip,pi0,pi0,pi0,pi0,kmi,k0}};
}

namespace {
  // pi+ n : Outgoing particle types of a given multiplicity
  static const G4int pipN2bfs[5][2] =
  {{pip,neu}, {pi0,pro}, {kpl,lam}, {kpl,s0}, {k0,sp}};

  static const G4int pipN3bfs[13][3] =
  {{pip,neu,pi0}, {pip,pro,pim}, {pi0,pro,pi0}, {pi0,kpl,lam},
   {pip,k0,lam},  {pim,kpl,sp},  {pi0,k0,sp},   {pip,kpl,sm},
   {pip,k0,s0},   {pi0,kpl,s0},  {k0b,neu,kpl}, {k0b,pro,k0},
   {kpl,pro,kmi}};

  static const G4int pipN4bfs[22][4] =
  {{pip,neu,pip,pim}, {pip,neu,pi0,pi0}, {pi0,pro,pip,pim},
   {pi0,pro,pi0,pi0}, {pip,pim,kpl,lam}, {pi0,pi0,kpl,lam},
   {pip,pi0,k0,lam},  {pip,pim,kpl,s0},  {pi0,pi0,kpl,s0},
   {pip,pi0,k0,s0},   {pip,pip,k0,sm},   {pip,pi0,kpl,sm},
   {pip,pim,k0,sp},   {pi0,pi0,k0,sp},   {pim,pi0,kpl,sp},
   {pip,neu,k0,k0b},  {pip,neu,kpl,kmi}, {pi0,neu,kpl,k0b},
   {pim,pro,kpl,k0b}, {pi0,pro,kpl,kmi}, {pi0,pro,k0,k0b},
   {pip,pro,kmi,k0}};

  static const G4int pipN5bfs[31][5] =
  {{pip,neu,pip,pim,pi0}, {pip,neu,pi0,pi0,pi0}, {pip,pro,pip,pim,pim},
   {pi0,pro,pip,pim,pi0}, {pi0,pro,pi0,pi0,pi0}, {pip,pim,pi0,kpl,lam},
   {pip,pi0,pi0,k0,lam},  {pip,pip,pim,k0,lam},  {pi0,pi0,pi0,kpl,lam},
   {pip,pip,pim,k0,s0},   {pip,pi0,pi0,k0,s0},   {pip,pim,pi0,kpl,s0},
   {pi0,pi0,pi0,kpl,s0},  {pip,pip,pim,kpl,sm},  {pip,pi0,pi0,kpl,sm},
   {pip,pip,pi0,k0,sm},   {pip,pim,pim,kpl,sp},  {pim,pi0,pi0,kpl,sp},
   {pip,pim,pi0,k0,sp},   {pi0,pi0,pi0,k0,sp},   {pi0,neu,pip,k0,k0b},
   {pip,neu,pi0,kpl,kmi}, {pip,neu,pim,kpl,k0b}, {pi0,neu,pi0,kpl,k0b},
   {pip,neu,pip,kmi,k0},  {pip,pro,pim,k0,k0b},  {pip,pro,pim,kpl,kmi},
   {pi0,pro,pim,kpl,k0b}, {pi0,pro,pip,kmi,k0},  {pi0,pro,pi0,kpl,kmi},
   {pi0,pro,pi0,k0,k0b}};

  static const G4int pipN6bfs[39][6] =
  {{pip,neu,pip,pip,pim,pim}, {pip,neu,pip,pim,pi0,pi0},
   {pi0,neu,pip,pi0,pi0,pi0}, {pi0,pro,pip,pip,pim,pim},
   {pip,pro,pim,pi0,pi0,pi0}, {pi0,pro,pi0,pi0,pi0,pi0},
   {pip,pip,pim,pi0,k0,lam},  {pip,pi0,pi0,pi0,k0,lam},
   {pip,pim,pi0,pi0,kpl,lam}, {pip,pip,pim,pim,kpl,lam},
   {pi0,pi0,pi0,pi0,kpl,lam}, {pip,pip,pim,pi0,k0,s0},
   {pip,pi0,pi0,pi0,k0,s0},   {pip,pim,pi0,pi0,kpl,s0},
   {pip,pip,pim,pim,kpl,s0},  {pi0,pi0,pi0,pi0,kpl,s0},
   {pip,pip,pi0,pi0,k0,sm},   {pip,pip,pip,pim,k0,sm},
   {pip,pi0,pi0,pi0,kpl,sm},  {pip,pip,pim,pi0,kpl,sm},
   {pip,pim,pi0,pi0,k0,sp},   {pip,pip,pim,pim,k0,sp},
   {pip,pim,pim,pi0,kpl,sp},  {pim,pi0,pi0,pi0,kpl,sp},
   {pip,neu,pim,pi0,kpl,k0b}, {pip,neu,pip,pi0,kmi,k0},
   {pip,neu,pi0,pi0,k0,k0b},  {pip,neu,pip,pim,k0,k0b},
   {pip,neu,pi0,pi0,kpl,kmi}, {pip,neu,pip,pim,kpl,kmi},
   {pi0,neu,pi0,pi0,kpl,k0b}, {pi0,pro,pip,pim,k0,k0b},
   {pi0,pro,pip,pim,kpl,kmi}, {pi0,pro,pip,pi0,kmi,k0},
   {pip,pro,pip,pim,kmi,k0},  {pi0,pro,pim,pi0,kpl,k0b},
   {pip,pro,pim,pim,kpl,k0b}, {pi0,pro,pi0,pi0,k0,k0b},
   {pi0,pro,pi0,pi0,kpl,kmi}};

  static const G4int pipN7bfs[46][7] =
  {{pip,neu,pip,pip,pim,pim,pi0}, {pip,neu,pip,pim,pi0,pi0,pi0},
   {pip,neu,pi0,pi0,pi0,pi0,pi0}, {pip,pro,pip,pip,pim,pim,pim},
   {pi0,pro,pip,pip,pim,pim,pi0}, {pi0,pro,pip,pim,pi0,pi0,pi0},
   {pi0,pro,pi0,pi0,pi0,pi0,pi0}, {pip,pip,pim,pi0,pi0,k0,lam},
   {pip,pip,pip,pim,pim,k0,lam},  {pip,pi0,pi0,pi0,pi0,k0,lam},
   {pip,pip,pim,pim,pi0,kpl,lam}, {pip,pim,pi0,pi0,pi0,kpl,lam},
   {pip,pip,pim,pi0,pi0,k0,s0},   {pip,pip,pip,pim,pim,k0,s0},
   {pip,pi0,pi0,pi0,pi0,k0,s0},   {pip,pip,pim,pim,pi0,kpl,s0},
   {pip,pim,pi0,pi0,pi0,kpl,s0},  {pip,pip,pim,pi0,pi0,kpl,sm},
   {pip,pip,pip,pim,pi0,k0,sm},   {pip,pip,pi0,pi0,pi0,k0,sm},
   {pip,pip,pip,pim,pim,kpl,sm},  {pip,pi0,pi0,pi0,pi0,kpl,sm},
   {pip,pip,pim,pim,pi0,k0,sp},   {pip,pim,pi0,pi0,pi0,k0,sp},
   {pip,pim,pim,pi0,pi0,kpl,sp},  {pip,pip,pim,pim,pim,kpl,sp},
   {pim,pi0,pi0,pi0,pi0,kpl,sp},  {pip,neu,pip,pim,pi0,k0,k0b},
   {pip,neu,pi0,pi0,pi0,k0,k0b},  {pip,neu,pip,pi0,pi0,kmi,k0},
   {pip,neu,pip,pip,pim,kmi,k0},  {pip,neu,pip,pim,pi0,kpl,kmi},
   {pip,neu,pi0,pi0,pi0,kpl,kmi}, {pip,neu,pim,pi0,pi0,kpl,k0b},
   {pip,neu,pip,pim,pim,kpl,k0b}, {pi0,neu,pi0,pi0,pi0,kpl,k0b},
   {pi0,pro,pip,pip,pim,kmi,k0},  {pi0,pro,pip,pi0,pi0,kmi,k0},
   {pi0,pro,pi0,pi0,pi0,k0,k0b},  {pi0,pro,pip,pim,pi0,k0,k0b},
   {pip,pro,pip,pim,pim,k0,k0b},  {pi0,pro,pip,pim,pi0,kpl,kmi},
   {pip,pro,pip,pim,pim,kpl,kmi}, {pi0,pro,pi0,pi0,pi0,kpl,kmi},
   {pi0,pro,pip,pim,pim,kpl,k0b}, {pi0,pro,pim,pi0,pi0,kpl,k0b}};

  static const G4int pipN8bfs[51][8] =
  {{pip,neu,pip,pip,pip,pim,pim,pim}, {pip,neu,pip,pip,pim,pim,pi0,pi0},
   {pip,neu,pip,pim,pi0,pi0,pi0,pi0}, {pip,neu,pi0,pi0,pi0,pi0,pi0,pi0},
   {pi0,pro,pi0,pi0,pi0,pi0,pi0,pi0}, {pi0,pro,pip,pim,pi0,pi0,pi0,pi0},
   {pi0,pro,pip,pip,pim,pim,pi0,pi0}, {pi0,pro,pip,pip,pip,pim,pim,pim},
   {pip,pip,pip,pim,pim,pi0,k0,lam},  {pip,pip,pim,pi0,pi0,pi0,k0,lam},
   {pip,pi0,pi0,pi0,pi0,pi0,k0,lam},  {pip,pip,pim,pim,pi0,pi0,kpl,lam},
   {pip,pim,pi0,pi0,pi0,pi0,kpl,lam}, {pip,pip,pip,pim,pim,pim,kpl,lam},
   {pip,pip,pip,pim,pim,pi0,k0,s0},   {pip,pip,pim,pi0,pi0,pi0,k0,s0},
   {pip,pip,pim,pim,pi0,pi0,kpl,s0},  {pip,pim,pi0,pi0,pi0,pi0,kpl,s0},
   {pip,pip,pip,pim,pim,pim,kpl,s0},  {pip,pip,pip,pim,pi0,pi0,k0,sm},
   {pip,pip,pi0,pi0,pi0,pi0,k0,sm},   {pip,pip,pip,pip,pim,pim,k0,sm},
   {pip,pip,pim,pi0,pi0,pi0,kpl,sm},  {pip,pip,pip,pim,pim,pi0,kpl,sm},
   {pip,pi0,pi0,pi0,pi0,pi0,kpl,sm},  {pip,pip,pim,pim,pi0,pi0,k0,sp},
   {pip,pim,pi0,pi0,pi0,pi0,k0,sp},   {pip,pip,pip,pim,pim,pim,k0,sp},
   {pip,pim,pim,pi0,pi0,pi0,kpl,sp},  {pip,pip,pim,pim,pim,pi0,kpl,sp},
   {pim,pi0,pi0,pi0,pi0,pi0,kpl,sp},  {pip,neu,pip,pim,pi0,pi0,k0,k0b},
   {pip,neu,pip,pip,pim,pim,k0,k0b},  {pip,neu,pi0,pi0,pi0,pi0,k0,k0b},
   {pip,neu,pip,pi0,pi0,pi0,kmi,k0},  {pip,neu,pip,pip,pim,pi0,kmi,k0},
   {pip,neu,pi0,pi0,pi0,pi0,kpl,kmi}, {pip,neu,pip,pim,pi0,pi0,kpl,kmi},
   {pip,neu,pip,pip,pim,pim,kpl,kmi}, {pip,neu,pip,pim,pim,pi0,kpl,k0b},
   {pip,neu,pim,pi0,pi0,pi0,kpl,k0b}, {pi0,pro,pip,pip,pim,pi0,kmi,k0},
   {pip,pro,pip,pip,pim,pim,kmi,k0},  {pi0,pro,pip,pi0,pi0,pi0,kmi,k0},
   {pi0,pro,pip,pip,pim,pim,k0,k0b},  {pi0,pro,pip,pim,pi0,pi0,k0,k0b},
   {pi0,pro,pip,pip,pim,pim,kpl,kmi}, {pi0,pro,pip,pim,pi0,pi0,kpl,kmi},
   {pi0,pro,pip,pim,pim,pi0,kpl,k0b}, {pip,pro,pip,pim,pim,pim,kpl,k0b},
   {pi0,pro,pim,pi0,pi0,pi0,kpl,k0b}};

  static const G4int pipN9bfs[58][9] =
  {{pip,neu,pip,pip,pip,pim,pim,pim,pi0}, {pip,neu,pip,pip,pim,pim,pi0,pi0,pi0},
   {pip,neu,pip,pim,pi0,pi0,pi0,pi0,pi0}, {pip,neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
   {pi0,pro,pi0,pi0,pi0,pi0,pi0,pi0,pi0}, {pi0,pro,pip,pim,pi0,pi0,pi0,pi0,pi0},
   {pi0,pro,pip,pip,pim,pim,pi0,pi0,pi0}, {pi0,pro,pip,pip,pip,pim,pim,pim,pi0},
   {pip,pro,pip,pip,pip,pim,pim,pim,pim}, {pip,pip,pip,pim,pim,pi0,pi0,k0,lam},
   {pip,pip,pim,pi0,pi0,pi0,pi0,k0,lam},  {pip,pip,pip,pip,pim,pim,pim,k0,lam},
   {pip,pi0,pi0,pi0,pi0,pi0,pi0,k0,lam},  {pip,pip,pim,pim,pi0,pi0,pi0,kpl,lam},
   {pip,pip,pip,pim,pim,pim,pi0,kpl,lam}, {pip,pim,pi0,pi0,pi0,pi0,pi0,kpl,lam},
   {pip,pip,pip,pim,pim,pi0,pi0,k0,s0},   {pip,pip,pim,pi0,pi0,pi0,pi0,k0,s0},
   {pip,pip,pip,pip,pim,pim,pim,k0,s0},   {pip,pip,pim,pim,pi0,pi0,pi0,kpl,s0},
   {pip,pip,pip,pim,pim,pim,pi0,kpl,s0},  {pip,pim,pi0,pi0,pi0,pi0,pi0,kpl,s0},
   {pip,pip,pip,pim,pi0,pi0,pi0,k0,sm},   {pip,pip,pip,pip,pim,pim,pi0,k0,sm},
   {pip,pip,pi0,pi0,pi0,pi0,pi0,k0,sm},   {pip,pip,pip,pim,pim,pi0,pi0,kpl,sm},
   {pip,pip,pim,pi0,pi0,pi0,pi0,kpl,sm},  {pip,pip,pip,pip,pim,pim,pim,kpl,sm},
   {pip,pip,pim,pim,pi0,pi0,pi0,k0,sp},   {pip,pip,pip,pim,pim,pim,pi0,k0,sp},
   {pip,pim,pi0,pi0,pi0,pi0,pi0,k0,sp},   {pip,pip,pim,pim,pim,pi0,pi0,kpl,sp},
   {pip,pim,pim,pi0,pi0,pi0,pi0,kpl,sp},  {pip,pip,pip,pim,pim,pim,pim,kpl,sp},
   {pip,neu,pip,pip,pim,pi0,pi0,kmi,k0},  {pip,neu,pip,pip,pip,pim,pim,kmi,k0},
   {pip,neu,pip,pi0,pi0,pi0,pi0,kmi,k0},  {pip,neu,pip,pim,pi0,pi0,pi0,k0,k0b},
   {pip,neu,pip,pip,pim,pim,pi0,k0,k0b},  {pip,neu,pi0,pi0,pi0,pi0,pi0,k0,k0b},
   {pip,neu,pip,pim,pi0,pi0,pi0,kpl,kmi}, {pip,neu,pip,pip,pim,pim,pi0,kpl,kmi},
   {pip,neu,pi0,pi0,pi0,pi0,pi0,kpl,kmi}, {pip,neu,pip,pim,pim,pi0,pi0,kpl,k0b},
   {pip,neu,pim,pi0,pi0,pi0,pi0,kpl,k0b}, {pip,neu,pip,pip,pim,pim,pim,kpl,k0b},
   {pi0,pro,pip,pip,pim,pi0,pi0,kmi,k0},  {pi0,pro,pip,pip,pip,pim,pim,kmi,k0},
   {pi0,pro,pip,pi0,pi0,pi0,pi0,kmi,k0},  {pi0,pro,pip,pip,pim,pim,pi0,k0,k0b},
   {pi0,pro,pip,pim,pi0,pi0,pi0,k0,k0b},  {pip,pro,pip,pip,pim,pim,pim,k0,k0b},
   {pi0,pro,pip,pip,pim,pim,pi0,kpl,kmi}, {pi0,pro,pip,pim,pi0,pi0,pi0,kpl,kmi},
   {pip,pro,pip,pip,pim,pim,pim,kpl,kmi}, {pi0,pro,pip,pim,pim,pi0,pi0,kpl,k0b},
   {pi0,pro,pip,pip,pim,pim,pim,kpl,k0b}, {pi0,pro,pim,pi0,pi0,pi0,pi0,kpl,k0b}};
}

namespace {
  // Total pi- p cross section as a function of kinetic energy
  static const G4double pimPtotXSec[30] = 
   { 6.13,  6.4,   6.67,  6.94,  7.22,  7.5,   8.3,   12.0,   14.4,   24.0,
    46.0,  72.04, 43.02, 27.19, 27.32, 42.72, 37.638, 51.546, 33.687, 34.971,
    32.57, 31.19, 30.24, 28.5,  27.1,  26.2,  25.48,  25.22,  24.9,   24.7};

  // pi- p cross sections as functions of kinetic energy and multiplicity
  const G4double pimPCrossSections[265][30] = {
  //
  // multiplicity 2 (5 channels)
  //
  //  pi- p (pi+ n)
   { 1.43, 1.5,   1.57,   1.64,  1.72,   1.8,   2.0,   3.0,   3.4,   7.0,
    14.0, 24.04, 14.905, 11.27, 11.004, 20.99, 13.97, 25.84, 12.01, 10.696,
     8.3,  7.1,   6.0,    5.4,   4.8,    4.4,   4.12,  3.95,  3.8,   3.6},

  //  n pi0  (p pi0)
   { 4.7,  4.9,  5.1,  5.3,   5.5,   5.7,   6.3,   9.0,  11.0,  17.0,
    32.0, 48.0, 28.0, 14.5,  11.04,  8.99,  4.79,  5.02,  2.08,  1.0,
     0.5, 0.25, 0.15,  0.095, 0.065, 0.046, 0.035, 0.026, 0.019, 0.015},

  //  L K0  (L K+)
   { 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.65,  0.29,  0.17,
     0.15, 0.08, 0.05, 0.032, 0.021, 0.015, 0.011, 0.008, 0.006, 0.005},

  //  S0 K0  (S0 K+)
   { 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,    0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.28,   0.19,  0.12,
     0.11, 0.06, 0.04, 0.026, 0.018, 0.013, 0.009, 0.006,  0.004, 0.003},

  //  S- K+  (S+ K0)
   { 0.0,  0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.20, 0.25, 0.09,
     0.04, 0.015, 0.007, 0.003, 0.001, 0.0,  0.0,  0.0,  0.0,  0.0},
  //
  // multiplicity 3 (13 channels)
  //
  //  p pi- pi0  (n pi+ pi0)
   { 0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,
     0.0,  0.0,  0.015, 0.15, 0.82, 3.4,  5.4,  5.5, 5.1,  5.0,
     3.3,  2.24, 1.65,  1.1,  0.8,  0.59, 0.42, 0.3, 0.22, 0.15},

  //  n pi+ pi-  (p pi+ pi-)
   { 0.0, 0.0, 0.0,   0.0,  0.0, 0.0, 0.0, 0.0,  0.0,  0.0,
     0.0, 0.0, 0.075, 0.95, 3.3, 6.5, 9.0, 9.2,  8.0,  8.0,
     5.0, 3.3, 2.4,   1.62, 1.1, 0.8, 0.6, 0.42, 0.33, 0.24},

  //  n pi0 pi0  (p pi0 pi0)
   { 0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.025, 0.32, 1.1,  2.2,   3.0,   2.3,   1.4,   0.7,
     0.39, 0.21, 0.11,  0.06, 0.03, 0.015, 0.008, 0.004, 0.002, 0.001},

  //  L K0 pi0  (L K+ pi0)
   { 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.043, 0.143,
     0.179, 0.143, 0.104, 0.074, 0.051, 0.035, 0.026, 0.018, 0.013, 0.009},
 
  //  L K+ pi-  (L K0 pi+)
   { 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.038, 0.11,
     0.135, 0.11, 0.08, 0.055, 0.039, 0.028, 0.02, 0.014, 0.01,  0.007},

  //  S- K0 pi+  (S+ K+ pi-)
   { 0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.08,
     0.13, 0.07, 0.036, 0.02, 0.01, 0.005, 0.003, 0.002, 0.001, 0.0},

  //  S- K+ pi0  (S+ K0 pi0)
  {  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0, 0.0, 0.0,
     0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0, 0.0, 0.04,
     0.07, 0.04, 0.02, 0.009, 0.004, 0.002, 0.001, 0.0, 0.0, 0.0},

  //  S+ K0 pi-  (S- K+ pi+)
  { 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0, 0.0, 0.0,
    0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0, 0.0, 0.04,
    0.06, 0.04, 0.02, 0.009, 0.004, 0.002, 0.001, 0.0, 0.0, 0.0},

  //  S0 K+ pi-  (S0 K0 pi+)
   { 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.04,
     0.08, 0.05, 0.03, 0.019, 0.012, 0.008, 0.005, 0.003, 0.002, 0.001},

  //  S0 K0 pi0  (S0 K+ pi0)
   { 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.03,  0.08,
     0.07, 0.05, 0.03, 0.019, 0.012, 0.008, 0.005, 0.003, 0.002, 0.001},

  //  p K0 K-  (n K+ K0bar)
   { 0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
     0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.03,
     0.07, 0.07, 0.055, 0.04, 0.03, 0.022, 0.018, 0.013, 0.01, 0.008},

  //  n K+ K-  (p K0 K0bar)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.03,
     0.12, 0.18, 0.18, 0.12, 0.07, 0.035, 0.02, 0.01, 0.005, 0.003},

  //  n K0 K0bar  (p K+ K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.03,
     0.10, 0.15, 0.16, 0.10, 0.05, 0.025, 0.014, 0.006, 0.003, 0.002},
  //
  // multiplicity 4 (22 channels)
  //
  //  p pi+ pi- pi-  (n pi+ pi+ pi-)
   { 0.0, 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.004, 0.035, 0.13, 0.39, 0.75, 1.4,
     1.8, 1.85, 1.65, 1.45, 1.3,   1.2,   1.05, 0.95, 0.83, 0.75},

  //  p pi- pi0 pi0  (n pi+ pi0 pi0)
   { 0.0, 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
     0.0, 0.0,  0.0,  0.0,  0.004, 0.035, 0.13, 0.39, 0.75, 1.4,
     1.8, 1.85, 1.65, 1.45, 1.3,   1.2,   1.05, 0.95, 0.83, 0.75},

  //  n pi+ pi- pi0  (p pi+ pi- pi0)
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0, 0.008, 0.07, 0.26, 0.78, 1.5,  2.8,
     3.6, 3.6, 3.3, 2.9, 2.6,   2.35, 2.1,  1.85, 1.66, 1.5},

  //  n pi0 pi0 pi0  (p pi0 pi0 pi0)
   { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0, 0.0, 0.0,  0.0,  0.04, 0.5,  0.95, 0.9,  0.8,  0.71,
     0.6, 0.5, 0.43, 0.35, 0.3,  0.25, 0.21, 0.17, 0.14, 0.115},

  //  L K0 pi+ pi-  (L K+ pi+ pi-)
   { 0.0,  0.0, 0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.012,
     0.06, 0.1, 0.09, 0.08, 0.065, 0.055, 0.047, 0.039, 0.035, 0.03},

  //  L K0 pi0 pi0  (L K+ pi0 pi0)
   { 0.0,  0.0, 0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0, 0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.012,
     0.06, 0.1, 0.09, 0.08, 0.065, 0.055, 0.047, 0.039, 0.035, 0.03},

  //  L K+ pi- pi0  (L K0 pi+ pi0)
   { 0.0,  0.0, 0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,  0.0,
     0.0,  0.0, 0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,  0.012,
     0.06, 0.1, 0.09, 0.075, 0.06, 0.05, 0.042, 0.035, 0.03, 0.025},

  //  S0 K0 pi+ pi-  (S0 K+ pi+ pi-)
   { 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0,   0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.009,
     0.035, 0.05, 0.05, 0.043, 0.036, 0.03, 0.025, 0.02, 0.017, 0.014},

  //  S0 K0 pi0 pi0  (S0 K+ pi0 pi0)
  {  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
     0.0,   0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.009,
     0.035, 0.05, 0.05, 0.043, 0.036, 0.03, 0.025, 0.02, 0.017, 0.014},

  //  S0 K+ pi- pi0  (S0 K0 pi+ pi0)
   { 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.007,
     0.03, 0.045, 0.04, 0.035, 0.03, 0.025, 0.022, 0.018, 0.015, 0.012},

  //  S+ K+ pi- pi-  (S- K0 pi+ pi+)
   { 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.007,
     0.03, 0.045, 0.04, 0.035, 0.03, 0.025, 0.022, 0.018, 0.015, 0.012},

  //  S+ K0 pi- pi0  (S- K+ pi+ pi0)
   { 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.007,
     0.03, 0.045, 0.04, 0.035, 0.03, 0.025, 0.022, 0.018, 0.015, 0.012},

  //  S- K+ pi+ pi-  (S+ K0 pi+ pi-)
   { 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.007,
     0.03, 0.045, 0.04, 0.035, 0.03, 0.025, 0.022, 0.018, 0.015, 0.012},

  //  S- K+ pi0 pi0  (S+ K0 pi0 pi0)
   { 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.007,
     0.03, 0.045, 0.04, 0.035, 0.03, 0.025, 0.022, 0.018, 0.015, 0.012},

  //  S- K0 pi+ pi0  (S+ K+ pi- pi0)
   { 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.007,
     0.03, 0.045, 0.04, 0.035, 0.03, 0.025, 0.022, 0.018, 0.015, 0.012},

  //  p pi- K+ K-  (n pi+ K0 K0bar)
   { 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.047, 0.09, 0.09, 0.083, 0.075, 0.068, 0.057, 0.05, 0.045},

  //  p pi- K0 K0bar  (n pi+ K+ K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.02, 0.05, 0.07, 0.07, 0.065, 0.06, 0.055, 0.049, 0.045, 0.04},

  //  p pi0 K0 K-  (n pi0 K+ K0bar)
   { 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.02, 0.05, 0.07, 0.07, 0.065, 0.06, 0.055, 0.049, 0.045, 0.04},

  //  n pi+ K0 K-  (p pi- K+ K0bar)
   { 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.02, 0.05, 0.07, 0.07, 0.065, 0.06, 0.055, 0.049, 0.045, 0.04},

  //  n pi0 K0 K0bar  (p pi0 K+ K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.02, 0.05, 0.07, 0.07, 0.065, 0.06, 0.055, 0.049, 0.045, 0.04},

  //  n pi0 K+ K-  (p pi0 K0 K0bar)
   { 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.02, 0.05, 0.07, 0.07, 0.065, 0.06, 0.055, 0.049, 0.045, 0.04},

  //  n pi- K+ K0bar  (p pi+ K0 K-)
   { 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.02, 0.05, 0.07, 0.07, 0.065, 0.06, 0.055, 0.049, 0.045, 0.04},
  //
  //  multiplicity 5 (31 channels)
  //
  //  p pi+ pi- pi- pi0  (n pi+ pi+ pi- pi0)
   { 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.0,   0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.002, 0.02, 0.10, 0.61,
     1.5, 1.9, 2.1, 1.85, 1.6, 1.35, 1.15,  0.91, 0.74, 0.62},

  //  p pi- pi0 pi0 pi0  (n pi+ pi0 pi0 pi0)
   {  0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  0.0,   0.0,  0.0,  0.0,
      0.0, 0.0,  0.0,  0.0,  0.0, 0.0,  0.001, 0.01, 0.05, 0.26,
      0.7, 0.98, 1.05, 0.93, 0.8, 0.67, 0.56,  0.45, 0.37, 0.3},

  //  n pi+ pi+ pi- pi-  (p pi+ pi+ pi- pi-)
   {  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0,   0.0,  0.0,  0.0,
      0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.003, 0.03, 0.12, 0.34,
      0.67, 0.93, 1.15, 1.2, 1.1, 0.94, 0.76,  0.59, 0.47, 0.37},

  //  n pi+ pi- pi0 pi0  (p pi+ pi- pi0 pi0)
   { 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.0,   0.0,  0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,  0.0, 0.0,  0.002, 0.02, 0.10, 0.52,
     1.4, 1.9, 2.1, 1.85, 1.6, 1.35, 1.15,  0.91, 0.74, 0.62},

  //  n pi0 pi0 pi0 pi0  (p pi0 pi0 pi0 pi0)
   {  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.003, 0.014, 0.073,
      0.198, 0.269, 0.292, 0.262, 0.227, 0.189, 0.162, 0.128, 0.102, 0.083},

  //  L K0 pi+ pi- pi0  (L K+ pi+ pi- pi0)
   {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
      0.01, 0.04, 0.07, 0.08, 0.08, 0.072, 0.065, 0.06, 0.055, 0.05},

  //  L K+ pi- pi0 pi0  (L K0 pi+ pi0 pi0)
   {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
      0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
      0.0, 0.01, 0.03, 0.04, 0.04, 0.035, 0.03, 0.024, 0.02, 0.017},

  //  L K+ pi+ pi- pi-  (L K0 pi+ pi+ pi-)
   {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
      0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
      0.0, 0.01, 0.03, 0.04, 0.04, 0.035, 0.03, 0.024, 0.02, 0.017},

  //  L K0 pi0 pi0 pi0  (L K+ pi0 pi0 pi0)
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.006, 0.019, 0.025, 0.022, 0.02, 0.018, 0.015, 0.013, 0.011},

  //  S0 K+ pi+ pi- pi-  (S0 K0 pi+ pi+ pi-)
   {  0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.005, 0.015, 0.02, 0.02, 0.017, 0.015, 0.012, 0.01, 0.008},

  //  S0 K+ pi- pi0 pi0  (S0 K0 pi+ pi0 pi0)
   {  0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.005, 0.015, 0.02, 0.02, 0.017, 0.015, 0.012, 0.01, 0.008},

  //  S0 K0 pi+ pi- pi0  (S0 K+ pi+ pi- pi0)
   {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
      0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
      0.0, 0.01, 0.03, 0.04, 0.04, 0.035, 0.03, 0.024, 0.02, 0.017},

  //  S0 K0 pi0 pi0 pi0  (S0 K+ pi0 pi0 pi0)
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.004, 0.006, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002},

  //  S+ K0 pi+ pi- pi-  (S- K+ pi+ pi+ pi-)
   {  0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.004, 0.01, 0.012, 0.011, 0.01, 0.009, 0.008, 0.007},

  //  S+ K0 pi- pi0 pi0  (S- K+ pi+ pi0 pi0)
   {  0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.004, 0.01, 0.012, 0.011, 0.01, 0.009, 0.008, 0.007},

  //  S+ K+ pi- pi- pi0  (S- K0 pi+ pi+ pi0)
   {  0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.004, 0.01, 0.012, 0.011, 0.01, 0.009, 0.008, 0.007},

  //  S- K0 pi+ pi+ pi-  (S+ K+ pi+ pi- pi-)
   {  0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.003, 0.01, 0.013, 0.012, 0.011, 0.01, 0.009, 0.008, 0.007},

  //  S- K0 pi+ pi0 pi0  (S+ K+ pi- pi0 pi0)
   {  0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.003, 0.01, 0.013, 0.012, 0.011, 0.01, 0.009, 0.008, 0.007},

  //  S- K+ pi+ pi- pi0  (S+ K0 pi+ pi- pi0)
   {  0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.006, 0.02, 0.027, 0.025, 0.023, 0.021, 0.019, 0.017, 0.015},

  //  S- K+ pi0 pi0 pi0  (S+ K0 pi0 pi0 pi0)
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.005, 0.005, 0.004, 0.003, 0.002, 0.001, 0.0},

  //  p pi- pi0 K+ K-  (n pi+ pi0 K0 K0bar)
   {  0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
      0.0, 0.008, 0.04, 0.065, 0.065, 0.06, 0.055, 0.05, 0.045, 0.04},

  //  p pi- pi0 K0 K0bar  (n pi+ pi0 K+ K-)
   {  0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
      0.0, 0.008, 0.04, 0.065, 0.065, 0.06, 0.055, 0.05, 0.045, 0.04},

  //  p pi+ pi- K0 K-  (n pi+ pi- K+ K0bar)
   {  0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,
      0.0, 0.008, 0.04, 0.065, 0.065, 0.06, 0.055, 0.05, 0.045, 0.04},

  //  p pi0 pi0 K0 K-  (n pi0 pi0 K+ K0bar)
   {  0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.004, 0.02, 0.032, 0.032, 0.03, 0.027, 0.025, 0.022, 0.02},

  //  p pi- pi- K+ K0bar  (n pi+ pi+ K0 K-)
   {  0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.004, 0.02, 0.032, 0.032, 0.03, 0.027, 0.025, 0.022, 0.02},

  //  n pi+ pi- K+ K-  (p pi+ pi- K0 K0bar)
   {  0.0, 0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
      0.0, 0.008, 0.04, 0.05, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025},

  //  n pi+ pi- K0 K0bar  (p pi+ pi- K+ K-)
   {  0.0, 0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
      0.0, 0.008, 0.04, 0.05, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025},

  //  n pi+ pi0 K0 K-  (p pi- pi0 K+ K0bar)
   {  0.0, 0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
      0.0, 0.008, 0.04, 0.05, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025},

  //  n pi- pi0 K+ K0bar  (p pi+ pi0 K0 K-)
   {  0.0, 0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,  0.0,
      0.0, 0.008, 0.04, 0.05, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025},

  //  n pi0 pi0 K0 K0bar  (p pi0 pi0 K+ K-)
   {  0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.004, 0.02, 0.025, 0.025, 0.022, 0.02, 0.017, 0.015, 0.012},

  //  n pi0 pi0 K+ K-  (p pi0 pi0 K0 K0bar)
   {  0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.004, 0.02, 0.025, 0.025, 0.022, 0.02, 0.017, 0.015, 0.012},
  //
  //  multiplicity 6 (39 channels)
  //
  //  p pi+ pi+ pi- pi- pi-  (n pi+ pi+ pi+ pi- pi-)
   {  0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.001, 0.007, 0.033,
      0.088, 0.167, 0.22, 0.27, 0.33, 0.37, 0.38, 0.39,  0.4,   0.4},

  //  p pi+ pi- pi- pi0 pi0  (n pi+ pi+ pi- pi0 pi0)
   {  0.0, 0.0,   0.0, 0.0, 0.0,  0.0,  0.0,  0.0,   0.0,  0.0,
      0.0, 0.0,   0.0, 0.0, 0.0,  0.0,  0.0,  0.004, 0.02, 0.095,
      0.2, 0.321, 0.4, 0.5, 0.55, 0.57, 0.58, 0.59,  0.6,  0.6},

  //  p pi- pi0 pi0 pi0 pi0  (n pi+ pi0 pi0 pi0 pi0)
   {  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.001, 0.004, 0.017,
      0.045, 0.083, 0.12, 0.145, 0.165, 0.185, 0.192, 0.195, 0.2,   0.2},

  //  n pi+ pi+ pi- pi- pi0  (p pi+ pi+ pi- pi- pi0)
   {  0.0, 0.0,   0.0, 0.0, 0.0,  0.0,  0.0,  0.0,   0.0,  0.0,
      0.0, 0.0,   0.0, 0.0, 0.0,  0.0,  0.0,  0.004, 0.02, 0.095,
      0.2, 0.321, 0.4, 0.5, 0.55, 0.57, 0.58, 0.59,  0.6,  0.6},

  //  n pi+ pi- pi0 pi0 pi0  (p pi+ pi- pi0 pi0 pi0)
   {  0.0,   0.0,   0.0,  0.0, 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,  0.0, 0.0,  0.0,  0.0,  0.003, 0.013, 0.067,
      0.167, 0.333, 0.43, 0.5, 0.57, 0.62, 0.63, 0.63,  0.63,  0.63},

  //  n pi0 pi0 pi0 pi0 pi0  (p pi0 pi0 pi0 pi0 pi0)
   {  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.001, 0.004,
      0.009, 0.017, 0.024, 0.03, 0.036, 0.04, 0.04, 0.04, 0.04,  0.04},

  //  L K+ pi+ pi- pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.008, 0.025, 0.033, 0.033, 0.03, 0.027, 0.024, 0.021, 0.019},

  //  L K+ pi- pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.002, 0.007, 0.011, 0.011, 0.01, 0.009, 0.008, 0.007, 0.006},

  //  L K0 pi+ pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.008, 0.025, 0.033, 0.033, 0.03, 0.027, 0.024, 0.021, 0.019},

  //  L K0 pi+ pi+ pi- pi-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.004, 0.012, 0.017, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011},

  //  L K0 pi0 pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003, 0.002, 0.002, 0.002},

  //  S0 K+ pi+ pi- pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.004, 0.012, 0.017, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011},

  //  S0 K+ pi- pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.005, 0.006, 0.006, 0.006, 0.005, 0.004, 0.003},

  //  S0 K0 pi+ pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.004, 0.012, 0.017, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011},

  //  S0 K0 pi+ pi+ pi- pi-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.002, 0.006, 0.009, 0.009, 0.009, 0.008, 0.007, 0.006, 0.006},

  //  S0 K0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001},

  //  S+ K+ pi- pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.002, 0.006, 0.008, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003},

  //  S+ K+ pi+ pi- pi- pi-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.005, 0.005, 0.005, 0.005, 0.004, 0.004, 0.004},

  //  S+ K0 pi- pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.005, 0.005, 0.005, 0.005, 0.004, 0.004, 0.004},

  //  S+ K0 pi+ pi- pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.004, 0.012, 0.016, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011},

  //  S- K+ pi+ pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.004, 0.012, 0.016, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011},

  //  S- K+ pi+ pi+ pi- pi-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.002, 0.006, 0.008, 0.008, 0.008, 0.007, 0.006, 0.005, 0.004},

  //  S- K+ pi0 pi0 pi0 pi0 negligible

  //  S- K0 pi+ pi+ pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.004, 0.012, 0.016, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011},

  //  S- K0 pi+ pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.004, 0.006, 0.006, 0.006, 0.005, 0.004, 0.004, 0.003},

  //  p pi+ pi- pi0 K0 K-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.005, 0.026, 0.046, 0.054, 0.054, 0.05, 0.044, 0.041, 0.036},

  //  p pi- pi- pi0 K+ K0b
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  //  p pi- pi0 pi0 K+ K-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  //  p pi+ pi- pi- K+ K-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  //  p pi- pi0 pi0 K0 K0b
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  //  p pi+ pi- pi- K0 K0b
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  //  p pi0 pi0 pi0 K0 K-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.005, 0.008, 0.008, 0.008, 0.007, 0.006, 0.005, 0.004},

  //  n pi+ pi- pi0 K+ K-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.005, 0.026, 0.046, 0.054, 0.054, 0.05, 0.044, 0.041, 0.036},

  //  n pi+ pi- pi0 K0 K0b
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.005, 0.026, 0.046, 0.054, 0.054, 0.05, 0.044, 0.041, 0.036},

  //  n pi- pi0 pi0 K+ K0b
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  //  n pi+ pi- pi- K+ K0b
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  //  n pi+ pi0 pi0 K0 K-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  //  n pi+ pi+ pi- K0 K-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.003, 0.013, 0.023, 0.027, 0.027, 0.025, 0.022, 0.02, 0.018},

  //  n pi0 pi0 pi0 K+ K-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.005, 0.008, 0.008, 0.008, 0.007, 0.006, 0.005, 0.004},

  //  n pi0 pi0 pi0 K0 K0b
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.005, 0.008, 0.008, 0.008, 0.007, 0.006, 0.005, 0.004},
  //
  //  multiplicity 7 (46 channels)
  //
  //  p pi+ pi+ pi- pi- pi- pi0
   {  0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.002, 0.012,
      0.045, 0.12, 0.213, 0.34, 0.45, 0.55, 0.641, 0.739, 0.784, 0.816},

  //  p pi+ pi- pi- pi0 pi0 pi0
   {  0.0,  0.0,   0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,  0.0,   0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.001, 0.006,
      0.02, 0.048, 0.1, 0.172, 0.25, 0.315, 0.377, 0.438, 0.451, 0.46},

  //  p pi- pi0 pi0 pi0 pi0 pi0
   {  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.001,
      0.004, 0.011, 0.025, 0.04, 0.057, 0.07, 0.08, 0.091, 0.095, 0.098},

  //  n pi+ pi+ pi+ pi- pi- pi-
   {  0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.001, 0.005,
      0.017, 0.04, 0.085, 0.145, 0.21, 0.265, 0.319, 0.365, 0.38,  0.38},

  //  n pi+ pi+ pi- pi- pi0 pi0 
   {  0.0,   0.0,  0.0,   0.0,  0.0, 0.0,  0.0,  0.0,   0.0,   0.0,
      0.0,   0.0,  0.0,   0.0,  0.0, 0.0,  0.0,  0.0,   0.002, 0.013,
      0.045, 0.12, 0.203, 0.31, 0.4, 0.48, 0.53, 0.572, 0.593, 0.604},

  //  n pi+ pi- pi0 pi0 pi0 pi0
   {  0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.001, 0.005,
      0.017, 0.04, 0.085, 0.145, 0.21, 0.265, 0.319, 0.365, 0.38,  0.38},

  //  n pi0 pi0 pi0 pi0 pi0 pi0
   {  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.001, 0.003, 0.006, 0.009, 0.012, 0.014, 0.016, 0.018, 0.019, 0.02},

  //  L K+ pi+ pi- pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.002, 0.006, 0.012, 0.018, 0.021, 0.021, 0.021, 0.021, 0.021},

  //  L K+ pi+ pi+ pi- pi- pi-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.002, 0.004, 0.006, 0.007, 0.007, 0.007, 0.007, 0.007},

  //  L K+ pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004, 0.004},

  //  L K0 pi+ pi+ pi- pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.002, 0.006, 0.012, 0.018, 0.021, 0.021, 0.021, 0.021, 0.021},

  //  L K0 pi+ pi- pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
      0.0, 0.001, 0.004, 0.008, 0.012, 0.014,  0.014,  0.014, 0.014, 0.014},

  //  L K0 pi0 pi0 pi0 pi0 pi0  negligible

  //  S0 K+ pi+ pi- pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
      0.0, 0.001, 0.003, 0.006, 0.009, 0.01,  0.01,  0.01, 0.01, 0.01},

  //  S0 K+ pi+ pi+ pi- pi- pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004, 0.004},

  //  S0 K+ pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.002,  0.002,  0.002, 0.002, 0.002},

  //  S0 K0 pi+ pi+ pi- pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
      0.0, 0.001, 0.003, 0.006, 0.009, 0.01,  0.01,  0.01, 0.01, 0.01},

  //  S0 K0 pi+ pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.004, 0.006, 0.007,  0.007,  0.007, 0.007, 0.007},

  //  S0 K0 pi0 pi0 pi0 pi0 pi0  negligible

  //  S+ K0 pi+ pi- pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,    0.0,    0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.007, 0.011, 0.012,  0.012,  0.012, 0.012, 0.012},

  //  S+ K+ pi+ pi- pi- pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.002, 0.004, 0.007, 0.008, 0.008, 0.008, 0.008, 0.008},

  //  S+ K+ pi- pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004, 0.004},

  //  S+ K0 pi+ pi+ pi- pi- pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004, 0.004},

  //  S+ K0 pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002},

  //  S- K+ pi+ pi+ pi- pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.007, 0.011, 0.012, 0.012, 0.012, 0.012, 0.012},

  //  S- K+ pi+ pi- pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.002, 0.004, 0.007, 0.008, 0.008, 0.008, 0.008, 0.008},

  //  S- K+ pi0 pi0 pi0 pi0 pi0  negligible

  //  S- K0 pi+ pi+ pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.007, 0.011, 0.012, 0.012, 0.012, 0.012, 0.012},

  //  S- K0 pi+ pi+ pi+ pi- pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004, 0.004},

  //  S- K0 pi+ pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002},

  //  p K+ K- pi+ pi- pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.003, 0.013, 0.023, 0.029, 0.029, 0.027, 0.025, 0.025, 0.025},

  //  p K+ K- pi- pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.004, 0.008, 0.01, 0.01, 0.009, 0.008, 0.007, 0.007},

  //  p K+ K0b pi- pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.002, 0.006, 0.012, 0.015, 0.015, 0.014, 0.013, 0.013, 0.013},

  //  p K+ K0b pi+ pi- pi- pi-
   {  0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.004, 0.008, 0.01, 0.01, 0.009, 0.008, 0.007, 0.007},

  //  p K0 K0b pi+ pi- pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.003, 0.013, 0.023, 0.029, 0.029, 0.027, 0.025, 0.025, 0.025},

  //  p K0 K0b pi- pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.004, 0.008, 0.01, 0.01, 0.009, 0.008, 0.007, 0.007},

  //  p K0 K- pi+ pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.003, 0.013, 0.023, 0.029, 0.029, 0.027, 0.025, 0.025, 0.025},

  //  p K0 K- pi+ pi+ pi- pi-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.002, 0.006, 0.012, 0.015, 0.015, 0.014, 0.013, 0.013, 0.013},

  //  p K0 K- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002},

  //  n K+ K0b pi+ pi- pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.003, 0.013, 0.023, 0.029, 0.029, 0.027, 0.025, 0.025, 0.025},

  //  n K+ K0b pi- pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.004, 0.008, 0.01, 0.01, 0.009, 0.008, 0.008, 0.008},

  //  n K+ K- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002},

  //  n K+ K- pi+ pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.003, 0.013, 0.023, 0.029, 0.029, 0.027, 0.025, 0.025, 0.025},

  //  n K+ K- pi+ pi+ pi- pi-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.002, 0.006, 0.012, 0.015, 0.015, 0.014, 0.014, 0.014, 0.014},

  //  n K0 K0b pi+ pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.003, 0.013, 0.023, 0.029, 0.029, 0.028, 0.027, 0.027, 0.027},

  //  n K0 K0b pi+ pi+ pi- pi-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.002, 0.006, 0.012, 0.015, 0.015, 0.014, 0.014, 0.014, 0.014},

  //  n K0 K0b pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002},

  //  n K0 K- pi+ pi+ pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
      0.0, 0.003, 0.013, 0.023, 0.029, 0.03, 0.03, 0.03, 0.03, 0.03},

  //  n K0 K- pi+ pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
      0.0, 0.001, 0.004, 0.008, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01},
  //
  //  multiplicity 8 (51 channels)
  //
  //  p pi+ pi+ pi+ pi- pi- pi- pi-  (n pi+ pi+ pi+ pi+ pi- pi- pi-)
   {  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,
      0.004, 0.018, 0.039, 0.06, 0.081, 0.11, 0.13, 0.152, 0.172, 0.189},

  //  p pi+ pi+ pi- pi- pi- pi0 pi0  (n pi+ pi+ pi+ pi- pi- pi0 pi0)
   {  0.0,  0.0,   0.0,  0.0,   0.0,  0.0,  0.0, 0.0,   0.0,   0.0,
      0.0,  0.0,   0.0,  0.0,   0.0,  0.0,  0.0, 0.0,   0.0,   0.003,
      0.02, 0.065, 0.14, 0.235, 0.34, 0.47, 0.6, 0.778, 0.902, 1.04},

  //  p pi+ pi- pi- pi0 pi0 pi0 pi0  (n pi+ pi+ pi- pi0 pi0 pi0 pi0)
   {  0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.001,
      0.009, 0.032, 0.07, 0.12, 0.18, 0.25, 0.31, 0.389, 0.451, 0.518},

  //  p pi- pi0 pi0 pi0 pi0 pi0 pi0  (n pi+ pi0 pi0 pi0 pi0 pi0 pi0)
   {  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.001, 0.003, 0.006, 0.01, 0.014, 0.019, 0.022, 0.028, 0.031, 0.034},

  //  n pi0 pi0 pi0 pi0 pi0 pi0 pi0  (p pi0 pi0 pi0 pi0 pi0 pi0 pi0)
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.006, 0.006},

  //  n pi+ pi- pi0 pi0 pi0 pi0 pi0  (p pi+ pi- pi0 pi0 pi0 pi0 pi0)
   {  0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.001,
      0.005, 0.015, 0.034, 0.055, 0.08, 0.11, 0.139, 0.176, 0.196, 0.213},

  //  n pi+ pi+ pi- pi- pi0 pi0 pi0  (p pi+ pi+ pi- pi- pi0 pi0 pi0)
   {  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
      0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.003,
      0.02, 0.065, 0.14, 0.22, 0.32, 0.45, 0.58, 0.778, 0.902, 1.04},

  //  n pi+ pi+ pi+ pi- pi- pi- pi0  (p pi+ pi+ pi+ pi- pi- pi- pi0)
   {  0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0, 0.0,  0.0,   0.0,
      0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0, 0.0,  0.0,   0.002,
      0.013, 0.05, 0.113, 0.17, 0.24, 0.32, 0.4, 0.51, 0.605, 0.69},

  //  L K+ pi+ pi+ pi- pi- pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.004, 0.007, 0.011, 0.014, 0.016, 0.019, 0.021, 0.024},

  //  L K+ pi+ pi- pi- pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.004, 0.007, 0.011, 0.014, 0.016, 0.019, 0.021, 0.024},

  //  L K+ pi- pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003, 0.003},

  //  L K0 pi+ pi+ pi- pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.002, 0.005, 0.011, 0.017, 0.021, 0.024, 0.028, 0.032, 0.036},

  //  L K0 pi+ pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.001, 0.003, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011},

  //  L K0 pi+ pi+ pi+ pi- pi- pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.005, 0.005, 0.005},

  //  L K0 pi0 pi0 pi0 pi0 pi0 pi0  (negligible)

  //  S0 K+ pi+ pi+ pi- pi- pi- pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0, 0.002, 0.004, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011},

  //  S0 K+ pi+ pi- pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0, 0.002, 0.004, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011},

  //  S0 K+ pi- pi0 pi0 pi0 pi0 pi0 (negligible)

  //  S0 K0 pi+ pi+ pi- pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.006, 0.009, 0.011, 0.013, 0.015, 0.017, 0.019},

  //  S0 K0 pi+ pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.005, 0.005, 0.005},

  //  S0 K0 pi+ pi+ pi+ pi- pi- pi-
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003, 0.003},

  //  S0 K0 pi0 pi0 pi0 pi0 pi0 pi0  negligible

  //  S+ K+ pi+ pi- pi- pi- pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.005, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018},

  //  S+ K+ pi- pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004},

  //  S+ K+ pi+ pi+ pi- pi- pi- pi-
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004},

  //  S+ K0 pi+ pi- pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.005, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018},

  //  S+ K0 pi+ pi+ pi- pi- pi- pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.005, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018},

  //  S+ K0 pi- pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002},

  //  S- K+ pi+ pi+ pi- pi- pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.003, 0.008, 0.012, 0.015, 0.018, 0.021, 0.024, 0.027},

  //  S- K+ pi+ pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009},

  //  S- K+ pi+ pi+ pi+ pi- pi- pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.006, 0.006},

  //  S- K+ pi0 pi0 pi0 pi0 pi0 pi0 (negligible)

  //  S- K0 pi+ pi+ pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.005, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018},

  //  S- K0 pi+ pi+ pi+ pi- pi- pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.005, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018},

  //  S- K0 pi+ pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002},

  //  p K+ K- pi+ pi- pi- pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.003, 0.01, 0.017, 0.021, 0.025, 0.029, 0.033, 0.037},

  //  p K+ K- pi+ pi+ pi- pi- pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.01, 0.011, 0.012},

  //  p K+ K- pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.002, 0.003, 0.004, 0.005, 0.006, 0.006, 0.006},

  //  p K+ K0b pi- pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.01, 0.011, 0.012},

  //  p K+ K0b pi+ pi- pi- pi- pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.002, 0.007, 0.011, 0.015, 0.018, 0.02, 0.022, 0.024},

  //  p K0 K0b pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.005, 0.005, 0.005},

  //  p K0 K0b pi+ pi- pi- pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.003, 0.01, 0.016, 0.021, 0.025, 0.029, 0.033, 0.037},

  //  p K0 K0b pi+ pi+ pi- pi- pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.01, 0.011, 0.012},

  //  p K0 K- pi+ pi+ pi- pi- pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.002, 0.006, 0.011, 0.015, 0.018, 0.02, 0.022, 0.024},

  //  p K0 K- pi+ pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.002, 0.006, 0.011, 0.015, 0.018, 0.02, 0.022, 0.024},

  //  p K0 K- pi0 pi0 pi0 pi0 pi0 (negligible)

  //  n K+ K0b pi+ pi- pi- pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.003, 0.009, 0.016, 0.021, 0.025, 0.029, 0.033, 0.037},

  //  n K+ K0b pi+ pi+ pi- pi- pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.012, 0.013},

  //  n K+ K0b pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.005, 0.005, 0.005},

  //  n K+ K- pi+ pi+ pi- pi- pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.003, 0.009, 0.016, 0.021, 0.025, 0.029, 0.033, 0.037},

  //  n K+ K- pi+ pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.002, 0.006, 0.011, 0.015, 0.018, 0.02, 0.022, 0.024},

  //  n K+ K- pi0 pi0 pi0 pi0 pi0 (negligible)

  //  n K0 K0b pi+ pi+ pi- pi- pi0
   {  0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.003, 0.01, 0.016, 0.021, 0.025, 0.029, 0.033, 0.037},

  //  n K0 K0b pi+ pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.002, 0.006, 0.011, 0.015, 0.018, 0.02, 0.022, 0.024},

  //  n K0 K0b pi0 pi0 pi0 pi0 pi0 (negligible)

  //  n K0 K- pi+ pi+ pi- pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.003, 0.01, 0.016, 0.021, 0.025, 0.029, 0.033, 0.037},

  //  n K0 K- pi+ pi+ pi+ pi- pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.012, 0.013},

  //  n K0 K- pi+ pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.005, 0.005},
  //
  //  multiplicity 9 (58 channels)
  //
  //  p pi+ pi+ pi+ pi- pi- pi- pi- pi0
   { 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
     0.001, 0.005, 0.018, 0.04, 0.075, 0.13, 0.192, 0.274, 0.332, 0.387},

  //  p pi+ pi+ pi- pi- pi- pi0 pi0 pi0
   { 0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
     0.002, 0.01, 0.036, 0.08, 0.15, 0.26, 0.385, 0.547, 0.665, 0.77},

  //  p pi+ pi- pi- pi0 pi0 pi0 pi0 pi0
   { 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.001, 0.003, 0.011, 0.024, 0.05, 0.087, 0.129, 0.182, 0.222, 0.26},

  //  p pi- pi0 pi0 pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.006, 0.007, 0.008, 0.009},

  //  n pi0 pi0 pi0 pi0 pi0 pi0 pi0 pi0 
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004},

  //  n pi+ pi- pi0 pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
     0.0, 0.001, 0.004, 0.008, 0.016, 0.026, 0.04, 0.055, 0.066, 0.077},

  //  n pi+ pi+ pi- pi- pi0 pi0 pi0 pi0
   { 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
     0.001, 0.007, 0.024, 0.06, 0.12, 0.195, 0.292, 0.413, 0.499, 0.575},

  //  n pi+ pi+ pi+ pi- pi- pi- pi0 pi0
   { 0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.002, 0.01, 0.036, 0.08, 0.16, 0.26, 0.39, 0.55, 0.66, 0.77},

  //  n pi+ pi+ pi+ pi+ pi- pi- pi- pi-
   { 0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0, 0.0,   0.0,   0.0,
     0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0, 0.0,   0.0,   0.0,
     0.001, 0.005, 0.018, 0.04, 0.08, 0.13, 0.2, 0.274, 0.332, 0.385},

  //  L K+ pi+ pi+ pi- pi- pi- pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.006, 0.011, 0.017, 0.022, 0.027, 0.032, 0.037},

  //  L K+ pi+ pi- pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.012, 0.014, 0.016, 0.018},

  //  L K+ pi+ pi+ pi+ pi- pi- pi- pi-
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.005, 0.005},

  //  L K+ pi- pi0 pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002},

  //  L K0 pi+ pi+ pi- pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.006, 0.011, 0.017, 0.022, 0.027, 0.032, 0.037},

  //  L K0 pi+ pi+ pi+ pi- pi- pi- pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.004, 0.008, 0.012, 0.015, 0.018, 0.021, 0.025},

  //  L K0 pi+ pi- pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.005, 0.005},

  //  L K0 pi0 pi0 pi0 pi0 pi0 pi0 pi0  (negligible)

  //  S0 K+ pi+ pi+ pi- pi- pi- pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.011, 0.014, 0.016, 0.019},

  //  S0 K+ pi+ pi- pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.006, 0.007, 0.008, 0.009},

  //  S0 K+ pi+ pi+ pi+ pi- pi- pi- pi-
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004},

  //  S0 K+ pi- pi0 pi0 pi0 pi0 pi0 pi0  (negligible)

  //  S0 K0 pi+ pi+ pi- pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.011, 0.014, 0.016, 0.019},

  //  S0 K0 pi+ pi+ pi+ pi- pi- pi- pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014},

  //  S0 K0 pi+ pi- pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004},

  //  S0 K0 pi0 pi0 pi0 pi0 pi0 pi0 pi0  (negligible)

  //  S+ K+ pi+ pi- pi- pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.011, 0.014, 0.016, 0.019},

  //  S+ K+ pi+ pi+ pi- pi- pi- pi- pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013, 0.015},

  //  S+ K+ pi- pi- pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  S+ K0 pi+ pi+ pi- pi- pi- pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.002, 0.006, 0.01, 0.014, 0.016, 0.02, 0.024, 0.028},

  //  S+ K0 pi+ pi- pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013, 0.015},

  //  S+ K0 pi+ pi+ pi+ pi- pi- pi- pi-
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006},

  //  S+ K0 pi- pi0 pi0 pi0 pi0 pi0 pi0 (negligible)

  //  S- K+ pi+ pi+ pi- pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.002, 0.005, 0.008, 0.012, 0.016, 0.02, 0.024, 0.028},

  //  S- K+ pi+ pi+ pi+ pi- pi- pi- pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.005, 0.008, 0.011, 0.014, 0.016, 0.019},

  //  S- K+ pi+ pi- pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.005, 0.005},

  //  S- K+ pi0 pi0 pi0 pi0 pi0 pi0 pi0 (negligible)

  //  S- K0 pi+ pi+ pi+ pi- pi- pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.002, 0.005, 0.008, 0.012, 0.016, 0.02, 0.024, 0.028},

  //  S- K0 pi+ pi+ pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014},

  //  S- K0 pi+ pi+ pi+ pi+ pi- pi- pi-
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.005, 0.005},

  //  S- K0 pi+ pi0 pi0 pi0 pi0 pi0 pi0 (negligible)

  //  p K+ K0b pi+ pi- pi- pi- pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.002, 0.006, 0.011, 0.014, 0.018, 0.022, 0.026},

  //  p K+ K0b pi+ pi+ pi- pi- pi- pi-
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007},

  //  p K+ K0b pi- pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007},

  //  p K+ K- pi+ pi- pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.002, 0.006, 0.011, 0.014, 0.018, 0.022, 0.026},

  //  p K+ K- pi+ pi+ pi- pi- pi- pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.002, 0.006, 0.011, 0.014, 0.018, 0.022, 0.026},

  //  p K+ K- pi- pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003},

  //  p K0 K0b pi+ pi- pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.002, 0.006, 0.011, 0.014, 0.018, 0.022, 0.026},

  //  p K0 K0b pi+ pi+ pi- pi- pi- pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.002, 0.006, 0.011, 0.014, 0.018, 0.022, 0.026},

  //  p K0 K0b pi- pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003},

  //  p K0 K- pi+ pi+ pi- pi- pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.003, 0.009, 0.016, 0.021, 0.027, 0.033, 0.039},

  //  p K0 K- pi+ pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013},

  //  p K0 K- pi+ pi+ pi+ pi- pi- pi-
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008},

  //  p K0 K- pi0 pi0 pi0 pi0 pi0 pi0  (negligible)

  //  n K+ K0b pi+ pi- pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.002, 0.006, 0.011, 0.014, 0.018, 0.022, 0.026},

  //  n K+ K0b pi+ pi+ pi- pi- pi- pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.002, 0.006, 0.011, 0.014, 0.018, 0.022, 0.026},

  //  n K+ K0b pi- pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003},

  //  n K+ K- pi+ pi+ pi- pi- pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.003, 0.009, 0.016, 0.021, 0.027, 0.033, 0.039},

  //  n K+ K- pi+ pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013},

  //  n K+ K- pi+ pi+ pi+ pi- pi- pi-
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008},

  //  n K+ K- pi0 pi0 pi0 pi0 pi0 pi0 (negligible)

  //  n K0 K0b pi+ pi+ pi- pi- pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.003, 0.01, 0.016, 0.021, 0.027, 0.033, 0.039},

  //  n K0 K0b pi+ pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013},

  //  n K0 K0b pi+ pi+ pi+ pi- pi- pi-
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008},

  //  n K0 K0b pi0 pi0 pi0 pi0 pi0 pi0 (negligible)

  //  n K0 K- pi+ pi+ pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.002, 0.007, 0.011, 0.014, 0.018, 0.022, 0.026},

  //  n K0 K- pi+ pi+ pi+ pi- pi- pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.002, 0.007, 0.011, 0.014, 0.018, 0.022, 0.026},

  //  n K0 K- pi+ pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003} };
 }

// Initialize both |T Tz> = |3/2 1/2> channels, using pipP cross-section table

const G4CascadePiMinusPChannelData::data_t
G4CascadePiMinusPChannelData::data(pimP2bfs, pimP3bfs, pimP4bfs,
				   pimP5bfs, pimP6bfs, pimP7bfs,
				   pimP8bfs, pimP9bfs, pimPCrossSections,
				   pimPtotXSec, pim*pro, "PiMinusP");

const G4CascadePiPlusNChannelData::data_t
G4CascadePiPlusNChannelData::data(pipN2bfs, pipN3bfs, pipN4bfs,
				  pipN5bfs, pipN6bfs, pipN7bfs,
				  pipN8bfs, pipN9bfs, pimPCrossSections,
				  pimPtotXSec, pip*neu, "PiPlusN");
