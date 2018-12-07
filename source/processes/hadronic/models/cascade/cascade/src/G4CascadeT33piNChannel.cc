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

#include "G4CascadePiPlusPChannel.hh"
#include "G4CascadePiMinusNChannel.hh"
#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;

namespace {
  // pi+ p : Outgoing particle types of a given multiplicity

  static const G4int pipP2bfs[2][2] =
   {{pip,pro}, {kpl,sp}};

  static const G4int pipP3bfs[7][3] =
   {{pip,pro,pi0}, {pip,neu,pip}, {pi0,sp,kpl}, {pip,sp,k0}, 
    {pip,s0,kpl},  {pip,lam,kpl}, {kpl,pro,k0b}};

  static const G4int pipP4bfs[15][4] =
   {{pip,pro,pip,pim}, {pip,pro,pi0,pi0}, {pip,neu,pip,pi0},
    {pip,sp,kpl,pim},  {pi0,sp,kpl,pi0},  {pip,sp,k0,pi0},
    {pip,s0,k0,pip},   {pip,s0,kpl,pi0},  {pip,lam,kpl,pi0},
    {pip,lam,k0,pip},  {pip,sm,kpl,pip},  {pip,pro,kpl,kmi},
    {pip,pro,k0,k0b},  {pi0,pro,kpl,k0b}, {pip,neu,kpl,k0b}};

  static const G4int pipP5bfs[24][5] =
   {{pip,pro,pip,pim,pi0}, {pip,pro,pi0,pi0,pi0}, {pip,neu,pip,pip,pim},
    {pip,neu,pip,pi0,pi0}, {pip,sp,kpl,pim,pi0},  {pi0,sp,kpl,pi0,pi0},
    {pip,sp,k0,pip,pim},   {pip,sp,k0,pi0,pi0},   {pip,lam,k0,pip,pi0},
    {pip,lam,kpl,pip,pim}, {pip,lam,kpl,pi0,pi0}, {pip,s0,kpl,pip,pim},
    {pip,s0,kpl,pi0,pi0},  {pip,s0,k0,pip,pi0},   {pip,sm,kpl,pip,pi0},
    {pip,sm,k0,pip,pip},   {pip,pro,pim,kpl,k0b}, {pip,pro,pip,k0,kmi},
    {pip,pro,pi0,kpl,kmi}, {pip,pro,pi0,k0,k0b},  {pi0,pro,pi0,kpl,k0b},
    {pip,neu,pip,kpl,kmi}, {pip,neu,pip,k0,k0b},  {pip,neu,pi0,kpl,k0b}};

  static const G4int pipP6bfs[33][6] =
   {{pip,pro,pip,pip,pim,pim}, {pip,pro,pip,pim,pi0,pi0},
    {pip,pro,pi0,pi0,pi0,pi0}, {pip,neu,pip,pi0,pi0,pi0},
    {pip,neu,pip,pip,pim,pi0}, {pip,pip,pim,pi0,kpl,lam},
    {pip,pi0,pi0,pi0,kpl,lam}, {pip,pip,pi0,pi0,k0,lam},
    {pip,pip,pip,pim,k0,lam},  {pip,pip,pim,pim,kpl,sp},
    {pip,pim,pi0,pi0,kpl,sp},  {pi0,pi0,pi0,pi0,kpl,sp},
    {pip,pip,pim,pi0,k0,sp},   {pip,pi0,pi0,pi0,k0,sp}, 
    {pip,pip,pim,pi0,kpl,s0},  {pip,pi0,pi0,pi0,kpl,s0},
    {pip,pip,pip,pim,k0,s0},   {pip,pip,pi0,pi0,k0,s0},
    {pip,pip,pi0,pi0,kpl,sm},  {pip,pip,pip,pim,kpl,sm},
    {pip,pip,pip,pi0,k0,sm},   {pip,pro,pim,pi0,kpl,k0b},
    {pi0,pro,pi0,pi0,kpl,k0b}, {pip,pro,pi0,pi0,kpl,kmi},
    {pip,pro,pip,pim,kpl,kmi}, {pip,pro,pi0,pi0,k0,k0b},
    {pip,pro,pip,pim,k0,k0b},  {pip,pro,pip,pi0,kmi,k0},
    {pip,neu,pi0,pi0,kpl,k0b}, {pip,neu,pip,pim,kpl,k0b},
    {pip,neu,pip,pi0,kpl,kmi}, {pip,neu,pip,pim,kpl,k0b},
    {pip,neu,pip,pip,k0,kmi}};

  static const G4int pipP7bfs[41][7] =
   {{pip,pro,pip,pip,pim,pim,pi0}, {pip,pro,pip,pim,pi0,pi0,pi0},
    {pip,pro,pi0,pi0,pi0,pi0,pi0}, {pip,neu,pip,pip,pip,pim,pim},
    {pip,neu,pip,pip,pim,pi0,pi0}, {pip,neu,pip,pi0,pi0,pi0,pi0},
    {pip,pip,pim,pi0,pi0,kpl,lam}, {pip,pip,pip,pim,pim,kpl,lam},
    {pip,pi0,pi0,pi0,pi0,kpl,lam}, {pip,pip,pip,pim,pi0,k0,lam},
    {pip,pip,pi0,pi0,pi0,k0,lam},  {pip,pip,pim,pim,pi0,kpl,sp},
    {pip,pim,pi0,pi0,pi0,kpl,sp},  {pip,pip,pim,pi0,pi0,k0,sp},
    {pip,pip,pip,pim,pim,k0,sp},   {pip,pi0,pi0,pi0,pi0,k0,sp},
    {pip,pip,pim,pi0,pi0,kpl,s0},  {pip,pip,pip,pim,pim,kpl,s0},
    {pip,pi0,pi0,pi0,pi0,kpl,s0},  {pip,pip,pip,pim,pi0,k0,s0},
    {pip,pip,pi0,pi0,pi0,k0,s0},   {pip,pip,pip,pim,pi0,kpl,sm},
    {pip,pip,pi0,pi0,pi0,kpl,sm},  {pip,pip,pip,pip,pim,k0,sm},
    {pip,pip,pip,pi0,pi0,k0,sm},   {pip,pro,pim,pi0,pi0,kpl,k0b},
    {pip,pro,pip,pim,pim,kpl,k0b}, {pi0,pro,pi0,pi0,pi0,kpl,k0b},
    {pip,pro,pip,pim,pi0,kpl,kmi}, {pip,pro,pi0,pi0,pi0,kpl,kmi},
    {pip,pro,pip,pim,pi0,k0,k0b},  {pip,pro,pi0,pi0,pi0,k0,k0b},
    {pip,pro,pip,pi0,pi0,kmi,k0},  {pip,pro,pip,pip,pim,kmi,k0},
    {pip,neu,pip,pim,pi0,kpl,k0b}, {pip,neu,pi0,pi0,pi0,kpl,k0b},
    {pip,neu,pip,pi0,pi0,kpl,kmi}, {pip,neu,pip,pip,pim,kpl,kmi},
    {pip,neu,pip,pi0,pi0,k0,k0b},  {pip,neu,pip,pip,pim,k0,k0b},
    {pip,neu,pip,pip,pi0,kmi,k0}};

  static const G4int pipP8bfs[47][8] =
   {{pip,pro,pip,pip,pip,pim,pim,pim}, {pip,pro,pip,pip,pim,pim,pi0,pi0},
    {pip,pro,pip,pim,pi0,pi0,pi0,pi0}, {pip,pro,pi0,pi0,pi0,pi0,pi0,pi0},
    {pip,neu,pip,pip,pip,pim,pim,pi0}, {pip,neu,pip,pip,pim,pi0,pi0,pi0},
    {pip,neu,pip,pi0,pi0,pi0,pi0,pi0}, {pip,pip,pip,pim,pim,pi0,kpl,lam},
    {pip,pi0,pi0,pi0,pi0,pi0,kpl,lam}, {pip,pip,pip,pim,pi0,pi0,k0,lam},
    {pip,pip,pi0,pi0,pi0,pi0,k0,lam},  {pip,pip,pip,pip,pim,pim,k0,lam},
    {pip,pip,pim,pim,pi0,pi0,kpl,sp},  {pip,pim,pi0,pi0,pi0,pi0,kpl,sp},
    {pip,pip,pip,pim,pim,pim,kpl,sp},  {pip,pip,pip,pim,pim,pi0,k0,sp}, 
    {pip,pip,pim,pi0,pi0,pi0,k0,sp},   {pip,pip,pip,pim,pim,pi0,kpl,s0},
    {pip,pi0,pi0,pi0,pi0,pi0,kpl,s0},  {pip,pip,pip,pim,pi0,pi0,k0,s0},
    {pip,pip,pi0,pi0,pi0,pi0,k0,s0},   {pip,pip,pip,pip,pim,pim,k0,s0},
    {pip,pip,pip,pim,pi0,pi0,kpl,sm},  {pip,pip,pi0,pi0,pi0,pi0,kpl,sm},
    {pip,pip,pip,pip,pim,pim,kpl,sm},  {pip,pip,pip,pip,pim,pi0,k0,sm}, 
    {pip,pip,pip,pi0,pi0,pi0,k0,sm},   {pip,pro,pip,pim,pim,pi0,kpl,k0b},
    {pip,pro,pim,pi0,pi0,pi0,kpl,k0b}, {pi0,pro,pi0,pi0,pi0,pi0,kpl,k0b},
    {pip,pro,pip,pim,pi0,pi0,kpl,kmi}, {pip,pro,pip,pip,pim,pim,kpl,kmi},
    {pip,pro,pi0,pi0,pi0,pi0,kpl,kmi}, {pip,pro,pip,pim,pi0,pi0,k0,k0b},
    {pip,pro,pip,pip,pim,pim,k0,k0b},  {pip,pro,pi0,pi0,pi0,pi0,k0,k0b},
    {pip,pro,pip,pip,pim,pi0,k0,kmi},  {pip,pro,pip,pi0,pi0,pi0,k0,kmi},
    {pip,neu,pip,pip,pim,pim,kpl,k0b}, {pip,neu,pip,pim,pi0,pi0,kpl,k0b},
    {pip,neu,pi0,pi0,pi0,pi0,kpl,k0b}, {pip,neu,pip,pip,pim,pi0,kpl,kmi},
    {pip,neu,pip,pi0,pi0,pi0,kpl,kmi}, {pip,neu,pip,pip,pim,pi0,k0,k0b},
    {pip,neu,pip,pi0,pi0,pi0,k0,k0b},  {pip,neu,pip,pip,pi0,pi0,k0,kmi},
    {pip,neu,pip,pip,pip,pim,k0,kmi}};

  static const G4int pipP9bfs[55][9] =
   {{pip,pro,pip,pip,pip,pim,pim,pim,pi0}, {pip,pro,pip,pip,pim,pim,pi0,pi0,pi0},
    {pip,pro,pip,pim,pi0,pi0,pi0,pi0,pi0}, {pip,pro,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
    {pip,neu,pip,pip,pip,pip,pim,pim,pim}, {pip,neu,pip,pip,pip,pim,pim,pi0,pi0},
    {pip,neu,pip,pip,pim,pi0,pi0,pi0,pi0}, {pip,neu,pip,pi0,pi0,pi0,pi0,pi0,pi0},
    {pip,pip,pip,pim,pim,pi0,pi0,kpl,lam}, {pip,pip,pim,pi0,pi0,pi0,pi0,kpl,lam},
    {pip,pip,pip,pip,pim,pim,pim,kpl,lam}, {pip,pi0,pi0,pi0,pi0,pi0,pi0,kpl,lam}, 
    {pip,pip,pip,pim,pi0,pi0,pi0,k0,lam},  {pip,pip,pip,pip,pim,pim,pi0,k0,lam},
    {pip,pip,pi0,pi0,pi0,pi0,pi0,k0,lam},  {pip,pip,pim,pim,pi0,pi0,pi0,kpl,sp},
    {pip,pip,pip,pim,pim,pim,pi0,kpl,sp},  {pip,pim,pi0,pi0,pi0,pi0,pi0,kpl,sp},
    {pip,pip,pip,pim,pim,pi0,pi0,k0,sp},   {pip,pip,pim,pi0,pi0,pi0,pi0,k0,sp},
    {pip,pip,pip,pip,pim,pim,pim,k0,sp},   {pip,pip,pip,pim,pim,pi0,pi0,kpl,s0},
    {pip,pip,pim,pi0,pi0,pi0,pi0,kpl,s0},  {pip,pip,pip,pip,pim,pim,pim,kpl,s0},
    {pip,pip,pip,pim,pi0,pi0,pi0,k0,s0},   {pip,pip,pip,pip,pim,pim,pi0,k0,s0},
    {pip,pip,pi0,pi0,pi0,pi0,pi0,k0,s0},   {pip,pip,pip,pim,pi0,pi0,pi0,kpl,sm},
    {pip,pip,pip,pip,pim,pim,pi0,kpl,sm},  {pip,pip,pi0,pi0,pi0,pi0,pi0,kpl,sm},
    {pip,pip,pip,pip,pim,pi0,pi0,k0,sm},   {pip,pip,pip,pi0,pi0,pi0,pi0,k0,sm},
    {pip,pro,pip,pim,pim,pi0,pi0,kpl,k0b}, {pip,pro,pim,pi0,pi0,pi0,pi0,kpl,k0b},
    {pip,pro,pip,pip,pim,pim,pim,kpl,k0b}, {pip,pro,pip,pim,pi0,pi0,pi0,kpl,kmi},
    {pip,pro,pip,pip,pim,pim,pi0,kpl,kmi}, {pip,pro,pi0,pi0,pi0,pi0,pi0,kpl,kmi}, 
    {pip,pro,pip,pim,pi0,pi0,pi0,k0,k0b},  {pip,pro,pip,pip,pim,pim,pi0,k0,k0b},
    {pip,pro,pi0,pi0,pi0,pi0,pi0,k0,k0b},  {pip,pro,pip,pip,pim,pi0,pi0,k0,kmi},
    {pip,pro,pip,pi0,pi0,pi0,pi0,k0,kmi},  {pip,pro,pip,pip,pip,pim,pim,k0,kmi},
    {pip,neu,pip,pim,pi0,pi0,pi0,kpl,k0b}, {pip,neu,pip,pip,pim,pim,pi0,kpl,k0b},
    {pip,neu,pi0,pi0,pi0,pi0,pi0,kpl,k0b}, {pip,neu,pip,pip,pim,pi0,pi0,kpl,kmi},
    {pip,neu,pip,pi0,pi0,pi0,pi0,kpl,kmi}, {pip,neu,pip,pip,pip,pim,pim,kpl,kmi},
    {pip,neu,pip,pip,pim,pi0,pi0,k0,k0b},  {pip,neu,pip,pi0,pi0,pi0,pi0,k0,k0b},
    {pip,neu,pip,pip,pip,pim,pim,k0,k0b},  {pip,neu,pip,pip,pip,pim,pi0,kmi,k0},
    {pip,neu,pip,pip,pi0,pi0,pi0,kmi,k0}};
}

namespace {
  // pi- n : Outgoing particle types of a given multiplicity

  static const G4int pimN2bfs[2][2] =
   {{pim,neu}, {k0,sm}};

  static const G4int pimN3bfs[7][3] =
   {{pim,neu,pi0}, {pim,pro,pim}, {pi0,sm,k0}, {pim,sm,kpl},
    {pim,s0,k0},   {pim,lam,k0},  {k0,neu,kmi}};

  static const G4int pimN4bfs[15][4] =
   {{pim,neu,pip,pim},{pim,neu,pi0,pi0},{pim,pro,pim,pi0},
    {pim,sm,k0,pip},  {pi0,sm,k0,pi0},  {pim,sm,kpl,pi0},
    {pim,s0,kpl,pim},  {pim,s0,k0,pi0},  {pim,lam,k0,pi0},
    {pim,lam,kpl,pim}, {pim,sp,k0,pim},  {pim,neu,k0,k0b},
    {pim,neu,kpl,kmi},  {pi0,neu,k0,kmi},  {pim,pro,k0,kmi}};

  static const G4int pimN5bfs[24][5] =
   {{pim,neu,pip,pim,pi0}, {pim,neu,pi0,pi0,pi0}, {pim,pro,pip,pim,pim},
    {pim,pro,pim,pi0,pi0}, {pim,sm,k0,pip,pi0},   {pi0,sm,k0,pi0,pi0},
    {pim,sm,kpl,pip,pim},  {pim,sm,kpl,pi0,pi0},  {pim,lam,kpl,pim,pi0},
    {pim,lam,k0,pip,pim},  {pim,lam,k0,pi0,pi0},  {pim,s0,k0,pip,pim},
    {pim,s0,k0,pi0,pi0},   {pim,s0,kpl,pim,pi0},  {pim,sp,k0,pim,pi0},
    {pim,sp,kpl,pim,pim},  {pim,neu,pip,k0,kmi},  {pim,neu,pim,kpl,k0b},
    {pim,neu,pi0,k0,k0b},  {pim,neu,pi0,kpl,kmi}, {pi0,neu,pi0,k0,kmi},
    {pim,pro,pim,k0,k0b},  {pim,pro,pim,kpl,kmi}, {pim,pro,pi0,k0,kmi}};

  static const G4int pimN6bfs[33][6] =
   {{pim,neu,pip,pip,pim,pim}, {pim,neu,pip,pim,pi0,pi0},
    {pim,neu,pi0,pi0,pi0,pi0}, {pim,pro,pim,pi0,pi0,pi0},
    {pim,pro,pip,pim,pim,pi0}, {pim,pip,pim,pi0,k0,lam},
    {pim,pi0,pi0,pi0,k0,lam},  {pim,pim,pi0,pi0,kpl,lam},
    {pim,pip,pim,pim,kpl,lam}, {pim,pip,pip,pim,k0,sm},
    {pim,pip,pi0,pi0,k0,sm},   {pi0,pi0,pi0,pi0,k0,sm},
    {pim,pip,pim,pi0,kpl,sm},  {pim,pi0,pi0,pi0,kpl,sm},
    {pim,pip,pim,pi0,k0,s0},   {pim,pi0,pi0,pi0,k0,s0},
    {pim,pip,pim,pim,kpl,s0},  {pim,pim,pi0,pi0,kpl,s0},
    {pim,pim,pi0,pi0,k0,sp},   {pim,pip,pim,pim,k0,sp},
    {pim,pim,pim,pi0,kpl,sp},  {pim,neu,pip,pi0,kmi,k0},
    {pi0,neu,pi0,pi0,kmi,k0},  {pim,neu,pi0,pi0,k0,k0b},
    {pim,neu,pip,pim,k0,k0b},  {pim,neu,pi0,pi0,kpl,kmi},
    {pim,neu,pip,pim,kpl,kmi}, {pim,neu,pim,pi0,kpl,k0b},
    {pim,pro,pi0,pi0,kmi,k0},  {pim,pro,pip,pim,kmi,k0},
    {pim,pro,pim,pi0,k0,k0b},  {pim,pro,pip,pim,kmi,k0},
    {pim,pro,pim,pim,kpl,k0b}};

  static const G4int pimN7bfs[41][7] =
   {{pim,neu,pip,pip,pim,pim,pi0}, {pim,neu,pip,pim,pi0,pi0,pi0},
    {pim,neu,pi0,pi0,pi0,pi0,pi0}, {pim,pro,pip,pip,pim,pim,pim},
    {pim,pro,pip,pim,pim,pi0,pi0}, {pim,pro,pim,pi0,pi0,pi0,pi0},
    {pim,pip,pim,pi0,pi0,k0,lam},  {pim,pip,pip,pim,pim,k0,lam},
    {pim,pi0,pi0,pi0,pi0,k0,lam},  {pim,pip,pim,pim,pi0,kpl,lam},
    {pim,pim,pi0,pi0,pi0,kpl,lam}, {pim,pip,pip,pim,pi0,k0,sm},
    {pim,pip,pi0,pi0,pi0,k0,sm},   {pim,pip,pim,pi0,pi0,kpl,sm},
    {pim,pip,pip,pim,pim,kpl,sm},  {pim,pi0,pi0,pi0,pi0,kpl,sm},
    {pim,pip,pim,pi0,pi0,k0,s0},   {pim,pip,pip,pim,pim,k0,s0},
    {pim,pi0,pi0,pi0,pi0,k0,s0},   {pim,pip,pim,pim,pi0,kpl,s0},
    {pim,pim,pi0,pi0,pi0,kpl,s0},  {pim,pip,pim,pim,pi0,k0,sp},
    {pim,pim,pi0,pi0,pi0,k0,sp},   {pim,pip,pim,pim,pim,kpl,sp},
    {pim,pim,pim,pi0,pi0,kpl,sp},  {pim,neu,pip,pi0,pi0,kmi,k0},
    {pim,neu,pip,pip,pim,kmi,k0},  {pi0,neu,pi0,pi0,pi0,kmi,k0},
    {pim,neu,pip,pim,pi0,k0,k0b},  {pim,neu,pi0,pi0,pi0,k0,k0b},
    {pim,neu,pip,pim,pi0,kpl,kmi}, {pim,neu,pi0,pi0,pi0,kpl,kmi},
    {pim,neu,pim,pi0,pi0,kpl,k0b}, {pim,neu,pip,pim,pim,kpl,k0b},
    {pim,pro,pip,pim,pi0,kmi,k0},  {pim,pro,pi0,pi0,pi0,kmi,k0},
    {pim,pro,pim,pi0,pi0,k0,k0b},  {pim,pro,pip,pim,pim,k0,k0b},
    {pim,pro,pim,pi0,pi0,kpl,kmi}, {pim,pro,pip,pim,pim,kpl,kmi},
    {pim,pro,pim,pim,pi0,kpl,k0b}};

  static const G4int pimN8bfs[47][8] =
   {{pim,neu,pip,pip,pip,pim,pim,pim}, {pim,neu,pip,pip,pim,pim,pi0,pi0},
    {pim,neu,pip,pim,pi0,pi0,pi0,pi0}, {pim,neu,pi0,pi0,pi0,pi0,pi0,pi0},
    {pim,pro,pip,pip,pim,pim,pim,pi0}, {pim,pro,pip,pim,pim,pi0,pi0,pi0},
    {pim,pro,pim,pi0,pi0,pi0,pi0,pi0}, {pim,pip,pip,pim,pim,pi0,k0,lam},
    {pim,pi0,pi0,pi0,pi0,pi0,k0,lam},  {pim,pip,pim,pim,pi0,pi0,kpl,lam},
    {pim,pim,pi0,pi0,pi0,pi0,kpl,lam}, {pim,pip,pip,pim,pim,pim,kpl,lam},
    {pim,pip,pip,pim,pi0,pi0,k0,sm},   {pim,pip,pi0,pi0,pi0,pi0,k0,sm},
    {pim,pip,pip,pip,pim,pim,k0,sm},   {pim,pip,pip,pim,pim,pi0,kpl,sm},
    {pim,pip,pim,pi0,pi0,pi0,kpl,sm},  {pim,pip,pip,pim,pim,pi0,k0,s0},
    {pim,pi0,pi0,pi0,pi0,pi0,k0,s0},   {pim,pip,pim,pim,pi0,pi0,kpl,s0},
    {pim,pim,pi0,pi0,pi0,pi0,kpl,s0},  {pim,pip,pip,pim,pim,pim,kpl,s0},
    {pim,pip,pim,pim,pi0,pi0,k0,sp},   {pim,pim,pi0,pi0,pi0,pi0,k0,sp},
    {pim,pip,pip,pim,pim,pim,k0,sp},   {pim,pip,pim,pim,pim,pi0,kpl,sp},
    {pim,pim,pim,pi0,pi0,pi0,kpl,sp},  {pim,neu,pip,pip,pim,pi0,kmi,k0},
    {pim,neu,pip,pi0,pi0,pi0,kmi,k0},  {pi0,neu,pi0,pi0,pi0,pi0,kmi,k0},
    {pim,neu,pip,pim,pi0,pi0,k0,k0b},  {pim,neu,pip,pip,pim,pim,k0,k0b},
    {pim,neu,pi0,pi0,pi0,pi0,k0,k0b},  {pim,neu,pip,pim,pi0,pi0,kpl,kmi},
    {pim,neu,pip,pip,pim,pim,kpl,kmi}, {pim,neu,pi0,pi0,pi0,pi0,kpl,kmi},
    {pim,neu,pip,pim,pim,pi0,kpl,k0b}, {pim,neu,pim,pi0,pi0,pi0,kpl,k0b},
    {pim,pro,pip,pip,pim,pim,kmi,k0},  {pim,pro,pip,pim,pi0,pi0,kmi,k0},
    {pim,pro,pi0,pi0,pi0,pi0,kmi,k0},  {pim,pro,pip,pim,pim,pi0,k0,k0b},
    {pim,pro,pim,pi0,pi0,pi0,k0,k0b},  {pim,pro,pip,pim,pim,pi0,kpl,kmi},
    {pim,pro,pim,pi0,pi0,pi0,kpl,kmi}, {pim,pro,pim,pim,pi0,pi0,kpl,k0b},
    {pim,pro,pip,pim,pim,pim,kpl,k0b}};

  static const G4int pimN9bfs[55][9] =
   {{pim,neu,pip,pip,pip,pim,pim,pim,pi0}, {pim,neu,pip,pip,pim,pim,pi0,pi0,pi0},
    {pim,neu,pip,pim,pi0,pi0,pi0,pi0,pi0}, {pim,neu,pi0,pi0,pi0,pi0,pi0,pi0,pi0},
    {pim,pro,pip,pip,pip,pim,pim,pim,pim}, {pim,pro,pip,pip,pim,pim,pim,pi0,pi0},
    {pim,pro,pip,pim,pim,pi0,pi0,pi0,pi0}, {pim,pro,pim,pi0,pi0,pi0,pi0,pi0,pi0},
    {pim,pip,pip,pim,pim,pi0,pi0,k0,lam},  {pim,pip,pim,pi0,pi0,pi0,pi0,k0,lam},
    {pim,pip,pip,pip,pim,pim,pim,k0,lam},  {pim,pi0,pi0,pi0,pi0,pi0,pi0,k0,lam},
    {pim,pip,pim,pim,pi0,pi0,pi0,kpl,lam}, {pim,pip,pip,pim,pim,pim,pi0,kpl,lam},
    {pim,pim,pi0,pi0,pi0,pi0,pi0,kpl,lam}, {pim,pip,pip,pim,pi0,pi0,pi0,k0,sm},
    {pim,pip,pip,pip,pim,pim,pi0,k0,sm},   {pim,pip,pi0,pi0,pi0,pi0,pi0,k0,sm},
    {pim,pip,pip,pim,pim,pi0,pi0,kpl,sm},  {pim,pip,pim,pi0,pi0,pi0,pi0,kpl,sm},
    {pim,pip,pip,pip,pim,pim,pim,kpl,sm},  {pim,pip,pip,pim,pim,pi0,pi0,k0,s0},
    {pim,pip,pim,pi0,pi0,pi0,pi0,k0,s0},   {pim,pip,pip,pip,pim,pim,pim,k0,s0},
    {pim,pip,pim,pim,pi0,pi0,pi0,kpl,s0},  {pim,pip,pip,pim,pim,pim,pi0,kpl,s0},
    {pim,pim,pi0,pi0,pi0,pi0,pi0,kpl,s0},  {pim,pip,pim,pim,pi0,pi0,pi0,k0,sp},
    {pim,pip,pip,pim,pim,pim,pi0,k0,sp},   {pim,pim,pi0,pi0,pi0,pi0,pi0,k0,sp},
    {pim,pip,pim,pim,pim,pi0,pi0,kpl,sp},  {pim,pim,pim,pi0,pi0,pi0,pi0,kpl,sp},
    {pim,neu,pip,pip,pim,pi0,pi0,kmi,k0},  {pim,neu,pip,pi0,pi0,pi0,pi0,kmi,k0},
    {pim,neu,pip,pip,pip,pim,pim,kmi,k0},  {pim,neu,pip,pim,pi0,pi0,pi0,k0,k0b},
    {pim,neu,pip,pip,pim,pim,pi0,k0,k0b},  {pim,neu,pi0,pi0,pi0,pi0,pi0,k0,k0b},
    {pim,neu,pip,pim,pi0,pi0,pi0,kpl,kmi}, {pim,neu,pip,pip,pim,pim,pi0,kpl,kmi},
    {pim,neu,pi0,pi0,pi0,pi0,pi0,kpl,kmi}, {pim,neu,pip,pim,pim,pi0,pi0,kpl,k0b},
    {pim,neu,pim,pi0,pi0,pi0,pi0,kpl,k0b}, {pim,neu,pip,pip,pim,pim,pim,kpl,k0b},
    {pim,pro,pip,pim,pi0,pi0,pi0,kmi,k0},  {pim,pro,pip,pip,pim,pim,pi0,kmi,k0},
    {pim,pro,pi0,pi0,pi0,pi0,pi0,kmi,k0},  {pim,pro,pip,pim,pim,pi0,pi0,k0,k0b},
    {pim,pro,pim,pi0,pi0,pi0,pi0,k0,k0b},  {pim,pro,pip,pip,pim,pim,pim,k0,k0b},
    {pim,pro,pip,pim,pim,pi0,pi0,kpl,kmi}, {pim,pro,pim,pi0,pi0,pi0,pi0,kpl,kmi},
    {pim,pro,pip,pip,pim,pim,pim,kpl,kmi}, {pim,pro,pip,pim,pim,pim,pi0,kpl,k0b},
    {pim,pro,pim,pim,pi0,pi0,pi0,kpl,k0b}}; 
}

namespace {
  // Total pi+ p cross sections as a function of kinetic energy
  // New cs after 7, 8, 9-body strange channel additions (27 Oct 17)
  static const G4double pipPtotXSec[30] = 
   {  0.0,   1.2,   2.5,   3.8,   5.0,  7.0,   9.0,  15.0,  30.0,  64.0,
    130.0, 190.0, 130.01, 55.85, 28.0, 17.14, 19.27, 27.58, 40.07, 32.51,
     30.46, 28.61, 27.26, 25.81, 25.0, 24.38, 24.0,  23.5,  23.0,  23.0};

  // pi+ p cross sections as functions of kinetic energy and multiplicity
  static const G4double pipPCrossSections[224][30] = {
  //
  // multiplicity 2 (2 channels)
  //
  //  p pi+ (n pi-)
   {  0.0,   1.2,   2.5,  3.8,  5.0,  7.00, 9.00, 15.0,  30.0,  64.0,
    130.0, 190.0, 130.0, 55.7, 27.2, 13.95, 8.38, 12.98, 18.53, 11.81,
      9.4, 7.7,     6.8,  5.9,  5.2,  4.7,  4.429, 4.275, 4.159, 4.109},

  //  S+ K+ (S- K0)
   {  0.0,  0.0, 0.0, 0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,  0.0, 0.0, 0.0,  0.0,   0.0,   0.0,   0.16,  0.62,  0.34,
      0.25, 0.2, 0.1, 0.06, 0.036, 0.025, 0.018, 0.014, 0.012, 0.011},
  //
  // multiplicity 3 (7 channels)
  //
  //  p pi+ pi0 (n pi- pi0)
   {  0.0, 0.0,  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0,  0.007, 0.11, 0.6,  2.39,  8.67, 9.99, 12.36,  6.66,
      4.4, 2.85, 1.85,  1.25, 0.83, 0.57, 0.432, 0.326, 0.251, 0.194},

  //  n pi+ pi+ (p pi- pi-)
   {  0.0, 0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.003, 0.04, 0.2,  0.6,  1.48,  2.5,   3.71,  3.22,
      2.7, 2.1, 1.47,  1.05, 0.73, 0.52, 0.391, 0.288, 0.219, 0.171},

  //  S+ K+ pi0 (S- K0 pi0)
   {  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.01,  0.13,
      0.15, 0.12, 0.09, 0.06, 0.042, 0.03, 0.023, 0.016, 0.013, 0.011},

  //  S+ K0 pi+
   {  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.01,  0.13,
      0.15, 0.12, 0.09, 0.065, 0.043, 0.03, 0.023, 0.016, 0.013, 0.011},

  //  S0 K+ pi+
   {  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.01,  0.13,
      0.15, 0.12, 0.09, 0.065, 0.043, 0.03, 0.023, 0.016, 0.013, 0.011},

  //  L K+ pi+
   {  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.01,  0.13,
      0.15, 0.12, 0.09, 0.065, 0.043, 0.03, 0.023, 0.016, 0.013, 0.011},

  //  p K+ K0bar
   {  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.04,
      0.07, 0.062, 0.045, 0.035, 0.026, 0.021, 0.018, 0.015, 0.013, 0.011},
  //
  // multiplicity 4 (15 channels)
  //
  //  p pi+ pi+ pi-
   {  0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,  0.0, 0.0,  0.0,  0.05, 0.206, 0.736, 2.26,  3.43,
      3.5, 3.05, 2.6, 2.25, 1.88, 1.6,  1.442, 1.23,  1.073, 0.936},

  //  p pi+ pi0 pi0
   {  0.0,   0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,   0.0, 0.0,   0.0,   0.0,  0.075, 0.265, 0.575, 1.045, 1.835,
      2.125, 1.7, 1.135, 0.675, 0.43, 0.27,  0.185, 0.125, 0.092, 0.072},

  //  n pi+ pi+ pi0
   {  0.0,   0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,   0.0, 0.0,   0.0,   0.0,  0.075, 0.265, 0.575, 1.045, 1.835,
      2.125, 1.7, 1.135, 0.675, 0.43, 0.27,  0.185, 0.125, 0.092, 0.072},

  //  S+ K+ pi+ pi-
   {  0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.03, 0.08, 0.09, 0.06, 0.037, 0.025, 0.016, 0.011, 0.008},

  //  S+ K+ pi0 pi0
   {  0.0, 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.03, 0.08, 0.09, 0.055, 0.035, 0.025, 0.016, 0.011, 0.008},

  //  S+ K0 pi+ pi0
   {  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.02, 0.06, 0.07, 0.053, 0.04, 0.031, 0.025, 0.018, 0.014, 0.011},

  //  S0 K0 pi+ pi+
   {  0.0, 0.0,   0.0,   0.0,  0.0,   0.0,    0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,  0.0,   0.0,    0.0,   0.0,   0.0,   0.0,
      0.0, 0.007, 0.018, 0.02, 0.018, 0.0155, 0.014, 0.012, 0.011, 0.01},

  //  S0 K+ pi+ pi0
   {  0.0, 0.0,   0.0,   0.0,  0.0,   0.0,    0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,  0.0,   0.0,    0.0,   0.0,   0.0,   0.0,
      0.0, 0.007, 0.018, 0.02, 0.018, 0.0155, 0.014, 0.012, 0.011, 0.01},

  //  L K+ pi+ pi0
   {  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.013, 0.04, 0.04, 0.031, 0.025, 0.021, 0.017, 0.014, 0.011},

  //  L K0 pi+ pi+
   {  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.013, 0.04, 0.04, 0.031, 0.025, 0.021, 0.017, 0.014, 0.011},

  //  S- K+ pi+ pi+
   {  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,
      0.0,  0.0,  0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,
      0.02, 0.16, 0.18, 0.14, 0.1, 0.065, 0.05, 0.04, 0.033, 0.026},

  //  p pi+ K+ K-
   {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.02, 0.13, 0.16, 0.12, 0.09, 0.073, 0.063, 0.056, 0.049, 0.043},

  //  p pi+ K0 K0bar
   {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,  0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.01, 0.06, 0.15, 0.14, 0.11, 0.085, 0.067, 0.053, 0.049, 0.046},

  //  p pi0 K+ K0bar
   {  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.01, 0.05, 0.09, 0.085, 0.08, 0.07, 0.064, 0.053, 0.049, 0.046},

  //  n pi+ K+ K0bar
   {  0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,  0.0,  0.0,  0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.01, 0.05, 0.07, 0.065, 0.06, 0.05, 0.043, 0.037, 0.033, 0.029},
  //
  // multiplicity 5 (24 channels)
  //
  //  p pi+ pi+ pi- pi0
   {  0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,   0.03,  0.31,  2.14,
      3.45, 4.0, 3.5, 2.75, 2.15, 1.73, 1.442, 1.175, 0.963, 0.8},

  //  p pi+ pi0 pi0 pi0
   {  0.0,   0.0,  0.0, 0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,  0.0, 0.0,  0.0,   0.0,  0.0,   0.004, 0.03,  0.2,
      0.391, 0.45, 0.4, 0.31, 0.234, 0.19, 0.154, 0.118, 0.098, 0.08},

  //  n pi+ pi+ pi+ pi-
   {  0.0, 0.0,   0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,  0.0,  0.0,  0.002, 0.015, 0.06,  0.2,
      0.4, 0.715, 0.99, 0.93, 0.75, 0.58, 0.453, 0.342, 0.262, 0.205},

  //  n pi+ pi+ pi0 pi0
   {  0.0, 0.0,   0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,  0.0,  0.0,  0.002, 0.015, 0.06,  0.2,
      0.4, 0.715, 0.99, 0.93, 0.75, 0.58, 0.453, 0.342, 0.262, 0.205},

  //  S+ K+ pi+ pi- pi0
   {  0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.008, 0.05, 0.09, 0.09, 0.075, 0.062, 0.053, 0.046, 0.04},

  //  S+ K+ pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.006, 0.025, 0.04, 0.04, 0.035, 0.031, 0.027, 0.023, 0.019},

  //  S+ K0 pi+ pi+ pi-
   {  0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.005, 0.025, 0.04, 0.04, 0.035, 0.031, 0.027, 0.023, 0.019},

  //  S+ K0 pi+ pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.005, 0.025, 0.04, 0.04, 0.035, 0.031, 0.027, 0.023, 0.019},

  //  L K0 pi+ pi+ pi0
   {  0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0, 0.006, 0.035, 0.06, 0.06, 0.055, 0.049, 0.045, 0.039, 0.034},

  //  L K+ pi+ pi+ pi-
   {  0.0, 0.0, 0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0, 0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
      0, 0.004, 0.02, 0.045, 0.05, 0.045, 0.041, 0.037, 0.033, 0.029},

  //  L K+ pi+ pi0 pi0
   {  0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.004, 0.02, 0.04, 0.05, 0.045, 0.041, 0.037, 0.033, 0.029},

  //  S0 K+ pi+ pi+ pi-
   {  0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.01, 0.02, 0.019, 0.017, 0.016, 0.013, 0.011},

  //  S0 K+ pi+ pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.01, 0.02, 0.019, 0.017, 0.016, 0.013, 0.011},

  //  S0 K0 pi+ pi+ pi0
   {  0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.01, 0.02, 0.019, 0.017, 0.016, 0.013, 0.011},

  //  S- K+ pi+ pi+ pi0
   {  0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.012, 0.03, 0.027, 0.023, 0.02, 0.016, 0.013, 0.011},

  //  S- K0 pi+ pi+ pi+
   {  0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.012, 0.03, 0.027, 0.023, 0.02, 0.016, 0.013, 0.011},

  //  p pi+ pi- K+ K0bar
   {  0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.015, 0.06, 0.08, 0.075, 0.067, 0.061, 0.055, 0.051},

  //  p pi+ pi+ K0 K-
   {  0.0, 0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.015, 0.04, 0.05, 0.05, 0.046, 0.043, 0.038, 0.034},

  //  p pi+ pi0 K+ K-
   {  0.0, 0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,  0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.02, 0.05, 0.06, 0.055, 0.05, 0.046, 0.041, 0.036, 0.034},
 
  //  p pi+ pi0 K0 K0bar
   {  0.0, 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.01, 0.039, 0.05, 0.055, 0.05, 0.046, 0.041, 0.036, 0.034},

  //  p pi0 pi0 K+ K0bar
   {  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.01, 0.03, 0.04, 0.04, 0.038, 0.032, 0.027, 0.023},

  //  n pi+ pi+ K+ K-
   {  0.0, 0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.02, 0.045, 0.055, 0.05, 0.045, 0.041, 0.037, 0.033, 0.029},

  //  n pi+ pi+ K0 K0bar
   {  0.0, 0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.01, 0.035, 0.035, 0.03, 0.026, 0.024, 0.02, 0.017, 0.015},

  //  n pi+ pi0 K+ K0bar
   {  0.0, 0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.01, 0.035, 0.035, 0.03, 0.026, 0.024, 0.02, 0.017, 0.015},
  //
  // multiplicity 6 (33 channels)
  //
  //  p pi+ pi+ pi+ pi- pi-
   {  0.0,   0.0,   0.0,   0.0,  0.0, 0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,  0.0, 0.0,  0.0,   0.0,   0.0,   0.02,
      0.091, 0.275, 0.397, 0.48, 0.5, 0.44, 0.411, 0.363, 0.317, 0.274},

  //  p pi+ pi+ pi- pi0 pi0
   {  0.0,   0.0,   0.0,   0.0,  0.0, 0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,  0.0, 0.0,  0.0,   0.0,   0.0,   0.02,
      0.091, 0.275, 0.397, 0.48, 0.5, 0.44, 0.411, 0.363, 0.317, 0.274},

  //  p pi+ pi0 pi0 pi0 pi0
   {  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.01, 0.04, 0.064, 0.085, 0.09, 0.085, 0.082, 0.075, 0.066, 0.057},

  //  n pi+ pi+ pi0 pi0 pi0
   {  0.0,   0.0,   0.0,   0.0,  0.0, 0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,  0.0, 0.0,  0.0,   0.0,   0.0,   0.02,
      0.091, 0.275, 0.397, 0.48, 0.5, 0.44, 0.411, 0.363, 0.317, 0.274},

  //  n pi+ pi+ pi+ pi- pi0
   {  0.0,   0.0,   0.0,   0.0,  0.0, 0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,  0.0, 0.0,  0.0,   0.0,   0.0,   0.02,
      0.091, 0.275, 0.397, 0.48, 0.5, 0.44, 0.411, 0.363, 0.317, 0.274},

  //  L K+ pi+ pi+ pi- pi0
   {  0.0,   0.0,   0.0,  0.0,  0.0,  0.0, 0.0,   0.0,   0.0,  0.0,
      0.0,   0.0,   0.0,  0.0,  0.0,  0.0, 0.0,   0.0,   0.0,  0.0,
      0.001, 0.008, 0.03, 0.06, 0.08, 0.1, 0.103, 0.096, 0.08, 0.068},

  //  L K+ pi+ pi0 pi0 pi0
   {  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.003, 0.01, 0.02, 0.027, 0.033, 0.034, 0.032, 0.026, 0.023},

  //  L K0 pi+ pi+ pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.004, 0.015, 0.03, 0.04, 0.05, 0.051, 0.048, 0.039, 0.034},

  //  L K0 pi+ pi+ pi+ pi-
   {  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.003, 0.01, 0.02, 0.027, 0.033, 0.034, 0.032, 0.026, 0.023},

  //  S+ K+ pi+ pi+ pi- pi-
   {  0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.001, 0.006, 0.021, 0.042, 0.06, 0.07, 0.072, 0.067, 0.056, 0.048},

  //  S+ K+ pi+ pi- pi0 pi0
   {  0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.002, 0.012, 0.042, 0.084, 0.12, 0.14, 0.144, 0.135, 0.112, 0.096},

  //  S+ K+ pi0 pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.005, 0.006, 0.006, 0.006, 0.005, 0.004, 0.003},

  //  S+ K0 pi+ pi+ pi- pi0
   {  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.003, 0.01, 0.02, 0.028, 0.033, 0.034, 0.032, 0.026, 0.023},

  //  S+ K0 pi+ pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.001, 0.003, 0.007, 0.01, 0.011, 0.011, 0.011, 0.01, 0.008},

  //  S0 K+ pi+ pi+ pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.001, 0.003, 0.006, 0.009, 0.011, 0.011, 0.011, 0.01, 0.009},

  //  S0 K+ pi+ pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004, 0.003},

  //  S0 K0 pi+ pi+ pi+ pi-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.002, 0.003, 0.005, 0.006, 0.006, 0.006, 0.005, 0.005},

  //  S0 K0 pi+ pi+ pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.002, 0.004, 0.007, 0.008, 0.008, 0.007, 0.007, 0.006},

  //  S- K+ pi+ pi+ pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.004, 0.009, 0.014, 0.015, 0.015, 0.014, 0.012, 0.01},

  //  S- K+ pi+ pi+ pi+ pi-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,
      0.0, 0.001, 0.003, 0.006, 0.008, 0.01, 0.01, 0.01, 0.009, 0.008},

  //  S- K0 pi+ pi+ pi+ pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,
      0.0, 0.001, 0.003, 0.006, 0.008, 0.01, 0.01, 0.01, 0.008, 0.007},

  //  p K+ K0b pi+ pi- pi0
   {  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.003, 0.01, 0.03, 0.075, 0.14, 0.165, 0.156, 0.146, 0.137},

  //  p K+ K0b pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.006, 0.014, 0.023, 0.028, 0.027, 0.024, 0.023},

  //  p K+ K- pi+ pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.005, 0.015, 0.037, 0.07, 0.082, 0.078, 0.073, 0.068},

  //  p K+ K- pi+ pi+ pi-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.005, 0.015, 0.037, 0.07, 0.082, 0.078, 0.073, 0.068},

  //  p K0 K0b pi+ pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.009, 0.02, 0.035, 0.041, 0.038, 0.036, 0.034},

  //  p K0 K0b pi+ pi+ pi-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.009, 0.02, 0.035, 0.041, 0.038, 0.036, 0.034},

  //  p K0 K- pi+ pi+ pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.005, 0.015, 0.037, 0.07, 0.082, 0.078, 0.073, 0.068},

  //  n K+ K0b pi+ pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.005, 0.02, 0.045, 0.07, 0.082, 0.078, 0.073, 0.068},

  //  n K+ K0b pi+ pi+ pi-
   {  0.0, 0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.005, 0.02, 0.045, 0.07, 0.082, 0.078, 0.073, 0.068},

  //  n K+ K- pi+ p1i+ pi0
   {  0.0, 0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.005, 0.02, 0.045, 0.07, 0.082, 0.078, 0.073, 0.068},

  //  n K0 K0b pi+ pi+ pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.004, 0.015, 0.035, 0.052, 0.062, 0.059, 0.055, 0.051},

  //  n K0 K- pi+ pi+ pi+
   {  0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.004, 0.01, 0.018, 0.021, 0.019, 0.018, 0.017},
  //
  // multiplicity 7 (41 channels)
  //
  //  p pi+ pi+ pi+ pi- pi- pi0
   {  0.0,   0.0,  0.0,  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,  0.0,  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
      0.041, 0.22, 0.55, 0.8, 0.9, 0.9, 0.927, 0.961, 0.986, 1.028},

  //  p pi+ pi+ pi- pi0 pi0 pi0
   {  0.0,   0.0,  0.0,  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,  0.0,  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,
      0.041, 0.22, 0.55, 0.8, 0.9, 0.9, 0.927, 0.961, 0.986, 1.028},

  //  p pi+ pi0 pi0 pi0 pi0 pi0
   {  0.0,   0.0,   0.0,  0.0,  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,  0.0,  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
      0.004, 0.025, 0.06, 0.09, 0.1, 0.095, 0.093, 0.096, 0.098, 0.103},

  //  n pi+ pi+ pi+ pi+ pi- pi-
   {  0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0,   0.0,  0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0,   0.0,  0.0,   0.0,
      0.004, 0.024, 0.075, 0.15, 0.25, 0.3, 0.309, 0.32, 0.328, 0.342},

  //  n pi+ pi+ pi+ pi- pi0 pi0
   {  0.0,  0.0,   0.0,  0.0,  0.0,  0.0, 0.0,   0.0,   0.0,   0.0,
      0.0,  0.0,   0.0,  0.0,  0.0,  0.0, 0.0,   0.0,   0.0,   0.0,
      0.02, 0.089, 0.24, 0.45, 0.75, 0.9, 0.927, 0.961, 0.986, 1.028},

  //  n pi+ pi+ pi0 pi0 pi0 pi0
   {  0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0,   0.0,  0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0,   0.0,  0.0,   0.0,
      0.004, 0.023, 0.075, 0.15, 0.25, 0.3, 0.309, 0.32, 0.328, 0.342},

  //  L K+ pi+ pi+ pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,   0.0,
      0.0, 0.002, 0.006, 0.014, 0.027, 0.039, 0.04, 0.04, 0.038, 0.038},

  //  L K+ pi+ pi+ pi+ pi- pi-
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.007, 0.011, 0.013, 0.014, 0.014, 0.013, 0.013},

  //  L K+ pi+ pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.007, 0.007, 0.008, 0.007},

  //  L K0 pi+ pi+ pi+ pi- pi0
   {  0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
      0.0, 0.002, 0.007, 0.02, 0.04, 0.07, 0.082, 0.08, 0.073, 0.068},

  //  L K0 pi+ pi+ pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.01, 0.023, 0.035, 0.041, 0.038, 0.036, 0.034},

  //  S+ K+ pi+ pi+ pi- pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.007, 0.012, 0.017, 0.017, 0.017, 0.016, 0.016},

  //  S+ K+ pi+ pi- pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.002, 0.004, 0.008, 0.011, 0.011, 0.012, 0.011, 0.01},

  //  S+ K+ pi0 pi0 pi0 pi0 pi0   negligible

  //  S+ K0 pi+ pi+ pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,
      0.0, 0.001, 0.002, 0.004, 0.007, 0.01, 0.01, 0.01, 0.009, 0.009},

  //  S+ K0 pi+ pi+ pi+ pi- pi-
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.005},

  //  S+ K0 pi+ pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002},

  //  S0 K+ pi+ pi+ pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.006, 0.009, 0.013, 0.013, 0.013, 0.013, 0.013},

  //  S0 K+ pi+ pi+ pi+ pi- pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.004, 0.004, 0.005},

  //  S0 K+ pi+ pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003, 0.003},

  //  S0 K0 pi+ pi+ pi+ pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.003, 0.007, 0.014, 0.023, 0.028, 0.027, 0.026, 0.024},

  //  S0 K0 pi+ pi+ pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.008, 0.012, 0.013, 0.013, 0.012, 0.011},

  //  S- K+ pi+ pi+ pi+ pi- pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.010, 0.012, 0.012, 0.013},

  //  S- K+ pi+ pi+ pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.006, 0.007, 0.008, 0.008},

  //  S- K0 pi+ pi+ pi+ pi+ pi-
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.005, 0.006, 0.007, 0.007},

  //  S- K0 pi+ pi+ pi+ pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.002, 0.004, 0.006, 0.007, 0.008, 0.009, 0.009},

  //  p K+ K0b pi+ pi- pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,  0.0,
      0.0, 0.001, 0.004, 0.015, 0.037, 0.05, 0.05, 0.05, 0.05, 0.051},

  //  p K+ K0b pi+ pi+ pi- pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.007, 0.018, 0.025, 0.026, 0.026, 0.025, 0.025},

  //  p K+ K0b pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.005, 0.005, 0.005, 0.006},

  //  p K+ K- pi+ pi+ pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.004, 0.017, 0.043, 0.058, 0.059, 0.059, 0.058, 0.058},

  //  p K+ K- pi+ pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
      0.0, 0.0, 0.002, 0.006, 0.014, 0.019, 0.02, 0.02, 0.02, 0.019},

  //  p K0 K0b pi+ pi+ pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.006, 0.025, 0.065, 0.088, 0.088, 0.088, 0.087, 0.086},

  //  p K0 K0b pi+ pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
      0.0, 0.0, 0.002, 0.008, 0.022, 0.029, 0.029, 0.03, 0.03, 0.03},

  //  p K0 K- pi+ pi+ pi0 pi0
   {  0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.002, 0.01, 0.033, 0.087, 0.117, 0.116, 0.118, 0.117, 0.114},

  //  p K0 K- pi+ pi+ pi+ pi-
   {  0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.003, 0.015, 0.05, 0.13, 0.175, 0.175, 0.176, 0.175, 0.172},

  //  n K+ K0b pi+ pi+ pi- pi0
   {  0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.001, 0.005, 0.018, 0.043, 0.058, 0.059, 0.059, 0.058, 0.058},

  //  n K+ K0b pi+ pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
      0.0, 0.0, 0.002, 0.006, 0.014, 0.019, 0.02, 0.02, 0.02, 0.019},

  //  n K+ K- pi+ pi+ pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
      0.0, 0.0, 0.002, 0.008, 0.022, 0.029, 0.029, 0.03, 0.03, 0.029},

  //  n K+ K- pi+ pi+ pi+ pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
      0.0, 0.0, 0.002, 0.006, 0.014, 0.019, 0.02, 0.02, 0.02, 0.019},

  //  n K0 K0b pi+ pi+ pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,
      0.0, 0.0, 0.002, 0.008, 0.022, 0.029, 0.029, 0.03, 0.03, 0.029},

  //  n K0 K0b pi+ pi+ pi+ pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
      0.0, 0.0, 0.002, 0.006, 0.014, 0.019, 0.02, 0.02, 0.02, 0.019},

  //  n K0 K- pi+ pi+ pi+ pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
      0.0, 0.0, 0.002, 0.006, 0.014, 0.019, 0.02, 0.02, 0.02, 0.019},
  //
  // multiplicity 8 (47 channels)
  //
  //  p pi+ pi+ pi+ pi+ pi- pi- pi-
   {  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.01, 0.028, 0.045, 0.07, 0.085, 0.095, 0.103, 0.112, 0.12, 0.131},

  //  p pi+ pi+ pi+ pi- pi- pi0 pi0
   {  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.01, 0.05, 0.15, 0.3, 0.45, 0.57, 0.618, 0.673, 0.723, 0.787},

  //  p pi+ pi+ pi- pi0 pi0 pi0 pi0
   {  0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.005, 0.025, 0.075, 0.15, 0.22, 0.28, 0.309, 0.337, 0.361, 0.394},

  //  p pi+ pi0 pi0 pi0 pi0 pi0 pi0
   {  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,
      0.001, 0.005, 0.01, 0.015, 0.018, 0.02, 0.021, 0.022, 0.024, 0.026},

  //  n pi+ pi+ pi+ pi+ pi- pi- pi0
   {  0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
      0.0,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,
      0.01, 0.06, 0.12, 0.2, 0.27, 0.33, 0.36, 0.395, 0.426, 0.468},

  //  n pi+ pi+ pi+ pi- pi0 pi0 pi0
   {  0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,  0.0,
      0.0,   0.0,  0.0,  0.0,   0.0,  0.0,  0.0,  0.0,   0.0,  0.0,
      0.013, 0.07, 0.16, 0.254, 0.36, 0.44, 0.48, 0.527, 0.57, 0.624},

  //  n pi+ pi+ pi0 pi0 pi0 pi0 pi0
   {  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.002, 0.014, 0.024, 0.04, 0.052, 0.066, 0.072, 0.079, 0.085, 0.094},

  //  L K+ pi+ pi+ pi+ pi- pi- pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.006, 0.014, 0.026, 0.041, 0.053, 0.059, 0.068},

  //  L K+ pi+ pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.006, 0.007},

  //  L K0 pi+ pi+ pi+ pi- pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0, 0.001, 0.003, 0.008, 0.017, 0.027, 0.036, 0.04, 0.046},

  //  L K0 pi+ pi+ pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.008, 0.015, 0.021, 0.027, 0.031, 0.034},

  //  L K0 pi+ pi+ pi+ pi+ pi- pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.008, 0.015, 0.021, 0.027, 0.031, 0.034},

  //  S+ K+ pi+ pi+ pi- pi- pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.005, 0.009, 0.012, 0.016, 0.017, 0.019},

  //  S+ K+ pi+ pi- pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.006, 0.007},

  //  S+ K+ pi+ pi+ pi+ pi- pi- pi-
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.006, 0.007, 0.008, 0.009},

  //  S+ K+ pi0 pi0 pi0 pi0 pi0 pi0 (negligible)

  //  S+ K0 pi+ pi+ pi+ pi- pi- pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.01, 0.011, 0.012},

  //  S+ K0 pi+ pi+ pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.01, 0.011, 0.012},

  //  S+ K0 pi+ pi0 pi0 pi0 pi0 pi0 (negligible)

  //  S0 K+ pi+ pi+ pi+ pi- pi- pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.013, 0.018, 0.021, 0.023},

  //  S0 K+ pi+ pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002},

  //  S0 K0 pi+ pi+ pi+ pi- pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.012, 0.013, 0.015},

  //  S0 K0 pi+ pi+ pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.009, 0.01, 0.011},

  //  S0 K0 pi+ pi+ pi+ pi+ pi- pi-
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.009, 0.01, 0.011},

  //  S- K+ pi+ pi+ pi+ pi- pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.012, 0.016, 0.017, 0.019},

  //  S- K+ pi+ pi+ pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.005},

  //  S- K+ pi+ pi+ pi+ pi+ pi- pi-
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.005},

  //  S- K0 pi+ pi+ pi+ pi+ pi- pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.005, 0.006, 0.007, 0.008, 0.009},

  //  S- K0 pi+ pi+ pi+ pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.005, 0.006, 0.007, 0.008},

  //  p K+ K0b pi+ pi+ pi- pi- pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.007, 0.019, 0.036, 0.051, 0.064, 0.077, 0.091},

  //  p K+ K0b pi+ pi- pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.004, 0.011, 0.022, 0.034, 0.043, 0.051, 0.06},

  //  p K+ K0b pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002},

  //  p K+ K- pi+ pi+ pi- pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.007, 0.019, 0.036, 0.051, 0.064, 0.077, 0.091},

  //  p K+ K- pi+ pi+ pi+ pi- pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.007, 0.012, 0.017, 0.024, 0.027, 0.031},

  //  p K+ K- pi+ pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.011, 0.013, 0.015},

  //  p K0 K0b pi+ pi+ pi- pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.002, 0.007, 0.019, 0.036, 0.051, 0.064, 0.077, 0.091},

  //  p K0 K0b pi+ pi+ pi+ pi- pi-
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.007, 0.012, 0.017, 0.024, 0.027, 0.031},

  //  p K0 K0b pi+ pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.011, 0.013, 0.015},

  //  p K0 K- pi+ pi+ pi+ pi- pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.004, 0.01, 0.023, 0.034, 0.043, 0.051, 0.06},

  //  p K0 K- pi+ pi+ pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.007, 0.012, 0.017, 0.024, 0.027, 0.031},

  //  n K+ K0b pi+ pi+ pi+ pi- pi-
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.006, 0.016, 0.022, 0.03, 0.035, 0.037},

  //  n K+ K0b pi+ pi+ pi- pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0, 0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,  0.0,   0.0,   0.0,  0.0, 0.0,
      0.0, 0.0, 0.0, 0.004, 0.02, 0.044, 0.065, 0.09, 0.1, 0.11},

  //  n K+ K0b pi+ pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.003, 0.008, 0.01, 0.015, 0.017, 0.018},

  //  n K+ K- pi+ pi+ pi+ pi- pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.001, 0.003, 0.007, 0.012, 0.017, 0.022, 0.026, 0.031},

  //  n K+ K- pi+ pi+ pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.011, 0.013, 0.015},

  //  n K0 K0b pi+ pi+ pi+ pi- pi0
   {  0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.002, 0.007, 0.016, 0.021, 0.026, 0.031, 0.037},

  //  n K0 K0b pi+ pi+ pi0 pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.003, 0.008, 0.01, 0.013, 0.015, 0.018},

  //  n K0 K- pi+ pi+ pi+ pi0 pi0
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.003, 0.008, 0.01, 0.013, 0.015, 0.018},

  //  n K0 K- pi+ pi+ pi+ pi+ pi-
   {  0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0, 0.0, 0.0, 0.001, 0.004, 0.005, 0.006, 0.008, 0.009},
  //
  // multiplicity 9 (55 channels)
  //
  //  p pi+ pi+ pi+ pi+ pi- pi- pi- pi0
   {  0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,
      0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0,   0.0,  0.0,  0.0,
      0.003, 0.013, 0.037, 0.08, 0.14, 0.2, 0.257, 0.31, 0.35, 0.399},

  //  p pi+ pi+ pi+ pi- pi- pi0 pi0 pi0
   {  0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0,   0.0,  0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,  0.0,  0.0, 0.0,   0.0,  0.0,   0.0,
      0.006, 0.032, 0.085, 0.16, 0.27, 0.4, 0.515, 0.62, 0.701, 0.799},

  //  p pi+ pi+ pi- pi0 pi0 pi0 pi0 pi0
   {  0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,  0.0,
      0.0,   0.0,   0.0,   0.0,  0.0,   0.0,  0.0,   0.0,   0.0,  0.0,
      0.002, 0.011, 0.027, 0.05, 0.085, 0.12, 0.154, 0.186, 0.21, 0.24},

  //  p pi+ pi0 pi0 pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
      0.0, 0.0, 0.001, 0.003, 0.005, 0.007, 0.008, 0.009, 0.01, 0.011},

  //  n pi+ pi+ pi+ pi+ pi+ pi- pi- pi-
   {  0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,  0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.003, 0.01, 0.02, 0.036, 0.048, 0.066, 0.082, 0.093, 0.103},

  //  n pi+ pi+ pi+ pi+ pi- pi- pi0 pi0
   {  0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.003, 0.022, 0.075, 0.15, 0.24, 0.36, 0.495, 0.617, 0.698, 0.77},

  //  n pi+ pi+ pi+ pi- pi0 pi0 pi0 pi0
   {  0.0,   0.0,   0.0,  0.0, 0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.0,   0.0,   0.0,  0.0, 0.0,  0.0,  0.0,   0.0,   0.0,   0.0,
      0.002, 0.015, 0.05, 0.1, 0.17, 0.25, 0.329, 0.412, 0.465, 0.513},

  //  n pi+ pi+ pi0 pi0 pi0 pi0 pi0 pi0
   {  0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
      0.0, 0.002, 0.005, 0.01, 0.017, 0.025, 0.033, 0.041, 0.046, 0.051},

  //  L K+ pi+ pi+ pi+ pi- pi- pi0 pi0
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.002, 0.006, 0.013, 0.026, 0.041, 0.053, 0.059, 0.068},

  //  L K+ pi+ pi+ pi- pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.003, 0.007, 0.014, 0.021, 0.027, 0.03,  0.034},

  //  L K+ pi+ pi+ pi+ pi+ pi- pi- pi-
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.008, 0.01, 0.011, 0.012},

  //  L K+ pi+ pi0 pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002},

  //  L K0 pi+ pi+ pi+ pi- pi0 pi0 pi0
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.001, 0.004, 0.009, 0.018, 0.028, 0.035, 0.039, 0.046},

  //  L K0 pi+ pi+ pi+ pi+ pi- pi- pi0
   { 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.001, 0.003, 0.007, 0.014, 0.021, 0.027, 0.03, 0.034},

  //  L K0 pi+ pi+ pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.003, 0.004, 0.005, 0.006, 0.007},

  // S+ K+ pi+ pi+ pi- pi- pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.005, 0.009, 0.013, 0.018, 0.021, 0.023},

  // S+ K+ pi+ pi+ pi+ pi- pi- pi- pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.013, 0.014, 0.015},

  // S+ K+ pi+ pi- pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.004, 0.004, 0.005},

  // S+ K+ pi0 pi0 pi0 pi0 pi0 pi0 pi0 (negligible)

  // S+ K0 pi+ pi+ pi+ pi- pi- pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.005, 0.009, 0.013, 0.018, 0.021, 0.023},

  // S+ K0 pi+ pi+ pi- pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.006, 0.009, 0.01, 0.011},

  // S+ K0 pi+ pi+ pi+ pi+ pi- pi- pi-
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003},

  // S+ K0 pi+ pi0 pi0 pi0 pi0 pi0 pi0 (negligible)

  // S0 K+ pi+ pi+ pi+ pi- pi- pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.005, 0.009, 0.013, 0.018, 0.021, 0.023},

  // S0 K+ pi+ pi+ pi- pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.006, 0.009, 0.01, 0.011},

  // S0 K+ pi+ pi+ pi+ pi+ pi- pi- pi-
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.003, 0.003, 0.003},

  // S0 K+ pi+ pi0 pi0 pi0 pi0 pi0 pi0 (negligible)

  // S0 K0 pi+ pi+ pi+ pi- pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.012, 0.013, 0.015},

  // S0 K0 pi+ pi+ pi+ pi+ pi- pi- pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.006, 0.009, 0.01, 0.011},

  // S0 K0 pi+ pi+ pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002},

  // S- K+ pi+ pi+ pi+ pi- pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.003, 0.006, 0.009, 0.013, 0.014, 0.015},

  // S- K+ pi+ pi+ pi+ pi+ pi- pi- pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.006, 0.01, 0.012, 0.013},

  // S- K+ pi+ pi+ pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.002, 0.002},

  // S- K0 pi+ pi+ pi+ pi+ pi- pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,  0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.002, 0.004, 0.006, 0.01, 0.012, 0.013},

  // S- K0 pi+ pi+ pi+ pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.001, 0.002, 0.002, 0.003, 0.003, 0.003},

  // S- K0 pi+ pi+ pi+ pi+ pi+ pi- pi-  (negligible)

  // p K+ K0b pi+ pi+ pi- pi- pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.005, 0.027, 0.056, 0.078, 0.098, 0.114, 0.132},

  // p K+ K0b pi+ pi- pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.003, 0.009, 0.018, 0.026, 0.033, 0.038, 0.044},

  // p K+ K0b pi+ pi+ pi+ pi- pi- pi-
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.006, 0.012, 0.018, 0.022, 0.026, 0.03},

  //  p K+ K0b pi0 pi0 pi0 pi0 pi0 pi0 (negligible)

  // p K+ K- pi+ pi+ pi- pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.004, 0.015, 0.035, 0.052, 0.065, 0.075, 0.088},

  // p K+ K- pi+ pi+ pi+ pi- pi- pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.004, 0.015, 0.035, 0.052, 0.065, 0.075, 0.088},

  // p K+ K- pi+ pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.004, 0.005, 0.006, 0.008, 0.009},

  // p K0 K0b pi+ pi+ pi- pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.003, 0.015, 0.035, 0.052, 0.065, 0.075, 0.088},

  // p K0 K0b pi+ pi+ pi+ pi- pi- pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.003, 0.015, 0.035, 0.052, 0.065, 0.075, 0.088},

  // p K0 K0b pi+ pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.004, 0.005, 0.006, 0.008, 0.009},

  // p K0 K- pi+ pi+ pi+ pi- pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.003, 0.015, 0.035, 0.052, 0.065, 0.075, 0.088},

  // p K0 K- pi+ pi+ pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.009, 0.013, 0.016, 0.019, 0.022},

  // p K0 K- pi+ pi+ pi+ pi+ pi- pi-
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.009, 0.013, 0.016, 0.019, 0.022},

  // n K+ K0b pi+ pi+ pi- pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.003, 0.015, 0.035, 0.052, 0.065, 0.075, 0.088},

  // n K+ K0b pi+ pi+ pi+ pi- pi- pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.003, 0.015, 0.035, 0.052, 0.065, 0.075, 0.088},

  // n K+ K0b pi+ pi0 pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0, 0.002, 0.004, 0.005, 0.006, 0.008, 0.009},

  // n K+ K- pi+ pi+ pi+ pi- pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.003, 0.015, 0.035, 0.052, 0.065, 0.075, 0.088},

  // n K+ K- pi+ pi+ pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.009, 0.013, 0.016, 0.019, 0.022},

  // n K+ K- pi+ pi+ pi+ pi+ pi- pi-
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.009, 0.013, 0.016, 0.019, 0.022},

  // n K0 K0b pi+ pi+ pi+ pi- pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.003, 0.015, 0.035, 0.052, 0.065, 0.075, 0.088},

  // n K0 K0b pi+ pi+ pi0 pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.009, 0.013, 0.016, 0.019, 0.022},

  // n K0 K0b pi+ pi+ pi+ pi+ pi- pi-
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.004, 0.009, 0.013, 0.016, 0.019, 0.022},

  // n K0 K- pi+ pi+ pi+ pi+ pi- pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.002, 0.008, 0.018, 0.026, 0.033, 0.038, 0.044},

  // n K0 K- pi+ pi+ pi+ pi0 pi0 pi0
   { 0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
     0.0, 0.0, 0.0, 0.001, 0.006, 0.012, 0.019, 0.022, 0.026, 0.03}};
}

// Initialize both |T Tz> = |3/2 3/2> channels, using pipP cross-section table

const G4CascadePiPlusPChannelData::data_t
G4CascadePiPlusPChannelData::data(pipP2bfs, pipP3bfs, pipP4bfs,
				  pipP5bfs, pipP6bfs, pipP7bfs,
				  pipP8bfs, pipP9bfs, pipPCrossSections,
				  pipPtotXSec, pip*pro, "PiPlusP");

const G4CascadePiMinusNChannelData::data_t
G4CascadePiMinusNChannelData::data(pimN2bfs, pimN3bfs, pimN4bfs,
				   pimN5bfs, pimN6bfs, pimN7bfs,
				   pimN8bfs, pimN9bfs, pipPCrossSections,
				   pipPtotXSec, pim*neu, "PiMinusN");
