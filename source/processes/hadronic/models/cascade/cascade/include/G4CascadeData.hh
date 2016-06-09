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
#ifndef G4_CASCADE_DATA_HH
#define G4_CASCADE_DATA_HH

#include "globals.hh"

template <int n2, int n3, int n4, int n5, int n6, int n7, int nxs>
struct G4CascadeData
{
  G4double *  tot;

  typedef G4double multiplicities_t[31];
  multiplicities_t * multiplicities;

  typedef G4int index_t[2];
  index_t const * index;

  typedef G4int x2bfs_t[2];
  x2bfs_t const * x2bfs;

  typedef G4int x3bfs_t[3];
  x3bfs_t const * x3bfs;

  typedef G4int x4bfs_t[4];
  x4bfs_t const * x4bfs;

  typedef G4int x5bfs_t[5];
  x5bfs_t const * x5bfs;

  typedef G4int x6bfs_t[6];
  x6bfs_t const * x6bfs;

  typedef G4int x7bfs_t[7];
  x7bfs_t const * x7bfs;

  typedef G4float crossSections_t[31];
  crossSections_t const * crossSections;

  void initialize();

//   G4double tot[31];
//   G4double multiplicities[6][31];

//   G4int index[6][2];
//   G4int x2bfs[n2][2];
//   G4int x3bfs[n3][3];
//   G4int x4bfs[n4][4];
//   G4int x5bfs[n5][5];
//   G4int x6bfs[n6][6];
//   G4int x7bfs[n7][7];

//   G4float crossSections[nxs][31];
};

template <int n2, int n3, int n4, int n5, int n6, int n7, int nxs>
inline
void
G4CascadeData<n2, n3, n4, n5, n6, n7, nxs>::initialize()
{
  // Initialize multiplicity array
  
  for (G4int m = 0; m < 6; m++) {
    G4int start = index[m][0];
    G4int stop = index[m][1];
    for (G4int k = 0; k < 31; k++) {
      multiplicities[m][k] = 0.0;
      for (G4int i = start; i < stop; i++) {
 	multiplicities[m][k] += crossSections[i][k];
      }
    }
  }
  
  // Initialize total cross section array
  
  for (G4int k = 0; k < 31; k++) {
    tot[k] = 0.0;
    for (G4int m = 0; m < 6; m++) {
      tot[k] += multiplicities[m][k];
    }
  }
}

#endif
