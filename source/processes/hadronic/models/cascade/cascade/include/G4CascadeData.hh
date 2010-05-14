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
// $Id: G4CascadeData.hh,v 1.2 2010-05-14 18:28:02 mkelsey Exp $
// GEANT4 tag: $Name: not supported by cvs2svn $
//
// 20100507  M. Kelsey -- Use template arguments to dimension const-refs
//		to arrays,for use in passing to functions as dimensioned.
//		Reduce offset array to one dimension.

#ifndef G4_CASCADE_DATA_HH
#define G4_CASCADE_DATA_HH

#include "globals.hh"
#include "G4CascadeSampler.hh"		/* To get number of energy bins */

template <int N2, int N3, int N4, int N5, int N6, int N7>
struct G4CascadeData
{
  // NOTE: Need access to N2 by value to initialize index array
  enum { N02=N2, N23=N2+N3, N24=N23+N4, N25=N24+N5, N26=N25+N6, N27=N26+N7 };
  enum { NM=6, NE=G4CascadeSampler::energyBins, NXS=N27 };

  G4double tot[NE];			// Summed cross-sections, to be filled
  G4double multiplicities[NM][NE];	// Multiplicity distributions
  static G4int index[NM+1];		// Start and stop indices to xsec's

  const G4int (&x2bfs)[N2][2];		// Initialized from file-scope inputs
  const G4int (&x3bfs)[N3][3];
  const G4int (&x4bfs)[N4][4];
  const G4int (&x5bfs)[N5][5];
  const G4int (&x6bfs)[N6][6];
  const G4int (&x7bfs)[N7][7];
  const G4double (&crossSections)[NXS][31];

  G4CascadeData(const G4int (&the2bfs)[N2][2], const G4int (&the3bfs)[N3][3],
		const G4int (&the4bfs)[N4][4], const G4int (&the5bfs)[N5][5],
		const G4int (&the6bfs)[N6][6], const G4int (&the7bfs)[N7][7],
		const G4double (&xsec)[NXS][31])
    : x2bfs(the2bfs), x3bfs(the3bfs), x4bfs(the4bfs), x5bfs(the5bfs),
      x6bfs(the6bfs), x7bfs(the7bfs), crossSections(xsec) { initialize(); }

  void initialize();			// Fill summed arrays from input
};

// Index ranges for crossSections table determined by template arguments
template <int N2, int N3, int N4, int N5, int N6, int N7>
G4int G4CascadeData<N2,N3,N4,N5,N6,N7>::index[NM+1] =
  { 0, 
    G4CascadeData<N2,N3,N4,N5,N6,N7>::N02, 
    G4CascadeData<N2,N3,N4,N5,N6,N7>::N23, 
    G4CascadeData<N2,N3,N4,N5,N6,N7>::N24, 
    G4CascadeData<N2,N3,N4,N5,N6,N7>::N25, 
    G4CascadeData<N2,N3,N4,N5,N6,N7>::N26,
    G4CascadeData<N2,N3,N4,N5,N6,N7>::N27 };


template <int N2, int N3, int N4, int N5, int N6, int N7> inline
void G4CascadeData<N2,N3,N4,N5,N6,N7>::initialize() {
  // Initialize multiplicity array
  
  for (G4int m = 0; m < NM; m++) {
    G4int start = index[m];
    G4int stop = index[m+1];
    for (G4int k = 0; k < NE; k++) {
      multiplicities[m][k] = 0.0;
      for (G4int i = start; i < stop; i++) {
 	multiplicities[m][k] += crossSections[i][k];
      }
    }
  }
  
  // Initialize total cross section array
  
  for (G4int k = 0; k < NE; k++) {
    tot[k] = 0.0;
    for (G4int m = 0; m < NM; m++) {
      tot[k] += multiplicities[m][k];
    }
  }
}

#endif
