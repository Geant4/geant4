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
// $Id: G4CascadeData.hh,v 1.8 2010-08-03 23:09:36 mkelsey Exp $
// GEANT4 tag: $Name: not supported by cvs2svn $
//
// 20100507  M. Kelsey -- Use template arguments to dimension const-refs
//		to arrays,for use in passing to functions as dimensioned.
//		Add two additional optional(!) template args for piN/NN.
//		Add new data member "sum" to separate summed xsec values
//		from measured inclusive (tot) cross-sections.  Add two
//		ctors to pass inclusive xsec array as input (for piN/NN).
// 20100611  M. Kelsey -- Work around Intel ICC compiler warning about
//		index[] subscripts out of range.  Dimension to full [9].
// 20100803  M. Kelsey -- Add printing function for debugging

#ifndef G4_CASCADE_DATA_HH
#define G4_CASCADE_DATA_HH

#include "globals.hh"
#include "G4CascadeSampler.hh"		/* To get number of energy bins */
#include <iomanip>


template <int NE,int N2,int N3,int N4,int N5,int N6,int N7,int N8=0,int N9=0>
struct G4CascadeData
{
  // NOTE: Need access to N2 by value to initialize index array
  enum { N02=N2, N23=N2+N3, N24=N23+N4, N25=N24+N5, N26=N25+N6, N27=N26+N7,
	 N28=N27+N8, N29=N28+N9 };

  enum { N8D=N8?N8:1, N9D=N9?N9:1 };	// SPECIAL: Can't dimension arrays [0]

  enum { NM=N9?8:N8?7:6, NXS=N29 };	// Multiplicity and cross-section bins

  G4int index[9];			// Start and stop indices to xsec's
  G4double multiplicities[NM][NE];	// Multiplicity distributions

  const G4int (&x2bfs)[N2][2];		// Initialized from file-scope inputs
  const G4int (&x3bfs)[N3][3];
  const G4int (&x4bfs)[N4][4];
  const G4int (&x5bfs)[N5][5];
  const G4int (&x6bfs)[N6][6];
  const G4int (&x7bfs)[N7][7];
  const G4int (&x8bfs)[N8D][8];		// These may not be used if mult==7
  const G4int (&x9bfs)[N9D][9];
  const G4double (&crossSections)[NXS][NE];

  G4double sum[NE];			// Summed cross-sections, computed
  const G4double (&tot)[NE];		// Inclusive cross-sections (from input)

  static const G4int empty8bfs[1][8];	// For multiplicity==7 case
  static const G4int empty9bfs[1][9];

  G4int maxMultiplicity() const { return NM+1; }  // Used by G4CascadeFunctions

  void print(G4int mult=-1) const;	// Dump multiplicty tables (-1 == all)
  void printXsec(const G4double (&xsec)[NE]) const;

  // Constructor for kaon/hyperon channels, with multiplicity <= 7
  G4CascadeData(const G4int (&the2bfs)[N2][2], const G4int (&the3bfs)[N3][3],
		const G4int (&the4bfs)[N4][4], const G4int (&the5bfs)[N5][5],
		const G4int (&the6bfs)[N6][6], const G4int (&the7bfs)[N7][7],
		const G4double (&xsec)[NXS][NE])
    : x2bfs(the2bfs), x3bfs(the3bfs), x4bfs(the4bfs), x5bfs(the5bfs),
      x6bfs(the6bfs), x7bfs(the7bfs), x8bfs(empty8bfs), x9bfs(empty9bfs),
      crossSections(xsec), tot(sum) { initialize(); }

  // Constructor for kaon/hyperon channels, with multiplicity <= 7 and inclusive
  G4CascadeData(const G4int (&the2bfs)[N2][2], const G4int (&the3bfs)[N3][3],
		const G4int (&the4bfs)[N4][4], const G4int (&the5bfs)[N5][5],
		const G4int (&the6bfs)[N6][6], const G4int (&the7bfs)[N7][7],
		const G4double (&xsec)[NXS][NE], const G4double (&theTot)[NE])
    : x2bfs(the2bfs), x3bfs(the3bfs), x4bfs(the4bfs), x5bfs(the5bfs),
      x6bfs(the6bfs), x7bfs(the7bfs), x8bfs(empty8bfs), x9bfs(empty9bfs),
      crossSections(xsec), tot(theTot) { initialize(); }

  // Constructor for pion/nuleon channels, with multiplicity > 7
  G4CascadeData(const G4int (&the2bfs)[N2][2], const G4int (&the3bfs)[N3][3],
		const G4int (&the4bfs)[N4][4], const G4int (&the5bfs)[N5][5],
		const G4int (&the6bfs)[N6][6], const G4int (&the7bfs)[N7][7],
		const G4int (&the8bfs)[N8D][8], const G4int (&the9bfs)[N9D][9],
		const G4double (&xsec)[NXS][NE])
    : x2bfs(the2bfs), x3bfs(the3bfs), x4bfs(the4bfs), x5bfs(the5bfs),
      x6bfs(the6bfs), x7bfs(the7bfs), x8bfs(the8bfs), x9bfs(the9bfs),
      crossSections(xsec), tot(sum) { initialize(); }

  // Constructor for pion/nuleon channels, with multiplicity > 7 and inclusive
  G4CascadeData(const G4int (&the2bfs)[N2][2], const G4int (&the3bfs)[N3][3],
		const G4int (&the4bfs)[N4][4], const G4int (&the5bfs)[N5][5],
		const G4int (&the6bfs)[N6][6], const G4int (&the7bfs)[N7][7],
		const G4int (&the8bfs)[N8D][8], const G4int (&the9bfs)[N9D][9],
		const G4double (&xsec)[NXS][NE], const G4double (&theTot)[NE])
    : x2bfs(the2bfs), x3bfs(the3bfs), x4bfs(the4bfs), x5bfs(the5bfs),
      x6bfs(the6bfs), x7bfs(the7bfs), x8bfs(the8bfs), x9bfs(the9bfs),
      crossSections(xsec), tot(theTot) { initialize(); }
  void initialize();			// Fill summed arrays from input
};

template <int NE,int N2,int N3,int N4,int N5,int N6,int N7,int N8,int N9> inline
void G4CascadeData<NE,N2,N3,N4,N5,N6,N7,N8,N9>::initialize() {
  // Initialize index offsets for cross-section array (can't do globally)
  index[0] = 0;   index[1] = N02; index[2] = N23; index[3] = N24;
  index[4] = N25; index[5] = N26; index[6] = N27; index[7] = N28;
  index[8] = N29;

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
    sum[k] = 0.0;
    for (G4int m = 0; m < NM; m++) {
      sum[k] += multiplicities[m][k];
    }
  }
}


// Dump individual cross-section table
template <int NE,int N2,int N3,int N4,int N5,int N6,int N7,int N8,int N9> inline
void G4CascadeData<NE,N2,N3,N4,N5,N6,N7,N8,N9>::
printXsec(const G4double (&xsec)[NE]) const {
  for (G4int k=0; k<NE; k++) {
    G4cout << " " << std::setw(5) << xsec[k];
    if ((k+1)%12 == 0) G4cout << G4endl;
  }
  G4cout << G4endl;
}

// Dump tables for specified multiplicity
template <int NE,int N2,int N3,int N4,int N5,int N6,int N7,int N8,int N9> inline
void G4CascadeData<NE,N2,N3,N4,N5,N6,N7,N8,N9>::print(G4int mult) const {
  if (mult < 0) {
    G4cout << " Total cross section (" << NE << " bins):" << G4endl;
    printXsec(tot);
    G4cout << " Summed cross section (" << NE << " bins):" << G4endl;
    printXsec(sum);
    G4cout << " Individual channel cross sections" << G4endl;
    for (int m=2; m<NM+2; m++) print(m);
    return;
  }

  G4int im = mult-2;		// Convert multiplicity to array index

  G4int start = index[im];
  G4int stop = index[im+1];
  G4cout << "\n Mulitplicity " << mult << " (indices " << start << " to "
	 << stop-1 << ") summed cross section:" << G4endl;

  printXsec(multiplicities[im]);
  
  for (G4int i=start; i<stop; i++) {
    G4int ichan=i-start;
    G4cout << "\n final state x" << mult << "bfs[" << ichan << "] : ";
    for (G4int fsi=0; fsi<mult; fsi++) {
      switch (mult) {
      case 2: G4cout << " " << x2bfs[ichan][fsi]; break;
      case 3: G4cout << " " << x3bfs[ichan][fsi]; break;
      case 4: G4cout << " " << x4bfs[ichan][fsi]; break;
      case 5: G4cout << " " << x5bfs[ichan][fsi]; break;
      case 6: G4cout << " " << x6bfs[ichan][fsi]; break;
      case 7: G4cout << " " << x7bfs[ichan][fsi]; break;
      case 8: G4cout << " " << x8bfs[ichan][fsi]; break;
      case 9: G4cout << " " << x9bfs[ichan][fsi]; break;
      default: ;
      }
    }
    G4cout << " -- cross section [" << i << "]:" << G4endl;
    printXsec(crossSections[i]);
  }
}


// Dummy arrays for use when optional template arguments are skipped
template <int NE,int N2,int N3,int N4,int N5,int N6,int N7,int N8,int N9>
const G4int G4CascadeData<NE,N2,N3,N4,N5,N6,N7,N8,N9>::empty8bfs[1][8] = {{0}};

template <int NE,int N2,int N3,int N4,int N5,int N6,int N7,int N8,int N9>
const G4int G4CascadeData<NE,N2,N3,N4,N5,N6,N7,N8,N9>::empty9bfs[1][9] = {{0}};

#endif
