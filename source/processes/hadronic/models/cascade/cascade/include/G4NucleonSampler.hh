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
// $Id: G4NucleonSampler.hh,v 1.1 2009/09/24 18:50:27 dennis Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
// Author: D. H. Wright
// Date:   10 September 2009
//

#ifndef G4NucleonSampler_h
#define G4NucleonSampler_h 1
 
// Class Description:
// Samples partial cross sections, multiplicities and final state particle
// types required for p and n scattering within a nucleus

#include "G4FinalStateSampler.hh"
 
class G4NucleonSampler : public G4FinalStateSampler
{
 public:
    
   G4NucleonSampler();
    
   ~G4NucleonSampler() 
   { }

   void printCrossSections() const;
    
   G4int GetMultiplicityT1(G4double KE) const;   // pp, nn
   G4int GetMultiplicityT0(G4double KE) const;   // pn

   std::vector<G4int>
   GetFSPartTypesForT1(G4int mult, G4double KE, G4int tindex) const;
   std::vector<G4int>
   GetFSPartTypesForT0(G4int mult, G4double KE) const;

   std::vector<G4int> GetFSPartTypesForPP(G4int mult, G4double KE) const
     {return GetFSPartTypesForT1(mult, KE, 0); }

   std::vector<G4int> GetFSPartTypesForNN(G4int mult, G4double KE) const
     {return GetFSPartTypesForT1(mult, KE, 1); }

   std::vector<G4int> GetFSPartTypesForNP(G4int mult, G4double KE) const
     {return GetFSPartTypesForT0(mult, KE); }

 protected:

   G4int PPindex[8][2];
   G4int NPindex[8][2];

   G4int T1_2bfs[2][1][2];
   G4int T1_3bfs[2][6][3];
   G4int T1_4bfs[2][18][4];
   G4int T1_5bfs[2][32][5];
   G4int T1_6bfs[2][7][6];
   G4int T1_7bfs[2][8][7];
   G4int T1_8bfs[2][10][8];
   G4int T1_9bfs[2][11][9];

   G4int T0_2bfs[1][2];
   G4int T0_3bfs[9][3];
   G4int T0_4bfs[22][4];
   G4int T0_5bfs[38][5];
   G4int T0_6bfs[7][6];
   G4int T0_7bfs[9][7];
   G4int T0_8bfs[10][8];
   G4int T0_9bfs[12][9];

   G4double PPsummed[30];
   G4double PPtot[30];
   G4double NPsummed[30];
   G4double NPtot[30];
   G4double t1_dSigma_dMult[8][30];
   G4double t0_dSigma_dMult[8][30];

   G4float PPCrossSections[93][30];
   G4float NPCrossSections[108][30];

private:

   void initCrossSections();

 };
#endif
