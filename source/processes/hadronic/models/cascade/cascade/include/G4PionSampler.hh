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
// $Id: G4PionSampler.hh,v 1.2 2009/09/17 18:10:44 dennis Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
// Author: D. H. Wright
// Date:   26 March 2009
//

#ifndef G4PionSampler_h
#define G4PionSampler_h 1
 
// Class Description:
// Samples partial cross sections, multiplicities and final state particle
// types required for pi+, pi- and pi0 scattering within a nucleus

#include "G4FinalStateSampler.hh"
 
class G4PionSampler : public G4FinalStateSampler
{
 public:
    
   G4PionSampler();
    
   ~G4PionSampler() 
   { }

   void printCrossSections() const;
    
   G4int GetMultiplicityT33(G4double KE) const;   // pi+p, pi-n
   G4int GetMultiplicityT31(G4double KE) const;   // pi-p, pi+n
   G4int GetMultiplicityT11(G4double KE) const;   // pi0p, pi0n

   std::vector<G4int>
   GetFSPartTypesForT33(G4int mult, G4double KE, G4int tindex) const;
   std::vector<G4int>
   GetFSPartTypesForT31(G4int mult, G4double KE, G4int tindex) const;
   std::vector<G4int>
   GetFSPartTypesForT11(G4int mult, G4double KE, G4int tindex) const;

   std::vector<G4int> GetFSPartTypesForPipP(G4int mult, G4double KE) const
     {return GetFSPartTypesForT33(mult, KE, 0); }

   std::vector<G4int> GetFSPartTypesForPimN(G4int mult, G4double KE) const
     {return GetFSPartTypesForT33(mult, KE, 1); }

   std::vector<G4int> GetFSPartTypesForPimP(G4int mult, G4double KE) const
     {return GetFSPartTypesForT31(mult, KE, 0); }

   std::vector<G4int> GetFSPartTypesForPipN(G4int mult, G4double KE) const
     {return GetFSPartTypesForT31(mult, KE, 1); }

   std::vector<G4int> GetFSPartTypesForPizP(G4int mult, G4double KE) const
     {return GetFSPartTypesForT11(mult, KE, 0); }

   std::vector<G4int> GetFSPartTypesForPizN(G4int mult, G4double KE) const
     {return GetFSPartTypesForT11(mult, KE, 1); }

 protected:

   G4int pipPindex[8][2];
   G4int pimPindex[8][2];
   G4int pizPindex[8][2];

   G4int T33_2bfs[2][2][2];
   G4int T33_3bfs[2][7][3];
   G4int T33_4bfs[2][15][4];
   G4int T33_5bfs[2][24][5];
   G4int T33_6bfs[2][5][6];
   G4int T33_7bfs[2][6][7];
   G4int T33_8bfs[2][7][8];
   G4int T33_9bfs[2][8][9];

   G4int T31_2bfs[2][5][2];
   G4int T31_3bfs[2][13][3];
   G4int T31_4bfs[2][22][4];
   G4int T31_5bfs[2][31][5];
   G4int T31_6bfs[2][6][6];
   G4int T31_7bfs[2][7][7];
   G4int T31_8bfs[2][8][8];
   G4int T31_9bfs[2][9][9];

   G4int T11_2bfs[2][5][2];
   G4int T11_3bfs[2][13][3];
   G4int T11_4bfs[2][21][4];
   G4int T11_5bfs[2][30][5];
   G4int T11_6bfs[2][6][6];
   G4int T11_7bfs[2][7][7];
   G4int T11_8bfs[2][8][8];
   G4int T11_9bfs[2][9][9];

   G4double pipPsummed[30];
   G4double pipPtot[30];
   G4double pimPsummed[30];
   G4double pimPtot[30];
   G4double pizPsummed[30];
   G4double pizPtot[30];
   G4double t31_dSigma_dMult[8][30];
   G4double t33_dSigma_dMult[8][30];
   G4double t11_dSigma_dMult[8][30];

   G4float pipPCrossSections[74][30];
   G4float pimPCrossSections[101][30];
   G4float pizPCrossSections[99][30];

 private:

   void initCrossSections();

 };
#endif
