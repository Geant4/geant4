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
// $Id: G4RPGPionInelastic.hh 94553 2015-11-24 09:05:06Z gcosmo $
//
// Author: D. H. Wright
// Date:   15 August 2007
//

#ifndef G4RPGPionInelastic_h
#define G4RPGPionInelastic_h 1
 
// Class Description:
// Partial cross sections, multiplicities and final state particle types 
// required for pi+ and pi- inelastic scattering in the re-parameterized 
// Gheisha model

#include "G4RPGInelastic.hh"
 
 class G4RPGPionInelastic : public G4RPGInelastic
 {
 public:
    
   G4RPGPionInelastic(const G4String& modelName = "RPGPionInelastic");
    
   ~G4RPGPionInelastic() { }

   //   void printCrossSections() const;
    
 protected:

   G4int GetMultiplicityT12(G4double KE) const;
   G4int GetMultiplicityT32(G4double KE) const;

   std::vector<G4int>
   GetFSPartTypesForT32(G4int mult, G4double KE, G4int tindex) const;
   std::vector<G4int>
   GetFSPartTypesForT12(G4int mult, G4double KE, G4int tindex) const;

   std::vector<G4int> GetFSPartTypesForPipP(G4int mult, G4double KE) const
     {return GetFSPartTypesForT32(mult, KE, 0); }

   std::vector<G4int> GetFSPartTypesForPimN(G4int mult, G4double KE) const
     {return GetFSPartTypesForT32(mult, KE, 1); }

   std::vector<G4int> GetFSPartTypesForPipN(G4int mult, G4double KE) const
     {return GetFSPartTypesForT12(mult, KE, 1); }

   std::vector<G4int> GetFSPartTypesForPimP(G4int mult, G4double KE) const
     {return GetFSPartTypesForT12(mult, KE, 0); }

   static const G4int pipPindex[8][2];
   static const G4int pimPindex[8][2];

   static const G4int T32_2bfs[2][2][2];
   static const G4int T32_3bfs[2][7][3];
   static const G4int T32_4bfs[2][15][4];
   static const G4int T32_5bfs[2][24][5];
   static const G4int T32_6bfs[2][5][6];
   static const G4int T32_7bfs[2][6][7];
   static const G4int T32_8bfs[2][7][8];
   static const G4int T32_9bfs[2][8][9];

   static const G4int T12_2bfs[2][5][2];
   static const G4int T12_3bfs[2][13][3];
   static const G4int T12_4bfs[2][22][4];
   static const G4int T12_5bfs[2][31][5];
   static const G4int T12_6bfs[2][6][6];
   static const G4int T12_7bfs[2][7][7];
   static const G4int T12_8bfs[2][8][8];
   static const G4int T12_9bfs[2][9][9];

   static G4ThreadLocal G4double pipPtot[30];
   static G4ThreadLocal G4double pimPtot[30];
   static G4ThreadLocal G4double t12_dSigma_dMult[8][30];
   static G4ThreadLocal G4double t32_dSigma_dMult[8][30];

   static const G4float pipPCrossSections[74][30];
   static const G4float pimPCrossSections[101][30];

 };
#endif
