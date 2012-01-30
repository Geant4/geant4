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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4QFreeScattering -- header file
// M.V. Kossov, ITEP(Moscow), 24-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Oct-2006
// ----------------------------------------------------------------------
// Short description: Provides percentage of quasi-free and quasi-elastic
// reactions in the inelastic reactions.
// ----------------------------------------------------------------------

#ifndef G4QFreeScattering_h
#define G4QFreeScattering_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4QThd.hh"
#include <vector>
#include "G4QPDGCode.hh"
#include "G4QHadronVector.hh"

class G4QFreeScattering
{
 protected:

  G4QFreeScattering();                 // Constructor

 public:

  ~G4QFreeScattering();                 // Destructor

  static G4QFreeScattering* GetPointer(); // Gives a pointer to this singletone

  // scatter (pPDG,p4M) on a virtual nucleon (NPDG,N4M), result: final pair(newN4M,newp4M)
  // if(newN4M.e()==0.) - below threshold, XS=0, no scattering of the projectile happened
  std::pair<G4LorentzVector,G4LorentzVector> Scatter(G4int NPDG, G4LorentzVector N4M,
                                                     G4int pPDG, G4LorentzVector p4M);
  // Inelastic pion production Pi+N+H by (pPDG,p4M) on a virtual nucleon (NPDG,N4M)
  // Normally there are 3 QHadrons in the output, if ptr=0 -> DoNothing (under threshold) 
  G4QHadronVector* InElF(G4int NPDG, G4LorentzVector N4M, G4int pPDG, G4LorentzVector p4M);

  // hN El & Tot XS (IU) mean over nucleons of (Z,N): on p -> (Z=1,N=0), on n -> (Z=0,N=1)
  std::pair<G4double,G4double> GetElTotMean(G4double pIU, G4int hPDG, G4int Z, G4int N);

  // For hadron PDG with momentum Mom (GeV/c) on F(p/n) calculate <sig_el,sig_tot> pair(mb)
  // F=true corresponds to the Neutron target, F=false corresponds to the Proton target
  std::pair<G4double,G4double> GetElTotXS(G4double Mom, G4int PDG, G4bool F);//<sigEl,sigT>
  std::pair<G4double,G4double> FetchElTot(G4double pGeV,G4int PDG,G4bool F);//<E,T>fromAMDB

  // This member function is in internal CHIPS units and must not be used externally
  std::pair<G4double,G4double> CalcElTot(G4double pGeV,G4int Index);//(sigEl,sigTot)(Index)

 // Body
 private:
  static std::vector<std::pair<G4double,G4double>*> vX; // Vector of ETPointers to LogTable
};
#endif
