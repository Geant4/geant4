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
//
//
// GEANT4 physics class: G4QuasiElRatios -- header file
// M.V. Kossov, ITEP(Moscow), 24-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 15-Oct-2006
// ----------------------------------------------------------------------
// This class has been extracted from the CHIPS model. 
// All the dependencies on CHIPS classes have been removed.
// Short description: Provides percentage of quasi-free and quasi-elastic
// reactions in the inelastic reactions.
// ----------------------------------------------------------------------

#ifndef G4QuasiElRatios_h
#define G4QuasiElRatios_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include <vector>
#include "G4LorentzVector.hh"

#include "G4ChipsProtonElasticXS.hh"
#include "G4ChipsNeutronElasticXS.hh"


class G4QuasiElRatios
{
    
public:
    G4QuasiElRatios();                 // Constructor
    
    ~G4QuasiElRatios();                 // Destructor
    
    // Pair(QuasiFree/Inelastic,QuasiElastic/QuasiFree)
    std::pair<G4double,G4double> GetRatios(G4double pIU, G4int prPDG, G4int tgZ, G4int tgN);
    // ChargeExchange/QuasiElastic factor pair<for protons (Z), for neutrons(N)>
    std::pair<G4double,G4double> GetChExFactor(G4double pIU, G4int pPDG, G4int Z, G4int N);
    // scatter (pPDG,p4M) on a virtual nucleon (NPDG,N4M), result: final pair(newN4M,newp4M)
    // if(newN4M.e()==0.) - below threshold, XS=0, no scattering of the projectile happened
    std::pair<G4LorentzVector,G4LorentzVector> Scatter(G4int NPDG, G4LorentzVector N4M,
                                                       G4int pPDG, G4LorentzVector p4M);
    // ChExer (pPDG,p4M) on a virtual nucleon (NPDG,N4M), result: final pair(newN4M,newp4M)
    // if(newN4M.e()==0.) - keep projectile, XS=0, no interaction of the progectile happened
    // User should himself change the charge (PDG) (e.g. pn->np, pi+n->pi0p, pi-p->pi0n etc.)
    // Recepy: change target n to p or taget p to n and conserve enrgy, changing projectile
    // Do not use for the nucleon, as it is already included in quasielastic, and for pi0.
    std::pair<G4LorentzVector,G4LorentzVector> ChExer(G4int NPDG, G4LorentzVector N4M,
                                                      G4int pPDG, G4LorentzVector p4M);
    // Mean hN El and Tot XS(IU) for the isotopic (Z,N): on p -> (Z=1,N=0), on n -> (Z=0,N=1)
    std::pair<G4double,G4double> GetElTot(G4double pIU, G4int hPDG, G4int Z, G4int N); //(IU)
    
    // Calculate ChEx/El ratio coefficient (p is in independent units, (Z,N) is a target)
    G4double ChExElCoef(G4double p, G4int Z, G4int N, G4int pPDG);
    
    // For hadron PDG with momentum Mom (GeV/c) on F(p/n) calculate <sig_el,sig_tot> pair(mb)
    // F=true corresponds to the Nroton target, F=false corresponds to the Proton target
    std::pair<G4double,G4double> GetElTotXS(G4double Mom, G4int PDG, G4bool F);//<sigEl,sigT>
    std::pair<G4double,G4double> FetchElTot(G4double pGeV,G4int PDG,G4bool F);//<E,T>fromAMDB
    
    G4bool RelDecayIn2(G4LorentzVector& theMomentum, G4LorentzVector& f4Mom, 
                       G4LorentzVector& s4Mom, G4LorentzVector& dir,
                       G4double maxCost = 1., G4double minCost = -1.);
    
private:
    G4ChipsProtonElasticXS* PCSmanager;
    G4ChipsNeutronElasticXS* NCSmanager;
    
    // These working member functions are in CHIPS units and must not be used externally
    G4double GetQF2IN_Ratio(G4double TotCS_mb, G4int A); // QuasiFree/Inelastic (fast)
    G4double CalcQF2IN_Ratio(G4double TCSmb, G4int A); // R=QuasuFree/Inelastic (sig_t in mb)
    std::pair<G4double,G4double> CalcElTot(G4double pGeV,G4int Index);//(sigEl,sigTot)(Index)
    
    // Body
private:
    std::vector<G4double*> *vT;     // Vector of pointers to LinTable
    std::vector<G4double*> *vL;     // Vector of pointers to LogTable
    std::vector<std::pair<G4double,G4double>*> *vX; // Vector of ETPointers to LogTable
    
    G4double lastSRatio;             // The last sigma value for which R was calculated
    G4double lastRRatio;             // The last ratio R which was calculated
    std::vector<G4int> vARatio;
    std::vector<G4double> vHRatio;     // Vector of max s initialized in the LinTable
    std::vector<G4int> vNRatio;    // Vector of topBin number initialized in LinTable
    std::vector<G4double> vMRatio;     // Vector of rel max ln(s) initialized in LogTable
    std::vector<G4int> vKRatio;     // Vector of topBin number initialized in LogTable
 
    G4int     lastARatio;             // theLast of calculated A
    G4double  lastHRatio;            // theLast of max s initialized in the LinTable
    G4int     lastNRatio;             // theLast of topBin number initialized in LinTable
    G4double  lastMRatio;            // theLast of rel max ln(s) initialized in LogTable
    G4int     lastKRatio;             // theLast of topBin number initialized in LogTable
    G4double* lastTRatio;             // theLast of pointer to LinTable in the C++ heap
    G4double* lastLRatio;             // theLast of pointer to LogTable in the C++ heap
    
    G4double lastPtot;              // The last momentum for which XS was calculated
    G4int    lastHtot;               // The last projPDG for which XS was calculated
    G4bool   lastFtot;            // The last nucleon for which XS was calculated
    std::pair<G4double,G4double> lastRtot; // The last result
    std::vector<G4int>  vItot;      // Vector of index for which XS was calculated
    std::vector<G4double> vMtot;      // Vector of rel max ln(p) initialized in LogTable
    std::vector<G4int> vKtot;      // Vector of topBin number initialized in LogTable
    G4int     lastItot;              // The Last index for which XS was calculated
    G4double  lastMtot;             // The Last rel max ln(p) initialized in LogTable
    G4int     lastKtot;             // The Last topBin number initialized in LogTable
    std::pair<G4double,G4double>* lastXtot; // The Last ETPointers to LogTable in heap


    
    
};
#endif
