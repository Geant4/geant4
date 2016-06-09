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
// $Id$
//
//      ---------------- G4QEnvironment ----------------
//             by Mikhail Kossov, August 2000.
//      header for Multy Quasmon Environment in the CHIPS Model
// ------------------------------------------------------------
// Short description: The G4QEnvironment class corresponds to the nuclear
// environment,  containing excited nucleons or nuclear clusters
// (G4Quasmons). In the constructer the nucleus (G4QNucleus) is clusterized
// and then the projectile hadron (or hadrons) can create one (or a few)
// Quasmons, which are fragmented by the Fragment member function. As a
// result a vector of G4QHadrons is created, which can include the residual
// nucleus in the ground state.
//---------------------------------------------------------------------

#ifndef G4QEnvironment_h
#define G4QEnvironment_h 1

#include "G4RandomDirection.hh"
#include "G4QuasmonVector.hh"
#include "G4QFreeScattering.hh"
#include "G4QThd.hh"

class G4QEnvironment 
{
public:
  G4QEnvironment(const G4QNucleus& theEnv); // Create Env and add Quasmons
                                            // later
  G4QEnvironment(const G4QHadronVector& projHadrons, const G4int targPDG);
  G4QEnvironment(const G4QEnvironment& right);  // copy QEnvironment by value
  G4QEnvironment(G4QEnvironment* right);        // copy QEnvironment by pointer
  ~G4QEnvironment();                                   // Public Destructor

  // Overloaded operators
  const G4QEnvironment& operator=(const G4QEnvironment& right);
  G4bool operator==(const G4QEnvironment &right) const;
  G4bool operator!=(const G4QEnvironment &right) const;

  //Selectors
  G4QNucleus       GetEnvironment() const;
  G4QuasmonVector* GetQuasmons();           // User is responsible for Destroy/Clear/Delete
  G4QHadronVector* GetQHadrons();           // User is responsible for Destroy/Clear/Delete
  G4QHadronVector* GetProjectiles();        // User is responsible for Destroy/Clear/Delete

  // Modifiers
  void AddQuasmon(G4Quasmon* Q);            // Add aQuasmon to theEnvironment
  G4QHadronVector* Fragment();              // User must clear and destroy the G4QHadronVec

  // External managment
  void DecayBaryon(G4QHadron* dB, G4QHadronVector* HV);     // Decay baryon
  void DecayMeson(G4QHadron* dB, G4QHadronVector* HV);      // Decay meson
  void DecayAntistrange(G4QHadron* aS, G4QHadronVector* HV);// Decay AntiStrNuc
  void CheckMassShell(G4QHadronVector* HV);                 // All hadrons on mass shell

  // Static functions
  static void SetParameters(G4double solAn=0.4,G4bool efFlag=false,G4double piThresh=141.4,
                            G4double mpisq=20000., G4double dinum=1880.);
  static void OpenElectromagneticDecays();
  static void CloseElectromagneticDecays();

protected:
  void CleanUpQHadrons();  // Only another G4QEnvironment issue can use it (can make mess)
  void FillQHadrons(G4QHadronVector* input); //Only another G4QEnvironment issue can use it

private:
  G4QHadronVector* FSInteraction();         // Final State Interaction after Hadronization
  G4QHadronVector  HadronizeQEnvironment(); // Main HadronizationFunction used in Fragment
  void             CopyAndDeleteHadronVector(G4QHadronVector* HV);//Copy HadrVect to Output
  void CreateQuasmon(const G4QContent& pQC, const G4LorentzVector& p4M, G4bool f=false);
  void             InitClustersVector(G4int maxC, G4int maxA);//Init.NucClust's for 1st int
  void             CleanUp();               // Makes theEnvironment=vacuum & kill Quasmons
  void             PrepareInteractionProbabilities(const G4QContent& projQC, G4double AP);
  void             EvaporateResidual(G4QHadron* h, G4bool f=true);// Evaporate NuclearFragm
  G4bool           CheckGroundState(G4Quasmon* quasm,G4bool corFlag=false);//as G4Q for QHV
  G4bool           DecayInEnvQ(G4Quasmon* quasm);  // Use befor evaporation in PANIC case

// Body
private:
  // Static Parameters
  static G4bool      WeakDecays;     // Flag for opening WeakDecays (notUsed: allAreClosed)
  static G4bool      ElMaDecays;     // Flag for opening ElectroMagDecays (true by default)
  static G4bool      EnergyFlux;     // Flag for Energy Flux use instead of Multy Quasmon
  static G4double    SolidAngle;     // Part of Solid Angle to capture secondaries(@@A-dep)
  static G4double    PiPrThresh;     // Pion Production Threshold for gammas
  static G4double    M2ShiftVir;     // Shift for M2=-Q2=m_pi^2 of the virtual gamma
  static G4double    DiNuclMass;     // Double Nucleon Mass for virtual normalization
  // Output hadrons
  G4QHadronVector    theQHadrons;    // Vector of generated secondary hadrons
  // Internal working objects
  G4QHadronVector    intQHadrons;    // Vector of postponed secondary hadrons
  G4QCHIPSWorld*     theWorld;       // the CHIPS World
  G4int              nBarClust;      // Maximum BarionNumber of clusters (To optimize calc)
  G4double           f2all;          // Ratio of freeNucleons to free+denseNucleons
  G4QuasmonVector    theQuasmons;    // Intermediate vectorOfQuasmons before fragmentation
  G4QCandidateVector theQCandidates; // Vector of possible candidates to clusters
  G4QNucleus theEnvironment; // InitialNucleus (later ResidualNuclearEnvironment)
  G4LorentzVector    tot4Mom;        // Total 4-momentum in the reaction
  G4int              totCharge;      // Total charge in the reaction (for current control)
  G4int              totBaryoN;      // Total baryon number in the reaction (for cur.cont)
  G4QHadronVector    theProjectiles; // Vector of projectiles in the interaction
  G4int              theTargetPDG;   // PDG of the target nucleus in the interaction
  G4QFreeScattering* theQFScat;      // Pointer to the CHIPS Quasi-Free Scatterer
};

// Inline functions
inline G4bool G4QEnvironment::operator==(const G4QEnvironment &rhs) const
                                                                     {return this == &rhs;}
inline G4bool G4QEnvironment::operator!=(const G4QEnvironment &rhs) const
                                                                     {return this != &rhs;}
inline G4QNucleus G4QEnvironment::GetEnvironment() const {return theEnvironment;}

#endif
