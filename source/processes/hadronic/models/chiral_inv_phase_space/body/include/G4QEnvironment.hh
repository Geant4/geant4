//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4QEnvironment.hh,v 1.9 2001-11-26 14:11:45 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QEnvironment ----------------
//             by Mikhail Kossov, August 2000.
//      header for Multy Quasmon Environment in the CHIPS Model
// ------------------------------------------------------------

#ifndef G4QEnvironment_h
#define G4QEnvironment_h 1

// Standard G4-headers
#include "G4QuasmonVector.hh"

class G4QEnvironment 
{
public:
  G4QEnvironment(const G4QHadronVector& projHadrons, const G4int targPDG);
  G4QEnvironment(const G4QEnvironment& right);         // copy QEnvironment by value
  G4QEnvironment(G4QEnvironment* right);               // copy QEnvironment by pointer
  ~G4QEnvironment();

  // Overloaded operators
  const G4QEnvironment& operator=(const G4QEnvironment& right);
  int operator==(const G4QEnvironment &right) const;
  int operator!=(const G4QEnvironment &right) const;

  //Selectors
  G4QNucleus       GetEnvironment() const;
  G4QuasmonVector* GetQuasmons();
  G4QHadronVector* GetQHadrons();
  G4QHadronVector* Fragment();                        // Unresp. wrapper for HadronizeQEnvironment()

  // Static functions
  static void SetParameters( G4double solAn=0.4, G4bool efFlag=false, G4double piThresh=141.4,
                            G4double mpisq=20000., G4double dinum=1880.);
  // General functions
  G4ThreeVector    RndmDir();                         // Randomize 3D direction (@@subst by libFunc)

private:  
  G4QHadronVector  HadronizeQEnvironment();           // Main HadronizationFunction used in Fragment
  void             CreateQuasmon(const G4QContent& projQC, const G4LorentzVector& proj4M);
  void             InitClustersVector(G4int maxC, G4int maxA); // Init.NucCclusters for 1st interact
  void             CleanUp();                         // Makes theEnvironment=vacuum & kill Quasmons
  void             PrepareInteractionProbabilities(const G4QContent& projQC, G4double AP);
  void             EvaporateResidual(G4QHadron* evap, G4bool corFlag = false);// Final Evaporation
  void             DecayDibaryon(G4QHadron* dB);      // Decay of any di-baryon (deuteron is kept)
  void             DecayThreeBaryon(G4QHadron* dB);   // Decay of ppp, nnn or LLL states
  void             DecayAntiStrange(G4QHadron* dB);   // Decay of the nucleus, containing K+/K0
  void             DecayAlphaBar(G4QHadron* dB);      // Decay of alpha+p or alpha+n states
  void             DecayAlphaDiN(G4QHadron* dB);      // Decay of alpha+p+p states
  void             DecayAlphaAlpha(G4QHadron* dB);    // Decay of alpha+alpha state
  G4bool           CheckGroundState(G4Quasmon* quasm, G4bool corFlag = false);// as G4Q for TotHadrV
  G4bool           DecayInEnvQ(G4Quasmon* quasm);     // Use befor evaporation in PANIC case

// Body
private:
  // Static Parameters
  static G4bool      EnergyFlux;      // Flag for Energy Flux use instead of Multy Quasmon
  static G4double    SolidAngle;      // Part of Solid Angle to capture secondaries (@@A-dep)
  static G4double    PiPrThresh;      // Pion Production Threshold for gammas
  static G4double    M2ShiftVir;      // Shift for M2=-Q2=m_pi^2 of the virtual gamma
  static G4double    DiNuclMass;      // Double Nucleon Mass for virtual normalization
  // Output hadrons
  G4QHadronVector    theQHadrons;     // Vector of generated secondary hadrons
  // Internal working parameters
  G4QCHIPSWorld      theWorld;        // the CHIPS World
  G4int              nBarClust;       // Maximum barion number of clusters (Calc. @ Interaction)
  G4double           f2all;           // Ratio of free nucleons to free+dense nucleons (do we need it in Quasmon?)
  G4QuasmonVector    theQuasmons;     // Intermediate vector of Quasmons before fragmentation (***delete***)
  G4QCandidateVector theQCandidates;  // Vector of possible candidates to clusters (***delete***)
  G4QNucleus         theEnvironment;  // Initial Nucleus & later Residual Nuclear Environment
  G4LorentzVector    tot4Mom;         // Total 4-momentum in the reaction
};

inline int G4QEnvironment::operator==(const G4QEnvironment &right) const {return this == &right;}
inline int G4QEnvironment::operator!=(const G4QEnvironment &right) const {return this != &right;}
inline G4QNucleus G4QEnvironment::GetEnvironment()                 const {return theEnvironment;}
#endif





