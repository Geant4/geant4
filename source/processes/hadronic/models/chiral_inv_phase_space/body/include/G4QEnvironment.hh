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
// $Id: G4QEnvironment.hh,v 1.4 2001-08-01 17:03:34 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4QEnvironment_h
#define G4QEnvironment_h 1


// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4QEnvironment ----------------
//             by Mikhail Kossov, August 2000.
//      header for Multy Quasmon Environment in the CHIPS Model
// ------------------------------------------------------------

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
  G4QHadronVector* Fragment();                   // GetQHadrons

  // Static functions
  static void SetParameters( G4double solAn, G4bool efFlag=false, G4double piThresh=141.4,
                            G4double mpisq=20000., G4double dinum=1880.);
private:  
  G4QHadronVector  HadronizeQEnvironment();
  void             CreateQuasmon(const G4QContent& projQC, const G4LorentzVector& proj4M);
  void             InitClustersVector(G4int maxClust);
  void             PrepareClusters();            // Prepare nuclear clusters for interaction
  void             CleanUp();                    // Nulefies theEnvironment & kills Quasmons
  void             PrepareInteractionProbabilities(const G4QContent& projQC);
  void             EvaporateResidual(G4QHadron* evap);
  void             DecayDibaryon(G4QHadron* dB);
  G4bool           DecayOutHadron(G4QHadron* qHadron);
  G4ThreeVector    RndmDir();

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
};

inline int G4QEnvironment::operator==(const G4QEnvironment &right) const {return this == &right;}
inline int G4QEnvironment::operator!=(const G4QEnvironment &right) const {return this != &right;}
inline G4QNucleus G4QEnvironment::GetEnvironment()                 const {return theEnvironment;}
#endif





