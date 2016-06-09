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
// $Id: G4QCaptureAtRest.hh,v 1.2 2010-06-25 09:46:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QCaptureAtRest header ----------------
//                 by Mikhail Kossov, December 2003.
// Header of G4QCaptureAtRest class of the CHIPS Simulation Branch in GEANT4
// -------------------------------------------------------------------------------
// At present (May 2009) only pi-, K- and antiNucleon capture are tested, which
// are the most crucial for the in matter simulation. The hyperon capture (Sigma-,
// Xi-, Omega-, antiSigma+) is implemented, but not tested and it is not clear how
// frequently this kind of interaction takes place in the simulation of the hadronic
// showers. The antiNeutron Capture At Rest is implemented by this G4QCaptureAtRest
// class, but it is not clear how the anti-neutrons are stopped in Geant4 tracking.
// It can be stopped only by interactions with electrons, as the annihilation cross
// section is huge and any interaction with nucleus results in annihilation. The
// mu- & tau- Capture At Rest (mu-,nu) & (mu-,nu) are weak processes, which must
// be simulated together with the reversed Betha decay (e-,nu). While mu- capture is
// similar to the pi- capture from the nuclear fragmentation point of view (the energy
// scale is shrinked because m_mu < m_pi and a part of the energy is lost because of
// the neutrino radiation), the time scale of the mu- capture process is not exact,
// but it is clear, that it is well delayed. By this reason the mu- capture can be
// excluded from the G4QCaptureAtRest and can be implemented in the "LongLivingDecay"
// branch of simulation, which includes excited states of nuclei and short living
// isotopes. On the "Fast Simulation" Level all radioactive isotopes, long living
// nuclear excitations, mu-atoms etc, which can be important for the background
// signals, must be collected in the continuous database and simulated separately.
// CHIPS is SU(3) event generator, so it does not include reactions with the heavy
// (c,b,t) quarks involved such as antiDs-, which can be simulated only by SU(6)
// QUIPS (QUark Invariant Phase Space) model. - May 2009, M.Kossov.-
// -------------------------------------------------------------------------------
// All algorithms are similar: the captured particle is absorbed by a nuclear cluster
// with the subsequent Quark Exchange nuclear fragmentation. The Anti-Proton (antiSigma+)
// Capture algorithm is more complicated: the anti-baryon annihilates with the quasyfree
// nucleons on the nuclear periphery. The peripheral interaction results in a number
// of mesons. A part of them misses the nucleus and comes directly to the output,
// while others create Multy Quasmon Excitation in the nucleus with the subsequent
// Quark Excange Fragmentation of the nucleus. At present the two step mechanism of
// the antiProton-Nucleus interaction is hardwired in the G4QEnvironment class, but
// with time the first step of the interaction can be moved to this G4QCaptureAtRest
// class, to make the G4QEnvirement class simpler and better defined. This is
// necessary because the G4QEnvironment class is going to loos the previlage of
// the CHIPS Head Class (as previously the G4Quasmon class lost it) and G4QCollision
// class is going to be the CHIPS Head Class, where a few Nuclear Environments can
// exist (e.g. the Nuclear Environment of the Projectile Nucleus and the Nuclear
// Environment of the Target Nucleus). By the way, the antiProton-H1 interaction At
// Rest (CHIPSI) can be still simulated with only the G4Quasmon class, as this
// reaction does not have any nuclear environment.- May 2009, Mikhail Kossov.-
// --------------------------------------------------------------------------------
// ****************************************************************************************
// This Header is a part of the CHIPS physics package (author: M. Kosov)
// ****************************************************************************************
// Short Description: This is a universal process for nuclear capture
// (including annihilation) of all negative particles (cold neutrons, negative
// hadrons, negative leptons: mu- & tau-). It can be used for the cold neutron
// capture, but somebody should decide what is the probability (defined
// by the capture cross-section and atomic material properties) to switch
// the cold neutron to the at-rest neutron. - M.K. 2009.
// ----------------------------------------------------------------------

#ifndef G4QCaptureAtRest_hh
#define G4QCaptureAtRest_hh

// GEANT4 Headers
#include "globals.hh"
#include "G4ios.hh"
#include "G4VRestProcess.hh"
#include "G4ParticleTypes.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4RandomDirection.hh"

// CHIPS Headers
#include "G4QEnvironment.hh"
#include "G4QIsotope.hh"
#include "G4QPDGToG4Particle.hh"

class G4QCaptureAtRest : public G4VRestProcess
{  
private:

  // Hide assignment operator as private 
  G4QCaptureAtRest& operator=(const G4QCaptureAtRest &right);

  // Copy constructor
  G4QCaptureAtRest(const G4QCaptureAtRest& );

public:

  // Constructor
  G4QCaptureAtRest(const G4String& processName ="CHIPSNuclearCaptureAtRest");

  // Destructor
  virtual ~G4QCaptureAtRest();

  virtual G4bool IsApplicable(const G4ParticleDefinition& particle);

  G4VParticleChange* AtRestDoIt(const G4Track& aTrack, const G4Step& aStep); 

  G4LorentzVector GetEnegryMomentumConservation();

  G4int GetNumberOfNeutronsInTarget();

  // Static functions
  static void SetManual();
  static void SetStandard();
  static void SetParameters(G4double temper=180., G4double ssin2g=.1, G4double etaetap=.3,
                            G4double fN=0., G4double fD=0., G4double cP=1., G4double mR=1.,
                            G4int npCHIPSWorld=234, G4double solAn=.5, G4bool efFlag=false,
                            G4double piTh=141.4,G4double mpi2=20000.,G4double dinum=1880.);

protected:                         

  // zero mean lifetime
  G4double GetMeanLifeTime(const G4Track& aTrack, G4ForceCondition* );
  G4double RandomizeDecayElectron(G4int Z); // Randomize energy of decay electron (in MeV)
private:

  G4bool RandomizeMuDecayOrCapture(G4int Z, G4int N); // true=MuCapture, false=MuDecay
  void CalculateEnergyDepositionOfMuCapture(G4int Z); // (2p->1s, MeV) @@ Now N-independent
  G4bool RandomizeTauDecayOrCapture(G4int Z, G4int N);// true=TauCapture, false=TauDecay
  void CalculateEnergyDepositionOfTauCapture(G4int Z);// (2p->1s, MeV) @@N-independ,Improve

// BODY
private:
  // Static Parameters
  static G4bool   manualFlag;  // If false then standard parameters are used
  static G4int    nPartCWorld; // The#of particles for hadronization (limit of A of fragm.)
  // -> Parameters of the G4Quasmon class:
  static G4double Temperature; // Quasmon Temperature
  static G4double SSin2Gluons; // Percent of ssbar sea in a constituen gluon
  static G4double EtaEtaprime; // Part of eta-prime in all etas
  // -> Parameters of the G4QNucleus class:
  static G4double freeNuc;     // probability of the quasi-free baryon on surface
  static G4double freeDib;     // probability of the quasi-free dibaryon on surface
  static G4double clustProb;   // clusterization probability in dense region
  static G4double mediRatio;   // relative vacuum hadronization probability
  // -> Parameters of the G4QEnvironment class:
  static G4bool   EnergyFlux;  // Flag for Energy Flux use instead of Multy Quasmon
  static G4double SolidAngle;  // Part of Solid Angle to capture secondaries(@@A-dep)
  static G4double PiPrThresh;  // Pion Production Threshold for gammas
  static G4double M2ShiftVir;  // Shift for M2=-Q2=m_pi^2 of the virtual gamma
  static G4double DiNuclMass;  // Double Nucleon Mass for virtual normalization
  //
  // Working parameters
  G4LorentzVector EnMomConservation;                  // Residual of Energy/Momentum Cons.
  G4int nOfNeutrons;                                  // #of neutrons in the target nucleus
  // Modifires for the reaction
  G4double Time;                                      // Time shift of the capture reaction
  G4double EnergyDeposition;                          // Energy deposited in the reaction

};
#endif
