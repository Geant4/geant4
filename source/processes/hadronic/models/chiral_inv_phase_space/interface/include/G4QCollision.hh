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
// $Id: G4QCollision.hh,v 1.7 2006/10/30 10:33:36 mkossov Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
//      ---------------- G4QCollision header ----------------
//                 by Mikhail Kossov, December 2003.
//  Header of G4QCollision class (mu-,pi-,K-) of the CHIPS Simulation Branch in GEANT4
// -------------------------------------------------------------------------------
// This is a unique CHIPS class for the Hadron-Nuclear Inelastic Interaction Prosesses.
// -------------------------------------------------------------------------------
// At present (Dec.04) only pi+/-, K+/- proton, neutron, antiproton and antineutron
// collisions with protons are implemented, which are fundamental for the in matter
// simulation of hadronic reactions. The interactions of the same particles with
// nuclei are planned only. The collisions of nuclei with nuclei are possible...
// The simulation is based on the G4QuasmonString class, which extends the CHIPS model
// to the highest energyes, implementing the Quasmon string with the
// String->Quasmons->Hadrons scenario of the quark-gluon string fragmentation
// --> CHIPS is a SU(3) event generator, so it does not include reactions with the
// heavy (c,b,t), which can be simulated only by the SU(6) QUIPS (QUark Invariant
// Phase Space) model which is an expantion of the CHIPS.-December 2003.M.Kossov.-
// -------------------------------------------------------------------------------
// Algorithms: the interactions in CHIPS are described by the quark exchange (QE) process.
// The first step is the low energy quark exchange. If as a result of the QE one or
// both secondary hadrons are below the pi0 threshold (roughly) they are pushed to the
// Ground State (GS) value(s). The excited (above the pi0 production threshold) hadronic
// state is considered as a Quasmon, which is filled in the G4QuasmonVector of the
// G4QuasmonString class. On the second step all G4Quasmons are decayed by the
// G4Quasmon class and fiill the G4QHadronVector output. If the exchange quark is too far
// in the rapidity space (a parameter of the G4QuasmonString class) from any of the quarks
// of the other hadron it creates a string with the nearest in the rapidity space quark.
// This string is converted into a Quasmon. This forces the coalescence of the residuals
// in the another Quasmon, while the possibility exist to create more residual Quasmons
// instead of one - one per each target-quark+projectile-antiquark(diquark) pair. This
// possibility is tuned by the Drell-Yan pair production process. If the target (or
// pojectile) are nuclei, then the Quasmons are created not only in vacuum, where they
// can be fragmented by the G4Quasmon class, but in nuclear matter of the residual target
// (or projectile). If the Quasmons are crated in nuclear matter, they are fragmented by
// the G4QEnvironment class with the subsequent Quark Exchange nuclear fragmentation.
// This is the planned scenario.- December 2004.Mikhail Kossov.-
// --------------------------------------------------------------------------------
// ****************************************************************************************
// ********* This HEADER is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************

#ifndef G4QCollision_hh
#define G4QCollision_hh

// GEANT4 Headers
#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh" 
#include "G4VDiscreteProcess.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleTypes.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"

// CHIPS Headers
#include "G4QEnvironment.hh"
#include "G4VQCrossSection.hh"
#include "G4QIsotope.hh"
#include "G4QProtonNuclearCrossSection.hh"
#include "G4QPhotonNuclearCrossSection.hh"
#include "G4QElectronNuclearCrossSection.hh"
#include "G4QMuonNuclearCrossSection.hh"
#include "G4QTauNuclearCrossSection.hh"
#include "G4QNuMuNuclearCrossSection.hh"
#include "G4QANuMuNuclearCrossSection.hh"
//#include "G4QuasmonString.hh"
#include "G4QPDGToG4Particle.hh"
#include <vector>

class G4QCollision : public G4VDiscreteProcess
{
public:

  // Constructor
  G4QCollision(const G4String& processName ="CHIPSNuclearCollision");

  // Destructor
  ~G4QCollision();

  G4bool IsApplicable(const G4ParticleDefinition& particle);

  G4double GetMeanFreePath(const G4Track& aTrack, G4double previousStepSize,
                           G4ForceCondition* condition);
  // It returns the MeanFreePath of the process for the current track :
  // (energy, material)
  // The previousStepSize and G4ForceCondition* are not used.
  // This function overloads a virtual function of the base class.		      
  // It is invoked by the ProcessManager of the Particle.
 

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep); 
  // It computes the final state of the process (at end of step),
  // returned as a ParticleChange object.			    
  // This function overloads a virtual function of the base class.
  // It is invoked by the ProcessManager of the Particle.

  // Fake void functions
  void SetPhysicsTableBining(G4double, G4double, G4int) {;}
  void BuildPhysicsTable(const G4ParticleDefinition&) {;}
  void PrintInfoDefinition() {;}

  // Internal Energy-Momentum Residual
  G4LorentzVector GetEnegryMomentumConservation();

  // Number of neutrons in the target nucleus (primary)
  G4int GetNumberOfNeutronsInTarget();

  // Static functions ---------------------------------------------------------------------
  static void SetManual();
  static void SetStandard();
  static void SetParameters(G4double temper=180., G4double ssin2g=.1, G4double etaetap=.3,
                            G4double fN=0., G4double fD=0., G4double cP=1., G4double mR=1.,
                            G4int npCHIPSWorld=234, G4double solAn=.5, G4bool efFlag=false,
                            G4double piTh=141.4,G4double mpi2=20000.,G4double dinum=1880.);
  static void SetPhotNucBias(G4double phnB=1.);
  static void SetWeakNucBias(G4double ccnB=1.);
  //--- End of static member functions ----------------------------------------------------

  G4double GetPhotNucBias(){return photNucBias;}
  G4double GetWeakNucBias(){return weakNucBias;}

private:

  // Hide assignment operator as private 
  G4QCollision& operator=(const G4QCollision &right);

  // Copy constructor
  G4QCollision(const G4QCollision&);

		// BODY
  // Static Parameters --------------------------------------------------------------------
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
  // -> Biasing parameters:
  static G4double photNucBias; // Biasing parameter for photo-($e,mu,tau)Nuclear reactions
  static G4double weakNucBias; // Biasing parameter for Charged Currents (nu,mu) reactions
  //--------------------------------- End of static parameters ---------------------------
  // Working parameters
  G4VQCrossSection* theCS;
  G4LorentzVector EnMomConservation;                  // Residual of Energy/Momentum Cons.
  G4int nOfNeutrons;                                  // #of neutrons in the target nucleus

  // Modifires for the reaction
  G4double Time;                                      // Time shift of the capture reaction
  G4double EnergyDeposition;                          // Energy deposited in the reaction
  static std::vector <G4int> ElementZ;                // Z of the element(i) in theLastCalc
  static std::vector <G4double> ElProbInMat;          // SumProbabilityElements in Material
  static std::vector <std::vector<G4int>*> ElIsoN;    // N of isotope(j) of Element(i)
  static std::vector <std::vector<G4double>*> IsoProbInEl;// SumProbabIsotopes in Element i
};
#endif

