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
// $Id: G4QGluonString.hh,v 1.2 2008/10/02 21:10:07 dennis Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//      ---------------- G4QGluonString header ----------------
//                 by Mikhail Kossov, October 2006.
// Header of G4QGluonString (N,pi,K,aN,Hyperons,aHyp) of CHIPS Simulation Branch
// -------------------------------------------------------------------------------
// -->> At present (17.10.2006) started only with proton-nuclear reactions <<--
// -------------------------------------------------------------------------------
// This is a QGS CHIPS class for the Hadron-Nuclear Inelastic Interaction Prosesses
// -------------------------------------------------------------------------------
// This class follows the structure of the G4QCollision class and its content can be
// eventually moved to G4QCollision class. At present the alternative realisation of
// the hadronic inelastic process (G4QGluonString) is kept for comparison with other
// G4QCollision algorithms. It can be attached as an inelastic process to all hadrons
// --------------------------------------------------------------------------------
// ****************************************************************************************
// ***** This HEADER is a property of the CHIPS hadronic package in Geant4 (M. Kosov) *****
// *********** DO NOT MAKE ANY CHANGE without approval of Mikhail.Kossov@cern.ch **********
// ****************************************************************************************

#ifndef G4QGluonString_hh
#define G4QGluonString_hh

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
#include "G4HadronicProcessType.hh"

// CHIPS Headers
#include "G4QEnvironment.hh"
#include "G4VQCrossSection.hh"
#include "G4QIsotope.hh"
#include "G4QProtonNuclearCrossSection.hh"
#include "G4QPDGToG4Particle.hh"
#include <vector>

class G4QGluonString : public G4VDiscreteProcess
{
public:

  // Constructor
  G4QGluonString(const G4String& processName ="CHIPS_QGS_Inelastic");

  // Destructor
  ~G4QGluonString();

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
  //--- End of static member functions ----------------------------------------------------

private:

  // Hide assignment operator as private 
  G4QGluonString& operator=(const G4QGluonString &right);

  // Copy constructor
  G4QGluonString(const G4QGluonString&);

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

