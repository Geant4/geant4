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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4QCollision.hh,v 1.2 2005/02/04 08:53:50 mkossov Exp $
// GEANT4 tag $Name: geant4-07-00-patch-01 $
//
//      ---------------- G4QCollision header ----------------
//                 by Mikhail Kossov, December 2003.
//  Header of G4QCollision class (mu-,pi-,K-) of the CHIPS Simulation Branch in GEANT4
// -------------------------------------------------------------------------------
// This is a unique CHIPS class for the Nuclear Capture At Rest Prosesses.
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
#include "G4QuasmonString.hh"
//<vector> is included in G4QIsotope.hh
//#include <vector>

class G4QCollision : public G4VDiscreteProcess
{  

public:

  // Constructor
  G4QCollision(const G4String& processName ="CHIPSNuclearCollision");

  // Destructor
  ~G4QCollision();

  G4bool IsApplicable(const G4ParticleDefinition& particle); // Now only for protons

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


  G4LorentzVector GetEnegryMomentumConservation();

  G4int GetNumberOfNeutronsInTarget();

protected:                         

private:

  // Hide assignment operator as private 
  G4QCollision& operator=(const G4QCollision &right);

  // Copy constructor
  G4QCollision(const G4QCollision&);

		// Body
  G4VQCrossSection* theCS;
  G4LorentzVector EnMomConservation;                  // Residual of Energy/Momentum Cons.
  G4int nOfNeutrons;                                  // #of neutrons in the target nucleus
};
#endif

