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
// $Id: G4QCoherentChargeExchange.hh,v 1.2 2010-01-14 11:24:36 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QCoherentChargeExchange header ----------------
//                 by Mikhail Kossov, December 2003.
//  Header of G4QCoherentChargeExchange class (hA) of the CHIPS Simulation Branch
// -------------------------------------------------------------------------------
// This is a unique CHIPS class for the Hadron-Nuclear Elastic Scattering Prosesses
// -------------------------------------------------------------------------------
// At present (Jan-06) only proton-to-neutron & neutron-to-proton scattering on nuclei
// are implemented. The scattering of mesons and nuclei on nuclei are possible...
// The simulation is based on the CHIPS approximation of total elastic and differential
// elastic cross sections from E=0 to the highest energyes.
// -------------------------------------------------------------------------------
// Short description: This class resolves an ambiguity in the definition of the
// "inelastic" cross section. As it was shown in Ph.D.Thesis (M.Kosov,ITEP,1979)
// it is more reasonable to subdivide the total cross-section in the coherent &
// incoherent parts, but the measuring method for the "inelastic" cross-sections
// consideres the lack of the projectile within the narrow forward solid angle
// with the consequent extrapolation of these partial cross-sections, corresponding
// to the particular solid angle, to the zero solid angle. The low angle region
// is shadowed by the elastic (coherent) scattering. BUT the coherent charge
// exchange (e.g. conversion p->n) is included by this procedure as a constant term
// in the extrapolation, so the "inelastic" cross-section differes from the
// incoherent cross-section by the value of the coherent charge exchange cross
// section. Fortunately, this cross-sectoion drops ruther fast with energy increasing.
// All Geant4 inelastic hadronic models (including CHIPS) simulate the incoherent
// reactions. So the incoherent (including quasielastic) cross-section must be used
// instead of the inelastic cross-section. For that the "inelastic" cross-section
// must be reduced by the value of the coherent charge-exchange cross-section, which
// is estimated (it must be tuned!) in this CHIPS class. The angular distribution
// is made (at present) identical to the corresponding coherent-elastic scattering 
// -----------------------------------------------------------------------------------

#ifndef G4QCoherentChargeExchange_hh
#define G4QCoherentChargeExchange_hh

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
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"

// CHIPS Headers
#include "G4QuasiFreeRatios.hh"
#include "G4QProtonElasticCrossSection.hh"
#include "G4QNeutronElasticCrossSection.hh"
#include "G4QIsotope.hh"
#include "G4QCHIPSWorld.hh"
#include "G4QHadron.hh"
#include <vector>

class G4QCoherentChargeExchange : public G4VDiscreteProcess
{
public:

  // Constructor
  G4QCoherentChargeExchange(const G4String& processName ="CHIPS_CoherChargeExScattering");

  // Destructor
  ~G4QCoherentChargeExchange();

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


  G4LorentzVector GetEnegryMomentumConservation();

  G4int GetNumberOfNeutronsInTarget();

private:

  // Hide assignment operator as private 
  G4QCoherentChargeExchange& operator=(const G4QCoherentChargeExchange &right);

  // Copy constructor
  G4QCoherentChargeExchange(const G4QCoherentChargeExchange&);

  // Calculate XS/t: oxs=true - only CS; xst=true - calculate XS, xst=false(oxs=f/t) - t/tm
  G4double CalculateXSt(G4bool oxs, G4bool xst, G4double p, G4int Z, G4int N, G4int pPDG);

  // BODY
  // Static Parameters --------------------------------------------------------------------
  static G4int    nPartCWorld; // The#of particles for hadronization (limit of A of fragm.)
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
