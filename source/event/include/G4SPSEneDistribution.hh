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
// G4SPSEneDistribution
//
// Class Description:
//
// To generate the energy of a primary vertex according to the
// defined distribution. This is a shared class between threads.
// Only one thread should use the set-methods here.
// Note that this is exactly what is achieved using UI commands.
// If you use the set methods to set defaults in your application take
// care that only one thread is executing them.
// In addition take care of calling these methods before the run is
// started. Do not use the setters during the event loop

// Author:       Fan Lei, QinetiQ ltd.
// Customer:     ESA/ESTEC
// History:
// - 05/02/2004, Fan Lei - Created.
//     Based on the G4GeneralParticleSource class.
// - 26/03/2014, Andrew Green.
//     Modification to use STL vectors instead of C-style arrays.
//     Also moved to dynamically allocated memory in the LinearInterpolation(),
//     ExpInterpolation() and LogInterpolation() functions.
// - 06/06/2014, Andrea Dotti.
//     For thread safety: this is a shared object.
//     Added mutex to control access to shared resources (data members).
//     in Getters and Setters, mutex is NOT used in GenerateOne() because it
//     is assumed that properties are not changed during event loop.
// - 24/11/2017, Fan Lei
//    Added cutoff power-law distribution option. Implementation is similar
//    to that of the BlackBody one.
// --------------------------------------------------------------------
#ifndef G4SPSEneDistribution_hh
#define G4SPSEneDistribution_hh 1

#include "G4PhysicsOrderedFreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"
#include "G4DataInterpolation.hh"
#include "G4Threading.hh"
#include "G4Cache.hh"
#include <vector>

#include "G4SPSRandomGenerator.hh"

class G4SPSEneDistribution
{
  public:

    G4SPSEneDistribution();
      // Constructor: initializes variables
   ~G4SPSEneDistribution();
      // Destructor

    void SetEnergyDisType(const G4String&);
      // Allows the user to choose the energy distribution type.
      // The arguments are: Mono (mono-energetic), Lin (linear),
      // Pow (power-law), Exp (exponential), Gauss (gaussian),
      // Brem (bremsstrahlung), BBody (black-body),
      // Cdg (cosmic diffuse gamma-ray), User (user-defined),
      // Arb (arbitrary point-wise), Epn (energy per nucleon)

    const G4String& GetEnergyDisType();

    void SetEmin(G4double);
      // Sets the minimum energy

    G4double GetEmin() const;
    G4double GetArbEmin();

    void SetEmax(G4double);
      // Sets the maximum energy

    G4double GetEmax() const;
    G4double GetArbEmax();

    void SetMonoEnergy(G4double);
      // Sets energy for mono-energetic distribution

    void SetAlpha(G4double);
      // Sets alpha for a power-law distribution

    void SetBiasAlpha(G4double);

    void SetTemp(G4double);
      // Sets Temperature for a Brem or BBody distributions

    void SetBeamSigmaInE(G4double);

    void SetEzero(G4double);
      // Sets Ezero for an exponential distribution

    void SetGradient(G4double);
      // Sets gradient for a linear distribution

    void SetInterCept(G4double);
      // Sets intercept for a linear distribution

    void UserEnergyHisto(const G4ThreeVector&);
      // Allows user to defined a histogram for the energy distribution

    void ArbEnergyHisto(const G4ThreeVector&);
      // Allows the user to define an Arbitrary set of points for the
      // energy distribution

    void ArbEnergyHistoFile(const G4String&);

    void EpnEnergyHisto(const G4ThreeVector&);
      // Allows the user to define an Energy per nucleon histogram

    void InputEnergySpectra(G4bool);
      // Allows the user to choose between momentum and energy histograms
      // for user-defined histograms and arbitrary point-wise spectra.
      // The default is true (energy)

    void InputDifferentialSpectra(G4bool);
      // Allows the user to choose between integral and differential 
      // distributions when using the arbitrary point-wise option

    void ArbInterpolate(const G4String&);
      // Allows the user to specify the type of function to
      // interpolate the Arbitrary points spectrum with

    const G4String& GetIntType();
    
    void Calculate();
      // Controls the calculation of Integral PDF for the Cdg and BBody
      // distributions

    void SetBiasRndm(G4SPSRandomGenerator* a);
      // Sets the biased random number generator

    void ReSetHist(const G4String&);
      // Resets the histogram for user defined distribution

    void SetVerbosity(G4int a);
      // Sets the verbosity level

    G4double GetWeight() const;

    G4double GetMonoEnergy();
      // Mono-energetic energy

    G4double GetSE();
      // Standard deviation for Gaussian distribution in energy

    G4double Getalpha() const;
      // Alpha (pow)

    G4double GetEzero() const;
      // E0 (exp)

    G4double GetTemp();
      // Temp (bbody,brem)

    G4double Getgrad() const;
      // Gradient and intercept for linear spectra

    G4double Getcept() const;

    G4PhysicsOrderedFreeVector GetUserDefinedEnergyHisto();

    G4PhysicsOrderedFreeVector GetArbEnergyHisto();

    G4double GenerateOne(G4ParticleDefinition*);
       // Generate one random energy for the specified particle

    G4double GetProbability (G4double);

    G4double GetArbEneWeight(G4double);

    inline void ApplyEnergyWeight(G4bool val) { applyEvergyWeight = val; }
    inline G4bool IfApplyEnergyWeight() const { return applyEvergyWeight; }

  private:

    void LinearInterpolation();
    void LogInterpolation();
    void ExpInterpolation();
    void SplineInterpolation();
    void CalculateCdgSpectrum();
    void CalculateBbodySpectrum();
    void CalculateCPowSpectrum();

    // The following methods generate energies according
    // to the spectral parameters defined above

    void GenerateMonoEnergetic();
    void GenerateBiasPowEnergies();
    void GenerateGaussEnergies();
    void GenerateBremEnergies();
    void GenerateBbodyEnergies();
    void GenerateCdgEnergies();
    void GenUserHistEnergies();
    void GenEpnHistEnergies();
    void GenArbPointEnergies(); // NOTE: REQUIRES UPDATE OF DATA MEMBERS
    void GenerateExpEnergies(G4bool);
    void GenerateLinearEnergies(G4bool);
    void GeneratePowEnergies(G4bool);
    void GenerateCPowEnergies();

    void ConvertEPNToEnergy();
      // Converts energy per nucleon to energy
    
    void BBInitHists();
    void CPInitHists();

  private:  // Non invariant data members become G4Cache

    G4String EnergyDisType; // energy dis type Variable  - Mono,Lin,Exp,etc
    G4double weight; // particle weight //// NOT INVARIANT
    G4double MonoEnergy; //Mono-energteic energy
    G4double SE; // Standard deviation for Gaussian distribution in energy

    G4double Emin, Emax; // emin and emax //// NOT INVARIANT
    G4double alpha, Ezero;// alpha (pow), E0 (exp) //// NOT INVARIANT
    G4double Temp; // Temp (bbody,brem)
    G4double biasalpha; // biased power index
    G4double grad, cept; // gradient and intercept for linear spectra //// NOT INVARIANT
    G4double prob_norm; // normalisation factor use in calculate the probability
    G4bool Biased = false; // biased to power-law
    G4bool EnergySpec = true; // energy spectra, false - momentum spectra
    G4bool DiffSpec = true; // differential spec, false integral spec

    G4PhysicsOrderedFreeVector UDefEnergyH; // energy hist data
    G4PhysicsOrderedFreeVector IPDFEnergyH;
    G4bool IPDFEnergyExist = false, IPDFArbExist = false, Epnflag = false;
    G4PhysicsOrderedFreeVector ArbEnergyH; // Arb x,y histogram
    G4PhysicsOrderedFreeVector IPDFArbEnergyH; // IPDF for Arb
    G4PhysicsOrderedFreeVector EpnEnergyH;
    G4double CDGhist[3]; // cumulative histo for cdg
    
    std::vector<G4double>* BBHist = nullptr;
    std::vector<G4double>* Bbody_x = nullptr;
    G4bool BBhistInit = false;
    G4bool BBhistCalcd = false;

    //  For cutoff power-law
    //
    std::vector<G4double>* CPHist = nullptr;
    std::vector<G4double>* CP_x = nullptr;
    G4bool CPhistInit = false;
    G4bool CPhistCalcd = false;

    G4String IntType; // Interpolation type
    G4double* Arb_grad = nullptr;
    G4double* Arb_cept = nullptr;
    G4bool Arb_grad_cept_flag = false;
    G4double* Arb_alpha = nullptr;
    G4double* Arb_Const = nullptr;
    G4bool Arb_alpha_Const_flag = false;
    G4double* Arb_ezero = nullptr;
    G4bool Arb_ezero_flag = false;

    G4bool applyEvergyWeight = false;

    G4double ArbEmin, ArbEmax;
      // Emin and Emax for the whole arb distribution used primarily for debug.

    G4double particle_energy;

    G4SPSRandomGenerator* eneRndm = nullptr;

    G4int verbosityLevel;

    G4PhysicsOrderedFreeVector ZeroPhysVector; // for re-set only

    std::vector<G4DataInterpolation*> SplineInt;
      // Holds Spline stuff required for sampling
    G4DataInterpolation* Splinetemp = nullptr;
      // Holds a temp Spline used for calculating area

    G4Mutex mutex; // protect access to shared resources

    // Thread local data (non-invariant during event loop).
    // These are copied from master one at the beginning of
    // generation of each event
    //
    struct threadLocal_t
    {
      G4double Emin;
      G4double Emax;
      G4double alpha;
      G4double Ezero;
      G4double grad;
      G4double cept;
      G4ParticleDefinition* particle_definition;
      G4double weight;
      G4double particle_energy;
    };
    G4Cache<threadLocal_t> threadLocalData;
};

#endif
