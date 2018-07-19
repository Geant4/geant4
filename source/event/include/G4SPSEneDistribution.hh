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
///////////////////////////////////////////////////////////////////////////////
//
// MODULE:        G4SPSEneDistribution.hh
//
// Version:      1.0
// Date:         5/02/04
// Author:       Fan Lei 
// Organisation: QinetiQ ltd.
// Customer:     ESA/ESTEC
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
//
// Version 1.0, 05/02/2004, Fan Lei, Created.
//    Based on the G4GeneralParticleSource class in Geant4 v6.0
//
//
//  26/03/2014, Andrew Green.
//      Modification to used STL vectors instead of C-style arrays. This should save some space,
//      particularly when the blackbody function is not used. Also moved to dynamically allocated
//      memory in the LinearInterpolation, ExpInterpolation and LogInterpolation functions. Again,
//      this will save space if these functions are unused.
//
// 06/06/2014  A Dotti
//    For thread safety: this is a shared object,
//    mutex has been added to control access to shared resources (data members).
//    in Getters and Setters, mutex is NOT used in GenerateOne because it is
//    assumed that properties are not changed during event loop.
//
// 24/11/2017  Fan Lei
//    Added cutoff power-law distribution option. Implementation is similar to that of the BlackBody one.
//
///////////////////////////////////////////////////////////////////////////////
//
//
// Class Description:
//
// To generate the energy of a primary vertex according to the defined distribution 
//
///////////////////////////////////////////////////////////////////////////////
//
// MEMBER FUNCTIONS
// ----------------
//
// G4SPSEneDistribution ()
//    Constructor: Initializes variables
//
// ~G4SPSEneDistribution ()
//    Destructor: 
//
// void SetEnergyDisType(G4String)
//    Allows the user to choose the energy distribution type. The arguments
//    are Mono (mono-energetic), Lin (linear), Pow (power-law), Exp 
//    (exponential), Gauss (gaussian), Brem (bremsstrahlung), BBody (black-body), Cdg
//    (cosmic diffuse gamma-ray), User (user-defined), Arb (arbitrary
//    point-wise), Epn (energy per nucleon).
//
// void SetEmin(G4double)
//    Sets the minimum energy.
//
// void SetEmax(G4double)
//    Sets the maximum energy.
//
// void SetMonoEnergy(G4double)
//    Sets energy for mono-energetic distribution.
//
// void SetAlpha(G4double)
//    Sets alpha for a power-law distribution.
//
// void SetTemp(G4double)
//    Sets Temperature for a Brem or BBody distributions.
//
// void SetEzero(G4double)
//    Sets Ezero for an exponential distribution.
//
// void SetGradient(G4double)
//    Sets gradient for a linear distribution.
//
// void SetInterCept(G4double)
//    Sets intercept for a linear distribution.
//
// void UserEnergyHisto(G4ThreeVector)
//    Allows user to defined a histogram for the energy distribution.
//
// void ArbEnergyHisto(G4ThreeVector)
//    Allows the user to define an Arbitrary set of points for the
//    energy distribution.
//
// void EpnEnergyHisto(G4ThreeVector)
//    Allows the user to define an Energy per nucleon histogram.
//
// void Calculate()
//    Controls the calculation of Integral PDF for the Cdg and BBody
//    distributions.
//
// void InputEnergySpectra(G4bool)
//    Allows the user to choose between momentum and energy histograms
//    for user-defined histograms and arbitrary point-wise spectr.
//    The default is true (energy).
//
// void InputDifferentialSpectra(G4bool)
//    Allows the user to choose between integral and differential 
//    distributions when using the arbitrary point-wise option.
//
// void ArbInterpolate(G4String)
//    ArbInterpolate allows the user to specify the type of function to
//    interpolate the Arbitrary points spectrum with.
//
//  void SetBiasRndm (G4SPSRandomGenerator* a)
//    Sets the biased random number generator
//
//  G4double GenerateOne(G4ParticleDefinition*);
//    Generate one random energy for the specified particle
//
//  void ReSetHist(G4String);
//    Re-sets the histogram for user defined distribution
//
// void SetVerbosity(G4int)
//    Sets the verbosity level.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef G4SPSEneDistribution_h
#define G4SPSEneDistribution_h 1

#include "G4PhysicsOrderedFreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"
#include "G4DataInterpolation.hh"
#include "G4Threading.hh"
#include "G4Cache.hh"
#include <vector>

#include "G4SPSRandomGenerator.hh"

/** Andrea Dotti Feb 2015
 * Important: This is a shared class between threads.
 * Only one thread should use the set-methods here.
 * Note that this is exactly what is achieved using UI commands.
 * If you use the set methods to set defaults in your
 * application take care that only one thread is executing them.
 * In addition take care of calling these methods before the run is started
 * Do not use these setters during the event loop
 */

class G4SPSEneDistribution {
public:
	G4SPSEneDistribution();//
	~G4SPSEneDistribution();//

	void SetEnergyDisType(G4String);//
    G4String GetEnergyDisType();//
	void SetEmin(G4double);//
    G4double GetEmin();//
    G4double GetArbEmin();//
	void SetEmax(G4double);//
    G4double GetEmax();//
    G4double GetArbEmax();//
	void SetMonoEnergy(G4double);//
	void SetAlpha(G4double);//
	void SetBiasAlpha(G4double);//
	void SetTemp(G4double);//
	void SetBeamSigmaInE(G4double);////
	void SetEzero(G4double);//
	void SetGradient(G4double);//
	void SetInterCept(G4double);//
	void UserEnergyHisto(G4ThreeVector);//
	void ArbEnergyHisto(G4ThreeVector);//
	void ArbEnergyHistoFile(G4String);//
	void EpnEnergyHisto(G4ThreeVector);//

	void InputEnergySpectra(G4bool);//
	void InputDifferentialSpectra(G4bool);//
	void ArbInterpolate(G4String);//
	G4String GetIntType();//
    
	void Calculate();//

	void SetBiasRndm(G4SPSRandomGenerator* a);//
	// method to re-set the histograms
	void ReSetHist(G4String);//
	// Set the verbosity level.
	void SetVerbosity(G4int a);//

	//x
	G4double GetWeight();

	G4double GetMonoEnergy(); //Mono-energteic energy
	G4double GetSE();// Standard deviation for Gaussion distrbution in energy
	G4double Getalpha(); // alpha (pow)
	G4double GetEzero(); // E0 (exp)
	G4double GetTemp(); // Temp (bbody,brem)
	G4double Getgrad(); // gradient and intercept for linear spectra
	G4double Getcept(); //

    G4PhysicsOrderedFreeVector GetUserDefinedEnergyHisto(); //
	G4PhysicsOrderedFreeVector GetArbEnergyHisto(); //

	G4double GenerateOne(G4ParticleDefinition*);
	G4double GetProbability (G4double);


private:
	void LinearInterpolation();//
	void LogInterpolation();//
	void ExpInterpolation();//
	void SplineInterpolation();//
	void CalculateCdgSpectrum();//
	void CalculateBbodySpectrum();//
	void CalculateCPowSpectrum();//

	// The following methods generate energies according to the spectral
	// parameters defined above.
	void GenerateMonoEnergetic();//G4double& outputEne);//
	void GenerateBiasPowEnergies();//G4double& outputEne,G4double& outputWeight);//
	void GenerateGaussEnergies();//
	void GenerateBremEnergies();//
	void GenerateBbodyEnergies();//
	void GenerateCdgEnergies();//
	void GenUserHistEnergies();//
	void GenEpnHistEnergies();//
	void GenArbPointEnergies();//<<<<<<<<<<< DOES NOT WORK, REQUIRES UPDATE OF DATA MEMBERS.
    void GenerateExpEnergies(G4bool);//
	void GenerateLinearEnergies(G4bool);//
	void GeneratePowEnergies(G4bool);//
	void GenerateCPowEnergies();//

    // converts energy per nucleon to energy.
	void ConvertEPNToEnergy();
    
    void BBInitHists();//
	void CPInitHists();//

private:

	G4String EnergyDisType; // energy dis type Variable  - Mono,Lin,Exp,etc
	G4double weight; // particle weight //// NOT INVARIANT
	G4double MonoEnergy; //Mono-energteic energy
	G4double SE; // Standard deviation for Gaussion distrbution in energy
    //Non invariant data members become G4Cache
	G4double Emin, Emax; // emin and emax                                         ////// NOT INVARIANT
	G4double alpha, Ezero;// alpha (pow), E0 (exp)                                ////// NOT INVARIANT
    G4double Temp; // Temp (bbody,brem)
	G4double biasalpha; // biased power index
	G4double grad, cept; // gradient and intercept for linear spectra             ////// NOT INVARIANT
    G4double prob_norm; // normalisation factor use in calculate the probability
    G4bool Biased; // true - biased to power-law
	G4bool EnergySpec; // true - energy spectra, false - momentum spectra
	G4bool DiffSpec; // true - differential spec, false integral spec
	//G4bool ApplyRig; // false no rigidity cutoff, true then apply one
	//G4double ERig; // energy of rigidity cutoff
	G4PhysicsOrderedFreeVector UDefEnergyH; // energy hist data
	G4PhysicsOrderedFreeVector IPDFEnergyH;
	G4bool IPDFEnergyExist, IPDFArbExist, Epnflag;
	G4PhysicsOrderedFreeVector ArbEnergyH; // Arb x,y histogram
	G4PhysicsOrderedFreeVector IPDFArbEnergyH; // IPDF for Arb
	G4PhysicsOrderedFreeVector EpnEnergyH;
	G4double CDGhist[3]; // cumulative histo for cdg
    
    //AG: Begin edit to use STL vectors.
//	G4double BBHist[10001], Bbody_x[10001];
    std::vector<G4double>* BBHist;
    std::vector<G4double>* Bbody_x;
    G4bool BBhistInit;
    G4bool BBhistCalcd;
//  For cutoff power-law	
    std::vector<G4double>* CPHist;
    std::vector<G4double>* CP_x;
    G4bool CPhistInit;
    G4bool CPhistCalcd;
	
    
    //AG: Edit here to use dynamic memory, will save space inless these functions are used.
	G4String IntType; // Interpolation type
//	G4double Arb_grad[1024], Arb_cept[1024]; // grad and cept for 1024 segments AG: Switched to DMA
    G4double* Arb_grad;
    G4double* Arb_cept;
    G4bool Arb_grad_cept_flag;
//	G4double Arb_alpha[1024], Arb_Const[1024]; // alpha and constants AG: Switched to DMA
    G4double* Arb_alpha;
    G4double* Arb_Const;
    G4bool Arb_alpha_Const_flag;
//	G4double Arb_ezero[1024]; // ezero AG: Switched to DMA
    G4double* Arb_ezero;
    G4bool Arb_ezero_flag;
	G4double ArbEmin, ArbEmax; // Emin and Emax for the whole arb distribution used primarily for debug.

	G4double particle_energy;

	G4SPSRandomGenerator* eneRndm;

	// Verbosity
	G4int verbosityLevel;

	G4PhysicsOrderedFreeVector ZeroPhysVector; // for re-set only

    std::vector<G4DataInterpolation*> SplineInt;//[1024]; // holds Spline stuff required for sampling
	G4DataInterpolation *Splinetemp; // holds a temp Spline used for calculating area

    G4Mutex mutex; // protect access to shared resources
    //Thread local data (non-invariant during event loop).
    //These are copied from master one at the beginning of generation
    //of each event
    struct threadLocal_t {
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

