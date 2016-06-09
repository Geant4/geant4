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

//
#include "G4SPSRandomGenerator.hh"

class G4SPSEneDistribution {
public:
	G4SPSEneDistribution();
	~G4SPSEneDistribution();

	void SetEnergyDisType(G4String);
	inline G4String GetEnergyDisType() {
		return EnergyDisType;
	}
	;
	void SetEmin(G4double);
	inline G4double GetEmin() {
		return Emin;
	}
	;
	inline G4double GetArbEmin() {
		return ArbEmin;
	}
	;
	void SetEmax(G4double);
	inline G4double GetEmax() {
		return Emax;
	}
	;
	inline G4double GetArbEmax() {
		return ArbEmax;
	}
	;
	void SetMonoEnergy(G4double);
	void SetAlpha(G4double);
	void SetBiasAlpha(G4double);
	void SetTemp(G4double);
	void SetBeamSigmaInE(G4double);
	void SetEzero(G4double);
	void SetGradient(G4double);
	void SetInterCept(G4double);
	void UserEnergyHisto(G4ThreeVector);
	void ArbEnergyHisto(G4ThreeVector);
	void ArbEnergyHistoFile(G4String);
	void EpnEnergyHisto(G4ThreeVector);

	void InputEnergySpectra(G4bool);
	void InputDifferentialSpectra(G4bool);
	void ArbInterpolate(G4String);
	inline G4String GetIntType() {
		return IntType;
	}
	;
	void Calculate();
	//
	void SetBiasRndm(G4SPSRandomGenerator* a) {
		eneRndm = a;
	}
	;
	// method to re-set the histograms
	void ReSetHist(G4String);
	// Set the verbosity level.
	void SetVerbosity(G4int a) {
		verbosityLevel = a;
	}
	;
	//x
	G4double GetWeight() {
		return weight;
	}

	G4double GetMonoEnergy() {
		return MonoEnergy;
	}
	; //Mono-energteic energy
	G4double GetSE() {
		return SE;
	}
	; // Standard deviation for Gaussion distrbution in energy
	G4double Getalpha() {
		return alpha;
	}
	; // alpha (pow)
	G4double GetEzero() {
		return Ezero;
	}
	; // E0 (exp)
	G4double GetTemp() {
		return Temp;
	}
	; // Temp (bbody,brem)
	G4double Getgrad() {
		return grad;
	}
	; // gradient and intercept for linear spectra
	G4double Getcept() {
		return cept;
	}
	;

	inline G4PhysicsOrderedFreeVector GetUserDefinedEnergyHisto() {
		return UDefEnergyH;
	}
	;
	inline G4PhysicsOrderedFreeVector GetArbEnergyHisto() {
		return ArbEnergyH;
	}
	;

	G4double GenerateOne(G4ParticleDefinition*);
	G4double GetProbability (G4double);


private:
	void LinearInterpolation();
	void LogInterpolation();
	void ExpInterpolation();
	void SplineInterpolation();
	void CalculateCdgSpectrum();
	void CalculateBbodySpectrum();

	// The following methods generate energies according to the spectral
	// parameters defined above.
	void GenerateMonoEnergetic();
	void GenerateLinearEnergies(G4bool);
	void GeneratePowEnergies(G4bool);
	void GenerateBiasPowEnergies();
	void GenerateExpEnergies(G4bool);
	void GenerateGaussEnergies();
	void GenerateBremEnergies();
	void GenerateBbodyEnergies();
	void GenerateCdgEnergies();
	void GenUserHistEnergies();
	void GenEpnHistEnergies();
	void GenArbPointEnergies();
	// converts energy per nucleon to energy.
	void ConvertEPNToEnergy();


private:

	G4String EnergyDisType; // energy dis type Variable  - Mono,Lin,Exp,etc
	G4double weight; // particle weight
	G4double MonoEnergy; //Mono-energteic energy
	G4double SE; // Standard deviation for Gaussion distrbution in energy
	G4double Emin, Emax; // emin and emax
	G4double alpha, Ezero, Temp; // alpha (pow), E0 (exp) and Temp (bbody,brem)
	G4double biasalpha; // biased power index
	G4double grad, cept; // gradient and intercept for linear spectra
        G4double prob_norm; // normalisation factor use in calculate the probability 
        G4bool Biased; // true - biased to power-law
	G4bool EnergySpec; // true - energy spectra, false - momentum spectra
	G4bool DiffSpec; // true - differential spec, false integral spec
	G4bool ApplyRig; // false no rigidity cutoff, true then apply one
	G4double ERig; // energy of rigidity cutoff
	G4PhysicsOrderedFreeVector UDefEnergyH; // energy hist data
	G4PhysicsOrderedFreeVector IPDFEnergyH;
	G4bool IPDFEnergyExist, IPDFArbExist, Epnflag;
	G4PhysicsOrderedFreeVector ArbEnergyH; // Arb x,y histogram
	G4PhysicsOrderedFreeVector IPDFArbEnergyH; // IPDF for Arb
	G4PhysicsOrderedFreeVector EpnEnergyH;
	G4double CDGhist[3]; // cumulative histo for cdg
	G4double BBHist[10001], Bbody_x[10001];
	G4String IntType; // Interpolation type
	G4double Arb_grad[1024], Arb_cept[1024]; // grad and cept for 1024 segments
	G4double Arb_alpha[1024], Arb_Const[1024]; // alpha and constants
	G4double Arb_ezero[1024]; // ezero
	G4double ArbEmin, ArbEmax; // Emin and Emax for the whole arb distribution used primarily for debug.

	G4double particle_energy;
	G4ParticleDefinition* particle_definition;

	G4SPSRandomGenerator* eneRndm;

	// Verbosity
	G4int verbosityLevel;

	G4PhysicsOrderedFreeVector ZeroPhysVector; // for re-set only

	G4DataInterpolation *SplineInt[1024]; // holds Spline stuff required for sampling
	G4DataInterpolation *Splinetemp; // holds a temp Spline used for calculating area

};

#endif

