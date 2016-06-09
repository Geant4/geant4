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
// MODULE:        G4SPSEneDistribution.cc
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
#include "Randomize.hh"
//#include <cmath>

#include "G4SPSEneDistribution.hh"

G4SPSEneDistribution::G4SPSEneDistribution() {
	//
	// Initialise all variables
	particle_energy = 1.0 * MeV;

	EnergyDisType = "Mono";
	weight = 1.;
	MonoEnergy = 1 * MeV;
	Emin = 0.;
	Emax = 1.e30;
	alpha = 0.;
	biasalpha = 0.;
        prob_norm = 1.0;
	Ezero = 0.;
	SE = 0.;
	Temp = 0.;
	grad = 0.;
	cept = 0.;
	Biased = false; // not biased
	EnergySpec = true; // true - energy spectra, false - momentum spectra
	DiffSpec = true; // true - differential spec, false integral spec
	IntType = "NULL"; // Interpolation type
	IPDFEnergyExist = false;
	IPDFArbExist = false;

	ArbEmin = 0.;
	ArbEmax = 1.e30;

	verbosityLevel = 0;

}

G4SPSEneDistribution::~G4SPSEneDistribution() {
}

void G4SPSEneDistribution::SetEnergyDisType(G4String DisType) {
	EnergyDisType = DisType;
	if (EnergyDisType == "User") {
		UDefEnergyH = IPDFEnergyH = ZeroPhysVector;
		IPDFEnergyExist = false;
	} else if (EnergyDisType == "Arb") {
		ArbEnergyH = IPDFArbEnergyH = ZeroPhysVector;
		IPDFArbExist = false;
	} else if (EnergyDisType == "Epn") {
		UDefEnergyH = IPDFEnergyH = ZeroPhysVector;
		IPDFEnergyExist = false;
		EpnEnergyH = ZeroPhysVector;
	}
}

void G4SPSEneDistribution::SetEmin(G4double emi) {
	Emin = emi;
}

void G4SPSEneDistribution::SetEmax(G4double ema) {
	Emax = ema;
}

void G4SPSEneDistribution::SetMonoEnergy(G4double menergy) {
	MonoEnergy = menergy;
}

void G4SPSEneDistribution::SetBeamSigmaInE(G4double e) {
	SE = e;
}
void G4SPSEneDistribution::SetAlpha(G4double alp) {
	alpha = alp;
}

void G4SPSEneDistribution::SetBiasAlpha(G4double alp) {
	biasalpha = alp;
	Biased = true;
}

void G4SPSEneDistribution::SetTemp(G4double tem) {
	Temp = tem;
}

void G4SPSEneDistribution::SetEzero(G4double eze) {
	Ezero = eze;
}

void G4SPSEneDistribution::SetGradient(G4double gr) {
	grad = gr;
}

void G4SPSEneDistribution::SetInterCept(G4double c) {
	cept = c;
}

void G4SPSEneDistribution::UserEnergyHisto(G4ThreeVector input) {
	G4double ehi, val;
	ehi = input.x();
	val = input.y();
	if (verbosityLevel > 1) {
		G4cout << "In UserEnergyHisto" << G4endl;
		G4cout << " " << ehi << " " << val << G4endl;
	}
	UDefEnergyH.InsertValues(ehi, val);
	Emax = ehi;
}

void G4SPSEneDistribution::ArbEnergyHisto(G4ThreeVector input) {
	G4double ehi, val;
	ehi = input.x();
	val = input.y();
	if (verbosityLevel > 1) {
		G4cout << "In ArbEnergyHisto" << G4endl;
		G4cout << " " << ehi << " " << val << G4endl;
	}
	ArbEnergyH.InsertValues(ehi, val);
}

void G4SPSEneDistribution::ArbEnergyHistoFile(G4String filename) {
	std::ifstream infile(filename, std::ios::in);
	if (!infile)
		G4Exception("Unable to open the histo ASCII file");
	G4double ehi, val;
	while (infile >> ehi >> val) {
		ArbEnergyH.InsertValues(ehi, val);
	}
}

void G4SPSEneDistribution::EpnEnergyHisto(G4ThreeVector input) {
	G4double ehi, val;
	ehi = input.x();
	val = input.y();
	if (verbosityLevel > 1) {
		G4cout << "In EpnEnergyHisto" << G4endl;
		G4cout << " " << ehi << " " << val << G4endl;
	}
	EpnEnergyH.InsertValues(ehi, val);
	Emax = ehi;
	Epnflag = true;
}

void G4SPSEneDistribution::Calculate() {
	if (EnergyDisType == "Cdg")
		CalculateCdgSpectrum();
	else if (EnergyDisType == "Bbody")
		CalculateBbodySpectrum();
}

void G4SPSEneDistribution::CalculateCdgSpectrum() {
	// This uses the spectrum from The INTEGRAL Mass Model (TIMM)
	// to generate a Cosmic Diffuse X/gamma ray spectrum.
	G4double pfact[2] = { 8.5, 112 };
	G4double spind[2] = { 1.4, 2.3 };
	G4double ene_line[3] = { 1. * keV, 18. * keV, 1E6 * keV };
	G4int n_par;

	ene_line[0] = Emin;
	if (Emin < 18 * keV) {
		n_par = 2;
		ene_line[2] = Emax;
		if (Emax < 18 * keV) {
			n_par = 1;
			ene_line[1] = Emax;
		}
	} else {
		n_par = 1;
		pfact[0] = 112.;
		spind[0] = 2.3;
		ene_line[1] = Emax;
	}

	// Create a cumulative histogram.
	CDGhist[0] = 0.;
	G4double omalpha;
	G4int i = 0;

	while (i < n_par) {
		omalpha = 1. - spind[i];
		CDGhist[i + 1] = CDGhist[i] + (pfact[i] / omalpha) * (std::pow(
				ene_line[i + 1] / keV, omalpha) - std::pow(ene_line[i] / keV,
				omalpha));
		i++;
	}

	// Normalise histo and divide by 1000 to make MeV.
	i = 0;
	while (i < n_par) {
		CDGhist[i + 1] = CDGhist[i + 1] / CDGhist[n_par];
		//      G4cout << CDGhist[i] << CDGhist[n_par] << G4endl;
		i++;
	}
}

void G4SPSEneDistribution::CalculateBbodySpectrum() {
	// create bbody spectrum
	// Proved very hard to integrate indefinitely, so different
	// method. User inputs emin, emax and T. These are used to
	// create a 10,000 bin histogram.
	// Use photon density spectrum = 2 nu**2/c**2 * (std::exp(h nu/kT)-1)
	// = 2 E**2/h**2c**2 times the exponential
	G4double erange = Emax - Emin;
	G4double steps = erange / 10000.;
	G4double Bbody_y[10000];
	G4double k = 8.6181e-11; //Boltzmann const in MeV/K
	G4double h = 4.1362e-21; // Plancks const in MeV s
	G4double c = 3e8; // Speed of light
	G4double h2 = h * h;
	G4double c2 = c * c;
	G4int count = 0;
	G4double sum = 0.;
	BBHist[0] = 0.;
	while (count < 10000) {
		Bbody_x[count] = Emin + G4double(count * steps);
		Bbody_y[count] = (2. * std::pow(Bbody_x[count], 2.)) / (h2 * c2
				* (std::exp(Bbody_x[count] / (k * Temp)) - 1.));
		sum = sum + Bbody_y[count];
		BBHist[count + 1] = BBHist[count] + Bbody_y[count];
		count++;
	}

	Bbody_x[10000] = Emax;
	// Normalise cumulative histo.
	count = 0;
	while (count < 10001) {
		BBHist[count] = BBHist[count] / sum;
		count++;
	}
}

void G4SPSEneDistribution::InputEnergySpectra(G4bool value) {
	// Allows user to specifiy spectrum is momentum
	EnergySpec = value; // false if momentum
	if (verbosityLevel > 1)
		G4cout << "EnergySpec has value " << EnergySpec << G4endl;
}

void G4SPSEneDistribution::InputDifferentialSpectra(G4bool value) {
	// Allows user to specify integral or differential spectra
	DiffSpec = value; // true = differential, false = integral
	if (verbosityLevel > 1)
		G4cout << "Diffspec has value " << DiffSpec << G4endl;
}

void G4SPSEneDistribution::ArbInterpolate(G4String IType) {
	if (EnergyDisType != "Arb")
		G4cout << "Error: this is for arbitrary distributions" << G4endl;
	IntType = IType;
	ArbEmax = ArbEnergyH.GetMaxLowEdgeEnergy();
	ArbEmin = ArbEnergyH.GetMinLowEdgeEnergy();

	// Now interpolate points
	if (IntType == "Lin")
		LinearInterpolation();
	if (IntType == "Log")
		LogInterpolation();
	if (IntType == "Exp")
		ExpInterpolation();
	if (IntType == "Spline")
		SplineInterpolation();
}

void G4SPSEneDistribution::LinearInterpolation() {
	// Method to do linear interpolation on the Arb points
	// Calculate equation of each line segment, max 1024.
	// Calculate Area under each segment
	// Create a cumulative array which is then normalised Arb_Cum_Area

	G4double Area_seg[1024]; // Stores area under each segment
	G4double sum = 0., Arb_x[1024], Arb_y[1024], Arb_Cum_Area[1024];
	G4int i, count;
	G4int maxi = ArbEnergyH.GetVectorLength();
	for (i = 0; i < maxi; i++) {
		Arb_x[i] = ArbEnergyH.GetLowEdgeEnergy(size_t(i));
		Arb_y[i] = ArbEnergyH(size_t(i));
	}
	// Points are now in x,y arrays. If the spectrum is integral it has to be
	// made differential and if momentum it has to be made energy.
	if (DiffSpec == false) {
		// Converts integral point-wise spectra to Differential
		for (count = 0; count < maxi - 1; count++) {
			Arb_y[count] = (Arb_y[count] - Arb_y[count + 1])
					/ (Arb_x[count + 1] - Arb_x[count]);
		}
		maxi--;
	}
	//
	if (EnergySpec == false) {
		// change currently stored values (emin etc) which are actually momenta
		// to energies.
		if (particle_definition == NULL)
			G4cout << "Error: particle not defined" << G4endl;
		else {
			// Apply Energy**2 = p**2c**2 + m0**2c**4
			// p should be entered as E/c i.e. without the division by c
			// being done - energy equivalent.
			G4double mass = particle_definition->GetPDGMass();
			// convert point to energy unit and its value to per energy unit
			G4double total_energy;
			for (count = 0; count < maxi; count++) {
				total_energy = std::sqrt((Arb_x[count] * Arb_x[count]) + (mass
						* mass)); // total energy

				Arb_y[count] = Arb_y[count] * Arb_x[count] / total_energy;
				Arb_x[count] = total_energy - mass; // kinetic energy
			}
		}
	}
	//
	i = 1;
	Arb_grad[0] = 0.;
	Arb_cept[0] = 0.;
	Area_seg[0] = 0.;
	Arb_Cum_Area[0] = 0.;
	while (i < maxi) {
		// calc gradient and intercept for each segment
		Arb_grad[i] = (Arb_y[i] - Arb_y[i - 1]) / (Arb_x[i] - Arb_x[i - 1]);
		if (verbosityLevel == 2)
			G4cout << Arb_grad[i] << G4endl;
		if (Arb_grad[i] > 0.) {
			if (verbosityLevel == 2)
				G4cout << "Arb_grad is positive" << G4endl;
			Arb_cept[i] = Arb_y[i] - (Arb_grad[i] * Arb_x[i]);
		} else if (Arb_grad[i] < 0.) {
			if (verbosityLevel == 2)
				G4cout << "Arb_grad is negative" << G4endl;
			Arb_cept[i] = Arb_y[i] + (-Arb_grad[i] * Arb_x[i]);
		} else {
			if (verbosityLevel == 2)
				G4cout << "Arb_grad is 0." << G4endl;
			Arb_cept[i] = Arb_y[i];
		}

		Area_seg[i] = ((Arb_grad[i] / 2) * (Arb_x[i] * Arb_x[i] - Arb_x[i - 1]
				* Arb_x[i - 1]) + Arb_cept[i] * (Arb_x[i] - Arb_x[i - 1]));
		Arb_Cum_Area[i] = Arb_Cum_Area[i - 1] + Area_seg[i];
		sum = sum + Area_seg[i];
		if (verbosityLevel == 2)
			G4cout << Arb_x[i] << Arb_y[i] << Area_seg[i] << sum << Arb_grad[i]
					<< G4endl;
		i++;
	}

	i = 0;
	while (i < maxi) {
		Arb_Cum_Area[i] = Arb_Cum_Area[i] / sum; // normalisation
		IPDFArbEnergyH.InsertValues(Arb_x[i], Arb_Cum_Area[i]);
		i++;
	}

	// now scale the ArbEnergyH, needed by Probability()
	ArbEnergyH.ScaleVector(1., 1./sum);

	if (verbosityLevel >= 1) {
		G4cout << "Leaving LinearInterpolation" << G4endl;
		ArbEnergyH.DumpValues();
		IPDFArbEnergyH.DumpValues();
	}
}

void G4SPSEneDistribution::LogInterpolation() {
	// Interpolation based on Logarithmic equations
	// Generate equations of line segments
	// y = Ax**alpha => log y = alpha*logx + logA
	// Find area under line segments
	// create normalised, cumulative array Arb_Cum_Area
	G4double Area_seg[1024]; // Stores area under each segment
	G4double sum = 0., Arb_x[1024], Arb_y[1024], Arb_Cum_Area[1024];
	G4int i, count;
	G4int maxi = ArbEnergyH.GetVectorLength();
	for (i = 0; i < maxi; i++) {
		Arb_x[i] = ArbEnergyH.GetLowEdgeEnergy(size_t(i));
		Arb_y[i] = ArbEnergyH(size_t(i));
	}
	// Points are now in x,y arrays. If the spectrum is integral it has to be
	// made differential and if momentum it has to be made energy.
	if (DiffSpec == false) {
		// Converts integral point-wise spectra to Differential
		for (count = 0; count < maxi - 1; count++) {
			Arb_y[count] = (Arb_y[count] - Arb_y[count + 1])
					/ (Arb_x[count + 1] - Arb_x[count]);
		}
		maxi--;
	}
	//
	if (EnergySpec == false) {
		// change currently stored values (emin etc) which are actually momenta
		// to energies.
		if (particle_definition == NULL)
			G4cout << "Error: particle not defined" << G4endl;
		else {
			// Apply Energy**2 = p**2c**2 + m0**2c**4
			// p should be entered as E/c i.e. without the division by c
			// being done - energy equivalent.
			G4double mass = particle_definition->GetPDGMass();
			// convert point to energy unit and its value to per energy unit
			G4double total_energy;
			for (count = 0; count < maxi; count++) {
				total_energy = std::sqrt((Arb_x[count] * Arb_x[count]) + (mass
						* mass)); // total energy

				Arb_y[count] = Arb_y[count] * Arb_x[count] / total_energy;
				Arb_x[count] = total_energy - mass; // kinetic energy
			}
		}
	}
	//
	i = 1;
	Arb_alpha[0] = 0.;
	Arb_Const[0] = 0.;
	Area_seg[0] = 0.;
	Arb_Cum_Area[0] = 0.;
	if (Arb_x[0] <= 0. || Arb_y[0] <= 0.) {
		G4cout << "You should not use log interpolation with points <= 0."
				<< G4endl;
		G4cout << "These will be changed to 1e-20, which may cause problems"
				<< G4endl;
		if (Arb_x[0] <= 0.)
			Arb_x[0] = 1e-20;
		if (Arb_y[0] <= 0.)
			Arb_y[0] = 1e-20;
	}

	G4double alp;
	while (i < maxi) {
		// Incase points are negative or zero
		if (Arb_x[i] <= 0. || Arb_y[i] <= 0.) {
			G4cout << "You should not use log interpolation with points <= 0."
					<< G4endl;
			G4cout
					<< "These will be changed to 1e-20, which may cause problems"
					<< G4endl;
			if (Arb_x[i] <= 0.)
				Arb_x[i] = 1e-20;
			if (Arb_y[i] <= 0.)
				Arb_y[i] = 1e-20;
		}

		Arb_alpha[i] = (std::log10(Arb_y[i]) - std::log10(Arb_y[i - 1]))
				/ (std::log10(Arb_x[i]) - std::log10(Arb_x[i - 1]));
		Arb_Const[i] = Arb_y[i] / (std::pow(Arb_x[i], Arb_alpha[i]));
		alp = Arb_alpha[i] + 1;
		if (alp == 0.) {
		  Area_seg[i] =	Arb_Const[i] * (std::log(Arb_x[i]) - std::log(Arb_x[i - 1])); 
		} else {
		  Area_seg[i] = (Arb_Const[i] / alp) * (std::pow(Arb_x[i], alp)
				- std::pow(Arb_x[i - 1], alp));
		}
		sum = sum + Area_seg[i];
		Arb_Cum_Area[i] = Arb_Cum_Area[i - 1] + Area_seg[i];
		if (verbosityLevel == 2)
			G4cout << Arb_alpha[i] << Arb_Const[i] << Area_seg[i] << G4endl;
		i++;
	}

	i = 0;
	while (i < maxi) {
		Arb_Cum_Area[i] = Arb_Cum_Area[i] / sum;
		IPDFArbEnergyH.InsertValues(Arb_x[i], Arb_Cum_Area[i]);
		i++;
	}

	// now scale the ArbEnergyH, needed by Probability()
	ArbEnergyH.ScaleVector(1., 1./sum);

	if (verbosityLevel >= 1)
		G4cout << "Leaving LogInterpolation " << G4endl;
}

void G4SPSEneDistribution::ExpInterpolation() {
	// Interpolation based on Exponential equations
	// Generate equations of line segments
	// y = Ae**-(x/e0) => ln y = -x/e0 + lnA
	// Find area under line segments
	// create normalised, cumulative array Arb_Cum_Area
	G4double Area_seg[1024]; // Stores area under each segment
	G4double sum = 0., Arb_x[1024], Arb_y[1024], Arb_Cum_Area[1024];
	G4int i, count;
	G4int maxi = ArbEnergyH.GetVectorLength();
	for (i = 0; i < maxi; i++) {
		Arb_x[i] = ArbEnergyH.GetLowEdgeEnergy(size_t(i));
		Arb_y[i] = ArbEnergyH(size_t(i));
	}
	// Points are now in x,y arrays. If the spectrum is integral it has to be
	// made differential and if momentum it has to be made energy.
	if (DiffSpec == false) {
		// Converts integral point-wise spectra to Differential
		for (count = 0; count < maxi - 1; count++) {
			Arb_y[count] = (Arb_y[count] - Arb_y[count + 1])
					/ (Arb_x[count + 1] - Arb_x[count]);
		}
		maxi--;
	}
	//
	if (EnergySpec == false) {
		// change currently stored values (emin etc) which are actually momenta
		// to energies.
		if (particle_definition == NULL)
			G4cout << "Error: particle not defined" << G4endl;
		else {
			// Apply Energy**2 = p**2c**2 + m0**2c**4
			// p should be entered as E/c i.e. without the division by c
			// being done - energy equivalent.
			G4double mass = particle_definition->GetPDGMass();
			// convert point to energy unit and its value to per energy unit
			G4double total_energy;
			for (count = 0; count < maxi; count++) {
				total_energy = std::sqrt((Arb_x[count] * Arb_x[count]) + (mass
						* mass)); // total energy

				Arb_y[count] = Arb_y[count] * Arb_x[count] / total_energy;
				Arb_x[count] = total_energy - mass; // kinetic energy
			}
		}
	}
	//
	i = 1;
	Arb_ezero[0] = 0.;
	Arb_Const[0] = 0.;
	Area_seg[0] = 0.;
	Arb_Cum_Area[0] = 0.;
	while (i < maxi) {
		G4double test = std::log(Arb_y[i]) - std::log(Arb_y[i - 1]);
		if (test > 0. || test < 0.) {
			Arb_ezero[i] = -(Arb_x[i] - Arb_x[i - 1]) / (std::log(Arb_y[i])
					- std::log(Arb_y[i - 1]));
			Arb_Const[i] = Arb_y[i] / (std::exp(-Arb_x[i] / Arb_ezero[i]));
			Area_seg[i] = -(Arb_Const[i] * Arb_ezero[i]) * (std::exp(-Arb_x[i]
					/ Arb_ezero[i]) - std::exp(-Arb_x[i - 1] / Arb_ezero[i]));
		} else {
			G4cout << "Flat line segment: problem" << G4endl;
			Arb_ezero[i] = 0.;
			Arb_Const[i] = 0.;
			Area_seg[i] = 0.;
		}
		sum = sum + Area_seg[i];
		Arb_Cum_Area[i] = Arb_Cum_Area[i - 1] + Area_seg[i];
		if (verbosityLevel == 2)
			G4cout << Arb_ezero[i] << Arb_Const[i] << Area_seg[i] << G4endl;
		i++;
	}

	i = 0;
	while (i < maxi) {
		Arb_Cum_Area[i] = Arb_Cum_Area[i] / sum;
		IPDFArbEnergyH.InsertValues(Arb_x[i], Arb_Cum_Area[i]);
		i++;
	}

	// now scale the ArbEnergyH, needed by Probability()
	ArbEnergyH.ScaleVector(1., 1./sum);

	if (verbosityLevel >= 1)
		G4cout << "Leaving ExpInterpolation " << G4endl;
}

void G4SPSEneDistribution::SplineInterpolation() {
	// Interpolation using Splines.
	// Create Normalised arrays, make x 0->1 and y hold
	// the function (Energy)
        // 
        // Current method based on the above will not work in all cases. 
        // new method is implemented below.
  
	G4double sum, Arb_x[1024], Arb_y[1024], Arb_Cum_Area[1024];
	G4int i, count;

	G4int maxi = ArbEnergyH.GetVectorLength();
	for (i = 0; i < maxi; i++) {
		Arb_x[i] = ArbEnergyH.GetLowEdgeEnergy(size_t(i));
		Arb_y[i] = ArbEnergyH(size_t(i));
	}
	// Points are now in x,y arrays. If the spectrum is integral it has to be
	// made differential and if momentum it has to be made energy.
	if (DiffSpec == false) {
		// Converts integral point-wise spectra to Differential
		for (count = 0; count < maxi - 1; count++) {
			Arb_y[count] = (Arb_y[count] - Arb_y[count + 1])
					/ (Arb_x[count + 1] - Arb_x[count]);
		}
		maxi--;
	}
	//
	if (EnergySpec == false) {
		// change currently stored values (emin etc) which are actually momenta
		// to energies.
		if (particle_definition == NULL)
			G4cout << "Error: particle not defined" << G4endl;
		else {
			// Apply Energy**2 = p**2c**2 + m0**2c**4
			// p should be entered as E/c i.e. without the division by c
			// being done - energy equivalent.
			G4double mass = particle_definition->GetPDGMass();
			// convert point to energy unit and its value to per energy unit
			G4double total_energy;
			for (count = 0; count < maxi; count++) {
				total_energy = std::sqrt((Arb_x[count] * Arb_x[count]) + (mass
						* mass)); // total energy

				Arb_y[count] = Arb_y[count] * Arb_x[count] / total_energy;
				Arb_x[count] = total_energy - mass; // kinetic energy
			}
		}
	}

	//
	i = 1;
	Arb_Cum_Area[0] = 0.;
	sum = 0.;
	Splinetemp = new G4DataInterpolation(Arb_x, Arb_y, maxi, 0., 0.);
	G4double ei[101],prob[101];
	while (i < maxi) {
	  // 100 step per segment for the integration of area
	  G4double de = (Arb_x[i] - Arb_x[i - 1])/100.;
	  G4double area = 0.;

	  for (count = 0; count < 101; count++) {
	    ei[count] = Arb_x[i - 1] + de*count ;
	    prob[count] =  Splinetemp->CubicSplineInterpolation(ei[count]);
	    if (prob[count] < 0.) { 
	      G4cout <<   "Warning: G4DataInterpolation returns value < 0  " << prob[count] <<" "<<ei[count]<< G4endl;
	      G4Exception("         Please use an alternative method, e.g. Lin, for interpolation");
	    }
	    area += prob[count]*de;
	  }
	  Arb_Cum_Area[i] = Arb_Cum_Area[i - 1] + area;
	  sum += area; 

	  prob[0] = prob[0]/(area/de);
	  for (count = 1; count < 100; count++)
	    prob[count] = prob[count-1] + prob[count]/(area/de);

	  SplineInt[i] = new G4DataInterpolation(prob, ei, 101, 0., 0.);
	  // note i start from 1!
	  i++;
	}
	i = 0;
	while (i < maxi) {
		Arb_Cum_Area[i] = Arb_Cum_Area[i] / sum; // normalisation
		IPDFArbEnergyH.InsertValues(Arb_x[i], Arb_Cum_Area[i]);
		i++;
	}
	// now scale the ArbEnergyH, needed by Probability()
	ArbEnergyH.ScaleVector(1., 1./sum);

	if (verbosityLevel > 0)
	  G4cout << "Leaving SplineInterpolation " << G4endl;
}

void G4SPSEneDistribution::GenerateMonoEnergetic() {
	// Method to generate MonoEnergetic particles.
	particle_energy = MonoEnergy;
}

void G4SPSEneDistribution::GenerateGaussEnergies() {
	// Method to generate Gaussian particles.
	particle_energy = G4RandGauss::shoot(MonoEnergy,SE);
	if (particle_energy < 0) particle_energy = 0.;
}

void G4SPSEneDistribution::GenerateLinearEnergies(G4bool bArb = false) {
	G4double rndm;
	G4double emaxsq = std::pow(Emax, 2.); //Emax squared
	G4double eminsq = std::pow(Emin, 2.); //Emin squared
	G4double intersq = std::pow(cept, 2.); //cept squared

	if (bArb)
		rndm = G4UniformRand();
	else
		rndm = eneRndm->GenRandEnergy();

	G4double bracket = ((grad / 2.) * (emaxsq - eminsq) + cept * (Emax - Emin));
	bracket = bracket * rndm;
	bracket = bracket + (grad / 2.) * eminsq + cept * Emin;
	// Now have a quad of form m/2 E**2 + cE - bracket = 0
	bracket = -bracket;
	//  G4cout << "BRACKET" << bracket << G4endl;
	if (grad != 0.) {
		G4double sqbrack = (intersq - 4 * (grad / 2.) * (bracket));
		//      G4cout << "SQBRACK" << sqbrack << G4endl;
		sqbrack = std::sqrt(sqbrack);
		G4double root1 = -cept + sqbrack;
		root1 = root1 / (2. * (grad / 2.));

		G4double root2 = -cept - sqbrack;
		root2 = root2 / (2. * (grad / 2.));

		//      G4cout << root1 << " roots " << root2 << G4endl;

		if (root1 > Emin && root1 < Emax)
			particle_energy = root1;
		if (root2 > Emin && root2 < Emax)
			particle_energy = root2;
	} else if (grad == 0.)
		// have equation of form cE - bracket =0
		particle_energy = bracket / cept;

	if (particle_energy < 0.)
		particle_energy = -particle_energy;

	if (verbosityLevel >= 1)
		G4cout << "Energy is " << particle_energy << G4endl;
}

void G4SPSEneDistribution::GeneratePowEnergies(G4bool bArb = false) {
	// Method to generate particle energies distributed as
	// a power-law

	G4double rndm;
	G4double emina, emaxa;

	emina = std::pow(Emin, alpha + 1);
	emaxa = std::pow(Emax, alpha + 1);

	if (bArb)
		rndm = G4UniformRand();
	else
		rndm = eneRndm->GenRandEnergy();

	if (alpha != -1.) {
		particle_energy = ((rndm * (emaxa - emina)) + emina);
		particle_energy = std::pow(particle_energy, (1. / (alpha + 1.)));
	} else {
		particle_energy = (std::log(Emin) + rndm * (std::log(Emax) - std::log(
				Emin)));
		particle_energy = std::exp(particle_energy);
	}
	if (verbosityLevel >= 1)
		G4cout << "Energy is " << particle_energy << G4endl;
}

void G4SPSEneDistribution::GenerateBiasPowEnergies() {
	// Method to generate particle energies distributed as
	// in biased power-law and calculate its weight

        G4double rndm;
	G4double emina, emaxa, emin, emax;

	G4double normal = 1. ;

	emin = Emin;
	emax = Emax;
	//	if (EnergyDisType == "Arb") { 
	//  emin = ArbEmin;
	//  emax = ArbEmax;
	//}

	rndm = eneRndm->GenRandEnergy();

	if (biasalpha != -1.) {
	        emina = std::pow(emin, biasalpha + 1);
	        emaxa = std::pow(emax, biasalpha + 1);
		particle_energy = ((rndm * (emaxa - emina)) + emina);
		particle_energy = std::pow(particle_energy, (1. / (biasalpha + 1.)));
		normal = 1./(1+biasalpha) * (emaxa - emina);
	} else {
		particle_energy = (std::log(emin) + rndm * (std::log(emax) - std::log(
				emin)));
		particle_energy = std::exp(particle_energy);
		normal = std::log(emax) - std::log(emin) ;
	}
	weight = GetProbability(particle_energy) / (std::pow(particle_energy,biasalpha)/normal);

	if (verbosityLevel >= 1)
		G4cout << "Energy is " << particle_energy << G4endl;
}

void G4SPSEneDistribution::GenerateExpEnergies(G4bool bArb = false) {
	// Method to generate particle energies distributed according
	// to an exponential curve.
	G4double rndm;

	if (bArb)
		rndm = G4UniformRand();
	else
		rndm = eneRndm->GenRandEnergy();

	particle_energy = -Ezero * (std::log(rndm * (std::exp(-Emax / Ezero)
			- std::exp(-Emin / Ezero)) + std::exp(-Emin / Ezero)));
	if (verbosityLevel >= 1)
		G4cout << "Energy is " << particle_energy << G4endl;
}

void G4SPSEneDistribution::GenerateBremEnergies() {
	// Method to generate particle energies distributed according
	// to a Bremstrahlung equation of
	// form I = const*((kT)**1/2)*E*(e**(-E/kT))

	G4double rndm;
	rndm = eneRndm->GenRandEnergy();
	G4double expmax, expmin, k;

	k = 8.6181e-11; // Boltzmann's const in MeV/K
	G4double ksq = std::pow(k, 2.); // k squared
	G4double Tsq = std::pow(Temp, 2.); // Temp squared

	expmax = std::exp(-Emax / (k * Temp));
	expmin = std::exp(-Emin / (k * Temp));

	// If either expmax or expmin are zero then this will cause problems
	// Most probably this will be because T is too low or E is too high

	if (expmax == 0.)
		G4cout << "*****EXPMAX=0. Choose different E's or Temp" << G4endl;
	if (expmin == 0.)
		G4cout << "*****EXPMIN=0. Choose different E's or Temp" << G4endl;

	G4double tempvar = rndm * ((-k) * Temp * (Emax * expmax - Emin * expmin)
			- (ksq * Tsq * (expmax - expmin)));

	G4double bigc = (tempvar - k * Temp * Emin * expmin - ksq * Tsq * expmin)
			/ (-k * Temp);

	// This gives an equation of form: Ee(-E/kT) + kTe(-E/kT) - C =0
	// Solve this iteratively, step from Emin to Emax in 1000 steps
	// and take the best solution.

	G4double erange = Emax - Emin;
	G4double steps = erange / 1000.;
	G4int i;
	G4double etest, diff, err;

	err = 100000.;

	for (i = 1; i < 1000; i++) {
		etest = Emin + (i - 1) * steps;

		diff = etest * (std::exp(-etest / (k * Temp))) + k * Temp * (std::exp(
				-etest / (k * Temp))) - bigc;

		if (diff < 0.)
			diff = -diff;

		if (diff < err) {
			err = diff;
			particle_energy = etest;
		}
	}
	if (verbosityLevel >= 1)
		G4cout << "Energy is " << particle_energy << G4endl;
}

void G4SPSEneDistribution::GenerateBbodyEnergies() {
	// BBody_x holds Energies, and BBHist holds the cumulative histo.
	// binary search to find correct bin then lin interpolation.
	// Use the earlier defined histogram + RandGeneral method to generate
	// random numbers following the histos distribution.
	G4double rndm;
	G4int nabove, nbelow = 0, middle;
	nabove = 10001;
	rndm = eneRndm->GenRandEnergy();

	// Binary search to find bin that rndm is in
	while (nabove - nbelow > 1) {
		middle = (nabove + nbelow) / 2;
		if (rndm == BBHist[middle])
			break;
		if (rndm < BBHist[middle])
			nabove = middle;
		else
			nbelow = middle;
	}

	// Now interpolate in that bin to find the correct output value.
	G4double x1, x2, y1, y2, m, q;
	x1 = Bbody_x[nbelow];
	x2 = Bbody_x[nbelow + 1];
	y1 = BBHist[nbelow];
	y2 = BBHist[nbelow + 1];
	m = (y2 - y1) / (x2 - x1);
	q = y1 - m * x1;

	particle_energy = (rndm - q) / m;

	if (verbosityLevel >= 1) {
		G4cout << "Energy is " << particle_energy << G4endl;
	}
}

void G4SPSEneDistribution::GenerateCdgEnergies() {
	// Gen random numbers, compare with values in cumhist
	// to find appropriate part of spectrum and then
	// generate energy in the usual inversion way.
	//  G4double pfact[2] = {8.5, 112};
	// G4double spind[2] = {1.4, 2.3};
	// G4double ene_line[3] = {1., 18., 1E6};
	G4double rndm, rndm2;
	G4double ene_line[3];
	G4double omalpha[2];
	if (Emin < 18 * keV && Emax < 18 * keV) {
		omalpha[0] = 1. - 1.4;
		ene_line[0] = Emin;
		ene_line[1] = Emax;
	}
	if (Emin < 18 * keV && Emax > 18 * keV) {
		omalpha[0] = 1. - 1.4;
		omalpha[1] = 1. - 2.3;
		ene_line[0] = Emin;
		ene_line[1] = 18. * keV;
		ene_line[2] = Emax;
	}
	if (Emin > 18 * keV) {
		omalpha[0] = 1. - 2.3;
		ene_line[0] = Emin;
		ene_line[1] = Emax;
	}
	rndm = eneRndm->GenRandEnergy();
	rndm2 = eneRndm->GenRandEnergy();

	G4int i = 0;
	while (rndm >= CDGhist[i]) {
		i++;
	}
	// Generate final energy.
	particle_energy = (std::pow(ene_line[i - 1], omalpha[i - 1]) + (std::pow(
			ene_line[i], omalpha[i - 1]) - std::pow(ene_line[i - 1], omalpha[i
			- 1])) * rndm2);
	particle_energy = std::pow(particle_energy, (1. / omalpha[i - 1]));

	if (verbosityLevel >= 1)
		G4cout << "Energy is " << particle_energy << G4endl;
}

void G4SPSEneDistribution::GenUserHistEnergies() {
	// Histograms are DIFFERENTIAL.
	//  G4cout << "In GenUserHistEnergies " << G4endl;
	if (IPDFEnergyExist == false) {
		G4int ii;
		G4int maxbin = G4int(UDefEnergyH.GetVectorLength());
		G4double bins[1024], vals[1024], sum;
		sum = 0.;

		if ((EnergySpec == false) && (particle_definition == NULL))
			G4cout << "Error: particle definition is NULL" << G4endl;

		if (maxbin > 1024) {
			G4cout << "Maxbin > 1024" << G4endl;
			G4cout << "Setting maxbin to 1024, other bins are lost" << G4endl;
		}

		if (DiffSpec == false)
			G4cout << "Histograms are Differential!!! " << G4endl;
		else {
			bins[0] = UDefEnergyH.GetLowEdgeEnergy(size_t(0));
			vals[0] = UDefEnergyH(size_t(0));
			sum = vals[0];
			for (ii = 1; ii < maxbin; ii++) {
				bins[ii] = UDefEnergyH.GetLowEdgeEnergy(size_t(ii));
				vals[ii] = UDefEnergyH(size_t(ii)) + vals[ii - 1];
				sum = sum + UDefEnergyH(size_t(ii));
			}
		}

		if (EnergySpec == false) {
			G4double mass = particle_definition->GetPDGMass();
			// multiply the function (vals) up by the bin width
			// to make the function counts/s (i.e. get rid of momentum
			// dependence).
			for (ii = 1; ii < maxbin; ii++) {
				vals[ii] = vals[ii] * (bins[ii] - bins[ii - 1]);
			}
			// Put energy bins into new histo, plus divide by energy bin width
			// to make evals counts/s/energy
			for (ii = 0; ii < maxbin; ii++) {
				bins[ii] = std::sqrt((bins[ii] * bins[ii]) + (mass * mass))
						- mass; //kinetic energy
			}
			for (ii = 1; ii < maxbin; ii++) {
				vals[ii] = vals[ii] / (bins[ii] - bins[ii - 1]);
			}
			sum = vals[maxbin - 1];
			vals[0] = 0.;
		}
		for (ii = 0; ii < maxbin; ii++) {
			vals[ii] = vals[ii] / sum;
			IPDFEnergyH.InsertValues(bins[ii], vals[ii]);
		}

		// Make IPDFEnergyExist = true
		IPDFEnergyExist = true;
		if (verbosityLevel > 1)
			IPDFEnergyH.DumpValues();
	}

	// IPDF has been create so carry on
	G4double rndm = eneRndm->GenRandEnergy();
	particle_energy = IPDFEnergyH.GetEnergy(rndm);

	if (verbosityLevel >= 1)
		G4cout << "Energy is " << particle_energy << G4endl;
}

void G4SPSEneDistribution::GenArbPointEnergies() {
	if (verbosityLevel > 0)
		G4cout << "In GenArbPointEnergies" << G4endl;
	G4double rndm;
	rndm = eneRndm->GenRandEnergy();
	//      IPDFArbEnergyH.DumpValues();
	// Find the Bin
	// have x, y, no of points, and cumulative area distribution
	G4int nabove, nbelow = 0, middle;
	nabove = IPDFArbEnergyH.GetVectorLength();
	//      G4cout << nabove << G4endl;
	// Binary search to find bin that rndm is in
	while (nabove - nbelow > 1) {
	  middle = (nabove + nbelow) / 2;
	  if (rndm == IPDFArbEnergyH(size_t(middle)))
	    break;
	  if (rndm < IPDFArbEnergyH(size_t(middle)))
	    nabove = middle;
	  else
	    nbelow = middle;
	}
	if (IntType == "Lin") {
	  Emax = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(nbelow + 1));
	  Emin = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(nbelow));
	  grad = Arb_grad[nbelow + 1];
	  cept = Arb_cept[nbelow + 1];
	  //	  G4cout << rndm << " " << Emax << " " << Emin << " " << grad << " " << cept << G4endl;
	  GenerateLinearEnergies(true);
	} else if (IntType == "Log") {
	  Emax = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(nbelow + 1));
	  Emin = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(nbelow));
	  alpha = Arb_alpha[nbelow + 1];
	  //	  G4cout << rndm << " " << Emax << " " << Emin << " " << alpha << G4endl;
	  GeneratePowEnergies(true);
	} else if (IntType == "Exp") {
	  Emax = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(nbelow + 1));
	  Emin = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(nbelow));
	  Ezero = Arb_ezero[nbelow + 1];
	  //	  G4cout << rndm << " " << Emax << " " << Emin << " " << Ezero << G4endl;
	  GenerateExpEnergies(true);
	} else if (IntType == "Spline") {
	  Emax = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(nbelow + 1));
	  Emin = IPDFArbEnergyH.GetLowEdgeEnergy(size_t(nbelow));
	  particle_energy = -1e100;
	  rndm = eneRndm->GenRandEnergy();
	  while (particle_energy < Emin || particle_energy > Emax) {
	    particle_energy = SplineInt[nbelow+1]->CubicSplineInterpolation(rndm);
	    rndm = eneRndm->GenRandEnergy();
	  }
	  if (verbosityLevel >= 1)
	    G4cout << "Energy is " << particle_energy << G4endl;
	} else
		G4cout << "Error: IntType unknown type" << G4endl;
}

void G4SPSEneDistribution::GenEpnHistEnergies() {
	//  G4cout << "In GenEpnHistEnergies " << Epnflag << G4endl;

	// Firstly convert to energy if not already done.
	if (Epnflag == true)
	// epnflag = true means spectrum is epn, false means e.
	{
		// convert to energy by multiplying by A number
		ConvertEPNToEnergy();
		// EpnEnergyH will be replace by UDefEnergyH.
		//      UDefEnergyH.DumpValues();
	}

	//  G4cout << "Creating IPDFEnergy if not already done so" << G4endl;
	if (IPDFEnergyExist == false) {
		// IPDF has not been created, so create it
		G4double bins[1024], vals[1024], sum;
		G4int ii;
		G4int maxbin = G4int(UDefEnergyH.GetVectorLength());
		bins[0] = UDefEnergyH.GetLowEdgeEnergy(size_t(0));
		vals[0] = UDefEnergyH(size_t(0));
		sum = vals[0];
		for (ii = 1; ii < maxbin; ii++) {
			bins[ii] = UDefEnergyH.GetLowEdgeEnergy(size_t(ii));
			vals[ii] = UDefEnergyH(size_t(ii)) + vals[ii - 1];
			sum = sum + UDefEnergyH(size_t(ii));
		}

		for (ii = 0; ii < maxbin; ii++) {
			vals[ii] = vals[ii] / sum;
			IPDFEnergyH.InsertValues(bins[ii], vals[ii]);
		}
		// Make IPDFEpnExist = true
		IPDFEnergyExist = true;
	}
	//  IPDFEnergyH.DumpValues();
	// IPDF has been create so carry on
	G4double rndm = eneRndm->GenRandEnergy();
	particle_energy = IPDFEnergyH.GetEnergy(rndm);

	if (verbosityLevel >= 1)
		G4cout << "Energy is " << particle_energy << G4endl;
}

void G4SPSEneDistribution::ConvertEPNToEnergy() {
	// Use this before particle generation to convert  the
	// currently stored histogram from energy/nucleon
	// to energy.
	//  G4cout << "In ConvertEpntoEnergy " << G4endl;
	if (particle_definition == NULL)
		G4cout << "Error: particle not defined" << G4endl;
	else {
		// Need to multiply histogram by the number of nucleons.
		// Baryon Number looks to hold the no. of nucleons.
		G4int Bary = particle_definition->GetBaryonNumber();
		//      G4cout << "Baryon No. " << Bary << G4endl;
		// Change values in histogram, Read it out, delete it, re-create it
		G4int count, maxcount;
		maxcount = G4int(EpnEnergyH.GetVectorLength());
		//      G4cout << maxcount << G4endl;
		G4double ebins[1024], evals[1024];
		if (maxcount > 1024) {
			G4cout << "Histogram contains more than 1024 bins!" << G4endl;
			G4cout << "Those above 1024 will be ignored" << G4endl;
			maxcount = 1024;
		}
		if (maxcount < 1) {
			G4cout << "Histogram contains less than 1 bin!" << G4endl;
			G4cout << "Redefine the histogram" << G4endl;
			return;
		}
		for (count = 0; count < maxcount; count++) {
			// Read out
			ebins[count] = EpnEnergyH.GetLowEdgeEnergy(size_t(count));
			evals[count] = EpnEnergyH(size_t(count));
		}

		// Multiply the channels by the nucleon number to give energies
		for (count = 0; count < maxcount; count++) {
			ebins[count] = ebins[count] * Bary;
		}

		// Set Emin and Emax
		Emin = ebins[0];
		if (maxcount > 1)
			Emax = ebins[maxcount - 1];
		else
			Emax = ebins[0];
		// Put energy bins into new histogram - UDefEnergyH.
		for (count = 0; count < maxcount; count++) {
			UDefEnergyH.InsertValues(ebins[count], evals[count]);
		}
		Epnflag = false; //so that you dont repeat this method.
	}
}

//
void G4SPSEneDistribution::ReSetHist(G4String atype) {
	if (atype == "energy") {
		UDefEnergyH = IPDFEnergyH = ZeroPhysVector;
		IPDFEnergyExist = false;
		Emin = 0.;
		Emax = 1e30;
	} else if (atype == "arb") {
		ArbEnergyH = IPDFArbEnergyH = ZeroPhysVector;
		IPDFArbExist = false;
	} else if (atype == "epn") {
		UDefEnergyH = IPDFEnergyH = ZeroPhysVector;
		IPDFEnergyExist = false;
		EpnEnergyH = ZeroPhysVector;
	} else {
		G4cout << "Error, histtype not accepted " << G4endl;
	}
}

G4double G4SPSEneDistribution::GenerateOne(G4ParticleDefinition* a) {
	particle_definition = a;
	particle_energy = -1.;

	while ((EnergyDisType == "Arb") ? (particle_energy < ArbEmin
			|| particle_energy > ArbEmax) : (particle_energy < Emin
			|| particle_energy > Emax)) {
		if (Biased) {
			GenerateBiasPowEnergies();
		} else {
			if (EnergyDisType == "Mono")
				GenerateMonoEnergetic();
			else if (EnergyDisType == "Lin")
				GenerateLinearEnergies();
			else if (EnergyDisType == "Pow")
				GeneratePowEnergies();
			else if (EnergyDisType == "Exp")
				GenerateExpEnergies();
			else if (EnergyDisType == "Gauss")
				GenerateGaussEnergies();
			else if (EnergyDisType == "Brem")
				GenerateBremEnergies();
			else if (EnergyDisType == "Bbody")
				GenerateBbodyEnergies();
			else if (EnergyDisType == "Cdg")
				GenerateCdgEnergies();
			else if (EnergyDisType == "User")
				GenUserHistEnergies();
			else if (EnergyDisType == "Arb")
				GenArbPointEnergies();
			else if (EnergyDisType == "Epn")
				GenEpnHistEnergies();
			else
				G4cout << "Error: EnergyDisType has unusual value" << G4endl;
		}
	}
	return particle_energy;
}

G4double G4SPSEneDistribution::GetProbability(G4double ene) {
	G4double prob = 1.;

	if (EnergyDisType == "Lin") {
	  if (prob_norm == 1.) {
	    prob_norm = 0.5*grad*Emax*Emax + cept*Emax - 0.5*grad*Emin*Emin - cept*Emin;
	  }
	  prob = cept + grad * ene;
	  prob /= prob_norm;
	}
	else if (EnergyDisType == "Pow") {
	  if (prob_norm == 1.) {
	    if (alpha != -1.) {
	      G4double emina = std::pow(Emin, alpha + 1);
	      G4double emaxa = std::pow(Emax, alpha + 1);
	      prob_norm = 1./(1.+alpha) * (emaxa - emina);
	    } else {
	      prob_norm = std::log(Emax) - std::log(Emin) ;
	    }
	  }
	  prob = std::pow(ene, alpha)/prob_norm;
	}
	else if (EnergyDisType == "Exp"){
	  if (prob_norm == 1.) {
	    prob_norm = -Ezero*(std::exp(-Emax/Ezero) - std::exp(Emin/Ezero));
	  }  
	  prob = std::exp(-ene / Ezero);
	  prob /= prob_norm;
	}
	else if (EnergyDisType == "Arb") {
	  prob = ArbEnergyH.Value(ene);
	  //  prob = ArbEInt->CubicSplineInterpolation(ene);
	  //G4double deltaY;
	  //prob = ArbEInt->PolynomInterpolation(ene, deltaY);
	  if (prob <= 0.) {
	    //G4cout << " Warning:G4SPSEneDistribution::GetProbability: prob<= 0. "<<prob <<" "<<ene << " " <<deltaY<< G4endl;
	    G4cout << " Warning:G4SPSEneDistribution::GetProbability: prob<= 0. "<<prob <<" "<<ene << G4endl;
	    prob = 1e-30;
	  }
	  // already normalised
	}
	else
		G4cout << "Error: EnergyDisType not supported" << G4endl;
       
	return prob;
}
