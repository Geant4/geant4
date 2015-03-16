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
// MODULE:        G4SPSRandomGenerator.cc
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
#include "G4PrimaryParticle.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include <cmath>
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4AutoLock.hh"

#include "G4SPSRandomGenerator.hh"

//G4SPSRandomGenerator* G4SPSRandomGenerator::instance = 0;

G4SPSRandomGenerator::bweights_t::bweights_t() {
  for ( int i = 0 ; i < 9 ; ++i ) w[i] = 1;
}

G4double&
G4SPSRandomGenerator::bweights_t::operator [](const int i) { return w[i]; }

G4SPSRandomGenerator::G4SPSRandomGenerator()
{
	// Initialise all variables

	// Bias variables
	XBias = false;
	IPDFXBias = false;
	YBias = false;
	IPDFYBias = false;
	ZBias = false;
	IPDFZBias = false;
	ThetaBias = false;
	IPDFThetaBias = false;
	PhiBias = false;
	IPDFPhiBias = false;
	EnergyBias = false;
	IPDFEnergyBias = false;
	PosThetaBias = false;
	IPDFPosThetaBias = false;
	PosPhiBias = false;
	IPDFPosPhiBias = false;
	verbosityLevel = 0;
    G4MUTEXINIT(mutex);
}

G4SPSRandomGenerator::~G4SPSRandomGenerator() {
    G4MUTEXDESTROY(mutex);
}

//G4SPSRandomGenerator* G4SPSRandomGenerator::getInstance ()
//{
//  if (instance == 0) instance = new G4SPSRandomGenerator();
//  return instance;
//}

// Biasing methods

void G4SPSRandomGenerator::SetXBias(G4ThreeVector input) {
    G4AutoLock l(&mutex);
	G4double ehi, val;
	ehi = input.x();
	val = input.y();
	XBiasH.InsertValues(ehi, val);
	XBias = true;
}

void G4SPSRandomGenerator::SetYBias(G4ThreeVector input) {
    G4AutoLock l(&mutex);
	G4double ehi, val;
	ehi = input.x();
	val = input.y();
	YBiasH.InsertValues(ehi, val);
	YBias = true;
}

void G4SPSRandomGenerator::SetZBias(G4ThreeVector input) {
    G4AutoLock l(&mutex);
	G4double ehi, val;
	ehi = input.x();
	val = input.y();
	ZBiasH.InsertValues(ehi, val);
	ZBias = true;
}

void G4SPSRandomGenerator::SetThetaBias(G4ThreeVector input) {
    G4AutoLock l(&mutex);
	G4double ehi, val;
	ehi = input.x();
	val = input.y();
	ThetaBiasH.InsertValues(ehi, val);
	ThetaBias = true;
}

void G4SPSRandomGenerator::SetPhiBias(G4ThreeVector input) {
    G4AutoLock l(&mutex);
	G4double ehi, val;
	ehi = input.x();
	val = input.y();
	PhiBiasH.InsertValues(ehi, val);
	PhiBias = true;
}

void G4SPSRandomGenerator::SetEnergyBias(G4ThreeVector input) {
    G4AutoLock l(&mutex);
	G4double ehi, val;
	ehi = input.x();
	val = input.y();
	EnergyBiasH.InsertValues(ehi, val);
	EnergyBias = true;
}

void G4SPSRandomGenerator::SetPosThetaBias(G4ThreeVector input) {
    G4AutoLock l(&mutex);
	G4double ehi, val;
	ehi = input.x();
	val = input.y();
	PosThetaBiasH.InsertValues(ehi, val);
	PosThetaBias = true;
}

void G4SPSRandomGenerator::SetPosPhiBias(G4ThreeVector input) {
    G4AutoLock l(&mutex);
	G4double ehi, val;
	ehi = input.x();
	val = input.y();
	PosPhiBiasH.InsertValues(ehi, val);
	PosPhiBias = true;
}

void G4SPSRandomGenerator::SetIntensityWeight(G4double weight) {
  bweights.Get()[8] = weight;
}

G4double G4SPSRandomGenerator::GetBiasWeight() {
  bweights_t& w = bweights.Get();
  return w[0] * w[1] * w[2] * w[3]
    * w[4] * w[5] * w[6] * w[7]
    * w[8];
}

void G4SPSRandomGenerator::SetVerbosity(G4int a) {
     G4AutoLock l(&mutex);
    verbosityLevel = a;
}


namespace {
  G4PhysicsOrderedFreeVector ZeroPhysVector; // for re-set only
}

void G4SPSRandomGenerator::ReSetHist(G4String atype) {
    G4AutoLock l(&mutex);
	if (atype == "biasx") {
		XBias = false;
		IPDFXBias = false;
		local_IPDFXBias.Get().val = false;
		XBiasH = IPDFXBiasH = ZeroPhysVector;
	} else if (atype == "biasy") {
		YBias = false;
		IPDFYBias = false;
                local_IPDFYBias.Get().val = false;
		YBiasH = IPDFYBiasH = ZeroPhysVector;
	} else if (atype == "biasz") {
		ZBias = false;
		IPDFZBias = false;
                local_IPDFZBias.Get().val = false;
		ZBiasH = IPDFZBiasH = ZeroPhysVector;
	} else if (atype == "biast") {
		ThetaBias = false;
		IPDFThetaBias = false;
                local_IPDFThetaBias.Get().val = false;
		ThetaBiasH = IPDFThetaBiasH = ZeroPhysVector;
	} else if (atype == "biasp") {
		PhiBias = false;
		IPDFPhiBias = false;
                local_IPDFPhiBias.Get().val = false;
		PhiBiasH = IPDFPhiBiasH = ZeroPhysVector;
	} else if (atype == "biase") {
		EnergyBias = false;
		IPDFEnergyBias = false;
                local_IPDFEnergyBias.Get().val = false;
		EnergyBiasH = IPDFEnergyBiasH = ZeroPhysVector;
	} else if (atype == "biaspt") {
		PosThetaBias = false;
		IPDFPosThetaBias = false;
		local_IPDFPosThetaBias.Get().val = false;
		PosThetaBiasH = IPDFPosThetaBiasH = ZeroPhysVector;
	} else if (atype == "biaspp") {
		PosPhiBias = false;
		IPDFPosPhiBias = false;
		local_IPDFPosPhiBias.Get().val = false;
		PosPhiBiasH = IPDFPosPhiBiasH = ZeroPhysVector;
	} else {
		G4cout << "Error, histtype not accepted " << G4endl;
	}
}

G4double G4SPSRandomGenerator::GenRandX() {
	if (verbosityLevel >= 1)
		G4cout << "In GenRandX" << G4endl;
	if (XBias == false) {
		// X is not biased
		G4double rndm = G4UniformRand();
		return (rndm);
	} else {
	    // X is biased
	    //This is shared among threads, and we need to initialize
	    //only once. Multiple instances of this class can exists
	    //so we rely on a class-private, thread-private variable
	    //to check if we need an initialiation. We do not use a lock here
	    //because the boolean on which we check is thread private
	    if ( local_IPDFXBias.Get().val == false ) {
	        //For time that this thread arrived, here
	        //Now two cases are possible: it is the first time
	         //ANY thread has ever initialized this.
	        //Now we need a lock. In any case, the thread local
	        //variable can now be set to true
	        local_IPDFXBias.Get().val = true;
	        G4AutoLock l(&mutex);
	        if (IPDFXBias == false) {
	            // IPDF has not been created, so create it
	            G4double bins[1024], vals[1024], sum;
	            G4int ii;
	            G4int maxbin = G4int(XBiasH.GetVectorLength());
	            bins[0] = XBiasH.GetLowEdgeEnergy(size_t(0));
	            vals[0] = XBiasH(size_t(0));
	            sum = vals[0];
	            for (ii = 1; ii < maxbin; ii++) {
	                bins[ii] = XBiasH.GetLowEdgeEnergy(size_t(ii));
	                vals[ii] = XBiasH(size_t(ii)) + vals[ii - 1];
	                sum = sum + XBiasH(size_t(ii));
	            }

	            for (ii = 0; ii < maxbin; ii++) {
	                vals[ii] = vals[ii] / sum;
	                IPDFXBiasH.InsertValues(bins[ii], vals[ii]);
	            }
	            // Make IPDFXBias = true
	            IPDFXBias = true;
	        }
	    }
          // IPDF has been create so carry on
          G4double rndm = G4UniformRand();

          // Calculate the weighting: Find the bin that the determined
          // rndm is in and the weigthing will be the difference in the
          // natural probability (from the x-axis) divided by the
          // difference in the biased probability (the area).
          size_t numberOfBin = IPDFXBiasH.GetVectorLength();
          G4int biasn1 = 0;
          G4int biasn2 = numberOfBin / 2;
          G4int biasn3 = numberOfBin - 1;
          while (biasn1 != biasn3 - 1) {
              if (rndm > IPDFXBiasH(biasn2))
                biasn1 = biasn2;
              else
                biasn3 = biasn2;
              biasn2 = biasn1 + (biasn3 - biasn1 + 1) / 2;
          }
          // retrieve the areas and then the x-axis values
          bweights_t& w = bweights.Get();
          w[0] = IPDFXBiasH(biasn2) - IPDFXBiasH(biasn2 - 1);
          G4double xaxisl = IPDFXBiasH.GetLowEdgeEnergy(size_t(biasn2 - 1));
          G4double xaxisu = IPDFXBiasH.GetLowEdgeEnergy(size_t(biasn2));
          G4double NatProb = xaxisu - xaxisl;
          //G4cout << "X Bin weight " << bweights[0] << " " << rndm << G4endl;
          //G4cout << "lower and upper xaxis vals "<<xaxisl<<" "<<xaxisu<<G4endl;
          w[0] = NatProb / w[0];
          if (verbosityLevel >= 1)
            G4cout << "X bin weight " << w[0] << " " << rndm << G4endl;
          return (IPDFXBiasH.GetEnergy(rndm));
	}
}

G4double G4SPSRandomGenerator::GenRandY() {
	if (verbosityLevel >= 1)
		G4cout << "In GenRandY" << G4endl;
	if (YBias == false) {
		// Y is not biased
		G4double rndm = G4UniformRand();
		return (rndm);
	} else {
	    // Y is biased
	    if ( local_IPDFYBias.Get().val == false ) {
	        local_IPDFYBias.Get().val = true;
	        G4AutoLock l(&mutex);
		if (IPDFYBias == false) {
			// IPDF has not been created, so create it
			G4double bins[1024], vals[1024], sum;
			G4int ii;
			G4int maxbin = G4int(YBiasH.GetVectorLength());
			bins[0] = YBiasH.GetLowEdgeEnergy(size_t(0));
			vals[0] = YBiasH(size_t(0));
			sum = vals[0];
			for (ii = 1; ii < maxbin; ii++) {
				bins[ii] = YBiasH.GetLowEdgeEnergy(size_t(ii));
				vals[ii] = YBiasH(size_t(ii)) + vals[ii - 1];
				sum = sum + YBiasH(size_t(ii));
			}

			for (ii = 0; ii < maxbin; ii++) {
				vals[ii] = vals[ii] / sum;
				IPDFYBiasH.InsertValues(bins[ii], vals[ii]);
			}
			// Make IPDFYBias = true
			IPDFYBias = true;
		}
	    }		// IPDF has been create so carry on
	    G4double rndm = G4UniformRand();
	    size_t numberOfBin = IPDFYBiasH.GetVectorLength();
	    G4int biasn1 = 0;
	    G4int biasn2 = numberOfBin / 2;
	    G4int biasn3 = numberOfBin - 1;
	    while (biasn1 != biasn3 - 1) {
	        if (rndm > IPDFYBiasH(biasn2))
	          biasn1 = biasn2;
	        else
	          biasn3 = biasn2;
	        biasn2 = biasn1 + (biasn3 - biasn1 + 1) / 2;
	    }
	    bweights_t& w = bweights.Get();
	    w[1] = IPDFYBiasH(biasn2) - IPDFYBiasH(biasn2 - 1);
	    G4double xaxisl = IPDFYBiasH.GetLowEdgeEnergy(size_t(biasn2 - 1));
	    G4double xaxisu = IPDFYBiasH.GetLowEdgeEnergy(size_t(biasn2));
	    G4double NatProb = xaxisu - xaxisl;
	    w[1] = NatProb / w[1];
	    if (verbosityLevel >= 1)
	      G4cout << "Y bin weight " << w[1] << " " << rndm << G4endl;
	    return (IPDFYBiasH.GetEnergy(rndm));
	}
}

G4double G4SPSRandomGenerator::GenRandZ() {
	if (verbosityLevel >= 1)
		G4cout << "In GenRandZ" << G4endl;
	if (ZBias == false) {
		// Z is not biased
		G4double rndm = G4UniformRand();
		return (rndm);
	} else {
		// Z is biased
	    if (local_IPDFZBias.Get().val == false ) {
	        local_IPDFZBias.Get().val = true;
	        G4AutoLock l(&mutex);
		if (IPDFZBias == false) {
		    // IPDF has not been created, so create it
		    G4double bins[1024], vals[1024], sum;
		    G4int ii;
		    G4int maxbin = G4int(ZBiasH.GetVectorLength());
		    bins[0] = ZBiasH.GetLowEdgeEnergy(size_t(0));
		    vals[0] = ZBiasH(size_t(0));
		    sum = vals[0];
		    for (ii = 1; ii < maxbin; ii++) {
		        bins[ii] = ZBiasH.GetLowEdgeEnergy(size_t(ii));
		        vals[ii] = ZBiasH(size_t(ii)) + vals[ii - 1];
		        sum = sum + ZBiasH(size_t(ii));
		    }

		    for (ii = 0; ii < maxbin; ii++) {
		        vals[ii] = vals[ii] / sum;
		        IPDFZBiasH.InsertValues(bins[ii], vals[ii]);
		    }
		    // Make IPDFZBias = true
		    IPDFZBias = true;
		}
	    }
	    // IPDF has been create so carry on
	    G4double rndm = G4UniformRand();
	    //      size_t weight_bin_no = IPDFZBiasH.FindValueBinLocation(rndm);
	    size_t numberOfBin = IPDFZBiasH.GetVectorLength();
	    G4int biasn1 = 0;
	    G4int biasn2 = numberOfBin / 2;
	    G4int biasn3 = numberOfBin - 1;
	    while (biasn1 != biasn3 - 1) {
	        if (rndm > IPDFZBiasH(biasn2))
	          biasn1 = biasn2;
	        else
	          biasn3 = biasn2;
	        biasn2 = biasn1 + (biasn3 - biasn1 + 1) / 2;
	    }
	    bweights_t& w = bweights.Get();
	    w[2] = IPDFZBiasH(biasn2) - IPDFZBiasH(biasn2 - 1);
	    G4double xaxisl = IPDFZBiasH.GetLowEdgeEnergy(size_t(biasn2 - 1));
	    G4double xaxisu = IPDFZBiasH.GetLowEdgeEnergy(size_t(biasn2));
	    G4double NatProb = xaxisu - xaxisl;
	    w[2] = NatProb / w[2];
	    if (verbosityLevel >= 1)
	      G4cout << "Z bin weight " << w[2] << " " << rndm << G4endl;
	    return (IPDFZBiasH.GetEnergy(rndm));
	}
}

G4double G4SPSRandomGenerator::GenRandTheta() {
	if (verbosityLevel >= 1) {
		G4cout << "In GenRandTheta" << G4endl;
		G4cout << "Verbosity " << verbosityLevel << G4endl;
	}
	if (ThetaBias == false) {
		// Theta is not biased
		G4double rndm = G4UniformRand();
		return (rndm);
	} else {
		// Theta is biased
	    if ( local_IPDFThetaBias.Get().val == false ) {
	        local_IPDFThetaBias.Get().val = true;
	        G4AutoLock l(&mutex);
		if (IPDFThetaBias == false) {
		    // IPDF has not been created, so create it
		    G4double bins[1024], vals[1024], sum;
		    G4int ii;
		    G4int maxbin = G4int(ThetaBiasH.GetVectorLength());
		    bins[0] = ThetaBiasH.GetLowEdgeEnergy(size_t(0));
		    vals[0] = ThetaBiasH(size_t(0));
		    sum = vals[0];
		    for (ii = 1; ii < maxbin; ii++) {
		        bins[ii] = ThetaBiasH.GetLowEdgeEnergy(size_t(ii));
		        vals[ii] = ThetaBiasH(size_t(ii)) + vals[ii - 1];
		        sum = sum + ThetaBiasH(size_t(ii));
		    }

		    for (ii = 0; ii < maxbin; ii++) {
		        vals[ii] = vals[ii] / sum;
		        IPDFThetaBiasH.InsertValues(bins[ii], vals[ii]);
		    }
		    // Make IPDFThetaBias = true
		    IPDFThetaBias = true;
		}
	    }
	    // IPDF has been create so carry on
	    G4double rndm = G4UniformRand();
	    //      size_t weight_bin_no = IPDFThetaBiasH.FindValueBinLocation(rndm);
	    size_t numberOfBin = IPDFThetaBiasH.GetVectorLength();
	    G4int biasn1 = 0;
	    G4int biasn2 = numberOfBin / 2;
	    G4int biasn3 = numberOfBin - 1;
	    while (biasn1 != biasn3 - 1) {
	        if (rndm > IPDFThetaBiasH(biasn2))
	          biasn1 = biasn2;
	        else
	          biasn3 = biasn2;
	        biasn2 = biasn1 + (biasn3 - biasn1 + 1) / 2;
	    }
	    bweights_t& w = bweights.Get();
	    w[3] = IPDFThetaBiasH(biasn2) - IPDFThetaBiasH(biasn2 - 1);
	    G4double xaxisl = IPDFThetaBiasH.GetLowEdgeEnergy(size_t(biasn2 - 1));
	    G4double xaxisu = IPDFThetaBiasH.GetLowEdgeEnergy(size_t(biasn2));
	    G4double NatProb = xaxisu - xaxisl;
	    w[3] = NatProb / w[3];
	    if (verbosityLevel >= 1)
	      G4cout << "Theta bin weight " << w[3] << " " << rndm
	      << G4endl;
	    return (IPDFThetaBiasH.GetEnergy(rndm));
	}
}

G4double G4SPSRandomGenerator::GenRandPhi() {
	if (verbosityLevel >= 1)
		G4cout << "In GenRandPhi" << G4endl;
	if (PhiBias == false) {
		// Phi is not biased
		G4double rndm = G4UniformRand();
		return (rndm);
	} else {
	    // Phi is biased
	    if ( local_IPDFPhiBias.Get().val == false ) {
	        local_IPDFPhiBias.Get().val = true;
	        G4AutoLock l(&mutex);
	        if (IPDFPhiBias == false) {
	            // IPDF has not been created, so create it
	            G4double bins[1024], vals[1024], sum;
	            G4int ii;
	            G4int maxbin = G4int(PhiBiasH.GetVectorLength());
	            bins[0] = PhiBiasH.GetLowEdgeEnergy(size_t(0));
	            vals[0] = PhiBiasH(size_t(0));
	            sum = vals[0];
	            for (ii = 1; ii < maxbin; ii++) {
	                bins[ii] = PhiBiasH.GetLowEdgeEnergy(size_t(ii));
	                vals[ii] = PhiBiasH(size_t(ii)) + vals[ii - 1];
	                sum = sum + PhiBiasH(size_t(ii));
	            }

	            for (ii = 0; ii < maxbin; ii++) {
	                vals[ii] = vals[ii] / sum;
	                IPDFPhiBiasH.InsertValues(bins[ii], vals[ii]);
	            }
	            // Make IPDFPhiBias = true
	            IPDFPhiBias = true;
		}
	    }
	    // IPDF has been create so carry on
	    G4double rndm = G4UniformRand();
	    //      size_t weight_bin_no = IPDFPhiBiasH.FindValueBinLocation(rndm);
	    size_t numberOfBin = IPDFPhiBiasH.GetVectorLength();
	    G4int biasn1 = 0;
	    G4int biasn2 = numberOfBin / 2;
	    G4int biasn3 = numberOfBin - 1;
	    while (biasn1 != biasn3 - 1) {
	        if (rndm > IPDFPhiBiasH(biasn2))
	          biasn1 = biasn2;
	        else
	          biasn3 = biasn2;
	        biasn2 = biasn1 + (biasn3 - biasn1 + 1) / 2;
	    }
	    bweights_t& w = bweights.Get();
	    w[4] = IPDFPhiBiasH(biasn2) - IPDFPhiBiasH(biasn2 - 1);
	    G4double xaxisl = IPDFPhiBiasH.GetLowEdgeEnergy(size_t(biasn2 - 1));
	    G4double xaxisu = IPDFPhiBiasH.GetLowEdgeEnergy(size_t(biasn2));
	    G4double NatProb = xaxisu - xaxisl;
	    w[4] = NatProb / w[4];
	    if (verbosityLevel >= 1)
	      G4cout << "Phi bin weight " << w[4] << " " << rndm << G4endl;
	    return (IPDFPhiBiasH.GetEnergy(rndm));
	}
}

G4double G4SPSRandomGenerator::GenRandEnergy() {
	if (verbosityLevel >= 1)
		G4cout << "In GenRandEnergy" << G4endl;
	if (EnergyBias == false) {
		// Energy is not biased
		G4double rndm = G4UniformRand();
		return (rndm);
	} else {
	    if ( local_IPDFEnergyBias.Get().val == false ) {
	        local_IPDFEnergyBias.Get().val = true;
	        // ENERGY is biased
	        G4AutoLock l(&mutex);
		if (IPDFEnergyBias == false) {
		    // IPDF has not been created, so create it
		    G4double bins[1024], vals[1024], sum;
		    G4int ii;
		    G4int maxbin = G4int(EnergyBiasH.GetVectorLength());
		    bins[0] = EnergyBiasH.GetLowEdgeEnergy(size_t(0));
		    vals[0] = EnergyBiasH(size_t(0));
		    sum = vals[0];
		    for (ii = 1; ii < maxbin; ii++) {
		        bins[ii] = EnergyBiasH.GetLowEdgeEnergy(size_t(ii));
		        vals[ii] = EnergyBiasH(size_t(ii)) + vals[ii - 1];
		        sum = sum + EnergyBiasH(size_t(ii));
		    }
		    IPDFEnergyBiasH = ZeroPhysVector;
		    for (ii = 0; ii < maxbin; ii++) {
		        vals[ii] = vals[ii] / sum;
		        IPDFEnergyBiasH.InsertValues(bins[ii], vals[ii]);
		    }
		    // Make IPDFEnergyBias = true
		    IPDFEnergyBias = true;
		}
	    }
	    // IPDF has been create so carry on
	    G4double rndm = G4UniformRand();
	    //  size_t weight_bin_no = IPDFEnergyBiasH.FindValueBinLocation(rndm);
	    size_t numberOfBin = IPDFEnergyBiasH.GetVectorLength();
	    G4int biasn1 = 0;
	    G4int biasn2 = numberOfBin / 2;
	    G4int biasn3 = numberOfBin - 1;
	    while (biasn1 != biasn3 - 1) {
	        if (rndm > IPDFEnergyBiasH(biasn2))
	          biasn1 = biasn2;
	        else
	          biasn3 = biasn2;
	        biasn2 = biasn1 + (biasn3 - biasn1 + 1) / 2;
	    }
	    bweights_t& w = bweights.Get();
	    w[5] = IPDFEnergyBiasH(biasn2) - IPDFEnergyBiasH(biasn2 - 1);
	    G4double xaxisl = IPDFEnergyBiasH.GetLowEdgeEnergy(size_t(biasn2 - 1));
	    G4double xaxisu = IPDFEnergyBiasH.GetLowEdgeEnergy(size_t(biasn2));
	    G4double NatProb = xaxisu - xaxisl;
	    w[5] = NatProb / w[5];
	    if (verbosityLevel >= 1)
	      G4cout << "Energy bin weight " << w[5] << " " << rndm
	      << G4endl;
	    return (IPDFEnergyBiasH.GetEnergy(rndm));
	}
}

G4double G4SPSRandomGenerator::GenRandPosTheta() {
	if (verbosityLevel >= 1) {
		G4cout << "In GenRandPosTheta" << G4endl;
		G4cout << "Verbosity " << verbosityLevel << G4endl;
	}
	if (PosThetaBias == false) {
		// Theta is not biased
		G4double rndm = G4UniformRand();
		return (rndm);
	} else {
		// Theta is biased
	    if ( local_IPDFPosThetaBias.Get().val == false ) {
	        local_IPDFPosThetaBias.Get().val = true;
	        G4AutoLock l(&mutex);
	        if (IPDFPosThetaBias == false) {
	            // IPDF has not been created, so create it
	            G4double bins[1024], vals[1024], sum;
	            G4int ii;
	            G4int maxbin = G4int(PosThetaBiasH.GetVectorLength());
	            bins[0] = PosThetaBiasH.GetLowEdgeEnergy(size_t(0));
	            vals[0] = PosThetaBiasH(size_t(0));
	            sum = vals[0];
	            for (ii = 1; ii < maxbin; ii++) {
	                bins[ii] = PosThetaBiasH.GetLowEdgeEnergy(size_t(ii));
	                vals[ii] = PosThetaBiasH(size_t(ii)) + vals[ii - 1];
	                sum = sum + PosThetaBiasH(size_t(ii));
	            }

	            for (ii = 0; ii < maxbin; ii++) {
	                vals[ii] = vals[ii] / sum;
	                IPDFPosThetaBiasH.InsertValues(bins[ii], vals[ii]);
	            }
	            // Make IPDFThetaBias = true
	            IPDFPosThetaBias = true;
		}
	    }
	    // IPDF has been create so carry on
	    G4double rndm = G4UniformRand();
	    //      size_t weight_bin_no = IPDFThetaBiasH.FindValueBinLocation(rndm);
	    size_t numberOfBin = IPDFPosThetaBiasH.GetVectorLength();
	    G4int biasn1 = 0;
	    G4int biasn2 = numberOfBin / 2;
	    G4int biasn3 = numberOfBin - 1;
	    while (biasn1 != biasn3 - 1) {
	        if (rndm > IPDFPosThetaBiasH(biasn2))
	          biasn1 = biasn2;
	        else
	          biasn3 = biasn2;
	        biasn2 = biasn1 + (biasn3 - biasn1 + 1) / 2;
	    }
	    bweights_t& w = bweights.Get();
	    w[6] = IPDFPosThetaBiasH(biasn2) - IPDFPosThetaBiasH(biasn2 - 1);
	      G4double xaxisl =
	          IPDFPosThetaBiasH.GetLowEdgeEnergy(size_t(biasn2 - 1));
	      G4double xaxisu = IPDFPosThetaBiasH.GetLowEdgeEnergy(size_t(biasn2));
	      G4double NatProb = xaxisu - xaxisl;
	      w[6] = NatProb / w[6];
	      if (verbosityLevel >= 1)
	        G4cout << "PosTheta bin weight " << w[6] << " " << rndm
	        << G4endl;
	      return (IPDFPosThetaBiasH.GetEnergy(rndm));
	}
}

G4double G4SPSRandomGenerator::GenRandPosPhi() {
	if (verbosityLevel >= 1)
		G4cout << "In GenRandPosPhi" << G4endl;
	if (PosPhiBias == false) {
		// PosPhi is not biased
		G4double rndm = G4UniformRand();
		return (rndm);
	} else {
		// PosPhi is biased
	    if (local_IPDFPosPhiBias.Get().val == false ) {
	        local_IPDFPosPhiBias.Get().val = true;
	        G4AutoLock l(&mutex);
		if (IPDFPosPhiBias == false) {
		    // IPDF has not been created, so create it
		    G4double bins[1024], vals[1024], sum;
		    G4int ii;
		    G4int maxbin = G4int(PosPhiBiasH.GetVectorLength());
		    bins[0] = PosPhiBiasH.GetLowEdgeEnergy(size_t(0));
		    vals[0] = PosPhiBiasH(size_t(0));
		    sum = vals[0];
		    for (ii = 1; ii < maxbin; ii++) {
		        bins[ii] = PosPhiBiasH.GetLowEdgeEnergy(size_t(ii));
		        vals[ii] = PosPhiBiasH(size_t(ii)) + vals[ii - 1];
		        sum = sum + PosPhiBiasH(size_t(ii));
		    }

		    for (ii = 0; ii < maxbin; ii++) {
		        vals[ii] = vals[ii] / sum;
		        IPDFPosPhiBiasH.InsertValues(bins[ii], vals[ii]);
		    }
		    // Make IPDFPosPhiBias = true
		    IPDFPosPhiBias = true;
		}
	    }
	    // IPDF has been create so carry on
	    G4double rndm = G4UniformRand();
	    //      size_t weight_bin_no = IPDFPosPhiBiasH.FindValueBinLocation(rndm);
	    size_t numberOfBin = IPDFPosPhiBiasH.GetVectorLength();
	    G4int biasn1 = 0;
	    G4int biasn2 = numberOfBin / 2;
	    G4int biasn3 = numberOfBin - 1;
	    while (biasn1 != biasn3 - 1) {
	        if (rndm > IPDFPosPhiBiasH(biasn2))
	          biasn1 = biasn2;
	        else
	          biasn3 = biasn2;
	        biasn2 = biasn1 + (biasn3 - biasn1 + 1) / 2;
	    }
	    bweights_t& w = bweights.Get();
	    w[7] = IPDFPosPhiBiasH(biasn2) - IPDFPosPhiBiasH(biasn2 - 1);
	    G4double xaxisl = IPDFPosPhiBiasH.GetLowEdgeEnergy(size_t(biasn2 - 1));
	    G4double xaxisu = IPDFPosPhiBiasH.GetLowEdgeEnergy(size_t(biasn2));
	    G4double NatProb = xaxisu - xaxisl;
	    w[7] = NatProb / w[7];
	    if (verbosityLevel >= 1)
	      G4cout << "PosPhi bin weight " << w[7] << " " << rndm
	      << G4endl;
	    return (IPDFPosPhiBiasH.GetEnergy(rndm));
	}
}
