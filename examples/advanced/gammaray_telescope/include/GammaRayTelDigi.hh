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
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelDigi ------
//           by F.Longo, R.Giannitrapani & G.Santin (24 oct 2001)
//
// ************************************************************

// This Class describe the digits 

#ifndef GammaRayTelDigi_h
#define GammaRayTelDigi_h 1

#include "G4VDigi.hh"
#include "G4TDigiCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class GammaRayTelDigi: public G4VDigi {
public:
	explicit GammaRayTelDigi();

	~GammaRayTelDigi() override;

	GammaRayTelDigi(const GammaRayTelDigi&);

	auto operator=(const GammaRayTelDigi& right) -> const GammaRayTelDigi&;

	auto operator==(const GammaRayTelDigi &right) const -> G4bool;

	inline auto operator new(size_t) -> void*;

	inline auto operator delete(void *digit) -> void;

	void Draw() override;

	void Print() override;

	inline void SetPlaneNumber(G4int value) {
		planeNumber = value;
	}

	inline void SetPlaneType(G4int value) {
		planeType = value;
	}

	inline void SetStripNumber(G4int value) {
		stripNumber = value;
	}

	inline void SetDigitType(G4int value) {
		digitType = value;
	}

	inline void SetEnergy(G4double value) {
		energy = value;
	}

	[[nodiscard]]
	inline auto GetPlaneNumber() const -> G4int {
		return planeNumber;
	}

	[[nodiscard]]
	inline auto GetPlaneType() const -> G4int {
		return planeType;
	}

	[[nodiscard]]
	inline auto GetStripNumber() const -> G4int {
		return stripNumber;
	}

	[[nodiscard]]
	inline auto GetDigitType() const -> G4int {
		return digitType;
	}

	[[nodiscard]]
	inline auto GetEnergy() const -> G4double {
		return energy;
	}

private:
    G4int planeType{0}; // (0: X plane, 1: Y plane)

    G4int planeNumber{0}; //  (active detector)

    G4int stripNumber{0}; // strip number

    G4int digitType{0}; // (0: TKR, 1: CAL, 2: ACD)

    G4double energy{0.}; // only for CAL
};

using GammaRayTelDigitsCollection = G4TDigiCollection<GammaRayTelDigi>;
extern G4ThreadLocal G4Allocator<GammaRayTelDigi> *digitAllocator;

inline auto GammaRayTelDigi::operator new(size_t) -> void* {
	if (digitAllocator == nullptr) {
	    digitAllocator = new G4Allocator<GammaRayTelDigi>;
	}
	return (void*) digitAllocator->MallocSingle();
}

inline auto GammaRayTelDigi::operator delete(void *digit) -> void {
    digitAllocator->FreeSingle((GammaRayTelDigi*) digit);
}
#endif
