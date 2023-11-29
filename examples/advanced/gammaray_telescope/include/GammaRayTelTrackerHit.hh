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
//      ------------ GammaRayTelTrackerHit  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************
// This Class describe the hits on the Tracker

#ifndef GammaRayTelTrackerHit_h
#define GammaRayTelTrackerHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"

class GammaRayTelTrackerHit: public G4VHit {
public:
	explicit GammaRayTelTrackerHit();

	~GammaRayTelTrackerHit() override;

	GammaRayTelTrackerHit(const GammaRayTelTrackerHit &right);

	auto operator=(const GammaRayTelTrackerHit &right) -> const GammaRayTelTrackerHit&;

	auto operator==(const GammaRayTelTrackerHit &right) const -> G4bool;

	inline auto operator new(size_t) -> void*;

	inline auto operator delete(void *hit) -> void;

	void Draw() override;

	void Print() override;

	inline void AddDepositedEnergy(G4double value) {
		siliconDepositedEnergy += value;
	}

	inline void SetStripNumber(const G4int &value) {
		stripNumber = value;
	}

	inline void SetSiliconPlaneNumber(const G4int &value) {
		siliconPlaneNumber = value;
	}

	inline void SetPlaneType(const G4int &value) {
		isXPlane = value;
	}

	inline void SetPosition(const G4ThreeVector &vector) {
		position = vector;
	}

	[[nodiscard]]
	inline auto GetDepositedEnergy() const -> G4double {
		return siliconDepositedEnergy;
	}

	[[nodiscard]]
	inline auto GetStripNumber() const -> G4int {
		return stripNumber;
	}

	[[nodiscard]]
	inline auto GetSiliconPlaneNumber() const -> G4int {
		return siliconPlaneNumber;
	}

	[[nodiscard]]
	inline auto GetPlaneType() const -> G4int {
		return isXPlane;
	}

	[[nodiscard]]
	inline auto GetPosition() const -> G4ThreeVector {
		return position;
	}

private:
    G4int stripNumber{0}; // Number of the strip

    G4int siliconPlaneNumber{0}; // Number of the plane

    G4int isXPlane{0}; // Type of the plane (0: Y, 1: X)

    G4double siliconDepositedEnergy{0.}; // Energy deposited on the silicon strip

    G4ThreeVector position{G4ThreeVector(0., 0., 0.)}; // Position of the hit
};

using GammaRayTelTrackerHitsCollection = G4THitsCollection<GammaRayTelTrackerHit>;
extern G4ThreadLocal G4Allocator<GammaRayTelTrackerHit> *trackerHitAllocator;

inline auto GammaRayTelTrackerHit::operator new(size_t) -> void* {
	if (trackerHitAllocator == nullptr) {
	    trackerHitAllocator = new G4Allocator<GammaRayTelTrackerHit>;
	}
	return (void*) trackerHitAllocator->MallocSingle();
}

inline void GammaRayTelTrackerHit::operator delete(void *hit) {
    trackerHitAllocator->FreeSingle((GammaRayTelTrackerHit*) hit);
}
#endif
