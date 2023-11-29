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
//      ------------ GammaRayTelCalorimeterHit  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************
// This Class describe the hits on the Calorimeter

#ifndef GammaRayTelCalorimeterHit_h
#define GammaRayTelCalorimeterHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"

class GammaRayTelCalorimeterHit: public G4VHit {
public:
	GammaRayTelCalorimeterHit();

	~GammaRayTelCalorimeterHit() override;

	GammaRayTelCalorimeterHit(const GammaRayTelCalorimeterHit &right);

	auto operator=(const GammaRayTelCalorimeterHit &right) -> const GammaRayTelCalorimeterHit&;

	auto operator==(const GammaRayTelCalorimeterHit &right) const -> G4bool;

	inline auto operator new(size_t) -> void*;

	inline auto operator delete(void* hit) -> void;

	void Draw() override;

	void Print() override;

	inline void AddEnergy(G4double value) {
		calDepositedEnergy += value;
	}

	inline void SetCALBarNumber(const G4int &value) {
		calBarNumber = value;
	}

	inline void SetCALPlaneNumber(const G4int &value) {
		calPlaneNumber = value;
	}

	inline void SetCALType(const G4int &value) {
		isCALPlane = value;
	}

	inline void SetPosition(const G4ThreeVector &vector) {
		position = vector;
	}

	[[nodiscard]]
	inline auto GetCALDepositedEnergy() const -> G4double {
		return calDepositedEnergy;
	}

	[[nodiscard]]
	inline auto GetCALBarNumber() const -> G4int {
		return calBarNumber;
	}

	[[nodiscard]]
	inline auto GetCALPlaneNumber() const -> G4int {
		return calPlaneNumber;
	}

	[[nodiscard]]
	inline auto GetCALType() const -> G4int {
		return isCALPlane;
	}

	[[nodiscard]]
	inline auto GetPosition() const -> G4ThreeVector {
		return position;
	}

private:
    G4int calBarNumber{0}; // Number of the calorimeter (CAL) tile

    G4int calPlaneNumber{0}; // Number of the calorimeter (CAL) plane

    G4int isCALPlane{0}; // Type of the plane (0: X plane, 1: Y plane)

    G4double calDepositedEnergy{0.}; // Energy deposited on the calorimeter (CAL) tile

    G4ThreeVector position{G4ThreeVector(0., 0., 0.)}; // Position of the hit
};

using GammaRayTelCalorimeterHitsCollection = G4THitsCollection<GammaRayTelCalorimeterHit>;
extern G4ThreadLocal G4Allocator<GammaRayTelCalorimeterHit> *calorimeterHitAllocator;

inline auto GammaRayTelCalorimeterHit::operator new(size_t) -> void* {
	if (calorimeterHitAllocator == nullptr) {
	    calorimeterHitAllocator = new G4Allocator<GammaRayTelCalorimeterHit>;
	}
	return (void*) calorimeterHitAllocator->MallocSingle();
}

inline auto GammaRayTelCalorimeterHit::operator delete(void* hit) -> void {
    calorimeterHitAllocator->FreeSingle((GammaRayTelCalorimeterHit*) hit);
}
#endif
