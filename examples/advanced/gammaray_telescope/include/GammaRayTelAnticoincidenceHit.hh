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
//      ------------ GammaRayTelAnticoincidenceHit  ------
//           by R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************
// This Class describe the hits on the Anticoincidence

#ifndef GammaRayTelAnticoincidenceHit_h
#define GammaRayTelAnticoincidenceHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"

class GammaRayTelAnticoincidenceHit: public G4VHit {
public:
	GammaRayTelAnticoincidenceHit();

	~GammaRayTelAnticoincidenceHit() override;

	GammaRayTelAnticoincidenceHit(const GammaRayTelAnticoincidenceHit &right);

	auto operator=(const GammaRayTelAnticoincidenceHit &right) -> const GammaRayTelAnticoincidenceHit&;

	auto operator==(const GammaRayTelAnticoincidenceHit &right) const -> G4bool;

	inline auto operator new(size_t) -> void*;

	inline auto operator delete(void* hit) -> void;

	void Draw() override;

	void Print() override;

	inline void AddEnergy(G4double value) {
		acdDepositedEnergy += value;
	}

	inline void SetACDTileNumber(const G4int &value) {
		acdTileNumber = value;
	}

	inline void SetACDType(const G4int &value) {
		isACDPlane = value;
	}

	inline void SetPosition(const G4ThreeVector &value) {
		position = value;
	}

	[[nodiscard]]
	inline auto GetEdepACD() const -> G4double {
		return acdDepositedEnergy;
	}

	[[nodiscard]]
	inline auto GetACDTileNumber() const -> G4int  {
		return acdTileNumber;
	}

	[[nodiscard]]
	inline auto GetACDType() const -> G4int {
		return isACDPlane;
	}

	[[nodiscard]]
	inline auto GetPosition() const -> G4ThreeVector {
		return position;
	}

private:
	G4int acdTileNumber{0}; // Number of the ACD tile

	G4int isACDPlane{0}; // Type of the plane (0: top, 1: L-R, 2: F-R)

    G4double acdDepositedEnergy{0.}; // Energy deposited on the ACD tile

    G4ThreeVector position{G4ThreeVector(0., 0., 0.)}; // Position of the hit
};

using GammaRayTelAnticoincidenceHitsCollection = G4THitsCollection<GammaRayTelAnticoincidenceHit>;
extern G4ThreadLocal G4Allocator<GammaRayTelAnticoincidenceHit> *anticoincidenceHitAllocator;

inline auto GammaRayTelAnticoincidenceHit::operator new(size_t) -> void* {
	if (anticoincidenceHitAllocator == nullptr) {
	    anticoincidenceHitAllocator = new G4Allocator<GammaRayTelAnticoincidenceHit>;
	}
	return (void*) anticoincidenceHitAllocator->MallocSingle();
}

inline auto GammaRayTelAnticoincidenceHit::operator delete(void *hit) -> void {
    anticoincidenceHitAllocator->FreeSingle((GammaRayTelAnticoincidenceHit*) hit);
}
#endif
