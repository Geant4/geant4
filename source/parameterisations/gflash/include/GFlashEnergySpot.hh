// Created by Joanna Weng, 9.11.04
#ifndef GFlashEnergySpot_h
#define GFlashEnergySpot_h
#include "G4ThreeVector.hh"

class GFlashEnergySpot
{
	public:
	GFlashEnergySpot();
	GFlashEnergySpot(const G4ThreeVector& point, G4double E);
	~GFlashEnergySpot();
	
	inline void SetEnergy(const G4double& E) {Energy = E;}
	inline G4double GetEnergy() const {return Energy;}
	
	inline void SetPosition(const G4ThreeVector& point) {Point = point;}
	inline G4ThreeVector GetPosition() const {return Point;}
		
	private:
	G4double Energy;  // energy deposition
	G4ThreeVector Point; // locus of energy deposition
};

#endif
