// Created by Joanna Weng, 9.11.04

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Polyline.hh"
#include "G4VVisManager.hh"
#include "G4Step.hh"
//GFlash
#include "GFlashEnergySpot.hh"

GFlashEnergySpot::GFlashEnergySpot() {}

GFlashEnergySpot::GFlashEnergySpot(const G4ThreeVector& point, G4double E)
{
	Point = point;
	Energy = E;
	// initialize shower start @@@@@
}

GFlashEnergySpot::~GFlashEnergySpot() {}
