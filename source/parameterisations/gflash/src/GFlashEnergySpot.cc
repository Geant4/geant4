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


void GFlashEnergySpot::Draw(G4Colour *color)
{
	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
	if (pVVisManager)
	{
		G4double localUnit = 1*cm; // possibly a user setable number @@@@
		G4Polyline polyline;
		G4Colour localColor(1.,.5,.5);
		if (color) localColor = *color;
		polyline.SetVisAttributes(localColor);
		G4ThreeVector pp(Point);
		pp.setZ(pp.z()+1*localUnit);
		polyline.push_back(pp);
		pp.setZ(pp.z()-2*localUnit);
		polyline.push_back(pp);
		pp = Point;
		polyline.push_back(pp);
		pp.setX(pp.x()+1*localUnit);
		polyline.push_back(pp);
		pp.setX(pp.x()-2*localUnit);
		polyline.push_back(pp);
		pp = Point;
		polyline.push_back(pp);
		pp.setY(pp.y()+1*localUnit);
		polyline.push_back(pp);
		pp.setY(pp.y()-2*localUnit);
		polyline.push_back(pp);
		pVVisManager -> Draw(polyline);
	}
}

void GFlashEnergySpot::Print()
{
	std::cout << " GFlashEnergySpot {E = " << Energy << "; Position = " 
	<< Point << " }"<< std::endl;
}

