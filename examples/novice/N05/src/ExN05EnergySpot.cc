// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05EnergySpot.cc,v 1.2 1999-12-15 14:49:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "ExN05EnergySpot.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Polyline.hh"
#include "G4VVisManager.hh"
#include "G4Step.hh"


ExN05EnergySpot::ExN05EnergySpot()
{;}

ExN05EnergySpot::ExN05EnergySpot(const G4ThreeVector& point, G4double E)
{
  Point = point;
  Energy = E;
}

ExN05EnergySpot::~ExN05EnergySpot()
{;}


void ExN05EnergySpot::Draw(G4Colour *color)
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (pVVisManager)
    {
      G4Polyline polyline;
      G4Colour colour(1.,.5,.5);
      if (color != 0) colour = *color;
      G4VisAttributes attribs(colour);
      polyline.SetVisAttributes(colour);
      G4ThreeVector pp(Point);
      // Draw a "home made" marker:
      // Will be better by using a real Marker:
      pp.setZ(pp.z()+1*cm);
      polyline.append(pp);
      pp.setZ(pp.z()-2*cm);
      polyline.append(pp);
      pp = Point;
      polyline.append(pp);
      pp.setX(pp.x()+1*cm);
      polyline.append(pp);
      pp.setX(pp.x()-2*cm);
      polyline.append(pp);
      pp = Point;
      polyline.append(pp);
      pp.setY(pp.y()+1*cm);
      polyline.append(pp);
      pp.setY(pp.y()-2*cm);
      polyline.append(pp);
      pVVisManager -> Draw(polyline);
    }
}

void ExN05EnergySpot::Print()
{
  G4cout << " ExN05EnergySpot {E = " << Energy << "; Position = " << Point << " }"<< G4endl;
}





