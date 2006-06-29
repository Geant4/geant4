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
// $Id: ExN05EnergySpot.hh,v 1.8 2006-06-29 17:52:28 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef ExN05EnergySpot_h
#define ExN05EnergySpot_h

#include "G4ThreeVector.hh"
class G4Colour;

class ExN05EnergySpot
{
public:
  ExN05EnergySpot();
  ExN05EnergySpot(const G4ThreeVector& point, G4double E);
  ~ExN05EnergySpot();

  inline void SetEnergy(const G4double& E) {Energy = E;}
  inline G4double GetEnergy() const {return Energy;}

  inline void SetPosition(const G4ThreeVector& point) {Point = point;}
  inline G4ThreeVector GetPosition() const {return Point;}

  G4int operator==(const ExN05EnergySpot& eSpot) const
  {
    return (Energy==eSpot.Energy && Point==eSpot.Point) ? 1 : 0;
  }

  // Draw:
  void Draw(G4Colour* color = 0);
  // Print:
  void Print();


private:
  G4double Energy;
  G4ThreeVector Point;
};

#endif
