// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05EnergySpot.hh,v 1.3 1999-12-15 14:49:28 gunter Exp $
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

  // Required by Rogue Wave explicitely,
  // when using G4RWTValOrderedVector<ExN05EnergySpot>
  int operator==(const ExN05EnergySpot& eSpot) const
  {
    return (Energy==eSpot.Energy && Point==eSpot.Point) ? 1 : 0;
  };

  // Draw:
  void Draw(G4Colour* color = NULL);
  // Print:
  void Print();


private:
  G4double Energy;
  G4ThreeVector Point;
};

#endif
