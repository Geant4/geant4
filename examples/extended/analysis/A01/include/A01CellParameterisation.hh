// $Id: A01CellParameterisation.hh,v 1.1 2002-11-13 07:17:35 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef A01CellParameterisation_H
#define A01CellParameterisation_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
class G4VPhysicalVolume;

class A01CellParameterisation : public G4VPVParameterisation
{
  public:
    A01CellParameterisation();
    virtual ~A01CellParameterisation();
    virtual void ComputeTransformation
                   (const G4int copyNo,G4VPhysicalVolume *physVol) const;

  private:
    G4double xCell[80];
    G4double yCell[80];
};

#endif


