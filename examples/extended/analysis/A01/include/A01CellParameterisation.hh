//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: A01CellParameterisation.hh,v 1.3 2002-12-13 11:34:27 gunter Exp $
// --------------------------------------------------------------
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


