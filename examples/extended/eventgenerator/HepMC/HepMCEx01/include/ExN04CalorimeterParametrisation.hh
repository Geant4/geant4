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

#ifndef ExN04CalorimeterParametrisation_H
#define ExN04CalorimeterParametrisation_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"

class G4VPhysicalVolume;
class G4Tubs;

class ExN04CalorimeterParametrisation : public G4VPVParameterisation
{ 
  public:
  
    ExN04CalorimeterParametrisation();
    virtual ~ExN04CalorimeterParametrisation();
    void ComputeTransformation(const G4int copyNo,
                                     G4VPhysicalVolume *physVol) const;
    void ComputeDimensions(      G4Tubs & calorimeterLayer,
                           const G4int copyNo,
                           const G4VPhysicalVolume * physVol) const;

  private:

#include "ExN04DetectorParameterDef.hh"

};

#endif


