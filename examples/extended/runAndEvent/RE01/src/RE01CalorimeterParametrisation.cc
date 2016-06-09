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
// $Id: RE01CalorimeterParametrisation.cc,v 1.1 2004/11/26 07:37:41 asaim Exp $
// GEANT4 tag $Name: geant4-07-01 $
//


#include "RE01CalorimeterParametrisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"

RE01CalorimeterParametrisation::RE01CalorimeterParametrisation()
{
#include "RE01DetectorParameterDef.icc"
}

RE01CalorimeterParametrisation::~RE01CalorimeterParametrisation()
{;}

void RE01CalorimeterParametrisation::ComputeTransformation
(const G4int,G4VPhysicalVolume *physVol) const
{
  G4ThreeVector origin;
  physVol->SetTranslation(origin);
}

void RE01CalorimeterParametrisation::ComputeDimensions
(G4Tubs & calorimeterLayer, const G4int copyNo, const G4VPhysicalVolume*) const
{
  G4double innerRad = caloTubs_rmin
              + copyNo*(absorber_thick+scinti_thick);
  calorimeterLayer.SetInnerRadius(innerRad);
  calorimeterLayer.SetOuterRadius(innerRad+absorber_thick);
  calorimeterLayer.SetZHalfLength(caloTubs_dz);
  calorimeterLayer.SetStartPhiAngle(caloTubs_sphi);
  calorimeterLayer.SetDeltaPhiAngle(caloTubs_dphi);
}
