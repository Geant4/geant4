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
/// \file electromagnetic/TestEm10/include/XTRTransparentRegRadModel.hh
/// \brief Definition of the XTRTransparentRegRadModel class
//
// 
///////////////////////////////////////////////////////////////////////////
// 
// Process describing a radiator of X-ray transition radiation.  
// Thicknesses of plates and gas gaps are fixed.
// We suppose that:
// formation zone ~ mean thickness << absorption length
// for each material and in the range 1-100 keV. This allows us to simplify
// interference effects in radiator stack (GetStackFactor method).
// 
// 
// History:
//
// 05.04.05 V. Grichine, first version 
//


#ifndef XTRTransparentRegRadModel_h
#define XTRTransparentRegRadModel_h 1

#include "G4VXTRenergyLoss.hh"

class XTRTransparentRegRadModel : public G4VXTRenergyLoss
{
public:

  XTRTransparentRegRadModel (G4LogicalVolume *anEnvelope,G4Material*,G4Material*,
                        G4double,G4double,G4int,
                        const G4String & processName = "XTRTransparentRegRadModel");
  ~XTRTransparentRegRadModel ();

  // reimplementation of base class function in analytical way

  G4double SpectralXTRdEdx(G4double energy);

  // Pure virtual function from base class

  G4double GetStackFactor( G4double energy, G4double gamma, G4double varAngle);
};

#endif
