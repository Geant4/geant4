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
//
//
// 
///////////////////////////////////////////////////////////////////////////
// 
// Model describing a regular radiator of X-ray transition radiation.
// It is Garibian like model for XTR after radiator  
// Thicknesses of plates and gas gaps are fixed.
// We suppose that:
// formation zone ~ mean thickness << absorption length
// for each material and in the range 1-100 keV. This allows us to simplify
// interference effects in radiator stack (GetStackFactor method).
// 
// 
// History:
// 10.10.05 V. Grichine, first version 
//


#ifndef G4XTRRegularRadModel_h
#define G4XTRRegularRadModel_h 1

#include "G4VXTRenergyLoss.hh"


class G4XTRRegularRadModel : public G4XTRenergyLoss
{
public:

   G4XTRRegularRadModel (G4LogicalVolume *anEnvelope,G4Material*,G4Material*,
                        G4double,G4double,G4int,
                        const G4String & processName = "XTRegularModel");
  ~G4XTRRegularRadModel ();

  // Pure virtual function from base class

  G4double GetStackFactor( G4double energy, G4double gamma, G4double varAngle);
};

#endif
