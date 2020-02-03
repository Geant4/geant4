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
// G4LineCurrentMagField
//
// Class description:
//
// Class describing line current magnetic field.
// 'fFieldConstant' determines the coefficient in the field law.
// The line current is directed along Z axis and crosses the XY
// plane in the origin point (0,0).

// Author: V.Grichine, 03.02.1997
// --------------------------------------------------------------------
#ifndef G4LINECURRENTMAGFIELD_HH
#define G4LINECURRENTMAGFIELD_HH

#include "G4MagneticField.hh"

class G4LineCurrentMagField : public G4MagneticField
{
  public:  // with description

    G4LineCurrentMagField(G4double pFieldConstant);
   ~G4LineCurrentMagField();

    void GetFieldValue(const G4double yTrack[],
                             G4double B[] ) const;
    G4Field* Clone() const;

  private:
  
    G4double fFieldConstant = 0.0;
};

#endif
