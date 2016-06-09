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
// $Id: G4LineCurrentMagField.hh,v 1.3 2003/10/31 14:35:51 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
//
// class G4LineCurrentMagField
//
// Class description:
//
// Class describing line current magnetic field.
// 'fFieldConstant' determines the coefficient in the field law.
// The line current is directed along Z axis and crosses the XY
// plane in the origin point (0,0).

// History:
// 3.2.97 - V. Grichine, created.
// --------------------------------------------------------------------

#ifndef G4LINECURRENTMAGFIELD_HH
#define G4LINECURRENTMAGFIELD_HH

#include "G4MagneticField.hh"

class G4LineCurrentMagField : public G4MagneticField
{
  public:  // with description

    G4LineCurrentMagField(G4double pFieldConstant);
   ~G4LineCurrentMagField();

    void GetFieldValue(const G4double yTrack[] ,
                             G4double B[]      ) const;
  private:
  
    G4double fFieldConstant;
};

#endif
