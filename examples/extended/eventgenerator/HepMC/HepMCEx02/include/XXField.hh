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

// ====================================================================
//
//   XXField.hh
//   $Id: XXField.hh,v 1.1 2002-04-29 20:44:35 asaim Exp $
//
// ====================================================================
#ifndef XX_FIELD_H
#define XX_FIELD_H

#include "globals.hh"
#include "G4MagneticField.hh"

class XXField : public G4MagneticField {
public:
  XXField() { }
  ~XXField() { }
  
  void GetFieldValue(const  G4double Point[3],  G4double* Bfield ) const;

};

#endif

