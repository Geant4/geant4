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
// $Id: G4QuadrupoleMagField.hh,v 1.3 2003/10/31 14:35:52 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
// class G4QuadrupoleMagField
//
// Class description:
//
// Class for creation of quadrupole magnetic field
//   fGradient - is the gradient value for quadrupole magnetic lense.
// Then the magnetic field components are:
//   Bx = B[0] = fGradient*X ,
//   By = B[1] = fGradient*Y ,
//   Bz = B[2] = 0 .
// Here X,Y,Z are the coordinates of a space point of interest.

// History:
// 3.2.97 - V.Grichine, created.
// -------------------------------------------------------------------

#ifndef G4QUADRUPOLEMAGFIELD_HH
#define G4QUADRUPOLEMAGFIELD_HH

#include "G4MagneticField.hh"

class G4QuadrupoleMagField : public G4MagneticField
{
  public: // with description

    G4QuadrupoleMagField(G4double pGradient);
   ~G4QuadrupoleMagField();

    void GetFieldValue(const G4double yTrack[],
                             G4double B[]     ) const;
  private:

    G4double fGradient;
};

#endif
