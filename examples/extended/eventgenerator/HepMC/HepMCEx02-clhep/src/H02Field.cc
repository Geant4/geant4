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
//   H02Field.hh
//   $Id: H02Field.cc,v 1.1 2002-11-19 10:36:20 murakami Exp $
//
// ====================================================================

#include "H02Field.hh"

////////////////////////////////////////////////////////////////////////////
void H02Field::GetFieldValue(const G4double Point[3], G4double* Bfield) const
////////////////////////////////////////////////////////////////////////////
{
  const G4double Bz= 3.0*tesla;
  const G4double rmax_sq = sqr(1.*m);
  const G4double zmax = 1.*m;

  Bfield[0]= 0.;
  Bfield[1] = 0.;
  if(abs(Point[2])<zmax && (sqr(Point[0])+sqr(Point[1]))<rmax_sq) {
    Bfield[2]= Bz;
  } else { 
    Bfield[2]= 0.; 
  }
}

