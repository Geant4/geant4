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
// $Id: G4FieldTrack.cc,v 1.7 2004-01-13 12:27:05 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------

#include "G4FieldTrack.hh"

std::ostream& operator<<( std::ostream& os, const G4FieldTrack& SixVec)
{
     const G4double *SixV = SixVec.SixVector;
     os << " ( ";
     os << " X= " << SixV[0] << " " << SixV[1] << " " << SixV[2] << " ";  // Position
     os << " V= " << SixV[3] << " " << SixV[4] << " " << SixV[5] << " ";  // Momentum
     os << " v2= " << G4ThreeVector(SixV[3], SixV[4], SixV[5]).mag();     // mom magnitude
     os << " mdm= " << SixVec.fMomentumDir.mag(); 
     os << " l= " << SixVec.GetCurveLength();
     os << " ) ";
     return os;
}
