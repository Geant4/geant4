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
// $Id: G4FieldTrack.cc,v 1.4 2001-07-11 09:59:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4FieldTrack.hh"

G4std::ostream& operator<<( G4std::ostream& os, G4FieldTrack& SixVec)
{
     G4double *SixV = SixVec.SixVector;
     os << " X= " << SixV[0] << " " << SixV[1] << " " << SixV[2] << " ";
     os << " V= " << SixV[3] << " " << SixV[4] << " " << SixV[5] << " ";
     os << " l= " << SixVec.GetCurveLength();
     return os;
}
