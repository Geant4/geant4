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
// $Id: G3SensVolVector.hh,v 1.6 2003/06/16 16:50:43 gunter Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// vector of logical volumes that were defined with
// tracking medium with ISVOL=1
//
// by I.Hrivnacova, 27 Sep 99

#ifndef G3SENSVOLVECTOR_HH
#define G3SENSVOLVECTOR_HH 1

#include <vector>
#include "G4LogicalVolume.hh"

typedef std::vector<G4LogicalVolume*> G3SensVolVector;

extern G3SensVolVector G3SensVol;
#endif
