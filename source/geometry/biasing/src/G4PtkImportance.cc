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
// $Id: G4PtkImportance.cc,v 1.4 2002-04-09 16:23:50 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4PtkImportance.cc
//
// ----------------------------------------------------------------------

#include "G4PtkImportance.hh"
#include "G4PStepStream.hh"

G4std::ostream& operator<<(G4std::ostream &out, const G4PtkImportance &ptki)
{
  for (G4PtkImportance::const_iterator it = ptki.begin();
       it != ptki.end(); it++) {
    out << (*it).first << ", importance = ";
    out << (*it).second << "\n";
  }
  return out;
}
