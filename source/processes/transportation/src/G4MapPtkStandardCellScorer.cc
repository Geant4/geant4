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
// $Id: G4MapPtkStandardCellScorer.cc,v 1.1 2002-07-10 15:51:04 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4MapPtkStandardCellScorer.cc
//
// ----------------------------------------------------------------------

#include "G4MapPtkStandardCellScorer.hh"
#include "G4PStepStream.hh"

G4std::ostream& operator << (G4std::ostream &out,
                             const G4MapPtkStandardCellScorer &Scorestore)
{
  for (G4MapPtkStandardCellScorer::const_iterator it = Scorestore.begin();
       it != Scorestore.end(); it++) {
    out << (*it).first << "\n";
    out << (*it).second;
  }
  return out;
}
