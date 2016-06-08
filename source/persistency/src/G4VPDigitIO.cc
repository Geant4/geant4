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
// File: G4VPDigitIO.cc
//
// History:
//   '01.08.10  Youhei Morita  Initial creation (with "fadsclass3")

#include "G4VPDigitIO.hh"

G4VPDigitIO* G4VPDigitIO::f_G4VPDigitIO = 0;

// Implementation of Constructor #1
G4VPDigitIO::G4VPDigitIO()
 : m_verbose(0)
{
  f_catalog = G4DCIOcatalog::GetDCIOcatalog();
}

// Implementation of SetVerboseLevel
void G4VPDigitIO::SetVerboseLevel(int v)
{
  m_verbose = v;

  // Loop through the registered Digit I/O managers
  for ( size_t i=0; i < f_catalog->NumberOfDCIOmanager(); i++ ) {
    G4VPDigitsCollectionIO* digitIOman = f_catalog->GetDCIOmanager(i);
    digitIOman->SetVerboseLevel(v);
  }
}

// End of G4VPDigitIO.cc

