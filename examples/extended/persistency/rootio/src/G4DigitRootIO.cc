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
// $Id: G4DigitRootIO.cc,v 1.4 2002-12-13 14:45:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4DigitRootIO.cc
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#include "G4DigitRootIO.hh"

// Implementation of Constructor #1
G4DigitRootIO::G4DigitRootIO()
{
  f_G4VPDigitIO = (G4VPDigitIO*) this;
}

// Implementation of Destructor #1
G4DigitRootIO::~G4DigitRootIO()
{}

// Implementation of GetDigitRootIO
G4DigitRootIO* G4DigitRootIO::GetDigitRootIO()
{
  G4DigitRootIO* dio=0;
  if ( f_G4VPDigitIO == 0 ) {
    dio = new G4DigitRootIO;
  }
  return dio;
}

// Implementation of Store
bool G4DigitRootIO::Store(const G4DCofThisEvent* dcevt)
{
  return true;
}

// Implementation of Retrieve
bool G4DigitRootIO::Retrieve(G4DCofThisEvent*& dcevt)
{
  return true;
}

// End of G4DigitRootIO.cc

