// $Id: pyG4Version.hh,v 1.1 2006-02-27 09:56:05 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Version.hh
//
//                                         2005 Q
// ====================================================================
#ifndef PYG4_VERSION_H
#define PYG4_VERSION_H

// Geant4 version
#if   G4VERSION_NUMBER == 700
#include "G4Version_7.0.hh"

#elif G4VERSION_NUMBER == 701
#include "G4Version_7.0.p01.hh"

#elif G4VERSION_NUMBER == 710
#include "G4Version_7.1.hh"

#elif G4VERSION_NUMBER == 711
#include "G4Version_7.1.p01.hh"

#else
#include "G4Version.hh"
#endif

#endif

