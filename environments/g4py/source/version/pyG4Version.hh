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
// $Id: pyG4Version.hh,v 1.3 2006-06-04 21:34:29 kmura Exp $
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

