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
// $Id: OlapVisManager.hh,v 1.1 2002-06-04 07:40:20 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// OlapVisManager
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#ifndef OLAPVISMANAGER_HH
#define OLAPVISMANAGER_HH

#ifdef G4VIS_USE
#include "G4VisManager.hh"

class OlapVisManager: public G4VisManager
{

public:

  OlapVisManager (G4int verboseLevel = 0);
  // Controls initial verbose level of VisManager and VisMessenger.
  // Can be changed by /vis/set/verbose.

private:

  virtual void RegisterGraphicsSystems ();

};

#endif
#endif
