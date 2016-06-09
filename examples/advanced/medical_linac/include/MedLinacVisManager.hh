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
// $Id: MedLinacVisManager.hh,v 1.2 2004/04/02 17:48:41 mpiergen Exp $
//
//
// Code developed by: M. Piergentili

//*********************************************************************

#ifndef MedLinacVisManager_h
#define MedLinacVisManager_h 1

#ifdef G4VIS_USE

#include "G4VisManager.hh"

//*********************************************************************

class MedLinacVisManager: public G4VisManager {

public:

  MedLinacVisManager ();

private:

  void RegisterGraphicsSystems ();

};

//*********************************************************************

#endif

#endif
