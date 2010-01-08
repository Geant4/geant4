#ifndef G4DINEUTRON_HH
#define G4DINEUTRON_HH
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
// $Id: G4Dineutron.hh,v 1.1 2010-01-08 23:19:41 mkelsey Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      Bertini Cascade dineutron class header file
//
//      History: first implementation, inspired by G4Proton
//      17 Nov 2009:  Michael Kelsey
// ----------------------------------------------------------------

#include "G4Ions.hh"

// ######################################################################
// ###                        DINEUTRON                                ###
// ######################################################################

class G4Dineutron : public G4Ions {
private:
  static G4Dineutron* theInstance;
  G4Dineutron() {}
  ~G4Dineutron() {}
  
public:
  static G4Dineutron* Definition();
  static G4Dineutron* DineutronDefinition();
  static G4Dineutron* Dineutron();
};

#endif	/* G4DINEUTRON_HH */
