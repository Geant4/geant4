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
//
// $Id: G4Diproton.hh 69389 2013-05-02 15:54:09Z mkelsey $
//
// ------------------------------------------------------------
//      Bertini Cascade diproton class header file
//
//      History: first implementation, inspired by G4Proton
//      17 Nov 2009:  Michael Kelsey
//	06 Apr 2010:  Reset theInstance in dtor, implement ctor in .cc.
//	13 Apr 2010:  Per Kurashige, inherit from G4VShortLivedParticle.
//	01 May 2013:  Remove G4ThreadLocal from static pointer.
// ----------------------------------------------------------------

#ifndef G4DIPROTON_HH
#define G4DIPROTON_HH

#include "G4VShortLivedParticle.hh"

// ######################################################################
// ###                        DIPROTON                                ###
// ######################################################################

class G4Diproton : public G4VShortLivedParticle {
private:
  static G4Diproton* theInstance;
  G4Diproton();
  ~G4Diproton() { theInstance = 0; }
  
public:
  static G4Diproton* Definition();
  static G4Diproton* DiprotonDefinition();
  static G4Diproton* Diproton();
};

#endif	/* G4DIPROTON_HH */
