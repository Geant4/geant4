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
//
// $Id: G4ErrorTrajErr.hh,v 1.2 2007-05-29 14:41:35 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Class Description:
//
// Trajectory error. Implemented for the moment as a CLHEP HepSymMatrix,
// until other possibly faster implementations are studied

// History:
// - Created:   P. Arce
// --------------------------------------------------------------------

#ifndef G4ErrorTrajErr_hh
#define G4ErrorTrajErr_hh

#include <CLHEP/Matrix/SymMatrix.h>

typedef CLHEP::HepSymMatrix G4ErrorTrajErr;

#endif
