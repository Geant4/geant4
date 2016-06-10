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
// $Id: G4VStoreNotifier.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
// class G4VStoreNotifier
//
// Class description:
//
// Simple abstract class allowing for implementation of user notifiers
// to be activated at registration/deregistration of objects in the
// volume, solid and region stores.
// See G4VNotifier for the details.

// Author:
// 01.09.04 G.Cosmo Initial version
// --------------------------------------------------------------------
#ifndef G4VSTORENOTIFIER_HH
#define G4VSTORENOTIFIER_HH

#include "G4VNotifier.hh"

typedef G4VNotifier G4VStoreNotifier;

#endif
