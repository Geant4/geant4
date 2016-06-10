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
// $Id: G4ClassificationOfNewTrack.hh 66892 2013-01-17 10:57:59Z gunter $
//
//

#ifndef G4ClassificationOfNewTrack_hh
#define G4ClassificationOfNewTrack_hh 1

// class description:
//
//  This header file contain an enumeration for the possible classifications
// for trackes newly pushed to the stack. G4UserStackingAction can set the
// classification.

enum G4ClassificationOfNewTrack
{ 
  fUrgent=0,     // put into the urgent stack
  fWaiting=1,    // put into the waiting stack
  fPostpone=-1,  // postpone to the next event
  fKill=-9,      // kill
  //----------------------------------------------------------------
  // following ENUM are available only if the user increases the
  // number of waiting stacks
  //----------------------------------------------------------------
  fWaiting_1=11, fWaiting_2=12, fWaiting_3=13, fWaiting_4=14, fWaiting_5=15,
  fWaiting_6=16, fWaiting_7=17, fWaiting_8=18, fWaiting_9=19, fWaiting_10=20
};

#endif

