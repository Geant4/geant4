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
// $Id: G4ClassificationOfNewTrack.hh,v 1.5 2001-07-19 00:14:15 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

