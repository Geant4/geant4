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
// $Id: G4DatabaseTypes.hh,v 1.5 2001/07/11 10:02:25 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
//

// file description:
//
//      This file defines types of database and transaction mode
//      supported by G4TransactionManager.

// History:
// 99.11.26 Y.Morita  Initial version

#ifndef G4DatabaseTypes_HH
#define G4DatabaseTypes_HH 1

enum ETypeOfDB { kRunDB, kEventDB, kGeomDB };
enum ETransactionMode { kUpdate, kRead };
enum ESustainedState {
               kAlreadySelected,
               kCannotOverride,
               kStartNewSustained,
               kCommitAndStartNonSustained,
               kStartNonSustained };

#endif

