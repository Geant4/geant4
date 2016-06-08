// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DatabaseTypes.hh,v 1.3 1999/11/28 21:54:19 morita Exp $
// GEANT4 tag $Name: geant4-02-00 $
//

// file description:
//
//      This file defines types of database and transaction mode
//      supported by G4TransactionManager.
//
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

