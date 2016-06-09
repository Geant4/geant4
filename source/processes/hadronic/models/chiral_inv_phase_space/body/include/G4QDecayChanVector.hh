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
// $Id: G4QDecayChanVector.hh,v 1.18 2005/06/04 13:08:23 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//      ---------------- G4QCandidateVector ----------------
//             by Mikhail Kossov, Sept 1999.
// Type defenition for Decay Channel Vector in CHIPS model
// --------------------------------------------------------------

#ifndef G4QDecayChanVector_h
#define G4QDecayChanVector_h 1

#include "G4QDecayChan.hh"
#include <vector>

typedef std::vector<G4QDecayChan *> G4QDecayChanVector;
struct DeleteQDecayChan
{
  void operator()(G4QDecayChan *aQ)
  {
    //G4cout<<"G4QDecayChanVector::DeleteQDecayChan: Before aQ="<<aQ<<G4endl; // TMP
    if(aQ) delete aQ;
    else G4cerr<<"***G4QDecayChanVector::DeleteQDecayChan: aQ="<<aQ<<G4endl;
  }
};

#endif
