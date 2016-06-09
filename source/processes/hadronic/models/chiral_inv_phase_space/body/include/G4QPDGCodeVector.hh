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
// $Id: G4QPDGCodeVector.hh,v 1.15 2003/12/09 15:38:09 gunter Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
//      ---------------- G4QCandidateVector ----------------
//             by Mikhail Kossov, Sept 1999.
// Type defenition for Hadron definition in CHIPS model
// ---------------------------------------------------------------

#ifndef G4QPDGCodeVector_h
#define G4QPDGCodeVector_h 1

#include "G4QPDGCode.hh"
#include <vector>

typedef std::vector<G4QPDGCode *> G4QPDGCodeVector;
struct DeleteQPDGCode
{
  void operator()(G4QPDGCode *aN)
  {
    //G4cout<<"G4QPDGCodeVector::DeleteQPDGCode: Before aN="<<aN<<G4endl; // TMP
    if(aN) delete aN; // void Destructor
    else G4cerr<<"***G4QPDGCodeVector::DeleteQPDGCode: aN="<<aN<<G4endl; // TMP
  }
};


#endif


