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
// $Id: G4QDecayChanVector.hh,v 1.20 2009-02-23 09:49:24 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QCandidateVector ----------------
//             by Mikhail Kossov, Sept 1999.
// Type defenition for Decay Channel Vector in CHIPS model
// --------------------------------------------------------------
// Short description: In the CHIPS World the particles (G4QParticle)
// are defined. For unstable particles there is a G4QDecayChannelVector
// which describes different channels of decay for the particular
// particle. So the G4QDecayChannel class is the class for the description
// of such a decay channel in two or three particles (the secondaries can
// be unstable too and have firther decay).
// -------------------------------------------------------------------

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
