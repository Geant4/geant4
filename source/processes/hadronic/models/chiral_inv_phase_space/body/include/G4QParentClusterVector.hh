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
// $Id: G4QParentClusterVector.hh,v 1.19 2009-02-23 09:49:24 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QCandidateVector ----------------
//             by Mikhail Kossov, Sept 1999.
// Type defenition of Parent nuclear cluster Vector in CHIPS model
// ---------------------------------------------------------------
// Short description: The parent cluster is the cluster, which can be
// used for the nuclear fragment production. Different clusters csn be
// used as the parent cluser for the particular G4QCandidate (nuclear
// fragment), e.g. t and He3 for the t-fragment production. So the
// G4QParentClusterVector is needed.
// -------------------------------------------------------------------

#ifndef G4QParentClusterVector_h
#define G4QParentClusterVector_h 1

#include "G4QParentCluster.hh"
#include <vector>

typedef std::vector<G4QParentCluster *> G4QParentClusterVector;
struct DeleteQParentCluster{ void operator()(G4QParentCluster *aN){delete aN;} };

#endif
