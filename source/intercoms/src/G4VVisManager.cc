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
// $Id: G4VVisManager.cc,v 1.5 2002-11-20 14:46:01 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Abstract interface for GEANT4 Visualization Manager.
// John Allison 19/Oct/1996.

#include "G4VVisManager.hh"

G4VVisManager::~G4VVisManager () {}

G4VVisManager* G4VVisManager::fpConcreteInstance = 0;

G4VVisManager* G4VVisManager::GetConcreteInstance ()
{
  return fpConcreteInstance;
}

void G4VVisManager::SetConcreteInstance (G4VVisManager* m)
{
  fpConcreteInstance = m;
}
