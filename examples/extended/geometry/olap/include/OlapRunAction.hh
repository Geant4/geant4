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
// $Id: OlapRunAction.hh,v 1.1 2002-06-04 07:40:19 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// OlapRunAction
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#ifndef OlapRunAction_h
#define OlapRunAction_h

#include "g4std/vector"
#include "g4std/fstream"
#include "G4UserRunAction.hh"
#include "OlapEventAction.hh"
#include "OlapLogManager.hh"

class OlapRunAction : public G4UserRunAction
{

public:

   OlapRunAction();
   ~OlapRunAction();
   
   void BeginOfRunAction(const G4Run* aRun);
   void EndOfRunAction(const G4Run* aRun);
   void DrawOlaps();
   
   G4std::vector<OlapInfo *> theOlaps;

private:

   G4std::ofstream FILE;
};
#endif
