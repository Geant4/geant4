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
// $Id: AnaEx01RunAction.hh,v 1.4 2001-07-11 09:57:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef AnaEx01RunAction_h
#define AnaEx01RunAction_h

#include "G4UserRunAction.hh"

class AnaEx01AnalysisManager;

class AnaEx01RunAction : public G4UserRunAction {
public:
  AnaEx01RunAction(AnaEx01AnalysisManager* = 0);
  ~AnaEx01RunAction();
public:
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);  
private:
  AnaEx01AnalysisManager* fAnalysisManager;
};

#endif

