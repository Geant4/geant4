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
// $Id: AnaEx01AnalysisManager.hh,v 1.5 2001-11-16 14:30:50 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef AnaEx01AnalysisManager_h
#define AnaEx01AnalysisManager_h 1

#ifdef G4ANALYSIS_USE

class G4Run;
class G4Event;
class G4Step;

class IAnalysisFactory;
class ITree;
class IHistogram1D;
class ITuple;

class AnaEx01AnalysisManager {
public:
  AnaEx01AnalysisManager();
  virtual ~AnaEx01AnalysisManager();
public:
  virtual void BeginOfRun(const G4Run*); 
  virtual void EndOfRun(const G4Run*); 
  virtual void BeginOfEvent(const G4Event*); 
  virtual void EndOfEvent(const G4Event*); 
  virtual void Step(const G4Step*);
private:
  int fCalorimeterCollID;                
  IAnalysisFactory* fAnalysisFactory;
  ITree* fTree;
  IHistogram1D* fEAbs;
  IHistogram1D* fLAbs;
  IHistogram1D* fEGap;
  IHistogram1D* fLGap;
  ITuple* fTuple;
};

#endif

#endif
