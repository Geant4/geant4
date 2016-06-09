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
// $Id: AnaEx01AnalysisManager.hh,v 1.6 2003/06/20 14:55:44 gbarrand Exp $
// GEANT4 tag $Name: geant4-05-02-patch-01 $
//
// 

#ifndef AnaEx01AnalysisManager_h
#define AnaEx01AnalysisManager_h 1

#ifdef G4ANALYSIS_USE

class G4Run;
class G4Event;
class G4Step;

namespace AIDA {
  class IAnalysisFactory;
  class ITree;
  class IHistogram1D;
  class ITuple;
}

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
  AIDA::IAnalysisFactory* fAnalysisFactory;
  AIDA::ITree* fTree;
  AIDA::IHistogram1D* fEAbs;
  AIDA::IHistogram1D* fLAbs;
  AIDA::IHistogram1D* fEGap;
  AIDA::IHistogram1D* fLGap;
  AIDA::ITuple* fTuple;
};

#endif

#endif
