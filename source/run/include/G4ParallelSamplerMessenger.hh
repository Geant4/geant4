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
// $Id: G4ParallelSamplerMessenger.hh,v 1.4 2002-07-29 16:03:15 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelSamplerMessenger
//
// Class description:
// This is a meesenger which provides the comands:
// - /imp/importancesample <particletypename>
//   to set up one G4ParallelImportnaceScoreSampler for the
//   givem <particletypename>
// - /imp/score <particletypename>
//   to create several G4ParallelScoreSamplers using the G4StandardScorer
//   for given <particletypename>s
// - /imp/writefile <outputfilename>
//   to write the score results to the given <outputfilename>
//   using a G4StandardScoreTable
//
// a particle type may be given only once in a command 
//
// The parallel geometry and importance store is obtained from
// an G4ImportanceGeometryConstructor object given
// to the constructor of this class.
// 
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4ParallelSamplerMessenger_hh
#define G4ParallelSamplerMessenger_hh G4ParallelSamplerMessenger_hh

#include "G4UImessenger.hh"
#include "g4std/map"


class G4VSampler;
class G4ImportanceGeometryConstructor;
class G4UIcmdWithAString;
class G4ParallelImportanceScoreSampler;
class G4ParallelScoreSampler;
class G4VPrintableScorer;
class G4VPrintableScorerFactory;

typedef G4std::map<G4String, G4VSampler *> G4MapNameSampler;
typedef G4std::map<G4String, G4VPrintableScorer *> G4MapNameScorer;

class G4ParallelSamplerMessenger : public G4UImessenger {
public:
  G4ParallelSamplerMessenger(G4ImportanceGeometryConstructor *igc,
			     const G4VPrintableScorerFactory &psf);
    // constructor needs G4ImportanceGeometryConstructor to get
    // a parallel geometry and the importance store

  ~G4ParallelSamplerMessenger(){}
  
  void SetNewValue(G4UIcommand * command, G4String newValue);
    //  from G4UImessenger

private:
  void CheckNewParticle(const G4String &particlename);
    // Error if a particle type is given more than once

  void Error(const G4String &m){
    G4Exception("Error: G4ParallelSamplerMessenger:" + m);
  }

  void ImpCmdAction(const G4String &particlename);
  void ScoreCmdAction(const G4String &particlename);
  void PrintCmdAction(const G4String &particlename);
  

  G4ImportanceGeometryConstructor *fImpGeoConst;
  const G4VPrintableScorerFactory &fPrintableScorerFactory;
  
  G4UIcmdWithAString *fImpCmd;
  G4UIcmdWithAString *fScoreCmd;
  G4UIcmdWithAString *fPrintCmd;

  G4ParallelImportanceScoreSampler *fImpSampler;
  G4VPrintableScorer *fImpScorer;
  G4String fImpPartilce;
  G4MapNameScorer fMapNameScorer;
  G4MapNameSampler fMapNameSampler;

};

#endif


