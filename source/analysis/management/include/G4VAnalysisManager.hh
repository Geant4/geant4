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
// $Id: G4VAnalysisManager.hh,v 1.5.4.1 2001/06/28 19:07:47 gunter Exp $
// GEANT4 tag $Name:  $
//
// Guy Barrand 20th Mai 2000
//
// Class description
//
//  Abstract base class for the analysis manager proxy mechanism.
//

#ifndef G4VANALYSIS_HH
#define G4VANALYSIS_HH

#if defined(G4ANALYSIS_BUILD) || defined(G4ANALYSIS_USE)

#include "globals.hh"

class IHistogram;
class IHistogramFactory;

#ifdef G4ANALYSIS_NON_AIDA
class ICloudFactory;
class ICloud;
class ITuple;
#endif

class G4VAnalysisSystem;

class G4VAnalysisManager {
public:
  virtual ~G4VAnalysisManager() {}
  virtual G4bool RegisterAnalysisSystem(G4VAnalysisSystem*) = 0;
  //
  // Get the factories :
  virtual IHistogramFactory* GetHistogramFactory(const G4String&) = 0;
  //
  // Access to facilites methods :
  virtual void Store(IHistogram* = 0,const G4String& = "") = 0;
  virtual void Plot(IHistogram*) = 0;
  //
#ifdef G4ANALYSIS_NON_AIDA
  virtual ICloudFactory* GetCloudFactory(const G4String&) = 0;
  virtual ITuple* CreateTuple(const G4String&,const G4String&,const G4String&) = 0;
  virtual void Plot(ICloud*) = 0;
#endif
};

#endif

#endif
