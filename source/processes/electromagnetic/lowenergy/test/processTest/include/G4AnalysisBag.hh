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
// $Id: G4AnalysisBag.hh,v 1.2 ????
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Authors: MGP
//
// History:
// -----------
//  

//
// -------------------------------------------------------------------

// Class description:

// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4ANALYSISBAG_HH
#define G4ANALYSISBAG_HH 1

#include "g4std/map"
#include "g4std/vector"
#include "globals.hh"

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

class IHistoManager;
class G4Ntuple;
class NTupleFactory;

class G4AnalysisBag {

public: 

  static G4AnalysisBag* getInstance();

  ~G4AnalysisBag();

  void init(const G4String& file);

  void addNtuple(G4Ntuple* ntuple, const G4String& name);

  void addHisto1D(IHistogram* histo, const G4String& name);

  const G4Ntuple* getNtuple(const G4String& name) const;

  const IHistogram* getHisto1D(const G4String& name) const;

  void store();

  const IHistogramManager* getManager() const { return histoManager; }

private:
 
  G4AnalysisBag();

  static G4AnalysisBag* instance;
  
  G4std::vector<G4Ntuple*> ntuples;
  IHistoManager* histoManager;
  NTupleFactory* factory;

};

#endif
