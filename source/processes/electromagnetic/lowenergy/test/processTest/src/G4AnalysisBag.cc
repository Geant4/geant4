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
// $Id: G4AnalysisBag.cc,v 1.2 ????
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Authors: MGP
//
// History:
// -----------
// 
//
// -------------------------------------------------------------------

#include "G4AnalysisBag.hh"
#include "G4Ntuple.hh"

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

G4AnalysisBag::G4AnalysisBag()
{
  hbookManager = new HBookFile("processTest.hbook", 39);
}

G4AnalysisBag::~G4AnalysisBag()
{ 
  delete hbookManager;

  G4int n = ntuples.size();
  for (G4int i=0; i<n; ++i)
    {
      delete ntuples[i];
      ntuples[i] = 0;
    } 
}
 
G4AnalysisBag* G4AnalysisBag::instance = 0;

G4AnalysisBag* G4AnalysisBag::getInstance()
{
  if (instance == 0)
    {
      instance = new G4AnalysisBag;
     
    }
  return instance;
}

void G4AnalysisBag::init(const G4String& file)
{
  hbookManager = new HBookFile(file, 39);
}

void G4AnalysisBag::addNtuple(G4Ntuple* ntuple)
{
  ntuples.push_back(ntuple);
}

const G4Ntuple* G4AnalysisBag::retrieveNtuple(G4int id) const
{
  G4Ntuple* nt = ntuples[id];
  return nt;
}

void G4AnalysisBag::write()
{
  hbookManager->write();
}


const HepTupleManager* G4AnalysisBag::getManager() const 
{
  return hbookManager;
}
