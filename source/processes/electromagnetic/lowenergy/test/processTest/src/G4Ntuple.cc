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
// $Id: G4Ntuple.cc,v 1.2 ????
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Authors: AP & MGP 
//
// History:
// -----------
// 
//
// -------------------------------------------------------------------

#include "G4Ntuple.hh"

// New Histogramming (from AIDA and Anaphe):
#include "Interfaces/IHistoManager.h"
#include "G4AnalysisBag.hh"
// For NtupleTag from Anaphe
#include "NtupleTag/LizardNTupleFactory.h"
using namespace Lizard;


G4Ntuple::G4Ntuple()
{
  ntuple = 0;
  name = "";
}

G4Ntuple::~G4Ntuple()
{ 
  //  delete ntuple;
}
 
void G4Ntuple::Book(const G4String& name)
{
  G4AnalysisBag* container = G4AnalysisBag::getInstance(); 
  NTupleFactory* factory = container->getFactory();  
  const char* name2(name);
  ntuple = factory->createC(name2);

  //NTupleFactory* factory = container->GetNtupleFactory();

  // Next create the nTuples using the factory and open it for writing
  // ntuple-name is composition of <fileName>:<dirName>:<ntupleID>
  //NTuple* ntuple = factory->createC( "testhisto1.hbook::1" );
  // Check if successful
  //  assert ( ntuple != 0 );}

  // ---- primary ntuple ------
  //-old  hbookManager = new HBookFile("comptontest.hbook", 58);
  //-old Tuple* ntuple1 = hbookManager->ntuple("Primary Ntuple");

  // NOTE:
  // Presently (Anaphe-3.6.3) each ntuple needs to be in a separate
  // file (and different from the histograms). This will be fixed
  // as of the next release.
}

// G4bool G4Ntuple::AddAndBind(const G4String& attribute, const Quantity& quantity)
G4bool G4Ntuple::AddAndBind(const G4String& attribute, float quantity)
{
  G4bool ok = true;
  // Old
  attributes.push_back(attribute);
  // New
  // ntuple->addAndBind(attribute,quantity);
  return ok;
}

void G4Ntuple::AddRow(const G4DataVector& row)
{
  ntuple->addRow();
}



