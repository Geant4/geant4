//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
// GEANT4 tag $Name: geant4-09-04-ref-00 $
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file/  History:
//    18 Sep. 2002, H.Kuirashige : Structure created based on object model
// --------------------------------------------------------------

#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include <iomanip>

const G4ParticleDefinition* G4ProductionCuts::gammaDef = 0;
const G4ParticleDefinition* G4ProductionCuts::electDef = 0;
const G4ParticleDefinition* G4ProductionCuts::positDef = 0;
const G4ParticleDefinition* G4ProductionCuts::protonDef = 0;

G4ProductionCuts::G4ProductionCuts() :
  isModified(true)
{
  for (G4int i=0; i< NumberOfG4CutIndex; i++) {
    fRangeCuts.push_back(0.0);
  }
}

G4ProductionCuts::G4ProductionCuts(const G4ProductionCuts& right) 
{
  *this = right;
}

G4ProductionCuts::~G4ProductionCuts()
{
  fRangeCuts.clear();
}

G4ProductionCuts & G4ProductionCuts::operator=(const G4ProductionCuts &right)
{
  if (&right==this) return *this;

  for (G4int i=0; i< NumberOfG4CutIndex; i++) {
    fRangeCuts[i] = right.fRangeCuts[i];
  }
  isModified = right.isModified;
  return *this;
}



G4int G4ProductionCuts::operator==(const G4ProductionCuts &right) const
{
  return (this == &right);
}


G4int G4ProductionCuts::operator!=(const G4ProductionCuts &right) const
{
  return (this !=  &right);
}


G4int  G4ProductionCuts::GetIndex(const G4String& name)
{
  static G4String gamma("gamma");
  static G4String electron("e-");
  static G4String positron("e+");
  static G4String proton("proton");
  
  G4int index;
  if       ( name == gamma )        { index =  0; }
  else  if ( name == electron )     { index =  1; }
  else  if ( name == positron )     { index =  2; }
  else  if ( name == proton )       { index =  3; }
  else                              { index = -1; }

  return index;
}


G4int  G4ProductionCuts::GetIndex(const G4ParticleDefinition* ptcl)
{ 
  if(!ptcl) return -1;
  // In the first call, pointers are set 
  if(gammaDef==0  && ptcl->GetParticleName()=="gamma")  { gammaDef = ptcl; }
  if(electDef==0  && ptcl->GetParticleName()=="e-")     { electDef = ptcl; }
  if(positDef==0  && ptcl->GetParticleName()=="e+")     { positDef = ptcl; }
  if(protonDef==0 && ptcl->GetParticleName()=="proton") { protonDef = ptcl; }

  G4int index;
  if(ptcl==gammaDef)       { index = 0;  }
  else if(ptcl==electDef)  { index = 1;  }
  else if(ptcl==positDef)  { index = 2;  }
  else if(ptcl==protonDef) { index = 3;  }
  else                     { index = -1; }

  return index;
}


void  G4ProductionCuts::SetProductionCuts(std::vector<G4double>& cut)
{  
  G4int vSize = cut.size();
  if (vSize != NumberOfG4CutIndex) {
#ifdef G4VERBOSE
    if ( G4ProductionCutsTable::GetProductionCutsTable()->GetVerboseLevel()>1) {
      G4cerr << "G4ProductionCuts::SetProductionCuts ";
      G4cerr << " The size of given cut value vector [=" << vSize << "]  "
	     << " is not consitent with number of CutIndex [="  
	     << NumberOfG4CutIndex << G4endl;
    }
#endif
    G4Exception( "G4ProductionCuts::SetProductionCuts ",
		 "ProcCuts108",
		 JustWarning, "Given vector size is inconsitent ");
    if (NumberOfG4CutIndex<vSize) vSize = NumberOfG4CutIndex;
  }
  for(G4int i = 0; (i<vSize ); i++) {
    fRangeCuts[i] = cut[i];
  }
  isModified = true;
}
