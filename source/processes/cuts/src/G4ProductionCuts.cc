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
// G4ProductionCuts class implementation
//
// Author: H.Kurashige, 17 September 2002 - First implementation
// --------------------------------------------------------------------

#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"

#include <sstream>

G4ProductionCuts::G4ProductionCuts()
{
  for (G4int i=0; i<NumberOfG4CutIndex; ++i)
  {
    fRangeCuts.push_back(0.0);
  }
}

G4ProductionCuts::G4ProductionCuts(const G4ProductionCuts& right) 
{
  for (G4int i=0; i<NumberOfG4CutIndex; ++i)
  {
    fRangeCuts.push_back(0.0);
  }
  *this = right;
}

G4ProductionCuts::~G4ProductionCuts()
{
  fRangeCuts.clear();
}

G4ProductionCuts& G4ProductionCuts::operator=(const G4ProductionCuts& right)
{
  if (&right==this) return *this;

  for (G4int i=0; i<NumberOfG4CutIndex; ++i)
  {
    fRangeCuts[i] = right.fRangeCuts[i];
  }
  isModified = right.isModified;
  return *this;
}

G4bool G4ProductionCuts::operator==(const G4ProductionCuts& right) const
{
  return (this == &right);
}

G4bool G4ProductionCuts::operator!=(const G4ProductionCuts& right) const
{
  return (this != &right);
}

void G4ProductionCuts::SetProductionCut(G4double cut, G4int index)
{
  if(index >= 0 && index < NumberOfG4CutIndex)
  {
    fRangeCuts[index] = cut;
    isModified = true;
  }
  else
  {
    std::ostringstream os;
    os << "Setting cuts for particles other than photon, e-, e+ or proton has "
          "no effect.";
    G4Exception("G4ProductionCuts::SetProductionCut", "ProcCuts110",
                JustWarning, os.str().c_str());
  }
}

void G4ProductionCuts::SetProductionCut(G4double cut)
{
  for(G4int i = 0; i < NumberOfG4CutIndex; ++i)
  {
    fRangeCuts[i] = cut;
  }
  isModified = true;
}

void G4ProductionCuts::SetProductionCut(G4double cut, G4ParticleDefinition* ptr)
{
  SetProductionCut(cut, GetIndex(ptr));
}

void G4ProductionCuts::SetProductionCut(G4double cut, const G4String& pName)
{
  SetProductionCut(cut, GetIndex(pName));
}

G4double G4ProductionCuts::GetProductionCut(G4int index) const
{
  G4double cut = -1.0;
  if (index>=0 && index<NumberOfG4CutIndex)
  {
    cut = fRangeCuts[index];
  }
  return cut;
}

G4double G4ProductionCuts::GetProductionCut(const G4String& name) const
{
  return GetProductionCut(GetIndex(name));
}


const std::vector<G4double>&   G4ProductionCuts::GetProductionCuts() const
{
  return fRangeCuts;
}

G4bool G4ProductionCuts::IsModified() const
{
  return isModified;
}

void G4ProductionCuts::PhysicsTableUpdated()
{
  isModified = false;
}

G4int G4ProductionCuts::GetIndex(const G4String& name)
{
  G4int index = -1;
  if ( name == "gamma" ) { index = 0; }
  else if ( name == "e-" ) { index = 1; }
  else if ( name == "e+" ) { index = 2; }
  else if ( name == "proton" ) { index = 3; }

  return index;
}


G4int G4ProductionCuts::GetIndex(const G4ParticleDefinition* ptr)
{
  G4int pdg = (nullptr == ptr) ? 0 : ptr->GetPDGEncoding();
  G4int index = -1;
  if (pdg == 22) { index = 0; }
  else if (pdg == 11) { index = 1; }
  else if (pdg == -11) { index = 2; }
  else if (pdg == 2212) { index = 3; }

  return index;
}

void G4ProductionCuts::SetProductionCuts(std::vector<G4double>& cut)
{  
  G4int vSize = (G4int)cut.size();
  if (vSize != NumberOfG4CutIndex)
  {
#ifdef G4VERBOSE
    if (G4ProductionCutsTable::GetProductionCutsTable()->GetVerboseLevel() > 1)
    {
      G4cout << "G4ProductionCuts::SetProductionCuts ";
      G4cout << " The size of given cut value vector [=" << vSize << "]  "
             << " is not consistent with number of CutIndex [="  
             << NumberOfG4CutIndex << G4endl;
    }
#endif
    G4Exception( "G4ProductionCuts::SetProductionCuts ",
                 "ProcCuts108",
                 JustWarning, "Given vector size is inconsistent ");
    if (NumberOfG4CutIndex<vSize)  { vSize = NumberOfG4CutIndex; }
  }
  for(G4int i = 0; i<vSize; ++i)
  {
    fRangeCuts[i] = cut[i];
  }
  isModified = true;
}
