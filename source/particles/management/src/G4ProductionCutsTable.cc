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
// $Id: G4ProductionCutsTable.cc,v 1.4 2003-01-14 22:26:50 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file/  History:
//    06/Oct. 2002, M.Asai : First implementation
// --------------------------------------------------------------

#include "G4ProductionCutsTable.hh"
#include "G4ProductionCuts.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4LogicalVolume.hh"
#include "G4RToEConvForAntiProton.hh"
#include "G4RToEConvForElectron.hh"
#include "G4RToEConvForGamma.hh"
#include "G4RToEConvForPositron.hh"
#include "G4RToEConvForProton.hh"
#include "G4RToEConvForNeutron.hh"
#include "G4RToEConvForAntiNeutron.hh"
#include "G4MaterialTable.hh"
#include "G4Material.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

G4ProductionCutsTable* G4ProductionCutsTable::fG4ProductionCutsTable = 0;

G4ProductionCutsTable* G4ProductionCutsTable::GetProductionCutsTable()
{ 
  if(!fG4ProductionCutsTable)
  { fG4ProductionCutsTable = new G4ProductionCutsTable(); }
  return fG4ProductionCutsTable;
}

G4ProductionCutsTable::G4ProductionCutsTable()
: firstUse(true)
{
  for(size_t i=0;i< NumberOfG4CutIndex;i++)
  {
    rangeCutTable.push_back(new G4CutVectorForAParticle);
    energyCutTable.push_back(new G4CutVectorForAParticle);
    rangeDoubleVector[i] = 0;
    energyDoubleVector[i] = 0;
    converters[i] = 0;
  }
  fG4RegionStore = G4RegionStore::GetInstance();
  defaultProductionCuts = new G4ProductionCuts();
}

G4ProductionCutsTable::G4ProductionCutsTable(const G4ProductionCutsTable& right)
{;}

G4ProductionCutsTable::~G4ProductionCutsTable()
{
  for(CoupleTableIterator itr=coupleTable.begin();itr!=coupleTable.end();itr++)
  { delete (*itr); }
  coupleTable.clear();
  for(size_t i=0;i< NumberOfG4CutIndex;i++)
  {
    delete rangeCutTable[i];
    delete energyCutTable[i];
    delete converters[i];
    if(rangeDoubleVector[i]!=0) delete [] rangeDoubleVector[i];
    if(energyDoubleVector[i]!=0) delete [] energyDoubleVector[i];
  }
}

void G4ProductionCutsTable::UpdateCoupleTable()
{
  if(firstUse)
  {
    if(G4ParticleTable::GetParticleTable()->FindParticle("gamma"))
    { converters[0] = new G4RToEConvForGamma(); }
    if(G4ParticleTable::GetParticleTable()->FindParticle("e-"))
    { converters[1] = new G4RToEConvForElectron(); }
    if(G4ParticleTable::GetParticleTable()->FindParticle("e+"))
    { converters[2] = new G4RToEConvForPositron(); }
    if(G4ParticleTable::GetParticleTable()->FindParticle("proton"))
    { converters[3] = new G4RToEConvForProton(); }
    if(G4ParticleTable::GetParticleTable()->FindParticle("anti_proton"))
    { converters[4] = new G4RToEConvForAntiProton(); }
    if(G4ParticleTable::GetParticleTable()->FindParticle("neutron"))
    { converters[5] = new G4RToEConvForNeutron(); }
    if(G4ParticleTable::GetParticleTable()->FindParticle("anti_neutron"))
    { converters[6] = new G4RToEConvForAntiNeutron(); }
    firstUse = false;
  }

  // Reset "used" flags of all couples
  for(CoupleTableIterator CoupleItr=coupleTable.begin();CoupleItr!=coupleTable.end();CoupleItr++)
  { (*CoupleItr)->SetUseFlag(false); }

  // Update Material-Cut-Couple
  typedef G4std::vector<G4Region*>::iterator regionIterator;
  G4Region* theWorldRegion = *(fG4RegionStore->begin());
  for(regionIterator rItr=fG4RegionStore->begin();rItr!=fG4RegionStore->end();rItr++)
  { 
    ///////////////////if(!((*rItr)->IsModified())) continue;
    G4ProductionCuts* fProductionCut = (*rItr)->GetProductionCuts();
    G4std::vector<G4Material*>::const_iterator mItr = (*rItr)->GetMaterialIterator();
    size_t nMaterial = (*rItr)->GetNumberOfMaterials();

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // The following part of the code should be removed once all EM processes
    // become "Region-aware"
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if((*rItr)==theWorldRegion)
    {
      mItr = G4Material::GetMaterialTable()->begin();
      nMaterial = G4Material::GetMaterialTable()->size();
    }
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // The previous part of the code should be removed once all EM processes 
    // become "Region-aware"
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for(size_t iMate=0;iMate<nMaterial;iMate++)
    {
      //check if this material cut couple has already been made
      G4bool coupleAlreadyDefined = false;
      G4MaterialCutsCouple* aCouple;
      for(CoupleTableIterator cItr=coupleTable.begin();cItr!=coupleTable.end();cItr++)
      {
        if((*cItr)->GetMaterial()==(*mItr) && (*cItr)->GetProductionCuts()==fProductionCut)
        { 
          coupleAlreadyDefined = true;
          aCouple = *cItr;
          break;
        }
      }
      
      //if this combination is new, cleate and register a couple
      if(!coupleAlreadyDefined)
      {
        aCouple = new G4MaterialCutsCouple((*mItr),fProductionCut);
        coupleTable.push_back(aCouple);
        aCouple->SetIndex(coupleTable.size()-1);
      }

      //Set the couple to the proper logical volumes in that region
      aCouple->SetUseFlag();
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // The following part of the code should be removed once all EM processes 
    // become "Region-aware"
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if((*rItr)==theWorldRegion)
      {
        aCouple->SetUseFlag(false);
        G4std::vector<G4Material*>::const_iterator mItr1 = (*rItr)->GetMaterialIterator();
        size_t nMaterial1 = (*rItr)->GetNumberOfMaterials();
        for(size_t iMate1=0;iMate1<nMaterial1;iMate1++)
        {
          if((*mItr1)==aCouple->GetMaterial()) aCouple->SetUseFlag();
          mItr1++;
        }
      }
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // The previous part of the code should be removed once all EM processes 
    // become "Region-aware"
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      G4std::vector<G4LogicalVolume*>::iterator rootLVItr
                         = (*rItr)->GetRootLogicalVolumeIterator();
      size_t nRootLV = (*rItr)->GetNumberOfRootVolumes();
      for(size_t iLV=0;iLV<nRootLV;iLV++)
      {
        //Set the couple to the proper logical volumes in that region
        G4LogicalVolume* aLV = *rootLVItr;
        G4Region* aR = *rItr;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // The following part of the code should be removed once all EM processes 
    // become "Region-aware"
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(aR==theWorldRegion) aR = 0;
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // The previous part of the code should be removed once all EM processes 
    // become "Region-aware"
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ScanAndSetCouple(aLV,aCouple,aR);

        //proceed to the next root logical volume in this region
        rootLVItr++;
      }

      //proceed to next material in this region
      mItr++;
    }
  }

  // Check if sizes of Range/Energy cuts tables are equal to the size of
  // the couple table
  // If new couples are made during the previous procedure, nCouple becomes
  // larger then nTable
  size_t nCouple = coupleTable.size();
  size_t nTable = energyCutTable[0]->size();
  G4bool newCoupleAppears = nCouple>nTable;
  if(newCoupleAppears)
  {
    for(size_t n=nCouple-nTable;n>0;n--)
    {
      for(size_t nn=0;nn< NumberOfG4CutIndex;nn++)
      {
        rangeCutTable[nn]->push_back(-1.);
        energyCutTable[nn]->push_back(-1.);
      }
    }
  }

  // Update RangeEnergy cuts tables
  size_t idx = 0;
  for(CoupleTableIterator cItr=coupleTable.begin();cItr!=coupleTable.end();cItr++)
  {
    G4ProductionCuts* aCut = (*cItr)->GetProductionCuts();
    const G4Material* aMat = (*cItr)->GetMaterial();
    if((*cItr)->IsRecalcNeeded())
    {
      for(size_t ptcl=0;ptcl< NumberOfG4CutIndex;ptcl++)
      {
        G4double rCut = aCut->GetProductionCut(ptcl);
        (*(rangeCutTable[ptcl]))[idx] = rCut;
        // if(converters[ptcl] && (*cItr)->IsUsed())
        if(converters[ptcl])
        { (*(energyCutTable[ptcl]))[idx] = converters[ptcl]->Convert(rCut,aMat); }
        else
        { (*(energyCutTable[ptcl]))[idx] = -1.; }
      }
    }
    idx++;  
  }

  // resize Range/Energy cuts double vectors if new couple is made
  if(newCoupleAppears)
  {
    for(size_t ix=0;ix<NumberOfG4CutIndex;ix++)
    {
      G4double* rangeVOld = rangeDoubleVector[ix];
      G4double* energyVOld = energyDoubleVector[ix];
      if(rangeVOld) delete [] rangeVOld;
      if(energyVOld) delete [] energyVOld;
      rangeDoubleVector[ix] = new G4double[(*(rangeCutTable[ix])).size()];
      energyDoubleVector[ix] = new G4double[(*(energyCutTable[ix])).size()];
    }
  }

  // Update Range/Energy cuts double vectors
  for(size_t ix=0;ix<NumberOfG4CutIndex;ix++)
  {
    for(size_t ixx=0;ixx<(*(rangeCutTable[ix])).size();ixx++)
    {
      rangeDoubleVector[ix][ixx] = (*(rangeCutTable[ix]))[ixx];
      energyDoubleVector[ix][ixx] = (*(energyCutTable[ix]))[ixx];
    }
  }
}

void G4ProductionCutsTable::SetEnergyRange(G4double lowedge, G4double highedge)
{
  G4VRangeToEnergyConverter::SetEnergyRange(lowedge,highedge);
}

G4double  G4ProductionCutsTable::GetLowEdgeEnergy() const
{
  return G4VRangeToEnergyConverter::GetLowEdgeEnergy();
}

G4double G4ProductionCutsTable::GetHighEdgeEnergy() const
{
  return G4VRangeToEnergyConverter::GetHighEdgeEnergy();
}
 

void G4ProductionCutsTable::ScanAndSetCouple(G4LogicalVolume* aLV,G4MaterialCutsCouple* aCouple,G4Region* aRegion)
{
  //Check whether or not this logical volume belongs to the same region
  if((aRegion!=0) && aLV->GetRegion()!=aRegion) return;

  //Check if this particular volume has a material matched to the couple
  if(aLV->GetMaterial()==aCouple->GetMaterial())
  {
    aLV->SetMaterialCutsCouple(aCouple);
  }

  size_t noDaughters = aLV->GetNoDaughters();
  if(noDaughters==0) return;

  //Loop over daughters with same region
  for(size_t i=0;i<noDaughters;i++)
  {
    G4LogicalVolume* daughterLVol = aLV->GetDaughter(i)->GetLogicalVolume();
    ScanAndSetCouple(daughterLVol,aCouple,aRegion);
  }
}

const G4MaterialCutsCouple* 
     G4ProductionCutsTable::GetMaterialCutsCouple(const G4Material* aMat, 
                           const G4ProductionCuts* aCut) const
{
  for(CoupleTableIterator cItr=coupleTable.begin();cItr!=coupleTable.end();cItr++)
  {
    if((*cItr)->GetMaterial()!=aMat) continue;
    if((*cItr)->GetProductionCuts()==aCut) return (*cItr);
  }
  return 0;
}

G4int G4ProductionCutsTable::GetCoupleIndex(const G4MaterialCutsCouple* aCouple) const
{
  G4int idx = 0;
  for(CoupleTableIterator cItr=coupleTable.begin();cItr!=coupleTable.end();cItr++)
  {
    if((*cItr)==aCouple) return idx;
    idx++;
  }
  return -1;
}

G4int G4ProductionCutsTable:: GetCoupleIndex(const G4Material* aMat,
                           const G4ProductionCuts* aCut) const
{
  const G4MaterialCutsCouple* aCouple = GetMaterialCutsCouple(aMat,aCut);
  return GetCoupleIndex(aCouple);
}

void G4ProductionCutsTable::DumpCouples() const
{
  G4cout << G4endl;
  G4cout << "========= Table of registored couples ==============================" << G4endl;
  for(CoupleTableIterator cItr=coupleTable.begin();cItr!=coupleTable.end();cItr++)
  {
    G4MaterialCutsCouple* aCouple = (*cItr);
    G4ProductionCuts* aCut = aCouple->GetProductionCuts();
    G4cout << G4endl;
    G4cout << "Index : " << aCouple->GetIndex() 
           << "     used in the geometry : ";
    if(aCouple->IsUsed()) G4cout << "Yes";
    else                  G4cout << "No ";
    G4cout << "     recalcuration needed : ";
    if(aCouple->IsRecalcNeeded()) G4cout << "Yes";
    else                          G4cout << "No ";
    G4cout << G4endl;
    G4cout << " Material : " << aCouple->GetMaterial()->GetName() << G4endl;
    G4cout << " Range cuts        : " 
           << " gamma " << G4BestUnit(aCut->GetProductionCut("gamma"),"Length")
           << "    e- " << G4BestUnit(aCut->GetProductionCut("e-"),"Length")
           << "    e+ " << G4BestUnit(aCut->GetProductionCut("e+"),"Length")
           << G4endl;
    G4cout << " Energy thresholds : " ;
    // if(!(aCouple->IsUsed()) || aCouple->IsRecalcNeeded()) G4cout << " is not ready to print";
    if(aCouple->IsRecalcNeeded()) G4cout << " is not ready to print";
    else
    G4cout << " gamma " << G4BestUnit((*(energyCutTable[0]))[aCouple->GetIndex()],"Energy")
           << "    e- " << G4BestUnit((*(energyCutTable[1]))[aCouple->GetIndex()],"Energy")
           << "    e+ " << G4BestUnit((*(energyCutTable[2]))[aCouple->GetIndex()],"Energy");
    G4cout << G4endl;
    if(aCouple->IsUsed())
    {
      G4cout << " Region(s) which use this couple : " << G4endl;
      typedef G4std::vector<G4Region*>::iterator regionIterator;
      for(regionIterator rItr=fG4RegionStore->begin();rItr!=fG4RegionStore->end();rItr++)
      {
        G4ProductionCuts* fProductionCut = (*rItr)->GetProductionCuts();
        G4std::vector<G4Material*>::const_iterator mItr = (*rItr)->GetMaterialIterator();
        size_t nMaterial = (*rItr)->GetNumberOfMaterials();
        for(size_t iMate=0;iMate<nMaterial;iMate++)
        {
          if(aCouple->GetMaterial()==(*mItr) &&
             aCouple->GetProductionCuts()==fProductionCut)
          {
            G4cout << "    " << (*rItr)->GetName() << G4endl;
            break;
          }
        }
      }
    }
  }
  G4cout << G4endl;
  G4cout << "====================================================================" << G4endl;
  G4cout << G4endl;
}
           

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//========================= Below has not yet implemented =====================
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// Store cuts and material information in files under the specified directory.
G4bool  G4ProductionCutsTable::StoreCutsTable(const G4String& dir, 
					      G4bool          ascii)
{
  if (!StoreMaterialInfo(dir, ascii)) return false;
  if (!StoreMaterialCutsCoupleInfo(dir, ascii)) return false;
  if (!StoreCutsInfo(dir,ascii)) return false;

  return true;
}
  
// Retrieve cuts values information in files under the specified directory.
G4bool  G4ProductionCutsTable::RetrieveCutsTable(const G4String& dir,
						 G4bool          ascii)
{
  if (!CheckForRetrieveCutsTable(dir, ascii)) return false;
  
  if (!RetrieveCutsInfo(dir,ascii)) return false;


  return true;
}

// check stored material and cut values are consistent with the current detector setup. 
G4bool G4ProductionCutsTable::CheckForRetrieveCutsTable(const G4String& directory, 
							G4bool          ascii)
{
  if (!CheckMaterialInfo(directory, ascii)) return false;
  if (!CheckMaterialCutsCoupleInfo(directory, ascii)) return false;
  return true;
}
  
// Store material information in files under the specified directory.
G4bool  G4ProductionCutsTable::StoreMaterialInfo(const G4String& directory, 
						 G4bool          ascii)
{
  // just dummy (should be implemented later)
  return true;
}

// check stored material is consistent with the current detector setup. 
G4bool  G4ProductionCutsTable::CheckMaterialInfo(const G4String& directory, 
						 G4bool          ascii)
{
  // just dummy (should be implemented later)
  return true;

}
  
// Store materialCutsCouple information in files under the specified directory.
G4bool  G4ProductionCutsTable::StoreMaterialCutsCoupleInfo(const G4String& directory, 
				    G4bool          ascii)
{
  // just dummy (should be implemented later)
  return true;
}

// check stored materialCutsCouple is consistent with the current detector setup. 
G4bool  G4ProductionCutsTable::CheckMaterialCutsCoupleInfo(const G4String& directory,
							   G4bool          ascii )
{
  // just dummy (should be implemented later)
  return true;
}
  
  
// Store cut values information in files under the specified directory.
G4bool  G4ProductionCutsTable::StoreCutsInfo(const G4String& directory, 
					     G4bool          ascii)
{
  // just dummy (should be implemented later)
  return true;
}
  
// Retrieve cut values information in files under the specified directory.
G4bool  G4ProductionCutsTable::RetrieveCutsInfo(const G4String& directory,
						G4bool          ascii)
{
  // just dummy (should be implemented later)
  return true;
}
  










