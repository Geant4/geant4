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
// G4ProductionCutsTable class implementation
//
// Author: M.Asai, 5 October 2002 - First implementation
// Modifications: H.Kurashige, 2004-2008
// --------------------------------------------------------------------

#include "G4ProductionCutsTable.hh"
#include "G4ProductionCuts.hh"
#include "G4MCCIndexConversionTable.hh"
#include "G4ProductionCutsTableMessenger.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4RegionStore.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RToEConvForElectron.hh"
#include "G4RToEConvForGamma.hh"
#include "G4RToEConvForPositron.hh"
#include "G4RToEConvForProton.hh"
#include "G4MaterialTable.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"

#include "G4Timer.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include <iomanip>                
#include <fstream>       

G4ProductionCutsTable* G4ProductionCutsTable::fProductionCutsTable = nullptr;

// --------------------------------------------------------------------
G4ProductionCutsTable* G4ProductionCutsTable::GetProductionCutsTable()
{ 
   static G4ProductionCutsTable theProductionCutsTable;
   if(fProductionCutsTable == nullptr)
   {
     fProductionCutsTable = &theProductionCutsTable;
   }
  return fProductionCutsTable;
}

// --------------------------------------------------------------------
G4ProductionCutsTable::G4ProductionCutsTable()
{
  for(std::size_t i=0; i< NumberOfG4CutIndex; ++i)
  {
    rangeCutTable.push_back(new std::vector<G4double>);
    energyCutTable.push_back(new std::vector<G4double>);
    rangeDoubleVector[i] = nullptr;
    energyDoubleVector[i] = nullptr;
    converters[i] = nullptr;
  }
  fG4RegionStore = G4RegionStore::GetInstance();
  defaultProductionCuts = new G4ProductionCuts();

  // add messenger for UI 
  fMessenger = new G4ProductionCutsTableMessenger(this);
}

// --------------------------------------------------------------------
G4ProductionCutsTable::~G4ProductionCutsTable()
{
  delete defaultProductionCuts;
  defaultProductionCuts = nullptr;

  for(auto itr=coupleTable.cbegin(); itr!=coupleTable.cend(); ++itr)
  {
    delete (*itr); 
  }
  coupleTable.clear();

  for(std::size_t i=0; i< NumberOfG4CutIndex; ++i)
  {
    delete rangeCutTable[i];
    delete energyCutTable[i];
    delete converters[i];
    if(rangeDoubleVector[i] != nullptr) delete [] rangeDoubleVector[i];
    if(energyDoubleVector[i] != nullptr) delete [] energyDoubleVector[i];
    rangeCutTable[i] = nullptr;
    energyCutTable[i] = nullptr;
    converters[i] = nullptr;
    rangeDoubleVector[i] = nullptr;
    energyDoubleVector[i] = nullptr;
  }
  fProductionCutsTable = nullptr;

  delete fMessenger;
  fMessenger = nullptr;
}

// --------------------------------------------------------------------
void G4ProductionCutsTable::UpdateCoupleTable(G4VPhysicalVolume* /*currWorld*/)
{
  if(firstUse)
  {
    if(G4ParticleTable::GetParticleTable()->FindParticle("gamma"))
    {
      converters[0] = new G4RToEConvForGamma(); 
      converters[0]->SetVerboseLevel(GetVerboseLevel());
    }
    if(G4ParticleTable::GetParticleTable()->FindParticle("e-"))
    {
      converters[1] = new G4RToEConvForElectron(); 
      converters[1]->SetVerboseLevel(GetVerboseLevel());
    }
    if(G4ParticleTable::GetParticleTable()->FindParticle("e+"))
    {
      converters[2] = new G4RToEConvForPositron(); 
      converters[2]->SetVerboseLevel(GetVerboseLevel());
    }
    if(G4ParticleTable::GetParticleTable()->FindParticle("proton"))
    {
      converters[3] = new G4RToEConvForProton(); 
      converters[3]->SetVerboseLevel(GetVerboseLevel());
    }
    firstUse = false;
  }

  // Reset "used" flags of all couples
  for(auto CoupleItr=coupleTable.cbegin();
           CoupleItr!=coupleTable.cend(); ++CoupleItr) 
  {
    (*CoupleItr)->SetUseFlag(false); 
  }

  // Update Material-Cut-Couple
  for(auto rItr=fG4RegionStore->cbegin(); rItr!=fG4RegionStore->cend(); ++rItr)
  {
    // Material scan is to be done only for the regions appear in the 
    // current tracking world.
    //    if((*rItr)->GetWorldPhysical()!=currentWorld) continue;

    if( (*rItr)->IsInMassGeometry() || (*rItr)->IsInParallelGeometry() )
    {
      G4ProductionCuts* fProductionCut = (*rItr)->GetProductionCuts();
      auto mItr = (*rItr)->GetMaterialIterator();
      std::size_t nMaterial = (*rItr)->GetNumberOfMaterials();
      (*rItr)->ClearMap();

      for(std::size_t iMate=0; iMate<nMaterial; ++iMate)
      {
        //check if this material cut couple has already been made
        G4bool coupleAlreadyDefined = false;
        G4MaterialCutsCouple* aCouple;
        for(auto cItr=coupleTable.cbegin(); cItr!=coupleTable.cend(); ++cItr)
        {
          if( (*cItr)->GetMaterial()==(*mItr)
           && (*cItr)->GetProductionCuts()==fProductionCut)
          { 
            coupleAlreadyDefined = true;
            aCouple = *cItr;
            break;
          }
        }
      
        // If this combination is new, cleate and register a couple
        if(!coupleAlreadyDefined)
        {
          aCouple = new G4MaterialCutsCouple((*mItr),fProductionCut);
          coupleTable.push_back(aCouple);
          aCouple->SetIndex(G4int(coupleTable.size()-1));
        }

        // Register this couple to the region
        (*rItr)->RegisterMaterialCouplePair((*mItr),aCouple);

        // Set the couple to the proper logical volumes in that region
        aCouple->SetUseFlag();

        auto rootLVItr = (*rItr)->GetRootLogicalVolumeIterator();
        std::size_t nRootLV = (*rItr)->GetNumberOfRootVolumes();
        for(std::size_t iLV=0; iLV<nRootLV; ++iLV)
        {
          // Set the couple to the proper logical volumes in that region
          G4LogicalVolume* aLV = *rootLVItr;
          G4Region* aR = *rItr;

          ScanAndSetCouple(aLV,aCouple,aR);

          // Proceed to the next root logical volume in this region
          ++rootLVItr;
        }

        // Proceed to next material in this region
        ++mItr;
      }
    }
  }

  // Check if sizes of Range/Energy cuts tables are equal to the size of
  // the couple table. If new couples are made during the previous procedure,
  // nCouple becomes larger then nTable

  std::size_t nCouple = coupleTable.size();
  std::size_t nTable = energyCutTable[0]->size();
  G4bool newCoupleAppears = nCouple>nTable;
  if(newCoupleAppears)
  {
    for(std::size_t n=nCouple-nTable; n>0; --n)
    {
      for(std::size_t nn=0; nn< NumberOfG4CutIndex; ++nn)
      {
        rangeCutTable[nn]->push_back(-1.);
        energyCutTable[nn]->push_back(-1.);
      }
    }
  }

  // Update RangeEnergy cuts tables
  std::size_t idx = 0;
  G4Timer timer;
  if (verboseLevel>2)
  {
    timer.Start();
  }
  for(auto cItr=coupleTable.cbegin(); cItr!=coupleTable.cend(); ++cItr)
  {
    G4ProductionCuts* aCut = (*cItr)->GetProductionCuts();
    const G4Material* aMat = (*cItr)->GetMaterial();
    if((*cItr)->IsRecalcNeeded())
    {
      for(std::size_t ptcl=0; ptcl< NumberOfG4CutIndex; ++ptcl)
      {
        G4double rCut = aCut->GetProductionCut((G4int)ptcl);
        (*(rangeCutTable[ptcl]))[idx] = rCut;
        // if(converters[ptcl] && (*cItr)->IsUsed())
        if(converters[ptcl])
        {
          (*(energyCutTable[ptcl]))[idx] = converters[ptcl]->Convert(rCut,aMat);
        }
        else
        {
          (*(energyCutTable[ptcl]))[idx] = -1.; 
        }
      }
    }
    ++idx;  
  }
  if (verboseLevel>2)
  {
    timer.Stop();
    G4cout << "G4ProductionCutsTable::UpdateCoupleTable() - "
           << "Elapsed time for calculation of energy cuts: " << G4endl;
    G4cout << timer << G4endl;
  }

  // resize Range/Energy cuts double vectors if new couple is made
  if(newCoupleAppears)
  {
    for(std::size_t ix=0; ix<NumberOfG4CutIndex; ++ix)
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
  for(std::size_t ix=0; ix<NumberOfG4CutIndex; ++ix)
  {
    for(std::size_t ixx=0; ixx<(*(rangeCutTable[ix])).size(); ++ixx)
    {
      rangeDoubleVector[ix][ixx] = (*(rangeCutTable[ix]))[ixx];
      energyDoubleVector[ix][ixx] = (*(energyCutTable[ix]))[ixx];
    }
  }
}

// --------------------------------------------------------------------
G4double G4ProductionCutsTable::ConvertRangeToEnergy(
                                  const G4ParticleDefinition* particle,
                                  const G4Material*           material, 
                                  G4double                    range )
{
  // This method gives energy corresponding to range value  

  // protection against premature call
  if(firstUse)
  {
#ifdef G4VERBOSE
    if(verboseLevel>0)
    {
      G4ExceptionDescription ed;
      ed << "Invoked prematurely before it is fully initialized.";
      G4Exception("G4ProductionCutsTable::ConvertRangeToEnergy()",
                  "CUTS0100", JustWarning, ed);
    }
#endif
    return -1.0;
  }

  // check material
  if (material == nullptr) return -1.0;

  // check range
  if (range == 0.0) return 0.0;
  if (range <0.0) return -1.0;

  // check particle
  G4int index = G4ProductionCuts::GetIndex(particle);

  if (index<0 || converters[index] == nullptr)
  {
#ifdef G4VERBOSE
    if(verboseLevel>0)
    {
      G4ExceptionDescription ed;
      ed << "Invoked ";
      if(particle != nullptr)
      { ed << "for particle <" << particle->GetParticleName() << ">."; }
      else
      { ed << "without valid particle pointer."; }
      G4Exception("G4ProductionCutsTable::ConvertRangeToEnergy()",
                  "CUTS0101", JustWarning, ed);
    }
#endif
    return -1.0;
  }

  return converters[index]->Convert(range, material);
}

// --------------------------------------------------------------------
void G4ProductionCutsTable::ResetConverters()
{}

// --------------------------------------------------------------------
void G4ProductionCutsTable::SetEnergyRange(G4double lowedge, G4double highedge)
{
  G4VRangeToEnergyConverter::SetEnergyRange(lowedge,highedge);
}

// --------------------------------------------------------------------
G4double  G4ProductionCutsTable::GetLowEdgeEnergy() const
{
  return G4VRangeToEnergyConverter::GetLowEdgeEnergy();
}

// --------------------------------------------------------------------
G4double G4ProductionCutsTable::GetHighEdgeEnergy() const
{
  return G4VRangeToEnergyConverter::GetHighEdgeEnergy();
}
 
// --------------------------------------------------------------------
void G4ProductionCutsTable::ScanAndSetCouple(G4LogicalVolume* aLV,
                                             G4MaterialCutsCouple* aCouple,
                                             G4Region* aRegion)
{
  // Check whether or not this logical volume belongs to the same region
  if((aRegion!=nullptr) && aLV->GetRegion()!=aRegion) return;

  // Check if this particular volume has a material matched to the couple
  if(aLV->GetMaterial()==aCouple->GetMaterial())
  {
    aLV->SetMaterialCutsCouple(aCouple);
  }

  std::size_t noDaughters = aLV->GetNoDaughters();
  if(noDaughters==0) return;

  // Loop over daughters with same region
  for(std::size_t i=0; i<noDaughters; ++i)
  {
    G4LogicalVolume* daughterLVol = aLV->GetDaughter(i)->GetLogicalVolume();
    ScanAndSetCouple(daughterLVol,aCouple,aRegion);
  }
}

// --------------------------------------------------------------------
void G4ProductionCutsTable::DumpCouples() const
{
  G4cout << G4endl;
  G4cout << "========= Table of registered couples ============================"
         << G4endl;
  for(auto cItr=coupleTable.cbegin(); cItr!=coupleTable.cend(); ++cItr)
  {
    G4MaterialCutsCouple* aCouple = (*cItr);
    G4ProductionCuts* aCut = aCouple->GetProductionCuts();
    G4cout << G4endl;
    G4cout << "Index : " << aCouple->GetIndex() 
           << "     used in the geometry : ";
    if(aCouple->IsUsed()) G4cout << "Yes";
    else                  G4cout << "No ";
////    G4cout << "     recalculation needed : ";
////    if(aCouple->IsRecalcNeeded()) G4cout << "Yes";
////    else                          G4cout << "No ";
    G4cout << G4endl;
    G4cout << " Material : " << aCouple->GetMaterial()->GetName() << G4endl;
    G4cout << " Range cuts        : " 
           << " gamma  " << G4BestUnit(aCut->GetProductionCut("gamma"),"Length")
           << "    e-  " << G4BestUnit(aCut->GetProductionCut("e-"),"Length")
           << "    e+  " << G4BestUnit(aCut->GetProductionCut("e+"),"Length")
           << " proton " << G4BestUnit(aCut->GetProductionCut("proton"),"Length")
           << G4endl;
    G4cout << " Energy thresholds : " ;
////    if(aCouple->IsRecalcNeeded()) {
////      G4cout << " is not ready to print";
////    } else {
    G4cout << " gamma  " << G4BestUnit((*(energyCutTable[0]))[aCouple->GetIndex()],"Energy")
           << "    e-  " << G4BestUnit((*(energyCutTable[1]))[aCouple->GetIndex()],"Energy")
           << "    e+  " << G4BestUnit((*(energyCutTable[2]))[aCouple->GetIndex()],"Energy") 
           << " proton " << G4BestUnit((*(energyCutTable[3]))[aCouple->GetIndex()],"Energy");
////    }
    G4cout << G4endl;

    if(aCouple->IsUsed())
    {
      G4cout << " Region(s) which use this couple : " << G4endl;
      for(auto rItr=fG4RegionStore->cbegin();
               rItr!=fG4RegionStore->cend(); ++rItr)
      {
        if (IsCoupleUsedInTheRegion(aCouple, *rItr) )
        {
          G4cout << "    " << (*rItr)->GetName() << G4endl;
        }
      }
    }
  }
  G4cout << G4endl;
  G4cout << "==================================================================" << G4endl;
  G4cout << G4endl;
}
           
// --------------------------------------------------------------------
G4bool  G4ProductionCutsTable::StoreCutsTable(const G4String& dir, 
                                              G4bool ascii)
{
  // Store cuts and material information in files under the specified directory

  if (!StoreMaterialInfo(dir, ascii)) return false;
  if (!StoreMaterialCutsCoupleInfo(dir, ascii)) return false;
  if (!StoreCutsInfo(dir, ascii)) return false;
  
#ifdef G4VERBOSE  
  if (verboseLevel >2)
  {
    G4cout << "G4ProductionCutsTable::StoreCutsTable()" << G4endl;
    G4cout << " Material/Cuts information have been successfully stored ";
    if (ascii)
    {
      G4cout << " in Ascii mode ";
    }
    else
    {
      G4cout << " in Binary mode ";
    }
    G4cout << " under " << dir << G4endl;  
  }  
#endif
  return true;
}
  
// --------------------------------------------------------------------
G4bool  G4ProductionCutsTable::RetrieveCutsTable(const G4String& dir,
                                                 G4bool ascii)
{
  if (!CheckForRetrieveCutsTable(dir, ascii)) return false;
  if (!RetrieveCutsInfo(dir, ascii)) return false;
#ifdef G4VERBOSE  
  if (verboseLevel >2)
  {
    G4cout << "G4ProductionCutsTable::RetrieveCutsTable()" << G4endl;
    G4cout << " Material/Cuts information have been successfully retrieved ";
    if (ascii)
    {
      G4cout << " in Ascii mode ";
    }
    else
    {
      G4cout << " in Binary mode ";
    }
    G4cout << " under " << dir << G4endl;  
  }  
#endif
  return true;
}

// --------------------------------------------------------------------
G4bool
G4ProductionCutsTable::CheckForRetrieveCutsTable(const G4String& directory, 
                                                 G4bool ascii)
{
  // check stored material and cut values are consistent
  // with the current detector setup

  G4cerr << "G4ProductionCutsTable::CheckForRetrieveCutsTable()"<< G4endl;
  //  isNeedForRestoreCoupleInfo = false;
  if (!CheckMaterialInfo(directory, ascii)) return false;
  if (verboseLevel >2)
  {
      G4cerr << "G4ProductionCutsTable::CheckMaterialInfo passed !!"<< G4endl;
  }
  if (!CheckMaterialCutsCoupleInfo(directory, ascii)) return false;
  if (verboseLevel >2)
  {
    G4cerr << "G4ProductionCutsTable::CheckMaterialCutsCoupleInfo passed !!"
           << G4endl;
  }
  return true;
}
  
// --------------------------------------------------------------------
G4bool G4ProductionCutsTable::StoreMaterialInfo(const G4String& directory, 
                                                G4bool ascii)
{
  // Store material information in files under the specified directory

  const G4String fileName = directory + "/" + "material.dat";
  const G4String key = "MATERIAL-V3.0";
  std::ofstream fOut;  

  // open output file
  if (!ascii )  fOut.open(fileName,std::ios::out|std::ios::binary);
  else          fOut.open(fileName,std::ios::out);

  // check if the file has been opened successfully 
  if (!fOut)
  {
#ifdef G4VERBOSE
    if (verboseLevel>0)
    {
      G4cerr << "G4ProductionCutsTable::StoreMaterialInfo() - ";
      G4cerr << "Cannot open file: " << fileName << G4endl;
    }
#endif
    G4Exception( "G4ProductionCutsTable::StoreMaterialInfo()",
                 "ProcCuts102", JustWarning, "Cannot open file!");
    return false;
  }
  
  const G4MaterialTable* matTable = G4Material::GetMaterialTable(); 
  // number of materials in the table
  G4int numberOfMaterial = (G4int)matTable->size();

  if (ascii)
  {
    /////////////// ASCII mode  /////////////////
    // key word
    fOut  << key << G4endl;
    
    // number of materials in the table
    fOut  << numberOfMaterial << G4endl;
    
    fOut.setf(std::ios::scientific);
  
    // material name and density
    for (std::size_t idx=0; static_cast<G4int>(idx)<numberOfMaterial; ++idx)
    {
      fOut << std::setw(FixedStringLengthForStore)
           << ((*matTable)[idx])->GetName();
      fOut << std::setw(FixedStringLengthForStore)
           << ((*matTable)[idx])->GetDensity()/(g/cm3) << G4endl;
    }
    
    fOut.unsetf(std::ios::scientific);

  }
  else
  {
    /////////////// Binary mode  /////////////////
    char temp[FixedStringLengthForStore];
    std::size_t i;

    // key word
    for (i=0; i<FixedStringLengthForStore; ++i)
    {
      temp[i] = '\0'; 
    }
    for (i=0; i<key.length() && i<FixedStringLengthForStore-1; ++i)
    {
      temp[i]=key[(G4int)i];
    }
    fOut.write(temp, FixedStringLengthForStore);

    // number of materials in the table
    fOut.write( (char*)(&numberOfMaterial), sizeof(G4int));
    
    // material name and density
    for (std::size_t imat=0; static_cast<G4int>(imat)<numberOfMaterial; ++imat)
    {
      G4String name =  ((*matTable)[imat])->GetName();
      G4double density = ((*matTable)[imat])->GetDensity();
      for (i=0; i<FixedStringLengthForStore; ++i)
        temp[i] = '\0'; 
      for (i=0; i<name.length() && i<FixedStringLengthForStore-1; ++i)
        temp[i]=name[(G4int)i];
      fOut.write(temp, FixedStringLengthForStore);
      fOut.write( (char*)(&density), sizeof(G4double));
    }    
  }    

  fOut.close();
  return true;
}

// --------------------------------------------------------------------
G4bool G4ProductionCutsTable::CheckMaterialInfo(const G4String& directory, 
                                                G4bool ascii)
{
  // Check stored material is consistent with the current detector setup

  const G4String fileName = directory + "/" + "material.dat";
  const G4String key = "MATERIAL-V3.0";
  std::ifstream fIn;  

  // open input file
  if (!ascii ) fIn.open(fileName,std::ios::in|std::ios::binary);
  else         fIn.open(fileName,std::ios::in);

  // check if the file has been opened successfully 
  if (!fIn)
  {
#ifdef G4VERBOSE
    if (verboseLevel >0)
    {
      G4cerr << "G4ProductionCutsTable::CheckMaterialInfo() - ";
      G4cerr << "Cannot open file: " << fileName << G4endl;
    }
#endif
    G4Exception( "G4ProductionCutsTable::CheckMaterialInfo()",
                 "ProcCuts102", JustWarning, "Cannot open file!");
    return false;
  }
  
  char temp[FixedStringLengthForStore];

  // key word
  G4String keyword;    
  if (ascii)
  {
    fIn >> keyword;
  }
  else
  {
    fIn.read(temp, FixedStringLengthForStore);
    keyword = (const char*)(temp);
  }
  if (key!=keyword)
  {
#ifdef G4VERBOSE
    if (verboseLevel >0)
    {
      G4cerr << "G4ProductionCutsTable::CheckMaterialInfo() - ";
      G4cerr << "Key word in " << fileName << "= " << keyword ;
      G4cerr <<"( should be   "<< key << ")" <<G4endl;
    }
#endif
    G4Exception( "G4ProductionCutsTable::CheckMaterialInfo()",
                 "ProcCuts103", JustWarning, "Bad Data Format");
    return false;
  }

  // number of materials in the table
  G4int nmat;
  if (ascii)
  {
    fIn >> nmat;
  }
  else
  {
    fIn.read( (char*)(&nmat), sizeof(G4int));
  }
  if ((nmat<=0) || (nmat >100000))
  {
    G4Exception( "G4ProductionCutsTable::CheckMaterialInfo()",
                 "ProcCuts108", JustWarning, 
                 "Number of materials is less than zero or too big");
    return false;
  }

  // list of material
  for (G4int idx=0; idx<nmat ; ++idx)
  {
    // check eof
    if(fIn.eof())
    {
#ifdef G4VERBOSE
      if (verboseLevel >0)
      {
        G4cout << "G4ProductionCutsTable::CheckMaterialInfo() - ";
        G4cout << "Encountered End of File " ;
        G4cout << " at " << idx+1 << "th  material "<< G4endl;
      }
#endif
      fIn.close();
      return false;
    }

    // check material name and density
    char name[FixedStringLengthForStore];
    G4double density;
    if (ascii)
    {
      fIn >> name >> density;
      density *= (g/cm3);
      
    }
    else
    {
      fIn.read(name, FixedStringLengthForStore);
      fIn.read((char*)(&density), sizeof(G4double));
    }
    if (fIn.fail())
    {
#ifdef G4VERBOSE
      if (verboseLevel >0)
      {
        G4cerr << "G4ProductionCutsTable::CheckMaterialInfo() - ";
        G4cerr << "Bad data format ";
        G4cerr << " at " << idx+1 << "th  material "<< G4endl;        
      }
#endif
      G4Exception( "G4ProductionCutsTable::CheckMaterialInfo()",
                   "ProcCuts103", JustWarning, "Bad Data Format");
      fIn.close();
      return false;
    }

    G4Material* aMaterial = G4Material::GetMaterial(name);
    if (aMaterial == nullptr ) continue;

    G4double ratio = std::fabs(density/aMaterial->GetDensity() );
    if ((0.999>ratio) || (ratio>1.001) )
    {
#ifdef G4VERBOSE
      if (verboseLevel >0)
      {
        G4cerr << "G4ProductionCutsTable::CheckMaterialInfo() - ";
        G4cerr << "Inconsistent material density" << G4endl;;
        G4cerr << " at " << idx+1 << "th  material "<< G4endl;        
        G4cerr << "Name:   " << name << G4endl;
        G4cerr << "Density:" << std::setiosflags(std::ios::scientific)
               << density / (g/cm3) ;
        G4cerr << "(should be " << aMaterial->GetDensity()/(g/cm3)<< ")"
               << " [g/cm3]"<< G4endl;      
        G4cerr << std::resetiosflags(std::ios::scientific);
      }
#endif
      G4Exception( "G4ProductionCutsTable::CheckMaterialInfo()",
                   "ProcCuts104", JustWarning, "Inconsistent material density");
      fIn.close();
      return false;
    }
  }

  fIn.close();
  return true;
}
  
// --------------------------------------------------------------------
G4bool
G4ProductionCutsTable::StoreMaterialCutsCoupleInfo(const G4String& directory, 
                                                   G4bool ascii)
{  
  // Store materialCutsCouple information in files under the specified directory

  const G4String fileName = directory + "/" + "couple.dat";
  const G4String key = "COUPLE-V3.0";
  std::ofstream fOut;  
  char temp[FixedStringLengthForStore];

  // open output file
  if (!ascii ) fOut.open(fileName,std::ios::out|std::ios::binary);
  else         fOut.open(fileName,std::ios::out);
  
  
  // check if the file has been opened successfully 
  if (!fOut)
  {
#ifdef G4VERBOSE
    if (verboseLevel >0)
    {
      G4cerr << "G4ProductionCutsTable::StoreMaterialCutsCoupleInfo() - ";
      G4cerr << "Cannot open file: " << fileName << G4endl;
    }
#endif
    G4Exception( "G4ProductionCutsTable::StoreMaterialCutsCoupleInfo()",
                 "ProcCuts102",
                 JustWarning, "Cannot open file!");
     return false;
  }
  G4int numberOfCouples = (G4int)coupleTable.size();
  if (ascii)
  {
    /////////////// ASCII mode  /////////////////
    // key word
    fOut << std::setw(FixedStringLengthForStore) <<  key << G4endl;
    
    // number of couples in the table
    fOut  << numberOfCouples << G4endl;
  }
  else
  {
    /////////////// Binary mode  /////////////////
    // key word
    std::size_t i;
    for (i=0; i<FixedStringLengthForStore; ++i)
      temp[i] = '\0'; 
    for (i=0; i<key.length() && i<FixedStringLengthForStore-1; ++i)
      temp[i]=key[(G4int)i];
    fOut.write(temp, FixedStringLengthForStore);
    
    // number of couples in the table   
    fOut.write( (char*)(&numberOfCouples), sizeof(G4int));
  }

  // Loop over all couples
  for (auto cItr=coupleTable.cbegin(); cItr!=coupleTable.cend(); ++cItr)
  {
    G4MaterialCutsCouple* aCouple = (*cItr);
    G4int index = aCouple->GetIndex(); 
    // cut value
    G4ProductionCuts* aCut = aCouple->GetProductionCuts();
    G4double cutValues[NumberOfG4CutIndex];
    for (std::size_t idx=0; idx <NumberOfG4CutIndex; ++idx)
    {
      cutValues[idx] = aCut->GetProductionCut((G4int)idx);
    }
    // material/region info
    G4String materialName = aCouple->GetMaterial()->GetName();
    G4String regionName = "NONE";
    if (aCouple->IsUsed())
    {
      for(auto rItr=fG4RegionStore->cbegin();
               rItr!=fG4RegionStore->cend(); ++rItr)
      {
        if (IsCoupleUsedInTheRegion(aCouple, *rItr))
        {
          regionName = (*rItr)->GetName();
          break;
        }
      }
    }

    if (ascii)
    {
      /////////////// ASCII mode  /////////////////
      // index number
      fOut  << index << G4endl; 
  
      // material name 
      fOut << std::setw(FixedStringLengthForStore) << materialName<< G4endl;
  
      // region name 
      fOut << std::setw(FixedStringLengthForStore) << regionName<< G4endl;

      fOut.setf(std::ios::scientific);
      // cut values
      for (std::size_t idx=0; idx< NumberOfG4CutIndex; ++idx)
      {
        fOut << std::setw(FixedStringLengthForStore) << cutValues[idx]/(mm)
             << G4endl;
      }
      fOut.unsetf(std::ios::scientific);

    }
    else
    {
      /////////////// Binary mode  /////////////////
      // index
      fOut.write( (char*)(&index), sizeof(G4int));
    
      // material name
      std::size_t i;
      for (i=0; i<FixedStringLengthForStore; ++i)
        temp[i] = '\0'; 
      for (i=0; i<materialName.length() && i<FixedStringLengthForStore-1; ++i)
        temp[i]=materialName[(G4int)i];
      fOut.write(temp, FixedStringLengthForStore);

      // region name
      for (i=0; i<FixedStringLengthForStore; ++i)
        temp[i] = '\0'; 
      for (i=0; i<regionName.length() && i<FixedStringLengthForStore-1; ++i)
        temp[i]=regionName[(G4int)i];
      fOut.write(temp, FixedStringLengthForStore);

      // cut values
      for (std::size_t idx=0; idx< NumberOfG4CutIndex; ++idx)
      {
         fOut.write( (char*)(&(cutValues[idx])), sizeof(G4double));
      }    
    }
  }
  fOut.close();
  return true;
}

// --------------------------------------------------------------------
G4bool
G4ProductionCutsTable::CheckMaterialCutsCoupleInfo(const G4String& directory,
                                                   G4bool ascii)
{
  // Check stored materialCutsCouple is consistent
  // with the current detector setup. 

  const G4String fileName = directory + "/" + "couple.dat";
  const G4String key = "COUPLE-V3.0";
  std::ifstream fIn;  

  // open input file
  if (!ascii ) fIn.open(fileName,std::ios::in|std::ios::binary);
  else         fIn.open(fileName,std::ios::in);

  // check if the file has been opened successfully 
  if (!fIn)
  {
#ifdef G4VERBOSE
    if (verboseLevel >0)
    {
      G4cerr << "G4ProductionCutTable::CheckMaterialCutsCoupleInfo() - ";
      G4cerr << "Cannot open file!" << fileName << G4endl;
    }
#endif
    G4Exception( "G4ProductionCutsTable::CheckMaterialCutsCoupleInfo()",
                 "ProcCuts102", JustWarning, "Cannot open file!");
    return false;
  }
  
  char temp[FixedStringLengthForStore];

   // key word
  G4String keyword;    
  if (ascii)
  {
    fIn >> keyword;
  }
  else
  {
    fIn.read(temp, FixedStringLengthForStore);
    keyword = (const char*)(temp);
  }
  if (key!=keyword)
  {
#ifdef G4VERBOSE
    if (verboseLevel >0)
    {
      G4cerr << "G4ProductionCutTable::CheckMaterialCutsCoupleInfo() - ";
      G4cerr << "Key word in " << fileName << "= " << keyword ;
      G4cerr <<"( should be   "<< key << ")" << G4endl;
    }
#endif
    G4Exception( "G4ProductionCutsTable::CheckMaterialCutsCoupleInfo()",
                 "ProcCuts103", JustWarning, "Bad Data Format");
    fIn.close();
    return false;
  }

  // numberOfCouples
  G4int numberOfCouples;    
  if (ascii)
  {
    fIn >> numberOfCouples;
  }
  else
  {
    fIn.read( (char*)(&numberOfCouples), sizeof(G4int));
  }

  // Reset MCCIndexConversionTable
  mccConversionTable.Reset(numberOfCouples);  

  // Read in couple information
  for (G4int idx=0; idx<numberOfCouples; ++idx)
  {
    // read in index
    G4int index; 
    if (ascii)
    {
      fIn >> index;
    }
    else
    {
      fIn.read( (char*)(&index), sizeof(G4int));
    }
    // read in index material name
    char mat_name[FixedStringLengthForStore];
    if (ascii)
    {
      fIn >> mat_name;
    }
    else
    {
      fIn.read(mat_name, FixedStringLengthForStore);
    }
    // read in index and region name
    char region_name[FixedStringLengthForStore];
    if (ascii)
    {
      fIn >> region_name;
    }
    else
    {
      fIn.read(region_name, FixedStringLengthForStore);
    }
    // cut value
    G4double cutValues[NumberOfG4CutIndex];
    for (std::size_t i=0; i< NumberOfG4CutIndex; ++i)
    {
      if (ascii)
      {
        fIn >>  cutValues[i];
        cutValues[i] *= (mm);
      }
      else
      {
        fIn.read( (char*)(&(cutValues[i])), sizeof(G4double));
      }
    }
 
    // Loop over all couples
    G4bool fOK = false;
    G4MaterialCutsCouple* aCouple = nullptr;
    for (auto cItr=coupleTable.cbegin(); cItr!=coupleTable.cend(); ++cItr)
    {
      aCouple = (*cItr);
      // check material name
      if ( mat_name !=  aCouple->GetMaterial()->GetName() ) continue;
      // check cut values
      G4ProductionCuts* aCut = aCouple->GetProductionCuts();
      G4bool fRatio = true;
      for (std::size_t j=0; j< NumberOfG4CutIndex; ++j)
      {
        // check ratio only if values are not the same
        if (cutValues[j] != aCut->GetProductionCut((G4int)j))
        {  
          G4double ratio =  cutValues[j]/aCut->GetProductionCut((G4int)j);
          fRatio = fRatio && (0.999<ratio) && (ratio<1.001) ;
        }
      } 
      if (!fRatio) continue; 
      // MCC matched 
      fOK = true;
      mccConversionTable.SetNewIndex(index, aCouple->GetIndex()); 
      break;
    }

    if (fOK)
    {
#ifdef G4VERBOSE
      // debug information 
      if (verboseLevel >1)
      {
        G4String regionname(region_name);
        G4Region* fRegion = nullptr;
        if ( regionname != "NONE" )
        {
          fRegion = fG4RegionStore->GetRegion(region_name);
          if (fRegion == nullptr)
          {
            G4cout << "G4ProductionCutTable::CheckMaterialCutsCoupleInfo() - ";
            G4cout << "Region " << regionname << " is not found ";            
            G4cout << index << ": in " << fileName  << G4endl;
          } 
        }
        if  (((regionname == "NONE") && (aCouple->IsUsed()))
         || ((fRegion!=nullptr) && !IsCoupleUsedInTheRegion(aCouple, fRegion)))
        {
          G4cout << "G4ProductionCutTable::CheckMaterialCutsCoupleInfo()"
                 << G4endl;
          G4cout << "A Couple is used different region in the current setup ";
          G4cout << index << ": in " << fileName  << G4endl;
          G4cout << " material: " << mat_name ;
          G4cout << " region: " << region_name << G4endl;
          for (std::size_t ii=0; ii< NumberOfG4CutIndex; ++ii)
          {
            G4cout << "cut[" << ii << "]=" << cutValues[ii]/mm;
            G4cout << " mm   :  ";
          } 
          G4cout << G4endl;
        }
        else if ( index !=  aCouple->GetIndex() )
        {
          G4cout << "G4ProductionCutTable::CheckMaterialCutsCoupleInfo() - ";
          G4cout << "Index of couples was modified "<< G4endl;
          G4cout << aCouple->GetIndex() << ":"
                 << aCouple->GetMaterial()->GetName();
          G4cout <<" is defined as " ;
          G4cout << index << ":"  << mat_name << " in " << fileName << G4endl;
        }
        else
        {
          G4cout << "G4ProductionCutTable::CheckMaterialCutsCoupleInfo() - ";
          G4cout << index << ":"  << mat_name << " in " << fileName ;
          G4cout << " is consistent with current setup" << G4endl;
        }
      }
#endif
    }
    else
    {
#ifdef G4VERBOSE
      if (verboseLevel >0)
      {
        G4cout << "G4ProductionCutTable::CheckMaterialCutsCoupleInfo()"
               << G4endl;
        G4cout << "Couples are not defined in the current detector setup ";
        G4cout << index << ": in " << fileName  << G4endl;
        G4cout << " material: " << mat_name ;
        G4cout << " region: " << region_name << G4endl;
        for (std::size_t ii=0; ii< NumberOfG4CutIndex; ++ii)
        {
          G4cout << "cut[" << ii << "]=" << cutValues[ii]/mm;
          G4cout << " mm   :  ";
        }
        G4cout << G4endl;
      }
#endif
    }
  }
  fIn.close();
  return true;
}

// --------------------------------------------------------------------
G4bool G4ProductionCutsTable::StoreCutsInfo(const G4String& directory, 
                                            G4bool ascii)
{
  // Store cut values information in files under the specified directory

  const G4String fileName = directory + "/" + "cut.dat";
  const G4String key = "CUT-V3.0";
  std::ofstream fOut;  
  char temp[FixedStringLengthForStore];
  
  // open output file
  if (!ascii ) fOut.open(fileName,std::ios::out|std::ios::binary);
  else         fOut.open(fileName,std::ios::out);

  // check if the file has been opened successfully 
  if (!fOut)
  {
    if(verboseLevel>0)
    {         
      G4cerr << "G4ProductionCutsTable::StoreCutsInfo() - ";
      G4cerr << "Cannot open file: " << fileName << G4endl;
    }
    G4Exception( "G4ProductionCutsTable::StoreCutsInfo()",
                 "ProcCuts102", JustWarning, "Cannot open file!");
    return false;
  }

  G4int numberOfCouples = (G4int)coupleTable.size();
  if (ascii)
  {
    /////////////// ASCII mode  /////////////////
    // key word
    fOut << key << G4endl;
    
    // number of couples in the table
    fOut << numberOfCouples << G4endl;
  }
  else
  {
    /////////////// Binary mode  /////////////////
    // key word
    std::size_t i;
    for (i=0; i<FixedStringLengthForStore; ++i)
      temp[i] = '\0'; 
    for (i=0; i<key.length() && i<FixedStringLengthForStore-1; ++i)
      temp[i]=key[(G4int)i];
    fOut.write(temp, FixedStringLengthForStore);
    
    // number of couples in the table   
    fOut.write( (char*)(&numberOfCouples), sizeof(G4int));
  }

  for (std::size_t idx=0; idx <NumberOfG4CutIndex; ++idx)
  {
    const std::vector<G4double>* fRange  = GetRangeCutsVector(idx);
    const std::vector<G4double>* fEnergy = GetEnergyCutsVector(idx);
    std::size_t i=0;
    // Loop over all couples
    for (auto cItr=coupleTable.cbegin();cItr!=coupleTable.cend(); ++cItr, ++i)
    {      
      if (ascii)
      {
        /////////////// ASCII mode  /////////////////
        fOut.setf(std::ios::scientific);
        fOut << std::setw(20) << (*fRange)[i]/mm  ;
        fOut << std::setw(20) << (*fEnergy)[i]/keV << G4endl;
        fOut.unsetf(std::ios::scientific);
      }
      else
      {
        /////////////// Binary mode  /////////////////
        G4double cut =  (*fRange)[i];
        fOut.write((char*)(&cut), sizeof(G4double));
        cut =  (*fEnergy)[i];
        fOut.write((char*)(&cut), sizeof(G4double));
      }
    }
  }
  fOut.close();
  return true;
}
  
// --------------------------------------------------------------------
G4bool G4ProductionCutsTable::RetrieveCutsInfo(const G4String& directory,
                                               G4bool ascii)
{
  // Retrieve cut values information in files under the specified directory

  const G4String fileName = directory + "/" + "cut.dat";
  const G4String key = "CUT-V3.0";
  std::ifstream fIn;  

  // open input file
  if (!ascii ) fIn.open(fileName,std::ios::in|std::ios::binary);
  else         fIn.open(fileName,std::ios::in);

  // check if the file has been opened successfully 
  if (!fIn)
  {
    if (verboseLevel >0)
    {
      G4cerr << "G4ProductionCutTable::RetrieveCutsInfo() - ";
      G4cerr << "Cannot open file: " << fileName << G4endl;
    }
    G4Exception( "G4ProductionCutsTable::RetrieveCutsInfo()",
                 "ProcCuts102", JustWarning, "Cannot open file!");
    return false;
  }
  
  char temp[FixedStringLengthForStore];

   // key word
  G4String keyword;    
  if (ascii)
  {
    fIn >> keyword;
  }
  else
  {
    fIn.read(temp, FixedStringLengthForStore);
    keyword = (const char*)(temp);
  }
  if (key!=keyword)
  {
    if (verboseLevel >0)
    {
      G4cerr << "G4ProductionCutTable::RetrieveCutsInfo() - ";
      G4cerr << "Key word in " << fileName << "= " << keyword ;
      G4cerr <<"( should be   "<< key << ")" << G4endl;
    }
    G4Exception( "G4ProductionCutsTable::RetrieveCutsInfo()",
                 "ProcCuts103", JustWarning, "Bad Data Format");
    return false;
  }

  // numberOfCouples
  G4int numberOfCouples;    
  if (ascii)
  {
    fIn >> numberOfCouples;
    if (fIn.fail())
    {
      G4Exception( "G4ProductionCutsTable::RetrieveCutsInfo()",
                   "ProcCuts103", JustWarning, "Bad Data Format");
      return false;
    }
  }
  else
  {
    fIn.read( (char*)(&numberOfCouples), sizeof(G4int));
  }

  if (numberOfCouples > static_cast<G4int>(mccConversionTable.size()) )
  {
    G4Exception( "G4ProductionCutsTable::RetrieveCutsInfo()",
                 "ProcCuts109", JustWarning, 
                 "Number of Couples in the file exceeds defined couples");
  }
  numberOfCouples = (G4int)mccConversionTable.size();
  
  for (std::size_t idx=0; static_cast<G4int>(idx) <NumberOfG4CutIndex; ++idx)
  {
    std::vector<G4double>* fRange  = rangeCutTable[idx];
    std::vector<G4double>* fEnergy = energyCutTable[idx];
    fRange->clear();
    fEnergy->clear();

    // Loop over all couples
    for (std::size_t i=0; static_cast<G4int>(i)< numberOfCouples; ++i)
    {      
      G4double rcut, ecut;
      if (ascii)
      {
        fIn >> rcut >> ecut;
        if (fIn.fail())
        {
          G4Exception( "G4ProductionCutsTable::RetrieveCutsInfo()",
                       "ProcCuts103", JustWarning, "Bad Data Format");
          return false;
        }
        rcut *= mm;
        ecut *= keV;
      }
      else
      {
        fIn.read((char*)(&rcut), sizeof(G4double));
        fIn.read((char*)(&ecut), sizeof(G4double));
      }
      if (!mccConversionTable.IsUsed(i)) continue;
      std::size_t new_index = mccConversionTable.GetIndex(i);
      (*fRange)[new_index]  = rcut;
      (*fEnergy)[new_index] = ecut;
    }
  }
  return true;
}

// --------------------------------------------------------------------
void G4ProductionCutsTable::SetVerboseLevel(G4int value)
{
  // Set same verbosity to all registered RangeToEnergyConverters  

  verboseLevel = value;
  for (G4int ip=0; ip< NumberOfG4CutIndex; ++ip)
  {
    if (converters[ip] != nullptr )
    {
      converters[ip]->SetVerboseLevel(value);
    }
  }
}

// --------------------------------------------------------------------
G4double G4ProductionCutsTable::GetMaxEnergyCut()
{
  return G4VRangeToEnergyConverter::GetMaxEnergyCut();
}

// --------------------------------------------------------------------
void G4ProductionCutsTable::SetMaxEnergyCut(G4double value)
{
  G4VRangeToEnergyConverter::SetMaxEnergyCut(value);
}
