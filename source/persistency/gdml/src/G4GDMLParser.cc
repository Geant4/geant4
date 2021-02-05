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
// G4GDMLParser implementation
//
// Author: Zoltan Torzsok, November 2007
// --------------------------------------------------------------------

#include "G4GDMLParser.hh"

#include "G4UnitsTable.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RegionStore.hh"
#include "G4UserLimits.hh"
#include "G4ProductionCuts.hh"
#include "G4ReflectionFactory.hh"
#include "G4Track.hh"

// --------------------------------------------------------------------
G4GDMLParser::G4GDMLParser()
  : strip(true)
{
  reader    = new G4GDMLReadStructure;
  writer    = new G4GDMLWriteStructure;
  messenger = new G4GDMLMessenger(this);

  xercesc::XMLPlatformUtils::Initialize();
}

// --------------------------------------------------------------------
G4GDMLParser::G4GDMLParser(G4GDMLReadStructure* extr)
  : urcode(true), strip(true)
{
  reader    = extr;
  writer    = new G4GDMLWriteStructure;
  messenger = new G4GDMLMessenger(this);

  xercesc::XMLPlatformUtils::Initialize();
}

// --------------------------------------------------------------------
G4GDMLParser::G4GDMLParser(G4GDMLReadStructure* extr,
                           G4GDMLWriteStructure* extw)
  : urcode(true), uwcode(true), strip(true)
{
  reader    = extr;
  writer    = extw;
  messenger = new G4GDMLMessenger(this);

  xercesc::XMLPlatformUtils::Initialize();
}

// --------------------------------------------------------------------
G4GDMLParser::~G4GDMLParser()
{
  xercesc::XMLPlatformUtils::Terminate();
  if(!urcode)
  {
    delete reader;
  }
  if(!uwcode)
  {
    delete writer;
  }
  delete ullist;
  delete rlist;

  delete messenger;
}

// --------------------------------------------------------------------
void G4GDMLParser::ImportRegions()
{
  G4ReflectionFactory* reflFactory     = G4ReflectionFactory::Instance();
  const G4GDMLAuxListType* auxInfoList = GetAuxList();
  for(auto iaux = auxInfoList->cbegin(); iaux != auxInfoList->cend(); ++iaux)
  {
    if(iaux->type != "Region")
      continue;

    G4String name = iaux->value;
    if(strip)
    {
      reader->StripName(name);
    }
    if(name.contains("DefaultRegionForTheWorld"))
      continue;

    if(!iaux->auxList)
    {
      G4Exception("G4GDMLParser::ImportRegions()", "ReadError", FatalException,
                  "Invalid definition of geometrical region!");
    }
    else  // Create region and loop over all region attributes
    {
      G4Region* aRegion       = new G4Region(name);
      G4ProductionCuts* pcuts = new G4ProductionCuts();
      aRegion->SetProductionCuts(pcuts);
      for(auto raux = iaux->auxList->cbegin();
               raux != iaux->auxList->cend(); ++raux)
      {
        const G4String& tag = raux->type;
        if(tag == "volume")
        {
          G4String volname = raux->value;
          if(strip)
          {
            reader->StripName(volname);
          }
          G4LogicalVolume* lvol =
            G4LogicalVolumeStore::GetInstance()->GetVolume(volname);
          aRegion->AddRootLogicalVolume(lvol);
          if(reflFactory->IsConstituent(lvol))
            aRegion->AddRootLogicalVolume(reflFactory->GetReflectedLV(lvol));
        }
        else if(tag == "pcut")
        {
          const G4String& cvalue = raux->value;
          const G4String& cunit  = raux->unit;
          if(G4UnitDefinition::GetCategory(cunit) != "Length")
          {
            G4Exception("G4GDMLParser::ImportRegions()", "InvalidRead",
                        FatalException, "Invalid unit for length!");
          }
          G4double cut =
            eval.Evaluate(cvalue) * G4UnitDefinition::GetValueOf(cunit);
          pcuts->SetProductionCut(cut, "proton");
        }
        else if(tag == "ecut")
        {
          const G4String& cvalue = raux->value;
          const G4String& cunit  = raux->unit;
          if(G4UnitDefinition::GetCategory(cunit) != "Length")
          {
            G4Exception("G4GDMLParser::ImportRegions()", "InvalidRead",
                        FatalException, "Invalid unit for length!");
          }
          G4double cut =
            eval.Evaluate(cvalue) * G4UnitDefinition::GetValueOf(cunit);
          pcuts->SetProductionCut(cut, "e-");
        }
        else if(tag == "poscut")
        {
          const G4String& cvalue = raux->value;
          const G4String& cunit  = raux->unit;
          if(G4UnitDefinition::GetCategory(cunit) != "Length")
          {
            G4Exception("G4GDMLParser::ImportRegions()", "InvalidRead",
                        FatalException, "Invalid unit for length!");
          }
          G4double cut =
            eval.Evaluate(cvalue) * G4UnitDefinition::GetValueOf(cunit);
          pcuts->SetProductionCut(cut, "e+");
        }
        else if(tag == "gamcut")
        {
          const G4String& cvalue = raux->value;
          const G4String& cunit  = raux->unit;
          if(G4UnitDefinition::GetCategory(cunit) != "Length")
          {
            G4Exception("G4GDMLParser::ImportRegions()", "InvalidRead",
                        FatalException, "Invalid unit for length!");
          }
          G4double cut =
            eval.Evaluate(cvalue) * G4UnitDefinition::GetValueOf(cunit);
          pcuts->SetProductionCut(cut, "gamma");
        }
        else if(tag == "ulimits")
        {
          G4double ustepMax = DBL_MAX, utrakMax = DBL_MAX, utimeMax = DBL_MAX;
          G4double uekinMin = 0., urangMin = 0.;
          const G4String& ulname = raux->value;
          for(auto uaux = raux->auxList->cbegin();
                   uaux != raux->auxList->cend(); ++uaux)
          {
            const G4String& ultag  = uaux->type;
            const G4String& uvalue = uaux->value;
            const G4String& uunit  = uaux->unit;
            G4double ulvalue = eval.Evaluate(uvalue) * eval.Evaluate(uunit);
            if(ultag == "ustepMax")
            {
              ustepMax = ulvalue;
            }
            else if(ultag == "utrakMax")
            {
              utrakMax = ulvalue;
            }
            else if(ultag == "utimeMax")
            {
              utimeMax = ulvalue;
            }
            else if(ultag == "uekinMin")
            {
              uekinMin = ulvalue;
            }
            else if(ultag == "urangMin")
            {
              urangMin = ulvalue;
            }
            else
            {
              G4Exception("G4GDMLParser::ImportRegions()", "ReadError",
                          FatalException, "Invalid definition of user-limits!");
            }
          }
          G4UserLimits* ulimits = new G4UserLimits(
            ulname, ustepMax, utrakMax, utimeMax, uekinMin, urangMin);
          aRegion->SetUserLimits(ulimits);
        }
        else
          continue;  // Ignore unknown tags
      }
    }
  }
}

// --------------------------------------------------------------------
void G4GDMLParser::ExportRegions(G4bool storeReferences)
{
  G4RegionStore* rstore            = G4RegionStore::GetInstance();
  G4ReflectionFactory* reflFactory = G4ReflectionFactory::Instance();
  for(std::size_t i = 0; i < rstore->size(); ++i)
     // Skip default regions associated to worlds
  {
    const G4String& tname = (*rstore)[i]->GetName();
    if(tname.contains("DefaultRegionForParallelWorld"))
      continue;
    const G4String& rname    = writer->GenerateName(tname, (*rstore)[i]);
    rlist                    = new G4GDMLAuxListType();
    G4GDMLAuxStructType raux = { "Region", rname, "", rlist };
    auto rlvol_iter = (*rstore)[i]->GetRootLogicalVolumeIterator();
    for(std::size_t j = 0; j < (*rstore)[i]->GetNumberOfRootVolumes(); ++j)
    {
      G4LogicalVolume* rlvol = *rlvol_iter;
      if(reflFactory->IsReflected(rlvol))
        continue;
      G4String vname = writer->GenerateName(rlvol->GetName(), rlvol);
      if(!storeReferences)
      {
        reader->StripName(vname);
      }
      G4GDMLAuxStructType rsubaux = { "volume", vname, "", 0 };
      rlist->push_back(rsubaux);
      ++rlvol_iter;
    }
    G4double gam_cut
      = (*rstore)[i]->GetProductionCuts()->GetProductionCut("gamma");
    G4GDMLAuxStructType caux1
      = { "gamcut", eval.ConvertToString(gam_cut), "mm", 0 };
    rlist->push_back(caux1);
    G4double e_cut = (*rstore)[i]->GetProductionCuts()->GetProductionCut("e-");
    G4GDMLAuxStructType caux2
      = { "ecut", eval.ConvertToString(e_cut), "mm", 0 };
    rlist->push_back(caux2);
    G4double pos_cut
      = (*rstore)[i]->GetProductionCuts()->GetProductionCut("e+");
    G4GDMLAuxStructType caux3
      = { "poscut", eval.ConvertToString(pos_cut), "mm", 0 };
    rlist->push_back(caux3);
    G4double p_cut
      = (*rstore)[i]->GetProductionCuts()->GetProductionCut("proton");
    G4GDMLAuxStructType caux4
      = { "pcut", eval.ConvertToString(p_cut), "mm", 0 };
    rlist->push_back(caux4);
    if((*rstore)[i]->GetUserLimits())
    {
      const G4Track fake_trk;
      ullist                   = new G4GDMLAuxListType();
      const G4String& utype    = (*rstore)[i]->GetUserLimits()->GetType();
      G4GDMLAuxStructType uaux = { "ulimits", utype, "mm", ullist };
      G4double max_step
        = (*rstore)[i]->GetUserLimits()->GetMaxAllowedStep(fake_trk);
      G4GDMLAuxStructType ulaux1
        = { "ustepMax", eval.ConvertToString(max_step), "mm", 0 };
      ullist->push_back(ulaux1);
      G4double max_trk
        = (*rstore)[i]->GetUserLimits()->GetUserMaxTrackLength(fake_trk);
      G4GDMLAuxStructType ulaux2
        = { "utrakMax", eval.ConvertToString(max_trk), "mm", 0 };
      ullist->push_back(ulaux2);
      G4double max_time
        = (*rstore)[i]->GetUserLimits()->GetUserMaxTime(fake_trk);
      G4GDMLAuxStructType ulaux3
        = { "utimeMax", eval.ConvertToString(max_time), "mm", 0 };
      ullist->push_back(ulaux3);
      G4double min_ekin
        = (*rstore)[i]->GetUserLimits()->GetUserMinEkine(fake_trk);
      G4GDMLAuxStructType ulaux4
        = { "uekinMin", eval.ConvertToString(min_ekin), "mm", 0 };
      ullist->push_back(ulaux4);
      G4double min_rng
        = (*rstore)[i]->GetUserLimits()->GetUserMinRange(fake_trk);
      G4GDMLAuxStructType ulaux5
        = { "urangMin", eval.ConvertToString(min_rng), "mm", 0 };
      ullist->push_back(ulaux5);
      rlist->push_back(uaux);
    }
    AddAuxiliary(raux);
  }
}
