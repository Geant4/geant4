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
// G4tgrVolumeMgr implementation
//
// Author: P.Arce, CIEMAT (November 2007)
// --------------------------------------------------------------------

#include "G4tgrVolumeMgr.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMaterialFactory.hh"
#include "G4tgrRotationMatrixFactory.hh"
#include "G4tgrFileReader.hh"
#include "G4tgrMessenger.hh"
#include "G4tgrSolid.hh"
#include "G4tgrSolidBoolean.hh"
#include "G4tgrSolidMultiUnion.hh"
#include "G4tgrSolidScaled.hh"

G4ThreadLocal G4tgrVolumeMgr* G4tgrVolumeMgr::theInstance = nullptr;

// --------------------------------------------------------------------
G4tgrVolumeMgr::G4tgrVolumeMgr()
{
}

// --------------------------------------------------------------------
G4tgrVolumeMgr::~G4tgrVolumeMgr()
{
  delete theInstance;
}

// --------------------------------------------------------------------
G4tgrVolumeMgr* G4tgrVolumeMgr::GetInstance()
{
  if(theInstance == nullptr)
  {
    theInstance = new G4tgrVolumeMgr;
  }
  return theInstance;
}

// --------------------------------------------------------------------
G4tgrSolid* G4tgrVolumeMgr::CreateSolid(const std::vector<G4String>& wl,
                                        G4bool bVOLUtag)
{
  G4tgrSolid* sol = FindSolid(wl[1]);
  if(sol != nullptr)
  {
    G4String ErrMessage = "Solid already exists... " + wl[1];
    G4Exception("G4tgrVolumeMgr::CreateSolid()", "InvalidSetup", FatalException,
                ErrMessage);
  }

  std::vector<G4String> wlc = wl;
  if(bVOLUtag)
  {
    wlc.pop_back();
  }

  G4String wl2 = wlc[2];
  for(G4int ii = 0; ii < (G4int)wl2.length(); ++ii)
  {
    wl2[ii] = (char)std::toupper(wl2[ii]);
  }
  if((wl2 == "UNION") || (wl2 == "SUBTRACTION") || (wl2 == "INTERSECTION"))
  {
    //---------- Boolean solid
    //---------- Create G4tgrSolidBoolean and fill the solid params
    sol = new G4tgrSolidBoolean(wlc);
  }
  else if(wl2 == "SCALED")
  {
    //---------- Create G4tgrSolidScaled and fill the solid params
    sol = new G4tgrSolidScaled(wlc);
  }
  else if(wl2 == "MULTIUNION")
  {
    //---------- Create G4tgrSolidMultiUnion and fill the solid params
    sol = new G4tgrSolidMultiUnion(wlc);
  }
  else
  {
    //---------- Create G4tgrSolidSimple and fill the solid params
    sol = new G4tgrSolid(wlc);
  }

  return sol;
}

// --------------------------------------------------------------------
void G4tgrVolumeMgr::RegisterMe(G4tgrSolid* sol)
{
  if(theG4tgrSolidMap.find(sol->GetName()) != theG4tgrSolidMap.cend())
  {
    G4String ErrMessage =
      "Cannot be two solids with the same name... " + sol->GetName();
    G4Exception("G4tgrVolumeMgr::RegisterMe()", "InvalidSetup", FatalException,
                ErrMessage);
  }
  theG4tgrSolidMap.insert(G4mapssol::value_type(sol->GetName(), sol));
}

// --------------------------------------------------------------------
void G4tgrVolumeMgr::UnRegisterMe(G4tgrSolid* sol)
{
  if(theG4tgrSolidMap.find(sol->GetName()) != theG4tgrSolidMap.cend())
  {
    G4String ErrMessage =
      "Cannot unregister a solid that is not registered... " + sol->GetName();
    G4Exception("G4tgrSolidMgr::unRegisterMe()", "InvalidSetup", FatalException,
                ErrMessage);
  }
  else
  {
    theG4tgrSolidMap.erase(theG4tgrSolidMap.find(sol->GetName()));
  }
}

// --------------------------------------------------------------------
void G4tgrVolumeMgr::RegisterMe(G4tgrVolume* vol)
{
  theG4tgrVolumeList.push_back(vol);
  if(theG4tgrVolumeMap.find(vol->GetName()) != theG4tgrVolumeMap.cend())
  {
    G4String ErrMessage =
      "Cannot be two volumes with the same name... " + vol->GetName();
    G4Exception("G4tgrVolumeMgr::RegisterMe()", "InvalidSetup", FatalException,
                ErrMessage);
  }
  theG4tgrVolumeMap.insert(G4mapsvol::value_type(vol->GetName(), vol));
}

// --------------------------------------------------------------------
void G4tgrVolumeMgr::UnRegisterMe(G4tgrVolume* vol)
{
  std::vector<G4tgrVolume*>::const_iterator ite;
  for(ite = theG4tgrVolumeList.cbegin();
      ite != theG4tgrVolumeList.cend(); ++ite)
  {
    if((*ite) == vol)
    {
      break;
    }
  }
  if(ite == theG4tgrVolumeList.cend())
  {
    G4String ErrMessage =
      "Cannot unregister a volume not registered... " + vol->GetName();
    G4Exception("G4tgrVolumeMgr::unRegisterMe()", "InvalidSetup",
                FatalException, ErrMessage);
  }
  else
  {
    theG4tgrVolumeList.erase(ite);
  }
  theG4tgrVolumeMap.erase(theG4tgrVolumeMap.find(vol->GetName()));
}

// --------------------------------------------------------------------
void G4tgrVolumeMgr::RegisterParentChild(const G4String& parentName,
                                         const G4tgrPlace* pl)
{
  theG4tgrVolumeTree.insert(G4mmapspl::value_type(parentName, pl));
}

// --------------------------------------------------------------------
G4tgrSolid* G4tgrVolumeMgr::FindSolid(const G4String& volname, G4bool exists)
{
  G4tgrSolid* vol = nullptr;

  G4mapssol::const_iterator svite = theG4tgrSolidMap.find(volname);
  if(svite == theG4tgrSolidMap.cend())
  {
    if(exists)
    {
      for(svite = theG4tgrSolidMap.cbegin();
          svite != theG4tgrSolidMap.cend(); ++svite)
      {
        G4cerr << " VOL:" << (*svite).first << G4endl;
      }
      G4String ErrMessage = "Solid not found... " + volname;
      G4Exception("G4tgrVolumeMgr::FindSolid()", "InvalidSetup", FatalException,
                  ErrMessage);
    }
  }
  else
  {
    vol = const_cast<G4tgrSolid*>((*svite).second);
  }

  return vol;
}

// --------------------------------------------------------------------
G4tgrVolume* G4tgrVolumeMgr::FindVolume(const G4String& volname, G4bool exists)
{
  G4tgrVolume* vol = nullptr;

  G4mapsvol::const_iterator svite = theG4tgrVolumeMap.find(volname);
  if(svite == theG4tgrVolumeMap.cend())
  {
    if(exists)
    {
      for(svite = theG4tgrVolumeMap.cbegin();
          svite != theG4tgrVolumeMap.cend(); ++svite)
      {
        G4cerr << " VOL:" << (*svite).first << G4endl;
      }
      G4String ErrMessage = "Volume not found... " + volname;
      G4Exception("G4tgrVolumeMgr::FindVolume()", "InvalidSetup",
                  FatalException, ErrMessage);
    }
    else
    {
      G4String WarMessage = "Volume does not exists... " + volname;
      G4Exception("G4tgrVolumeMgr::FindVolume()", "SearchFailed", JustWarning,
                  WarMessage);
    }
  }
  else
  {
    vol = const_cast<G4tgrVolume*>((*svite).second);
  }

  return vol;
}

// --------------------------------------------------------------------
std::vector<G4tgrVolume*> G4tgrVolumeMgr::FindVolumes(const G4String& volname,
                                                      G4bool exists)
{
  std::vector<G4tgrVolume*> vols;

  G4mapsvol::const_iterator svite;
  for(svite = theG4tgrVolumeMap.cbegin();
      svite != theG4tgrVolumeMap.cend(); ++svite)
  {
    if(G4tgrUtils::AreWordsEquivalent(volname, (*svite).second->GetName()))
    {
      vols.push_back(const_cast<G4tgrVolume*>((*svite).second));
    }
  }

  if(vols.size() == 0)
  {
    if(exists)
    {
      for(svite = theG4tgrVolumeMap.cbegin();
          svite != theG4tgrVolumeMap.cend(); ++svite)
      {
        G4cerr << " VOL:" << (*svite).first << G4endl;
      }
      G4String ErrMessage = "Volume not found... " + volname;
      G4Exception("G4tgrVolumeMgr::FindVolumes()", "InvalidSetup",
                  FatalException, ErrMessage);
    }
    else
    {
      G4String WarMessage = "Volume does not exists... " + volname;
      G4Exception("G4tgrVolumeMgr::FindVolumes()", "SearchFailed", JustWarning,
                  WarMessage);
    }
  }

  return vols;
}

// --------------------------------------------------------------------
const G4tgrVolume* G4tgrVolumeMgr::GetTopVolume()
{
  //--- Start from any G4tgrVolume and go upwards until you get to the top.
  //    Check that indeed all volumes drive to the same top volume

  const G4tgrVolume* topVol = nullptr;
  for(auto itetv = theG4tgrVolumeMap.cbegin();
           itetv != theG4tgrVolumeMap.cend(); ++itetv)
  {
    const G4tgrVolume* vol = (*itetv).second;
#ifdef G4VERBOSE
    if(G4tgrMessenger::GetVerboseLevel() >= 3)
    {
      G4cout << " G4tgrVolumeMgr::GetTopVolume() - Vol: " << vol->GetName()
             << " no place = " << vol->GetPlacements().size() << G4endl;
    }
#endif

    while(vol->GetPlacements().size() != 0)
    {
      vol = FindVolume((*(vol->GetPlacements()).cbegin())->GetParentName(), 1);
#ifdef G4VERBOSE
      if(G4tgrMessenger::GetVerboseLevel() >= 3)
      {
        G4cout << " G4tgrVolumeMgr::GetTopVolume() - Vol: " << vol->GetName()
               << " N place = " << vol->GetPlacements().size() << G4endl;
      }
#endif
    }
    if((topVol != nullptr) && (topVol != vol) &&
       (topVol->GetType() != "VOLDivision") &&
       (vol->GetType() != "VOLDivision"))
    {
      G4Exception("G4tgrVolumeMgr::GetTopVolume()",
                  "Two world volumes found, second will be taken", JustWarning,
                  (G4String("Both volumes are at the top of a hierarchy: ") +
                   topVol->GetName() + " & " + vol->GetName())
                    .c_str());
    }
    topVol = vol;
  }

  return topVol;
}

// --------------------------------------------------------------------
std::pair<G4mmapspl::iterator, G4mmapspl::iterator>
G4tgrVolumeMgr::GetChildren(const G4String& name)
{
  std::pair<G4mmapspl::iterator, G4mmapspl::iterator> dite;
  dite = theG4tgrVolumeTree.equal_range(name);
  return dite;
}

// --------------------------------------------------------------------
void G4tgrVolumeMgr::DumpVolumeTree()
{
  G4cout << " @@@@@@@@@@@@@@@@ DUMPING G4tgrVolume's Tree  " << G4endl;

  const G4tgrVolume* vol = GetTopVolume();

  DumpVolumeLeaf(vol, 0, 0);
}

// --------------------------------------------------------------------
void G4tgrVolumeMgr::DumpVolumeLeaf(const G4tgrVolume* vol, unsigned int copyNo,
                                    unsigned int leafDepth)
{
  for(std::size_t ii = 0; ii < leafDepth; ++ii)
  {
    G4cout << "  ";
  }
  G4cout << " VOL:(" << leafDepth << ")" << vol->GetName() << "   copy No "
         << copyNo << G4endl;

  //---------- construct the children of this VOL
  std::pair<G4mmapspl::iterator, G4mmapspl::iterator> children =
    GetChildren(vol->GetName());
  G4mmapspl::const_iterator cite;

  ++leafDepth;
  for(cite = children.first; cite != children.second; ++cite)
  {
    //---- find G4tgrVolume pointed by G4tgrPlace
    const G4tgrPlace* pla       = (*cite).second;
    const G4tgrVolume* volchild = pla->GetVolume();
    //--- find copyNo
    unsigned int cn = pla->GetCopyNo();
    DumpVolumeLeaf(volchild, cn, leafDepth);
  }
}

// --------------------------------------------------------------------
void G4tgrVolumeMgr::DumpSummary()
{
  //---------- Dump number of objects of each class
  G4cout << " @@@@@@@@@@@@@@@@@@ Dumping Detector Summary " << G4endl;
  G4cout << " @@@ Geometry built inside world volume: "
         << GetTopVolume()->GetName() << G4endl;
  G4cout << " Number of G4tgrVolume's: " << theG4tgrVolumeMap.size() << G4endl;
  unsigned int nPlace = 0;
  for(auto cite = theG4tgrVolumeMap.cbegin();
      cite != theG4tgrVolumeMap.cend(); ++cite)
  {
    nPlace += ((*cite).second)->GetPlacements().size();
  }
  G4cout << " Number of G4tgrPlace's: " << nPlace << G4endl;

  G4tgrMaterialFactory* matef = G4tgrMaterialFactory::GetInstance();
  G4cout << " Number of G4tgrIsotope's: " << matef->GetIsotopeList().size()
         << G4endl;
  G4cout << " Number of G4tgrElement's: " << matef->GetElementList().size()
         << G4endl;
  G4cout << " Number of G4tgrMaterial's: " << matef->GetMaterialList().size()
         << G4endl;

  G4tgrRotationMatrixFactory* rotmf = G4tgrRotationMatrixFactory::GetInstance();
  G4cout << " Number of G4tgrRotationMatrix's: "
         << rotmf->GetRotMatList().size() << G4endl;

  //---------- Dump detail list of objects of each class
  DumpVolumeTree();

  matef->DumpIsotopeList();
  matef->DumpElementList();
  matef->DumpMaterialList();
  rotmf->DumpRotmList();
}
