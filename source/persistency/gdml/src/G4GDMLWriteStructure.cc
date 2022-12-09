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
// G4GDMLWriteStructure implementation
//
// Author: Zoltan Torzsok, November 2007
// --------------------------------------------------------------------

#include "G4GDMLWriteStructure.hh"
#include "G4GDMLEvaluator.hh"

#include "G4Material.hh"
#include "G4ReflectedSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ReflectionFactory.hh"
#include "G4PVDivision.hh"
#include "G4PVReplica.hh"
#include "G4Region.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

#include "G4VSensitiveDetector.hh"
#include "G4AssemblyStore.hh"
#include "G4AssemblyVolume.hh"

G4int G4GDMLWriteStructure::levelNo = 0;  // Counter for level being exported

// --------------------------------------------------------------------
G4GDMLWriteStructure::G4GDMLWriteStructure()
  : G4GDMLWriteParamvol()
  , maxLevel(INT_MAX)
{
  reflFactory = G4ReflectionFactory::Instance();
}

// --------------------------------------------------------------------
G4GDMLWriteStructure::~G4GDMLWriteStructure()
{
}

// --------------------------------------------------------------------
void G4GDMLWriteStructure::DivisionvolWrite(
  xercesc::DOMElement* volumeElement, const G4PVDivision* const divisionvol)
{
  EAxis axis       = kUndefined;
  G4int number     = 0;
  G4double width   = 0.0;
  G4double offset  = 0.0;
  G4bool consuming = false;

  divisionvol->GetReplicationData(axis, number, width, offset, consuming);
  axis = divisionvol->GetDivisionAxis();

  G4String unitString("mm");
  G4String axisString("kUndefined");
  if(axis == kXAxis)
  {
    axisString = "kXAxis";
  }
  else if(axis == kYAxis)
  {
    axisString = "kYAxis";
  }
  else if(axis == kZAxis)
  {
    axisString = "kZAxis";
  }
  else if(axis == kRho)
  {
    axisString = "kRho";
  }
  else if(axis == kPhi)
  {
    axisString = "kPhi";
    unitString = "rad";
  }

  const G4String name = GenerateName(divisionvol->GetName(), divisionvol);
  const G4String volumeref =
    GenerateName(divisionvol->GetLogicalVolume()->GetName(),
                 divisionvol->GetLogicalVolume());

  xercesc::DOMElement* divisionvolElement = NewElement("divisionvol");
  divisionvolElement->setAttributeNode(NewAttribute("axis", axisString));
  divisionvolElement->setAttributeNode(NewAttribute("number", number));
  divisionvolElement->setAttributeNode(NewAttribute("width", width));
  divisionvolElement->setAttributeNode(NewAttribute("offset", offset));
  divisionvolElement->setAttributeNode(NewAttribute("unit", unitString));
  xercesc::DOMElement* volumerefElement = NewElement("volumeref");
  volumerefElement->setAttributeNode(NewAttribute("ref", volumeref));
  divisionvolElement->appendChild(volumerefElement);
  volumeElement->appendChild(divisionvolElement);
}

// --------------------------------------------------------------------
void G4GDMLWriteStructure::PhysvolWrite(xercesc::DOMElement* volumeElement,
                                        const G4VPhysicalVolume* const physvol,
                                        const G4Transform3D& T,
                                        const G4String& ModuleName)
{
  HepGeom::Scale3D scale;
  HepGeom::Rotate3D rotate;
  HepGeom::Translate3D translate;

  T.getDecomposition(scale, rotate, translate);

  const G4ThreeVector scl(scale(0, 0), scale(1, 1), scale(2, 2));
  const G4ThreeVector rot = GetAngles(rotate.getRotation());
  const G4ThreeVector pos = T.getTranslation();

  const G4String name    = GenerateName(physvol->GetName(), physvol);
  const G4int copynumber = physvol->GetCopyNo();

  xercesc::DOMElement* physvolElement = NewElement("physvol");
  physvolElement->setAttributeNode(NewAttribute("name", name));
  if(copynumber)
  {
    physvolElement->setAttributeNode(NewAttribute("copynumber", copynumber));
  }

  volumeElement->appendChild(physvolElement);

  G4LogicalVolume* lv;
  // Is it reflected?
  if(reflFactory->IsReflected(physvol->GetLogicalVolume()))
  {
    lv = reflFactory->GetConstituentLV(physvol->GetLogicalVolume());
  }
  else
  {
    lv = physvol->GetLogicalVolume();
  }

  const G4String volumeref = GenerateName(lv->GetName(), lv);

  if(ModuleName.empty())
  {
    xercesc::DOMElement* volumerefElement = NewElement("volumeref");
    volumerefElement->setAttributeNode(NewAttribute("ref", volumeref));
    physvolElement->appendChild(volumerefElement);
  }
  else
  {
    xercesc::DOMElement* fileElement = NewElement("file");
    fileElement->setAttributeNode(NewAttribute("name", ModuleName));
    fileElement->setAttributeNode(NewAttribute("volname", volumeref));
    physvolElement->appendChild(fileElement);
  }

  if(std::fabs(pos.x()) > kLinearPrecision ||
     std::fabs(pos.y()) > kLinearPrecision ||
     std::fabs(pos.z()) > kLinearPrecision)
  {
    PositionWrite(physvolElement, name + "_pos", pos);
  }
  if(std::fabs(rot.x()) > kAngularPrecision ||
     std::fabs(rot.y()) > kAngularPrecision ||
     std::fabs(rot.z()) > kAngularPrecision)
  {
    RotationWrite(physvolElement, name + "_rot", rot);
  }
  if(std::fabs(scl.x() - 1.0) > kRelativePrecision ||
     std::fabs(scl.y() - 1.0) > kRelativePrecision ||
     std::fabs(scl.z() - 1.0) > kRelativePrecision)
  {
    ScaleWrite(physvolElement, name + "_scl", scl);
  }
}

// --------------------------------------------------------------------
void G4GDMLWriteStructure::ReplicavolWrite(
  xercesc::DOMElement* volumeElement, const G4VPhysicalVolume* const replicavol)
{
  EAxis axis       = kUndefined;
  G4int number     = 0;
  G4double width   = 0.0;
  G4double offset  = 0.0;
  G4bool consuming = false;
  G4String unitString("mm");

  replicavol->GetReplicationData(axis, number, width, offset, consuming);

  const G4String volumeref = GenerateName(
    replicavol->GetLogicalVolume()->GetName(), replicavol->GetLogicalVolume());

  xercesc::DOMElement* replicavolElement = NewElement("replicavol");
  replicavolElement->setAttributeNode(NewAttribute("number", number));
  xercesc::DOMElement* volumerefElement = NewElement("volumeref");
  volumerefElement->setAttributeNode(NewAttribute("ref", volumeref));
  replicavolElement->appendChild(volumerefElement);
  xercesc::DOMElement* replicateElement = NewElement("replicate_along_axis");
  replicavolElement->appendChild(replicateElement);

  xercesc::DOMElement* dirElement = NewElement("direction");
  if(axis == kXAxis)
  {
    dirElement->setAttributeNode(NewAttribute("x", "1"));
  }
  else if(axis == kYAxis)
  {
    dirElement->setAttributeNode(NewAttribute("y", "1"));
  }
  else if(axis == kZAxis)
  {
    dirElement->setAttributeNode(NewAttribute("z", "1"));
  }
  else if(axis == kRho)
  {
    dirElement->setAttributeNode(NewAttribute("rho", "1"));
  }
  else if(axis == kPhi)
  {
    dirElement->setAttributeNode(NewAttribute("phi", "1"));
    unitString = "rad";
  }
  replicateElement->appendChild(dirElement);

  xercesc::DOMElement* widthElement = NewElement("width");
  widthElement->setAttributeNode(NewAttribute("value", width));
  widthElement->setAttributeNode(NewAttribute("unit", unitString));
  replicateElement->appendChild(widthElement);

  xercesc::DOMElement* offsetElement = NewElement("offset");
  offsetElement->setAttributeNode(NewAttribute("value", offset));
  offsetElement->setAttributeNode(NewAttribute("unit", unitString));
  replicateElement->appendChild(offsetElement);

  volumeElement->appendChild(replicavolElement);
}

// --------------------------------------------------------------------
void G4GDMLWriteStructure::AssemblyWrite(xercesc::DOMElement* volumeElement,
                                         const G4int assemblyID)
{
  G4AssemblyStore* assemblies  = G4AssemblyStore::GetInstance();
  G4AssemblyVolume* myassembly = assemblies->GetAssembly(assemblyID);

  xercesc::DOMElement* assemblyElement = NewElement("assembly");
  G4String name = "Assembly_" + std::to_string(assemblyID);

  assemblyElement->setAttributeNode(NewAttribute("name", name));

  auto vit = myassembly->GetTripletsIterator();

  G4int depth = 0;
  const G4String ModuleName;

  for(std::size_t i5 = 0; i5 < myassembly->TotalTriplets(); ++i5)
  {
    G4LogicalVolume* lvol = (*vit).GetVolume();
    if (lvol == nullptr)
    {
      G4String message = "Nested assemblies not yet supported for exporting. Sorry!";
      G4Exception("G4GDMLWriteStructure::AssemblyWrite()", "InvalidSetup",
                  FatalException, message);
      return;
    }
    TraverseVolumeTree(lvol, depth + 1);

    const G4ThreeVector rot = GetAngles((*vit).GetRotation()->inverse());
    const G4ThreeVector pos = (*vit).GetTranslation();

    const G4String pname =
      GenerateName((*vit).GetVolume()->GetName() + "_pv", &(*vit));

    xercesc::DOMElement* physvolElement = NewElement("physvol");
    physvolElement->setAttributeNode(NewAttribute("name", pname));

    assemblyElement->appendChild(physvolElement);

    const G4String volumeref =
      GenerateName((*vit).GetVolume()->GetName(), (*vit).GetVolume());

    xercesc::DOMElement* volumerefElement = NewElement("volumeref");
    volumerefElement->setAttributeNode(NewAttribute("ref", volumeref));
    physvolElement->appendChild(volumerefElement);

    if(std::fabs(pos.x()) > kLinearPrecision ||
       std::fabs(pos.y()) > kLinearPrecision ||
       std::fabs(pos.z()) > kLinearPrecision)
    {
      PositionWrite(physvolElement,name+"_position_" + std::to_string(i5), pos);
    }

    if(std::fabs(rot.x()) > kAngularPrecision ||
       std::fabs(rot.y()) > kAngularPrecision ||
       std::fabs(rot.z()) > kAngularPrecision)
    {
      RotationWrite(physvolElement,name+"_rotation_" + std::to_string(i5), rot);
    }
    ++vit;
  }

  volumeElement->appendChild(assemblyElement);
}

// --------------------------------------------------------------------
void G4GDMLWriteStructure::BorderSurfaceCache(
  const G4LogicalBorderSurface* const bsurf)
{
  if(bsurf == nullptr)
  {
    return;
  }

  const G4SurfaceProperty* psurf = bsurf->GetSurfaceProperty();

  // Generate the new element for border-surface
  //
  const G4String& bsname             = GenerateName(bsurf->GetName(), bsurf);
  const G4String& psname             = GenerateName(psurf->GetName(), psurf);
  xercesc::DOMElement* borderElement = NewElement("bordersurface");
  borderElement->setAttributeNode(NewAttribute("name", bsname));
  borderElement->setAttributeNode(NewAttribute("surfaceproperty", psname));

  const G4String volumeref1 =
    GenerateName(bsurf->GetVolume1()->GetName(), bsurf->GetVolume1());
  const G4String volumeref2 =
    GenerateName(bsurf->GetVolume2()->GetName(), bsurf->GetVolume2());
  xercesc::DOMElement* volumerefElement1 = NewElement("physvolref");
  xercesc::DOMElement* volumerefElement2 = NewElement("physvolref");
  volumerefElement1->setAttributeNode(NewAttribute("ref", volumeref1));
  volumerefElement2->setAttributeNode(NewAttribute("ref", volumeref2));
  borderElement->appendChild(volumerefElement1);
  borderElement->appendChild(volumerefElement2);

  if(FindOpticalSurface(psurf))
  {
    const G4OpticalSurface* opsurf =
      dynamic_cast<const G4OpticalSurface*>(psurf);
    if(opsurf == nullptr)
    {
      G4Exception("G4GDMLWriteStructure::BorderSurfaceCache()", "InvalidSetup",
                  FatalException, "No optical surface found!");
      return;
    }
    OpticalSurfaceWrite(solidsElement, opsurf);
  }

  borderElementVec.push_back(borderElement);
}

// --------------------------------------------------------------------
void G4GDMLWriteStructure::SkinSurfaceCache(
  const G4LogicalSkinSurface* const ssurf)
{
  if(ssurf == nullptr)
  {
    return;
  }

  const G4SurfaceProperty* psurf = ssurf->GetSurfaceProperty();

  // Generate the new element for border-surface
  //
  const G4String& ssname           = GenerateName(ssurf->GetName(), ssurf);
  const G4String& psname           = GenerateName(psurf->GetName(), psurf);
  xercesc::DOMElement* skinElement = NewElement("skinsurface");
  skinElement->setAttributeNode(NewAttribute("name", ssname));
  skinElement->setAttributeNode(NewAttribute("surfaceproperty", psname));

  const G4String volumeref = GenerateName(ssurf->GetLogicalVolume()->GetName(),
                                          ssurf->GetLogicalVolume());
  xercesc::DOMElement* volumerefElement = NewElement("volumeref");
  volumerefElement->setAttributeNode(NewAttribute("ref", volumeref));
  skinElement->appendChild(volumerefElement);

  if(FindOpticalSurface(psurf))
  {
    const G4OpticalSurface* opsurf =
      dynamic_cast<const G4OpticalSurface*>(psurf);
    if(opsurf == nullptr)
    {
      G4Exception("G4GDMLWriteStructure::SkinSurfaceCache()", "InvalidSetup",
                  FatalException, "No optical surface found!");
      return;
    }
    OpticalSurfaceWrite(solidsElement, opsurf);
  }

  skinElementVec.push_back(skinElement);
}

// --------------------------------------------------------------------
G4bool G4GDMLWriteStructure::FindOpticalSurface(const G4SurfaceProperty* psurf)
{
  const G4OpticalSurface* osurf = dynamic_cast<const G4OpticalSurface*>(psurf);
  auto pos = std::find(opt_vec.cbegin(), opt_vec.cend(), osurf);
  if(pos != opt_vec.cend())
  {
    return false;
  }  // item already created!

  opt_vec.push_back(osurf);  // cache it for future reference
  return true;
}

// --------------------------------------------------------------------
const G4LogicalSkinSurface* G4GDMLWriteStructure::GetSkinSurface(
  const G4LogicalVolume* const lvol)
{
  G4LogicalSkinSurface* surf = 0;
  std::size_t nsurf = G4LogicalSkinSurface::GetNumberOfSkinSurfaces();
  if(nsurf)
  {
    const G4LogicalSkinSurfaceTable* stable =
      G4LogicalSkinSurface::GetSurfaceTable();
    for(auto pos = stable->cbegin(); pos != stable->cend(); ++pos)
    {
      if(lvol == (*pos)->GetLogicalVolume())
      {
        surf = *pos;
        break;
      }
    }
  }
  return surf;
}

// --------------------------------------------------------------------
const G4LogicalBorderSurface* G4GDMLWriteStructure::GetBorderSurface(
  const G4VPhysicalVolume* const pvol)
{
  G4LogicalBorderSurface* surf = nullptr;
  std::size_t nsurf = G4LogicalBorderSurface::GetNumberOfBorderSurfaces();
  if(nsurf)
  {
    const G4LogicalBorderSurfaceTable* btable =
      G4LogicalBorderSurface::GetSurfaceTable();
    for(auto pos = btable->cbegin(); pos != btable->cend(); ++pos)
    {
      if(pvol == pos->first.first)  // just the first in the couple
      {                             // could be enough?
        surf = pos->second;         // break;
        BorderSurfaceCache(surf);
      }
    }
  }
  return surf;
}

// --------------------------------------------------------------------
void G4GDMLWriteStructure::SurfacesWrite()
{
#ifdef G4VERBOSE
  G4cout << "G4GDML: Writing surfaces..." << G4endl;
#endif
  for(auto pos = skinElementVec.cbegin();
           pos != skinElementVec.cend(); ++pos)
  {
    structureElement->appendChild(*pos);
  }
  for(auto pos = borderElementVec.cbegin();
           pos != borderElementVec.cend(); ++pos)
  {
    structureElement->appendChild(*pos);
  }
}

// --------------------------------------------------------------------
void G4GDMLWriteStructure::StructureWrite(xercesc::DOMElement* gdmlElement)
{
#ifdef G4VERBOSE
  G4cout << "G4GDML: Writing structure..." << G4endl;
#endif

  // filling the list of phys volumes that are parts of assemblies

  G4AssemblyStore* assemblies = G4AssemblyStore::GetInstance();

  for(auto it = assemblies->cbegin(); it != assemblies->cend(); ++it)
  {
    auto vit = (*it)->GetVolumesIterator();

    for(std::size_t i5 = 0; i5 < (*it)->TotalImprintedVolumes(); ++i5)
    {
      G4String pvname = (*vit)->GetName();
      std::size_t pos = pvname.find("_impr_") + 6;
      G4String impID  = pvname.substr(pos);

      pos   = impID.find("_");
      impID = impID.substr(0, pos);

      assemblyVolMap[*vit] = (*it)->GetAssemblyID();
      imprintsMap[*vit]    = std::atoi(impID.c_str());
      ++vit;
    }
  }

  structureElement = NewElement("structure");
  gdmlElement->appendChild(structureElement);
}

// --------------------------------------------------------------------
G4Transform3D G4GDMLWriteStructure::TraverseVolumeTree(
  const G4LogicalVolume* const volumePtr, const G4int depth)
{
  if(VolumeMap().find(volumePtr) != VolumeMap().cend())
  {
    return VolumeMap()[volumePtr];  // Volume is already processed
  }

  G4VSolid* solidPtr = volumePtr->GetSolid();
  G4Transform3D R, invR;
  G4int trans = 0;

  std::map<const G4LogicalVolume*, G4GDMLAuxListType>::iterator auxiter;

  ++levelNo;

  while(true)  // Solve possible displacement/reflection
  {            // of the referenced solid!
    if(trans > maxTransforms)
    {
      G4String ErrorMessage = "Referenced solid in volume '" +
                              volumePtr->GetName() +
                              "' was displaced/reflected too many times!";
      G4Exception("G4GDMLWriteStructure::TraverseVolumeTree()", "InvalidSetup",
                  FatalException, ErrorMessage);
    }

    if(G4ReflectedSolid* refl = dynamic_cast<G4ReflectedSolid*>(solidPtr))
    {
      R        = R * refl->GetTransform3D();
      solidPtr = refl->GetConstituentMovedSolid();
      ++trans;
      continue;
    }

    if(G4DisplacedSolid* disp = dynamic_cast<G4DisplacedSolid*>(solidPtr))
    {
      R        = R * G4Transform3D(disp->GetObjectRotation(),
                            disp->GetObjectTranslation());
      solidPtr = disp->GetConstituentMovedSolid();
      ++trans;
      continue;
    }

    break;
  }

  // check if it is a reflected volume
  G4LogicalVolume* tmplv = const_cast<G4LogicalVolume*>(volumePtr);

  if(reflFactory->IsReflected(tmplv))
  {
    tmplv = reflFactory->GetConstituentLV(tmplv);
    if(VolumeMap().find(tmplv) != VolumeMap().cend())
    {
      return R;  // Volume is already processed
    }
  }

  // Only compute the inverse when necessary!
  //
  if(trans > 0)
  {
    invR = R.inverse();
  }

  const G4String name = GenerateName(tmplv->GetName(), tmplv);

  G4String materialref = "NULL";

  if(volumePtr->GetMaterial())
  {
    materialref = GenerateName(volumePtr->GetMaterial()->GetName(),
                               volumePtr->GetMaterial());
  }

  const G4String solidref = GenerateName(solidPtr->GetName(), solidPtr);

  xercesc::DOMElement* volumeElement = NewElement("volume");
  volumeElement->setAttributeNode(NewAttribute("name", name));
  xercesc::DOMElement* materialrefElement = NewElement("materialref");
  materialrefElement->setAttributeNode(NewAttribute("ref", materialref));
  volumeElement->appendChild(materialrefElement);
  xercesc::DOMElement* solidrefElement = NewElement("solidref");
  solidrefElement->setAttributeNode(NewAttribute("ref", solidref));
  volumeElement->appendChild(solidrefElement);

  std::size_t daughterCount = volumePtr->GetNoDaughters();

  if(levelNo == maxLevel)  // Stop exporting if reached levels limit
  {
    daughterCount = 0;
  }

  std::map<G4int, std::vector<G4int> > assemblyIDToAddedImprints;

  for(std::size_t i = 0; i < daughterCount; ++i)  // Traverse all the children!
  {
    const G4VPhysicalVolume* const physvol = volumePtr->GetDaughter(i);
    const G4String ModuleName              = Modularize(physvol, depth);

    G4Transform3D daughterR;

    if(ModuleName.empty())  // Check if subtree requested to be
    {                       // a separate module!
      daughterR = TraverseVolumeTree(physvol->GetLogicalVolume(), depth + 1);
    }
    else
    {
      G4GDMLWriteStructure writer;
      daughterR = writer.Write(ModuleName, physvol->GetLogicalVolume(),
                               SchemaLocation, depth + 1);
    }

    if(const G4PVDivision* const divisionvol =
         dynamic_cast<const G4PVDivision*>(physvol))  // Is it division?
    {
      if(!G4Transform3D::Identity.isNear(invR * daughterR, kRelativePrecision))
      {
        G4String ErrorMessage = "Division volume in '" + name +
                                "' can not be related to reflected solid!";
        G4Exception("G4GDMLWriteStructure::TraverseVolumeTree()",
                    "InvalidSetup", FatalException, ErrorMessage);
      }
      DivisionvolWrite(volumeElement, divisionvol);
    }
    else
    {
      if(physvol->IsParameterised())  // Is it a paramvol?
      {
        if(!G4Transform3D::Identity.isNear(invR * daughterR,
                                           kRelativePrecision))
        {
          G4String ErrorMessage = "Parameterised volume in '" + name +
                                  "' can not be related to reflected solid!";
          G4Exception("G4GDMLWriteStructure::TraverseVolumeTree()",
                      "InvalidSetup", FatalException, ErrorMessage);
        }
        ParamvolWrite(volumeElement, physvol);
      }
      else
      {
        if(physvol->IsReplicated())  // Is it a replicavol?
        {
          if(!G4Transform3D::Identity.isNear(invR * daughterR,
                                             kRelativePrecision))
          {
            G4String ErrorMessage = "Replica volume in '" + name +
                                    "' can not be related to reflected solid!";
            G4Exception("G4GDMLWriteStructure::TraverseVolumeTree()",
                        "InvalidSetup", FatalException, ErrorMessage);
          }
          ReplicavolWrite(volumeElement, physvol);
        }
        else  // Is it a physvol or an assembly?
        {
          if(assemblyVolMap.find(physvol) != assemblyVolMap.cend())
          {
            G4int assemblyID = assemblyVolMap[physvol];

            G4String assemblyref = "Assembly_" + std::to_string(assemblyID);

            // here I need to retrieve the imprint ID

            G4int imprintID = imprintsMap[physvol];

            // there are 2 steps:
            //
            // 1) add assembly to the structure if that has not yet been done
            // (but after the constituents volumes have been added)
            //

            if(std::find(addedAssemblies.cbegin(), addedAssemblies.cend(),
                         assemblyID) == addedAssemblies.cend())
            {
              AssemblyWrite(structureElement, assemblyID);
              addedAssemblies.push_back(assemblyID);
              assemblyIDToAddedImprints[assemblyID] = std::vector<G4int>();
            }

            // 2) add the assembly (as physical volume) to the mother volume
            // (but only once), using it's original position and rotation.
            //

            // here I need a check if assembly has been already added to the
            // mother volume
            std::vector<G4int>& addedImprints = assemblyIDToAddedImprints[assemblyID];
            if(std::find(addedImprints.cbegin(), addedImprints.cend(),
                         imprintID) == addedImprints.cend())
            {
              G4String imprintname = "Imprint_" + std::to_string(imprintID) + "_";
              imprintname          = GenerateName(imprintname, physvol);

              // I need to get those two from the  container of imprints from
              // the assembly I have the imprint ID, I need to get pos and rot
              //
              G4Transform3D& transf = G4AssemblyStore::GetInstance()
                                        ->GetAssembly(assemblyID)
                                        ->GetImprintTransformation(imprintID);

              HepGeom::Scale3D scale;
              HepGeom::Rotate3D rotate;
              HepGeom::Translate3D translate;

              transf.getDecomposition(scale, rotate, translate);

              const G4ThreeVector scl(scale(0, 0), scale(1, 1), scale(2, 2));
              const G4ThreeVector rot = GetAngles(rotate.getRotation().inverse());
              const G4ThreeVector pos = transf.getTranslation();

              // here I need a normal physvol referencing to my assemblyref

              xercesc::DOMElement* physvolElement = NewElement("physvol");
              physvolElement->setAttributeNode(NewAttribute("name", imprintname));

              xercesc::DOMElement* volumerefElement = NewElement("volumeref");
              volumerefElement->setAttributeNode(NewAttribute("ref", assemblyref));
              physvolElement->appendChild(volumerefElement);

              if(std::fabs(pos.x()) > kLinearPrecision ||
                 std::fabs(pos.y()) > kLinearPrecision ||
                 std::fabs(pos.z()) > kLinearPrecision)
              {
                PositionWrite(physvolElement, imprintname + "_pos", pos);
              }
              if(std::fabs(rot.x()) > kAngularPrecision ||
                 std::fabs(rot.y()) > kAngularPrecision ||
                 std::fabs(rot.z()) > kAngularPrecision)
              {
                RotationWrite(physvolElement, imprintname + "_rot", rot);
              }
              if(std::fabs(scl.x() - 1.0) > kRelativePrecision ||
                 std::fabs(scl.y() - 1.0) > kRelativePrecision ||
                 std::fabs(scl.z() - 1.0) > kRelativePrecision)
              {
                ScaleWrite(physvolElement, name + "_scl", scl);
              }

              volumeElement->appendChild(physvolElement);
              //
              addedImprints.push_back(imprintID);
            }
          }
          else  // not part of assembly, so a normal physical volume
          {
            G4RotationMatrix rot;
            if(physvol->GetFrameRotation() != nullptr)
            {
              rot = *(physvol->GetFrameRotation());
            }
            G4Transform3D P(rot, physvol->GetObjectTranslation());

            PhysvolWrite(volumeElement, physvol, invR * P * daughterR,
                         ModuleName);
          }
        }
      }
    }
    // BorderSurfaceCache(GetBorderSurface(physvol));
    GetBorderSurface(physvol);
  }

  if(cexport)
  {
    ExportEnergyCuts(volumePtr);
  }
  // Add optional energy cuts

  if(sdexport)
  {
    ExportSD(volumePtr);
  }
  // Add optional SDs

  // Here write the auxiliary info
  //
  auxiter = auxmap.find(volumePtr);
  if(auxiter != auxmap.cend())
  {
    AddAuxInfo(&(auxiter->second), volumeElement);
  }

  structureElement->appendChild(volumeElement);
  // Append the volume AFTER traversing the children so that
  // the order of volumes will be correct!

  VolumeMap()[tmplv] = R;

  AddExtension(volumeElement, volumePtr);
  // Add any possible user defined extension attached to a volume

  AddMaterial(volumePtr->GetMaterial());
  // Add the involved materials and solids!

  AddSolid(solidPtr);

  SkinSurfaceCache(GetSkinSurface(volumePtr));

  return R;
}

// --------------------------------------------------------------------
void G4GDMLWriteStructure::AddVolumeAuxiliary(G4GDMLAuxStructType myaux,
                                              const G4LogicalVolume* const lvol)
{
  auto pos = auxmap.find(lvol);

  if(pos == auxmap.cend())
  {
    auxmap[lvol] = G4GDMLAuxListType();
  }

  auxmap[lvol].push_back(myaux);
}

// --------------------------------------------------------------------
void G4GDMLWriteStructure::SetEnergyCutsExport(G4bool fcuts)
{
  cexport = fcuts;
}

// --------------------------------------------------------------------
void G4GDMLWriteStructure::ExportEnergyCuts(const G4LogicalVolume* const lvol)
{
  G4GDMLEvaluator eval;
  G4ProductionCuts* pcuts     = lvol->GetRegion()->GetProductionCuts();
  G4ProductionCutsTable* ctab = G4ProductionCutsTable::GetProductionCutsTable();
  G4Gamma* gamma              = G4Gamma::Gamma();
  G4Electron* eminus          = G4Electron::Electron();
  G4Positron* eplus           = G4Positron::Positron();
  G4Proton* proton            = G4Proton::Proton();

  G4double gamma_cut = ctab->ConvertRangeToEnergy(
    gamma, lvol->GetMaterial(), pcuts->GetProductionCut("gamma"));
  G4double eminus_cut = ctab->ConvertRangeToEnergy(
    eminus, lvol->GetMaterial(), pcuts->GetProductionCut("e-"));
  G4double eplus_cut = ctab->ConvertRangeToEnergy(
    eplus, lvol->GetMaterial(), pcuts->GetProductionCut("e+"));
  G4double proton_cut = ctab->ConvertRangeToEnergy(
    proton, lvol->GetMaterial(), pcuts->GetProductionCut("proton"));

  G4GDMLAuxStructType gammainfo  = { "gammaECut",
                                    eval.ConvertToString(gamma_cut), "MeV", 0 };
  G4GDMLAuxStructType eminusinfo = { "electronECut",
                                     eval.ConvertToString(eminus_cut), "MeV",
                                     0 };
  G4GDMLAuxStructType eplusinfo  = { "positronECut",
                                    eval.ConvertToString(eplus_cut), "MeV", 0 };
  G4GDMLAuxStructType protinfo   = { "protonECut",
                                   eval.ConvertToString(proton_cut), "MeV", 0 };

  AddVolumeAuxiliary(gammainfo, lvol);
  AddVolumeAuxiliary(eminusinfo, lvol);
  AddVolumeAuxiliary(eplusinfo, lvol);
  AddVolumeAuxiliary(protinfo, lvol);
}

// --------------------------------------------------------------------
void G4GDMLWriteStructure::SetSDExport(G4bool fsd)
{
  sdexport = fsd;
}

// --------------------------------------------------------------------
void G4GDMLWriteStructure::ExportSD(const G4LogicalVolume* const lvol)
{
  G4VSensitiveDetector* sd = lvol->GetSensitiveDetector();

  if(sd != nullptr)
  {
    G4String SDname = sd->GetName();

    G4GDMLAuxStructType SDinfo = { "SensDet", SDname, "", 0 };
    AddVolumeAuxiliary(SDinfo, lvol);
  }
}

// --------------------------------------------------------------------
G4int G4GDMLWriteStructure::GetMaxExportLevel() const
{
  return maxLevel;
}

// --------------------------------------------------------------------
void G4GDMLWriteStructure::SetMaxExportLevel(G4int level)
{
  if(level <= 0)
  {
    G4Exception("G4GDMLWriteStructure::TraverseVolumeTree()", "InvalidSetup",
                FatalException, "Levels to export must be greater than zero!");
    return;
  }
  maxLevel = level;
  levelNo  = 0;
}
