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
// G4tgbGeometryDumper implementation
//
// Author: P.Arce, CIEMAT (November 2007)
// --------------------------------------------------------------------

#include "G4tgbGeometryDumper.hh"

#include "G4tgrMessenger.hh"

#include "G4UIcommand.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4Trd.hh"
#include "G4Para.hh"
#include "G4Torus.hh"
#include "G4Hype.hh"
#include "G4Polycone.hh"
#include "G4GenericPolycone.hh"
#include "G4Polyhedra.hh"
#include "G4EllipticalTube.hh"
#include "G4Ellipsoid.hh"
#include "G4EllipticalCone.hh"
#include "G4Hype.hh"
#include "G4Tet.hh"
#include "G4TwistedBox.hh"
#include "G4TwistedTrap.hh"
#include "G4TwistedTrd.hh"
#include "G4TwistedTubs.hh"
#include "G4MultiUnion.hh"
#include "G4ScaledSolid.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4BooleanSolid.hh"
#include "G4ReflectionFactory.hh"
#include "G4ReflectedSolid.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GeometryTolerance.hh"
#include "G4VPVParameterisation.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

// --------------------------------------------------------------------
G4ThreadLocal G4tgbGeometryDumper* G4tgbGeometryDumper::theInstance = nullptr;

// --------------------------------------------------------------------
G4tgbGeometryDumper::G4tgbGeometryDumper()
{
}

// --------------------------------------------------------------------
G4tgbGeometryDumper* G4tgbGeometryDumper::GetInstance()
{
  if(theInstance == nullptr)
  {
    theInstance = new G4tgbGeometryDumper;
  }
  return theInstance;
}

// --------------------------------------------------------------------
void G4tgbGeometryDumper::DumpGeometry(const G4String& fname)
{
  theFile = new std::ofstream(fname);

  G4VPhysicalVolume* pv = GetTopPhysVol();
  DumpPhysVol(pv);  // dump volume and recursively it will dump all hierarchy
}

// --------------------------------------------------------------------
G4VPhysicalVolume* G4tgbGeometryDumper::GetTopPhysVol()
{
  G4PhysicalVolumeStore* pvstore = G4PhysicalVolumeStore::GetInstance();
  G4VPhysicalVolume* pv = *(pvstore->cbegin());
  for(;;)
  {
    G4LogicalVolume* lv = pv->GetMotherLogical();
    if(lv == 0)
    {
      break;
    }

    //----- look for one PV of this LV
    for(auto ite = pvstore->cbegin(); ite != pvstore->cend(); ++ite)
    {
      pv = (*ite);
      if(pv->GetLogicalVolume() == lv)
      {
        break;
      }
    }
  }

  return pv;
}

// --------------------------------------------------------------------
G4tgbGeometryDumper::~G4tgbGeometryDumper()
{
}

// --------------------------------------------------------------------
void G4tgbGeometryDumper::DumpPhysVol(G4VPhysicalVolume* pv)
{
  //--- Dump logical volume first
  G4LogicalVolume* lv = pv->GetLogicalVolume();

  G4ReflectionFactory* reffact = G4ReflectionFactory::Instance();

  //--- It is not needed to dump _refl volumes created when parent is reflected
  // !!WARNING : it must be avoided to reflect a volume hierarchy if children
  //             has also been reflected, as both will have same name

  if(reffact->IsReflected(lv) && reffact->IsReflected(pv->GetMotherLogical()))
  {
    return;
  }

  G4bool bVolExists = CheckIfLogVolExists(lv->GetName(), lv);

  //---- Construct this PV
  if(pv->GetMotherLogical() != nullptr)  // not WORLD volume
  {
    if(!pv->IsReplicated())
    {
      G4String lvName = lv->GetName();
      if(!bVolExists)
      {
        lvName = DumpLogVol(lv);
      }
      DumpPVPlacement(pv, lvName);
    }
    else if(pv->IsParameterised())
    {
      G4PVParameterised* pvparam = (G4PVParameterised*) (pv);
      DumpPVParameterised(pvparam);
    }
    else
    {
      G4String lvName = lv->GetName();
      if(!bVolExists)
      {
        lvName = DumpLogVol(lv);
      }
      G4PVReplica* pvrepl = (G4PVReplica*) (pv);
      DumpPVReplica(pvrepl, lvName);
    }
  }
  else
  {
    DumpLogVol(lv);
  }

  if(!bVolExists)
  {
    //---- Construct PV's who has this LV as mother
    std::vector<G4VPhysicalVolume*> pvChildren = GetPVChildren(lv);
    for(auto ite = pvChildren.cbegin(); ite != pvChildren.cend(); ++ite)
    {
      DumpPhysVol(*ite);
    }
  }
}

// --------------------------------------------------------------------
void G4tgbGeometryDumper::DumpPVPlacement(G4VPhysicalVolume* pv,
                                          const G4String& lvName, G4int copyNo)
{
  G4String pvName = pv->GetName();

  G4RotationMatrix* rotMat = pv->GetRotation();
  if(rotMat == nullptr)
    rotMat = new G4RotationMatrix();

  //---- Check if it is reflected
  G4ReflectionFactory* reffact = G4ReflectionFactory::Instance();
  G4LogicalVolume* lv          = pv->GetLogicalVolume();
  if(reffact->IsReflected(lv))
  {
#ifdef G4VERBOSE
    if(G4tgrMessenger::GetVerboseLevel() >= 1)
    {
      G4cout << " G4tgbGeometryDumper::DumpPVPlacement() - Reflected volume: "
             << pv->GetName() << G4endl;
    }
#endif
    G4ThreeVector colx = rotMat->colX();
    G4ThreeVector coly = rotMat->colY();
    G4ThreeVector colz = rotMat->colZ();
    // apply a Z reflection (reflection matrix is decomposed in new
    // reflection-free rotation + z-reflection)
    colz *= -1.;
    G4Rep3x3 rottemp(colx.x(), coly.x(), colz.x(), colx.y(), coly.y(), colz.y(),
                     colx.z(), coly.z(), colz.z());
    // matrix representation (inverted)
    *rotMat = G4RotationMatrix(rottemp);
    *rotMat = (*rotMat).inverse();
    pvName += "_refl";
  }
  G4String rotName  = DumpRotationMatrix(rotMat);
  G4ThreeVector pos = pv->GetTranslation();

  if(copyNo == -999)  // for parameterisations copy number is provided
  {
    copyNo = pv->GetCopyNo();
  }

  G4String fullname = pvName + "#" + G4UIcommand::ConvertToString(copyNo) +
                      "/" + pv->GetMotherLogical()->GetName();

  if(!CheckIfPhysVolExists(fullname, pv))
  {
    (*theFile) << ":PLACE " << SubstituteRefl(AddQuotes(lvName)) << " "
               << copyNo << " "
               << SubstituteRefl(AddQuotes(pv->GetMotherLogical()->GetName()))
               << " " << AddQuotes(rotName) << " " << pos.x() << " " << pos.y()
               << " " << pos.z() << G4endl;

    thePhysVols[fullname] = pv;
  }
}

// --------------------------------------------------------------------
void G4tgbGeometryDumper::DumpPVParameterised(G4PVParameterised* pv)
{
  G4String pvName = pv->GetName();

  EAxis axis;
  G4int nReplicas;
  G4double width;
  G4double offset;
  G4bool consuming;
  pv->GetReplicationData(axis, nReplicas, width, offset, consuming);

  G4VPVParameterisation* param = pv->GetParameterisation();

  G4LogicalVolume* lv             = pv->GetLogicalVolume();
  G4VSolid* solid1st              = param->ComputeSolid(0, pv);
  G4Material* mate1st             = param->ComputeMaterial(0, pv);
  std::vector<G4double> params1st = GetSolidParams(solid1st);
  std::vector<G4double> newParams;
  G4VSolid* newSolid = solid1st;
  G4String lvName;

  for(G4int ii = 0; ii < nReplicas; ++ii)
  {
    G4Material* newMate = param->ComputeMaterial(ii, pv);
    if(solid1st->GetEntityType() == "G4Box")
    {
      G4Box* box = (G4Box*) (solid1st);
      param->ComputeDimensions(*box, ii, pv);
      newParams = GetSolidParams(box);
      newSolid  = (G4VSolid*) box;
    }
    else if(solid1st->GetEntityType() == "G4Tubs")
    {
      G4Tubs* tubs = (G4Tubs*) (solid1st);
      param->ComputeDimensions(*tubs, ii, pv);
      newParams = GetSolidParams(tubs);
      newSolid  = (G4VSolid*) tubs;
    }
    else if(solid1st->GetEntityType() == "G4Trd")
    {
      G4Trd* trd = (G4Trd*) (solid1st);
      param->ComputeDimensions(*trd, ii, pv);
      newParams = GetSolidParams(trd);
      newSolid  = (G4VSolid*) trd;
    }
    else if(solid1st->GetEntityType() == "G4Trap")
    {
      G4Trap* trap = (G4Trap*) (solid1st);
      param->ComputeDimensions(*trap, ii, pv);
      newParams = GetSolidParams(trap);
      newSolid  = (G4VSolid*) trap;
    }
    else if(solid1st->GetEntityType() == "G4Cons")
    {
      G4Cons* cons = (G4Cons*) (solid1st);
      param->ComputeDimensions(*cons, ii, pv);
      newParams = GetSolidParams(cons);
      newSolid  = (G4VSolid*) cons;
    }
    else if(solid1st->GetEntityType() == "G4Sphere")
    {
      G4Sphere* sphere = (G4Sphere*) (solid1st);
      param->ComputeDimensions(*sphere, ii, pv);
      newParams = GetSolidParams(sphere);
      newSolid  = (G4VSolid*) sphere;
    }
    else if(solid1st->GetEntityType() == "G4Orb")
    {
      G4Orb* orb = (G4Orb*) (solid1st);
      param->ComputeDimensions(*orb, ii, pv);
      newParams = GetSolidParams(orb);
      newSolid  = (G4VSolid*) orb;
    }
    else if(solid1st->GetEntityType() == "G4Torus")
    {
      G4Torus* torus = (G4Torus*) (solid1st);
      param->ComputeDimensions(*torus, ii, pv);
      newParams = GetSolidParams(torus);
      newSolid  = (G4VSolid*) torus;
    }
    else if(solid1st->GetEntityType() == "G4Para")
    {
      G4Para* para = (G4Para*) (solid1st);
      param->ComputeDimensions(*para, ii, pv);
      newParams = GetSolidParams(para);
      newSolid  = (G4VSolid*) para;
    }
    else if(solid1st->GetEntityType() == "G4Polycone")
    {
      G4Polycone* polycone = (G4Polycone*) (solid1st);
      param->ComputeDimensions(*polycone, ii, pv);
      newParams = GetSolidParams(polycone);
      newSolid  = (G4VSolid*) polycone;
    }
    else if(solid1st->GetEntityType() == "G4Polyhedra")
    {
      G4Polyhedra* polyhedra = (G4Polyhedra*) (solid1st);
      param->ComputeDimensions(*polyhedra, ii, pv);
      newParams = GetSolidParams(polyhedra);
      newSolid  = (G4VSolid*) polyhedra;
    }
    else if(solid1st->GetEntityType() == "G4Hype")
    {
      G4Hype* hype = (G4Hype*) (solid1st);
      param->ComputeDimensions(*hype, ii, pv);
      newParams = GetSolidParams(hype);
      newSolid  = (G4VSolid*) hype;
    }
    if(ii == 0 || mate1st != newMate || params1st[0] != newParams[0])
    {
      G4String extraName = "";
      if(ii != 0)
      {
        extraName = "#" + G4UIcommand::ConvertToString(ii) + "/" +
                    pv->GetMotherLogical()->GetName();
      }
      lvName = DumpLogVol(lv, extraName, newSolid, newMate);
    }

    param->ComputeTransformation(ii, pv);
    DumpPVPlacement(pv, lvName, ii);
  }
}

// --------------------------------------------------------------------
void G4tgbGeometryDumper::DumpPVReplica(G4PVReplica* pv, const G4String& lvName)
{
  EAxis axis;
  G4int nReplicas;
  G4double width;
  G4double offset;
  G4bool consuming;
  pv->GetReplicationData(axis, nReplicas, width, offset, consuming);
  G4String axisName;
  switch(axis)
  {
    case kXAxis:
      axisName = "X";
      break;
    case kYAxis:
      axisName = "Y";
      break;
    case kZAxis:
      axisName = "Z";
      break;
    case kRho:
      axisName = "R";
      break;
    case kPhi:
      axisName = "PHI";
      break;
    case kRadial3D:
    case kUndefined:
      G4String ErrMessage =
        "Unknown axis of replication for volume" + pv->GetName();
      G4Exception("G4tgbGeometryDumper::DumpPVReplica", "Wrong axis ",
                  FatalException, ErrMessage);
      break;
  }

  G4String fullname = lvName + "/" + pv->GetMotherLogical()->GetName();

  if(!CheckIfPhysVolExists(fullname, pv))
  {
    (*theFile) << ":REPL " << SubstituteRefl(AddQuotes(lvName)) << " "
               << SubstituteRefl(AddQuotes(pv->GetMotherLogical()->GetName()))
               << " " << axisName << " " << nReplicas;
    if(axis != kPhi)
    {
      (*theFile) << " " << width << " " << offset << G4endl;
    }
    else
    {
      (*theFile) << " " << width / deg << "*deg"
                 << " " << offset / deg << "*deg" << G4endl;
    }

    thePhysVols[fullname] = pv;
  }
}

// --------------------------------------------------------------------
G4String G4tgbGeometryDumper::DumpLogVol(G4LogicalVolume* lv,
                                         const G4String& extraName,
                                         G4VSolid* solid,
                                         G4Material* mate)
{
  G4String lvName;

  if(extraName == "")  //--- take out the '_refl' in the name
  {
    lvName = GetObjectName(lv, theLogVols);
  }
  else
  {
    lvName = lv->GetName() + extraName;
  }

  if(theLogVols.find(lvName) != theLogVols.cend())  // alredy dumped
  {
    return lvName;
  }

  if(solid == nullptr)
  {
    solid = lv->GetSolid();
  }

  //---- Dump solid
  G4String solidName = DumpSolid(solid, extraName);

  //---- Dump material
  if(mate == nullptr)
  {
    mate = lv->GetMaterial();
  }
  G4String mateName = DumpMaterial(mate);

  //---- Dump logical volume (solid + material)
  (*theFile) << ":VOLU " << SubstituteRefl(AddQuotes(lvName)) << " "
             << SupressRefl(AddQuotes(solidName)) << " " << AddQuotes(mateName)
             << G4endl;

  theLogVols[lvName] = lv;

  return lvName;
}

// --------------------------------------------------------------------
G4String G4tgbGeometryDumper::DumpMaterial(G4Material* mat)
{
  G4String mateName = GetObjectName(mat, theMaterials);
  if(theMaterials.find(mateName) != theMaterials.cend())  // alredy dumped
  {
    return mateName;
  }

  std::size_t numElements = mat->GetNumberOfElements();
  G4double density   = mat->GetDensity() / g * cm3;

  // start tag
  //
  if(numElements == 1)
  {
    (*theFile) << ":MATE " << AddQuotes(mateName) << " " << mat->GetZ() << " "
               << mat->GetA() / (g / mole) << " " << density << G4endl;
  }
  else
  {
    const G4ElementVector* elems = mat->GetElementVector();
    const G4double* fractions    = mat->GetFractionVector();
    for(std::size_t ii = 0; ii < numElements; ++ii)
    {
      DumpElement(const_cast<G4Element*>((*elems)[ii]));
    }

    (*theFile) << ":MIXT " << AddQuotes(mateName) << " " << density << " "
               << numElements << G4endl;
    // close start element tag and get ready to do composit "parts"
    for(std::size_t ii = 0; ii < numElements; ++ii)
    {
      (*theFile) << "   " << AddQuotes(GetObjectName(const_cast<G4Element*>((*elems)[ii]), theElements))
                 << " " << fractions[ii] << G4endl;
    }
  }

  (*theFile) << ":MATE_MEE " << AddQuotes(mateName) << " "
             << mat->GetIonisation()->GetMeanExcitationEnergy() / eV << "*eV"
             << G4endl;

  (*theFile) << ":MATE_TEMPERATURE " << AddQuotes(mateName) << " "
             << mat->GetTemperature() / kelvin << "*kelvin" << G4endl;

  (*theFile) << ":MATE_PRESSURE " << AddQuotes(mateName) << " "
             << mat->GetPressure() / atmosphere << "*atmosphere" << G4endl;

  G4State state = mat->GetState();
  G4String stateStr;
  switch(state)
  {
    case kStateUndefined:
      stateStr = "Undefined";
      break;
    case kStateSolid:
      stateStr = "Solid";
      break;
    case kStateLiquid:
      stateStr = "Liquid";
      break;
    case kStateGas:
      stateStr = "Gas";
      break;
  }

  (*theFile) << ":MATE_STATE " << AddQuotes(mateName) << " " << stateStr
             << G4endl;

  theMaterials[mateName] = mat;

  return mateName;
}

// --------------------------------------------------------------------
void G4tgbGeometryDumper::DumpElement(G4Element* ele)
{
  G4String elemName = GetObjectName(ele, theElements);

  if(theElements.find(elemName) != theElements.cend())  // alredy dumped
  {
    return;
  }

  //--- Add symbol name: Material mixtures store the components as elements
  //    (even if the input are materials), but without symbol
  //
  G4String symbol = ele->GetSymbol();
  if(symbol == "" || symbol == " ")
  {
    symbol = elemName;
  }

  if(ele->GetNumberOfIsotopes() == 0)
  {
    (*theFile) << ":ELEM " << AddQuotes(elemName) << " " << AddQuotes(symbol)
               << " " << ele->GetZ() << " " << ele->GetA() / (g / mole) << " "
               << G4endl;
  }
  else
  {
    const G4IsotopeVector* isots = ele->GetIsotopeVector();
    for(std::size_t ii = 0; ii < ele->GetNumberOfIsotopes(); ++ii)
    {
      DumpIsotope((*isots)[ii]);
    }

    (*theFile) << ":ELEM_FROM_ISOT " << AddQuotes(elemName) << " "
               << AddQuotes(symbol) << " " << ele->GetNumberOfIsotopes()
               << G4endl;
    const G4double* fractions = ele->GetRelativeAbundanceVector();
    for(std::size_t ii = 0; ii < ele->GetNumberOfIsotopes(); ++ii)
    {
      (*theFile) << "   " << AddQuotes(GetObjectName((*isots)[ii], theIsotopes))
                 << " " << fractions[ii] << G4endl;
    }
  }
  theElements[elemName] = ele;
}

// --------------------------------------------------------------------
void G4tgbGeometryDumper::DumpIsotope(G4Isotope* isot)
{
  G4String isotName = GetObjectName(isot, theIsotopes);
  if(theIsotopes.find(isotName) != theIsotopes.cend())  // alredy dumped
  {
    return;
  }

  (*theFile) << ":ISOT " << AddQuotes(isotName) << " " << isot->GetZ() << " "
             << isot->GetN() << " " << isot->GetA() / (g / mole) << " "
             << G4endl;

  theIsotopes[isotName] = isot;
}

// --------------------------------------------------------------------
G4String G4tgbGeometryDumper::DumpSolid(G4VSolid* solid,
                                        const G4String& extraName)
{
  G4String solidName;
  if(extraName == "")
  {
    solidName = GetObjectName(solid, theSolids);
  }
  else
  {
    solidName = solid->GetName() + extraName;
  }

  if(theSolids.find(solidName) != theSolids.cend())  // alredy dumped
  {
    return solidName;
  }

  G4String solidType = solid->GetEntityType();
  solidType          = GetTGSolidType(solidType);
  
  if(solidType == "UNIONSOLID")
  {
    DumpBooleanVolume("UNION", solid);
  }
  else if(solidType == "SUBTRACTIONSOLID")
  {
    DumpBooleanVolume("SUBTRACTION", solid);
  }
  else if(solidType == "INTERSECTIONSOLID")
  {
    DumpBooleanVolume("INTERSECTION", solid);
  }
  else if(solidType == "REFLECTEDSOLID")
  {
    G4ReflectedSolid* solidrefl = dynamic_cast<G4ReflectedSolid*>(solid);
    if(solidrefl == nullptr)
    {
      G4Exception("G4tgbGeometryDumper::DumpSolid()", "InvalidType",
                  FatalException, "Invalid reflected solid!");
      return solidName;
    }
    G4VSolid* solidori = solidrefl->GetConstituentMovedSolid();
    DumpSolid(solidori);
  }
  else if(solidType == "MULTIUNION")
  {
    DumpMultiUnionVolume(solid);
  }
  else if(solidType == "SCALEDSOLID")
  {
    DumpScaledVolume(solid);
  }
  else
  {
    (*theFile) << ":SOLID " << AddQuotes(solidName) << " ";
    (*theFile) << AddQuotes(solidType) << " ";
    DumpSolidParams( solid );
    theSolids[solidName] = solid;
  }

  return solidName;
}

// --------------------------------------------------------------------
void G4tgbGeometryDumper::DumpBooleanVolume(const G4String& solidType,
                                            G4VSolid* so)
{
  G4BooleanSolid* bso = dynamic_cast<G4BooleanSolid*>(so);
  if(bso == nullptr)
  {
    return;
  }
  G4VSolid* solid0             = bso->GetConstituentSolid(0);
  G4VSolid* solid1             = bso->GetConstituentSolid(1);
  G4DisplacedSolid* solid1Disp = nullptr;
  G4bool displaced             = dynamic_cast<G4DisplacedSolid*>(solid1);
  if(displaced)
  {
    solid1Disp = dynamic_cast<G4DisplacedSolid*>(solid1);
    if(solid1Disp != nullptr)
    {
      solid1 = solid1Disp->GetConstituentMovedSolid();
    }
    else
    {
      return;
    }
  }
  DumpSolid(solid0);
  DumpSolid(solid1);

  G4String rotName;
  G4ThreeVector pos;
  if(displaced)
  {
    pos = solid1Disp->GetObjectTranslation();  // translation is of mother frame
    rotName = DumpRotationMatrix(new G4RotationMatrix(
      (solid1Disp->GetTransform().NetRotation()).inverse()));
  }
  else  // no displacement
  {
    rotName = DumpRotationMatrix(new G4RotationMatrix);
    pos     = G4ThreeVector();
  }

  G4String bsoName = GetObjectName(so, theSolids);
  if(theSolids.find(bsoName) != theSolids.cend())
    return;  // alredy dumped
  G4String solid0Name = FindSolidName(solid0);
  G4String solid1Name = FindSolidName(solid1);

  (*theFile) << ":SOLID " << AddQuotes(bsoName) << " " << AddQuotes(solidType)
             << " " << AddQuotes(solid0Name) << " " << AddQuotes(solid1Name)
             << " " << AddQuotes(rotName) << " " << approxTo0(pos.x()) << " "
             << approxTo0(pos.y()) << " " << approxTo0(pos.z()) << " "
             << G4endl;

  theSolids[bsoName] = bso;
}

// --------------------------------------------------------------------
void G4tgbGeometryDumper::DumpMultiUnionVolume( G4VSolid* so)
{
  const G4MultiUnion* muun = dynamic_cast<const G4MultiUnion*>(so);
  if(muun != nullptr)
    {
      G4int nSolids = muun->GetNumberOfSolids();
      std::vector<G4String> rotList;
      for( G4int iso = 0; iso < nSolids; iso++ ) {
	G4Transform3D trans = muun->GetTransformation(iso);
	G4String rotName = DumpRotationMatrix( new G4RotationMatrix(trans.getRotation()));
	rotList.push_back(rotName);
	G4VSolid* solN = muun->GetSolid(iso);
	DumpSolid(solN);
      }
      G4String bsoName = GetObjectName(const_cast<G4VSolid*>(so), theSolids);
      (*theFile) << ":SOLID " << AddQuotes(bsoName) << " MULTIUNION "
		 << nSolids;
      
      for( G4int iso = 0; iso < nSolids; iso++ ) {
	G4VSolid* solN = muun->GetSolid(iso);
	G4Transform3D trans = muun->GetTransformation(iso);
	G4ThreeVector pos = trans.getTranslation();  // translation is of mother frame
	(*theFile) << " " <<  solN->GetName()
		   << " " << " " << rotList[iso]
		   << " " << approxTo0(pos.x())
		   << " " << approxTo0(pos.y())
		   << " " << approxTo0(pos.z());
      }
      (*theFile) << G4endl;
      
    }
}

// --------------------------------------------------------------------
void G4tgbGeometryDumper::DumpScaledVolume( G4VSolid* so)
{
  const G4ScaledSolid* ssol = dynamic_cast<const G4ScaledSolid*>(so);
  if(ssol != nullptr)
    {
      G4VSolid* unscaledSolid = ssol->GetUnscaledSolid();
      G4Scale3D scaleTransf = ssol->GetScaleTransform();
      G4String bsoName = GetObjectName(const_cast<G4VSolid*>(so), theSolids);
      (*theFile) << ":SOLID " << AddQuotes(bsoName) << " SCALED "
		 << unscaledSolid->GetName() << " "
		 << scaleTransf.xx() << " "
		 << scaleTransf.yy() << " "
		 << scaleTransf.zz() << G4endl;
    }
}

// --------------------------------------------------------------------
void G4tgbGeometryDumper::DumpSolidParams(G4VSolid* so)
{
  std::vector<G4double> params = GetSolidParams(so);
  for(std::size_t ii = 0; ii < params.size(); ++ii)
  {
    (*theFile) << params[ii] << " ";
  }
  (*theFile) << G4endl;
}

// --------------------------------------------------------------------
std::vector<G4double> G4tgbGeometryDumper::GetSolidParams(const G4VSolid* so)
{
  std::vector<G4double> params;

  G4String solidType = so->GetEntityType();
  solidType          = GetTGSolidType(solidType);

  if(solidType == "BOX")
  {
    const G4Box* sb = dynamic_cast<const G4Box*>(so);
    if(sb != nullptr)
    {
      params.push_back(sb->GetXHalfLength());
      params.push_back(sb->GetYHalfLength());
      params.push_back(sb->GetZHalfLength());
    }
  }
  else if(solidType == "TUBS")
  {
    const G4Tubs* tu = dynamic_cast<const G4Tubs*>(so);
    if(tu != nullptr)
    {
      params.push_back(tu->GetInnerRadius());
      params.push_back(tu->GetOuterRadius());
      params.push_back(tu->GetZHalfLength());
      params.push_back(tu->GetStartPhiAngle() / deg);
      params.push_back(tu->GetDeltaPhiAngle() / deg);
    }
  }
  else if(solidType == "TRAP")
  {
    const G4Trap* trp = dynamic_cast<const G4Trap*>(so);
    if(trp != nullptr)
    {
      G4ThreeVector symAxis(trp->GetSymAxis());
      params.push_back(trp->GetZHalfLength());
      params.push_back(symAxis.theta() / deg);
      params.push_back(symAxis.phi() / deg);
      params.push_back(trp->GetYHalfLength1());
      params.push_back(trp->GetXHalfLength1());
      params.push_back(trp->GetXHalfLength2());
      params.push_back(std::atan(trp->GetTanAlpha1()) / deg);
      params.push_back(trp->GetYHalfLength2());
      params.push_back(trp->GetXHalfLength3());
      params.push_back(trp->GetXHalfLength4());
      params.push_back(std::atan(trp->GetTanAlpha2()) / deg);
    }
  }
  else if(solidType == "TRD")
  {
    const G4Trd* tr = dynamic_cast<const G4Trd*>(so);
    if(tr != nullptr)
    {
      params.push_back(tr->GetXHalfLength1());
      params.push_back(tr->GetXHalfLength2());
      params.push_back(tr->GetYHalfLength1());
      params.push_back(tr->GetYHalfLength2());
      params.push_back(tr->GetZHalfLength());
    }
  }
  else if(solidType == "PARA")
  {
    const G4Para* para = dynamic_cast<const G4Para*>(so);
    if(para != nullptr)
    {
      G4ThreeVector symAxis(para->GetSymAxis());
      params.push_back(para->GetXHalfLength());
      params.push_back(para->GetYHalfLength());
      params.push_back(para->GetZHalfLength());
      params.push_back(std::atan(para->GetTanAlpha()) / deg);
      params.push_back(symAxis.theta() / deg);
      params.push_back(symAxis.phi() / deg);
    }
  }
  else if(solidType == "CONS")
  {
    const G4Cons* cn = dynamic_cast<const G4Cons*>(so);
    if(cn != nullptr)
    {
      params.push_back(cn->GetInnerRadiusMinusZ());
      params.push_back(cn->GetOuterRadiusMinusZ());
      params.push_back(cn->GetInnerRadiusPlusZ());
      params.push_back(cn->GetOuterRadiusPlusZ());
      params.push_back(cn->GetZHalfLength());
      params.push_back(cn->GetStartPhiAngle() / deg);
      params.push_back(cn->GetDeltaPhiAngle() / deg);
    }
  }
  else if(solidType == "SPHERE")
  {
    const G4Sphere* sphere = dynamic_cast<const G4Sphere*>(so);
    if(sphere != nullptr)
    {
      params.push_back(sphere->GetInnerRadius());
      params.push_back(sphere->GetOuterRadius());
      params.push_back(sphere->GetStartPhiAngle() / deg);
      params.push_back(sphere->GetDeltaPhiAngle() / deg);
      params.push_back(sphere->GetStartThetaAngle() / deg);
      params.push_back(sphere->GetDeltaThetaAngle() / deg);
    }
  }
  else if(solidType == "ORB")
  {
    const G4Orb* orb = dynamic_cast<const G4Orb*>(so);
    if(orb != nullptr)
    {
      params.push_back(orb->GetRadius());
    }
  }
  else if(solidType == "TORUS")
  {
    const G4Torus* torus = dynamic_cast<const G4Torus*>(so);
    if(torus != nullptr)
    {
      params.push_back(torus->GetRmin());
      params.push_back(torus->GetRmax());
      params.push_back(torus->GetRtor());
      params.push_back(torus->GetSPhi() / deg);
      params.push_back(torus->GetDPhi() / deg);
    }
  }
  else if(solidType == "POLYCONE")
  {
    //--- Dump RZ corners, as original parameters will not be present
    //    if it was build from RZ corners
    const G4Polycone* plc = dynamic_cast<const G4Polycone*>(so);
    if(plc != nullptr)
    {
      G4double angphi = plc->GetStartPhi() / deg;
      if(angphi > 180 * deg)
      {
        angphi -= 360 * deg;
      }
      G4double endphi = plc->GetEndPhi() / deg;
      if(endphi > 180 * deg)
      {
        endphi -= 360 * deg;
      }
      params.push_back(angphi);
      params.push_back(endphi - angphi);
      //      params.push_back(plc->GetOriginalParameters()->Opening_angle / deg);
      G4int ncor = plc->GetNumRZCorner();
      params.push_back(ncor);

      for(G4int ii = 0; ii < ncor; ++ii)
      {
        params.push_back(plc->GetCorner(ii).r);
        params.push_back(plc->GetCorner(ii).z);
      }
    }
  }
  else if(solidType == "GENERICPOLYCONE")
  {
    //--- Dump RZ corners
    const G4GenericPolycone* plc = dynamic_cast<const G4GenericPolycone*>(so);
    if(plc != nullptr)
    {
      G4double angphi = plc->GetStartPhi() / deg;
      if(angphi > 180 * deg)
      {
        angphi -= 360 * deg;
      }
      G4double endphi = plc->GetEndPhi() / deg;
      if(endphi > 180 * deg)
      {
        endphi -= 360 * deg;
      }
      params.push_back(angphi);
      params.push_back(endphi - angphi);
      G4int ncor = plc->GetNumRZCorner();
      params.push_back(ncor);

      for(G4int ii = 0; ii < ncor; ++ii)
      {
        params.push_back(plc->GetCorner(ii).r);
        params.push_back(plc->GetCorner(ii).z);
      }
    }
  }
  else if(solidType == "POLYHEDRA")
  {
    //--- Dump RZ corners, as original parameters will not be present
    //    if it was build from RZ corners
    const G4Polyhedra* ph = (dynamic_cast<const G4Polyhedra*>(so));
    if(ph != nullptr)
    {
      G4double angphi = ph->GetStartPhi() / deg;
      if(angphi > 180 * deg)
        angphi -= 360 * deg;

      G4int ncor = ph->GetNumRZCorner();

      params.push_back(angphi);
      params.push_back(ph->GetOriginalParameters()->Opening_angle / deg);
      params.push_back(ph->GetNumSide());
      params.push_back(ncor);

      for(G4int ii = 0; ii < ncor; ++ii)
      {
        params.push_back(ph->GetCorner(ii).r);
        params.push_back(ph->GetCorner(ii).z);
      }
    }
  }
  else if(solidType == "ELLIPTICALTUBE")
  {
    const G4EllipticalTube* eltu = dynamic_cast<const G4EllipticalTube*>(so);
    if(eltu != nullptr)
    {
      params.push_back(eltu->GetDx());
      params.push_back(eltu->GetDy());
      params.push_back(eltu->GetDz());
    }
  }
  else if(solidType == "ELLIPSOID")
  {
    const G4Ellipsoid* dso = dynamic_cast<const G4Ellipsoid*>(so);
    if(dso != nullptr)
    {
      params.push_back(dso->GetSemiAxisMax(0));
      params.push_back(dso->GetSemiAxisMax(1));
      params.push_back(dso->GetSemiAxisMax(2));
      params.push_back(dso->GetZBottomCut());
      params.push_back(dso->GetZTopCut());
    }
  }
  else if(solidType == "ELLIPTICAL_CONE")
  {
    const G4EllipticalCone* elco = dynamic_cast<const G4EllipticalCone*>(so);
    if(elco != nullptr)
    {
      params.push_back(elco->GetSemiAxisX());
      params.push_back(elco->GetSemiAxisY());
      params.push_back(elco->GetZMax());
      params.push_back(elco->GetZTopCut());
    }
  }
  else if(solidType == "HYPE")
  {
    const G4Hype* hype = dynamic_cast<const G4Hype*>(so);
    if(hype != nullptr)
    {
      params.push_back(hype->GetInnerRadius());
      params.push_back(hype->GetOuterRadius());
      params.push_back(hype->GetInnerStereo() / deg);
      params.push_back(hype->GetOuterStereo() / deg);
      params.push_back(2 * hype->GetZHalfLength());
    }
    //  } else if( solidType == "TET" ) {
  }
  else if(solidType == "TWISTEDBOX")
  {
    const G4TwistedBox* tbox = dynamic_cast<const G4TwistedBox*>(so);
    if(tbox != nullptr)
    {
      params.push_back(tbox->GetPhiTwist() / deg);
      params.push_back(tbox->GetXHalfLength());
      params.push_back(tbox->GetYHalfLength());
      params.push_back(tbox->GetZHalfLength());
    }
  }
  else if(solidType == "TWISTEDTRAP")
  {
    const G4TwistedTrap* ttrap = dynamic_cast<const G4TwistedTrap*>(so);
    if(ttrap != nullptr)
    {
      params.push_back(ttrap->GetPhiTwist() / deg);
      params.push_back(ttrap->GetZHalfLength());
      params.push_back(ttrap->GetPolarAngleTheta() / deg);
      params.push_back(ttrap->GetAzimuthalAnglePhi() / deg);
      params.push_back(ttrap->GetY1HalfLength());
      params.push_back(ttrap->GetX1HalfLength());
      params.push_back(ttrap->GetX2HalfLength());
      params.push_back(ttrap->GetY2HalfLength());
      params.push_back(ttrap->GetX3HalfLength());
      params.push_back(ttrap->GetX4HalfLength());
      params.push_back(ttrap->GetTiltAngleAlpha() / deg);
    }
  }
  else if(solidType == "TWISTEDTRD")
  {
    const G4TwistedTrd* ttrd = dynamic_cast<const G4TwistedTrd*>(so);
    if(ttrd != nullptr)
    {
      params.push_back(ttrd->GetX1HalfLength());
      params.push_back(ttrd->GetX2HalfLength());
      params.push_back(ttrd->GetY1HalfLength());
      params.push_back(ttrd->GetY2HalfLength());
      params.push_back(ttrd->GetZHalfLength());
      params.push_back(ttrd->GetPhiTwist() / deg);
    }
  }
  else if(solidType == "TWISTEDTUBS")
  {
    const G4TwistedTubs* ttub = dynamic_cast<const G4TwistedTubs*>(so);
    if(ttub != nullptr)
    {
      params.push_back(ttub->GetInnerRadius());
      params.push_back(ttub->GetOuterRadius());
      params.push_back(ttub->GetZHalfLength());
      params.push_back(ttub->GetDPhi() / deg);
      params.push_back(ttub->GetPhiTwist() / deg);
    }
  }
  else
  {
    G4String ErrMessage = "Solid type not supported, sorry... " + solidType;
    G4Exception("G4tgbGeometryDumper::DumpSolidParams()", "NotImplemented",
                FatalException, ErrMessage);
  }

  return params;
}

// --------------------------------------------------------------------
G4String G4tgbGeometryDumper::DumpRotationMatrix(G4RotationMatrix* rotm)
{
  if(rotm == nullptr)
  {
    rotm = new G4RotationMatrix();
  }

  G4double de      = MatDeterminant(rotm);
  G4String rotName = LookForExistingRotation(rotm);
  if(rotName != "")
  {
    return rotName;
  }

  G4ThreeVector v(1., 1., 1.);
  if(de < -0.9)  // a reflection ....
  {
    (*theFile) << ":ROTM ";
    rotName = "RRM";
    rotName += G4UIcommand::ConvertToString(theRotationNumber++);

    (*theFile) << AddQuotes(rotName) << std::setprecision(9) << " "
               << approxTo0(rotm->xx()) << " " << approxTo0(rotm->yx()) << " "
               << approxTo0(rotm->zx()) << " " << approxTo0(rotm->xy()) << " "
               << approxTo0(rotm->yy()) << " " << approxTo0(rotm->zy()) << " "
               << approxTo0(rotm->xz()) << " " << approxTo0(rotm->yz()) << " "
               << approxTo0(rotm->zz()) << G4endl;
  }
  else if(de > 0.9)  // a rotation ....
  {
    (*theFile) << ":ROTM ";
    rotName = "RM";
    rotName += G4UIcommand::ConvertToString(theRotationNumber++);

    (*theFile) << AddQuotes(rotName) << " " << approxTo0(rotm->thetaX() / deg)
               << " " << approxTo0(rotm->phiX() / deg) << " "
               << approxTo0(rotm->thetaY() / deg) << " "
               << approxTo0(rotm->phiY() / deg) << " "
               << approxTo0(rotm->thetaZ() / deg) << " "
               << approxTo0(rotm->phiZ() / deg) << G4endl;
  }

  theRotMats[rotName] = rotm;

  return rotName;
}

// --------------------------------------------------------------------
std::vector<G4VPhysicalVolume*>
G4tgbGeometryDumper::GetPVChildren(G4LogicalVolume* lv)
{
  G4PhysicalVolumeStore* pvstore = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4VPhysicalVolume*> children;
  for(auto ite = pvstore->cbegin(); ite != pvstore->cend(); ++ite)
  {
    if((*ite)->GetMotherLogical() == lv)
    {
      children.push_back(*ite);
#ifdef G4VERBOSE
      if(G4tgrMessenger::GetVerboseLevel() >= 1)
      {
        G4cout << " G4tgbGeometryDumper::GetPVChildren() - adding children: "
               << (*ite)->GetName() << " of " << lv->GetName() << G4endl;
      }
#endif
    }
  }

  return children;
}

// --------------------------------------------------------------------
G4String G4tgbGeometryDumper::GetTGSolidType(const G4String& solidType)
{
  G4String newsolidType = solidType.substr(2, solidType.length());
  for(G4int ii = 0; ii < (G4int)newsolidType.length(); ++ii)
  {
    newsolidType[ii] = (char)std::toupper(newsolidType[ii]);
  }
  return newsolidType;
}

// --------------------------------------------------------------------
G4double G4tgbGeometryDumper::MatDeterminant(G4RotationMatrix* ro)
{
  G4Rep3x3 r = ro->rep3x3();
  return r.xx_ * (r.yy_ * r.zz_ - r.zy_ * r.yz_) -
         r.yx_ * (r.xy_ * r.zz_ - r.zy_ * r.xz_) +
         r.zx_ * (r.xy_ * r.yz_ - r.yy_ * r.xz_);
}

// --------------------------------------------------------------------
G4double G4tgbGeometryDumper::approxTo0(G4double val)
{
  G4double precision =
    G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  if(std::fabs(val) < precision)
  {
    val = 0.0;
  }
  return val;
}

// --------------------------------------------------------------------
G4String G4tgbGeometryDumper::AddQuotes(const G4String& str)
{
  //--- look if there is a separating blank

  G4bool bBlank = FALSE;
  std::size_t siz = str.length();
  for(std::size_t ii = 0; ii < siz; ++ii)
  {
    if(str.substr(ii, 1) == " ")
    {
      bBlank = TRUE;
      break;
    }
  }
  G4String str2 = str;
  if(bBlank)
  {
    str2 = G4String("\"") + str2 + G4String("\"");
  }
  return str2;
}

// --------------------------------------------------------------------
G4String G4tgbGeometryDumper::SupressRefl(G4String name)
{
  std::size_t irefl = name.rfind("_refl");
  if(irefl != G4String::npos)
  {
    name = name.substr(0, irefl);
  }
  return name;
}

// --------------------------------------------------------------------
G4String G4tgbGeometryDumper::SubstituteRefl(G4String name)
{
  std::size_t irefl = name.rfind("_refl");
  if(irefl != G4String::npos)
  {
    name = name.substr(0, irefl) + "_REFL";
  }
  return name;
}

// --------------------------------------------------------------------
G4String G4tgbGeometryDumper::GetIsotopeName(G4Isotope* isot)
{
  G4String isotName = isot->GetName();
  // first look if this is isotope is already dumped,
  // with original isotope name or new one
  //
  std::map<G4String, G4Isotope*>::const_iterator ite;
  for(ite = theIsotopes.cbegin(); ite != theIsotopes.cend(); ++ite)
  {
    if(isot == (*ite).second)
    {
      return (*ite).first;
    }
  }

  // Now look if there is another isotope dumped with same name,
  // and if found add _N to the name
  //
  ite = theIsotopes.find(isotName);
  if(ite != theIsotopes.cend())  // Isotope found with same name
  {
    G4Isotope* isotold = (*ite).second;
    if(isot != isotold)  // new isotope it is not the really
    {                    // the same one as isotope found
      if(!Same2G4Isotopes(isot, isotold))
      {                // if the two have same data, use the old one
        G4int ii = 2;  // G4Nist does names isotopes of same element
                       // with same name
        for(;; ++ii)
        {
          G4String newIsotName =
            isotName + "_" + G4UIcommand::ConvertToString(ii);
          std::map<G4String, G4Isotope*>::const_iterator ite2 =
            theIsotopes.find(newIsotName);
          if(ite2 == theIsotopes.cend())
          {
            isotName = newIsotName;
            break;
          }
          else
          {
            if(Same2G4Isotopes(isot, (*ite2).second))
            {
              isotName = newIsotName;
              break;
            }
          }
        }
      }
    }
  }
  return isotName;
}

// --------------------------------------------------------------------
template <class TYP>
G4String G4tgbGeometryDumper::GetObjectName(
  TYP* obj, std::map<G4String, TYP*> objectsDumped)
{
  G4String objName = obj->GetName();

  // first look if this is objecy is already dumped,
  // with original object name or new one
  //
  typename std::map<G4String, TYP*>::const_iterator ite;
  for(ite = objectsDumped.cbegin(); ite != objectsDumped.cend(); ++ite)
  {
    if(obj == (*ite).second)
    {
      return (*ite).first;
    }
  }

  // Now look if there is another object dumped with same name,
  // and if found add _N to the name
  //
  ite = objectsDumped.find(objName);

  if(ite != objectsDumped.cend())  // Object found with same name
  {
    TYP* objold = (*ite).second;
    if(obj != objold)  // new object it is not the really
    {                  // the same one as object found
      G4int ii = 2;
      for(;; ++ii)
      {
        G4String newObjName = objName + "_" + G4UIcommand::ConvertToString(ii);
        typename std::map<G4String, TYP*>::const_iterator ite2 =
          objectsDumped.find(newObjName);
        if(ite2 == objectsDumped.cend())
        {
          objName = newObjName;
          break;
        }
      }
    }
  }
  return objName;
}

// --------------------------------------------------------------------
G4bool G4tgbGeometryDumper::CheckIfLogVolExists(const G4String& name,
                                                G4LogicalVolume* pt)
{
  if(theLogVols.find(name) != theLogVols.cend())
  {
    G4LogicalVolume* lvnew = (*(theLogVols.find(name))).second;
    if(lvnew != pt)
    {
      /*
      //---- Reflected volumes are repeated

      G4ReflectionFactory* reffact = G4ReflectionFactory::Instance();
      if( !reffact->IsReflected( pt ) && !reffact->IsReflected( lvnew ) )
      {
        G4String ErrMessage = "LogVol found but not same as before: " + name;
        G4Exception("G4tgbGeometryDumper::CheckIfLogVolExists()",
                    "InvalidSetup", FatalException, ErrMessage);
      }
      */
    }
    return true;
  }
  else
  {
    return false;
  }
}

// --------------------------------------------------------------------
G4bool G4tgbGeometryDumper::CheckIfPhysVolExists(const G4String& name,
                                                 G4VPhysicalVolume* pt)
{
#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 1)
  {
    G4cout << " G4tgbGeometryDumper::CheckIfPhysVolExists() - " << name
           << G4endl;
  }
#endif
  if(thePhysVols.find(name) != thePhysVols.cend())
  {
    if((*(thePhysVols.find(name))).second != pt)
    {
      // G4String ErrMessage = "Placement found but not same as before: "
      //                     + name;
      // G4Exception("G4tgbGeometryDumper::CheckIfPhysVolExists()",
      //             "InvalidSetup", FatalException, ErrMessage);
      G4cerr << " G4tgbGeometryDumper::CheckIfPhysVolExists () -"
             << " Placement found but not same as before : " << name << G4endl;
    }
    return true;
  }
  else
  {
    return false;
  }
}

// --------------------------------------------------------------------
G4String
G4tgbGeometryDumper::LookForExistingRotation(const G4RotationMatrix* rotm)
{
  G4String rmName = "";

  for(auto ite = theRotMats.cbegin(); ite != theRotMats.cend(); ++ite)
  {
    if((*ite).second->isNear(*rotm))
    {
      rmName = (*ite).first;
      break;
    }
  }
  return rmName;
}

// --------------------------------------------------------------------
G4bool G4tgbGeometryDumper::Same2G4Isotopes(G4Isotope* isot1, G4Isotope* isot2)
{
  if((isot1->GetZ() != isot2->GetZ()) || (isot1->GetN() != isot2->GetN()) ||
     (isot1->GetA() != isot2->GetA()))
  {
    return false;
  }
  else
  {
    return true;
  }
}

// --------------------------------------------------------------------
const G4String& G4tgbGeometryDumper::FindSolidName(G4VSolid* solid)
{
  std::map<G4String, G4VSolid*>::const_iterator ite;
  for(ite = theSolids.cbegin(); ite != theSolids.cend(); ++ite)
  {
    if(solid == (*ite).second)
    {
      return (*ite).first;
    }
  }

  if(ite == theSolids.cend())
  {
    G4Exception("G4tgbGeometryDumper::FindSolidName()", "ReadError",
                FatalException, "Programming error.");
  }
  return (*ite).first;
}
