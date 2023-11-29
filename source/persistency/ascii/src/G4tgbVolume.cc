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
// G4tgbVolume implementation
//
// Author: P.Arce, CIEMAT (November 2007)
// --------------------------------------------------------------------

#include "G4tgbVolume.hh"

#include "G4tgbVolumeMgr.hh"
#include "G4tgbMaterialMgr.hh"
#include "G4tgbRotationMatrixMgr.hh"
#include "G4tgbPlaceParamLinear.hh"
#include "G4tgbPlaceParamSquare.hh"
#include "G4tgbPlaceParamCircle.hh"

#include "G4tgrSolid.hh"
#include "G4tgrSolidBoolean.hh"
#include "G4tgrSolidMultiUnion.hh"
#include "G4tgrSolidScaled.hh"
#include "G4tgrVolume.hh"
#include "G4tgrVolumeDivision.hh"
#include "G4tgrVolumeAssembly.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgrPlace.hh"
#include "G4tgrPlaceSimple.hh"
#include "G4tgrPlaceDivRep.hh"
#include "G4tgrPlaceParameterisation.hh"
#include "G4tgrUtils.hh"

#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4MultiUnion.hh"
#include "G4ScaledSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVDivision.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
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
#include "G4AssemblyVolume.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"
#include "G4ExtrudedSolid.hh"

#include "G4VisExtent.hh"
#include "G4Material.hh"
#include "G4RotationMatrix.hh"
#include "G4ReflectionFactory.hh"

#include "G4VisAttributes.hh"
#include "G4RegionStore.hh"
#include "G4tgrMessenger.hh"
#include "G4UIcommand.hh"
#include "G4GeometryTolerance.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// --------------------------------------------------------------------
G4tgbVolume::G4tgbVolume()
{
}

// --------------------------------------------------------------------
G4tgbVolume::~G4tgbVolume()
{
}

// --------------------------------------------------------------------
G4tgbVolume::G4tgbVolume(G4tgrVolume* vol)
{
  theTgrVolume = vol;
}

// --------------------------------------------------------------------
void G4tgbVolume::ConstructG4Volumes(const G4tgrPlace* place,
                                     const G4LogicalVolume* parentLV)
{
#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 2)
  {
    G4cout << G4endl << "@@@ G4tgbVolume::ConstructG4Volumes - " << GetName()
           << G4endl;
    if(place && parentLV)
      G4cout << "   place in LV " << parentLV->GetName() << G4endl;
  }
#endif
  G4tgbVolumeMgr* g4vmgr     = G4tgbVolumeMgr::GetInstance();
  G4LogicalVolume* logvol    = g4vmgr->FindG4LogVol(GetName());
  G4bool bFirstCopy          = false;
  G4VPhysicalVolume* physvol = nullptr;
  if(logvol == nullptr)
  {
    bFirstCopy = true;
    if(theTgrVolume->GetType() != "VOLDivision")
    {
      //--- If first time build solid and LogVol
      G4VSolid* solid = FindOrConstructG4Solid(theTgrVolume->GetSolid());
      if(solid != nullptr)  // for G4AssemblyVolume it is nullptr
      {
        g4vmgr->RegisterMe(solid);
        logvol = ConstructG4LogVol(solid);
        g4vmgr->RegisterMe(logvol);
        g4vmgr->RegisterChildParentLVs(logvol, parentLV);
      }
    }
    else
    {
      return;
    }
  }
  //--- Construct PhysVol
  physvol = ConstructG4PhysVol(place, logvol, parentLV);

  if(physvol != nullptr)  // nullptr for G4AssemblyVolumes
  {
    g4vmgr->RegisterMe(physvol);

    if(logvol == nullptr)  // case of divisions
    {
      logvol = physvol->GetLogicalVolume();
    }
  }
  else
  {
    return;
  }

  //--- If first copy build children placements in this LogVol
  if(bFirstCopy)
  {
    std::pair<G4mmapspl::iterator, G4mmapspl::iterator> children =
      G4tgrVolumeMgr::GetInstance()->GetChildren(GetName());
    for(auto cite = children.first; cite != children.second; ++cite)
    {
      //----- Call G4tgrPlace ->constructG4Volumes
      //---- find G4tgbVolume corresponding to the G4tgrVolume
      //     pointed by G4tgrPlace
      G4tgrPlace* pl    = const_cast<G4tgrPlace*>((*cite).second);
      G4tgbVolume* svol = g4vmgr->FindVolume(pl->GetVolume()->GetName());
      //--- find copyNo
#ifdef G4VERBOSE
      if(G4tgrMessenger::GetVerboseLevel() >= 2)
      {
        G4cout << " G4tgbVolume::ConstructG4Volumes - construct daughter "
               << pl->GetVolume()->GetName() << " # " << pl->GetCopyNo()
               << G4endl;
      }
#endif
      svol->ConstructG4Volumes(pl, logvol);
    }
  }
}

// --------------------------------------------------------------------
G4VSolid* G4tgbVolume::FindOrConstructG4Solid(const G4tgrSolid* sol)
{
  G4double angularTolerance =
    G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  if(sol == nullptr)
  {
    return nullptr;
  }

#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 2)
  {
    G4cout << " G4tgbVolume::FindOrConstructG4Solid():" << G4endl
           << "   SOLID = " << sol << G4endl << "   " << sol->GetName()
           << " of type " << sol->GetType() << G4endl;
  }
#endif

  //----- Check if solid exists already
  G4VSolid* solid = G4tgbVolumeMgr::GetInstance()->FindG4Solid(sol->GetName());
  if(solid != nullptr)
  {
    return solid;
  }

  // Give 'sol' as Boolean solids needs to call this method twice

#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 2)
  {
    G4cout << " G4tgbVolume::FindOrConstructG4Solid() - "
           << sol->GetSolidParams().size() << G4endl;
  }
#endif

  std::vector<G4double> solParam;

  // In case of BOOLEAN solids, solidParams are taken from components

  if(sol->GetSolidParams().size() == 1)
  {
    solParam = *sol->GetSolidParams()[0];
  }

  //----------- instantiate the appropiate G4VSolid type
  G4String stype = sol->GetType();
  G4String sname = sol->GetName();

  if(stype == "BOX")
  {
    CheckNoSolidParams(stype, 3, (G4int)solParam.size());
    solid = new G4Box(sname, solParam[0], solParam[1], solParam[2]);
  }
  else if(stype == "TUBE")
  {
    CheckNoSolidParams(stype, 3, (G4int)solParam.size());
    solid = new G4Tubs(sname, solParam[0], solParam[1], solParam[2], 0. * deg,
                       360. * deg);
  }
  else if(stype == "TUBS")
  {
    CheckNoSolidParams(stype, 5, (G4int)solParam.size());
    G4double phiDelta = solParam[4];
    if(std::fabs(phiDelta - twopi) < angularTolerance)
    {
      phiDelta = twopi;
    }
    solid = new G4Tubs(sname, solParam[0], solParam[1], solParam[2],
                       solParam[3], phiDelta);
  }
  else if(stype == "TRAP")
  {
    if(solParam.size() == 11)
    {
      solid = new G4Trap(sname, solParam[0], solParam[1], solParam[2],
                         solParam[3], solParam[4], solParam[5], solParam[6],
                         solParam[7], solParam[8], solParam[9], solParam[10]);
    }
    else if(solParam.size() == 4)
    {
      solid = new G4Trap(sname, solParam[0], solParam[1] / deg,
                         solParam[2] / deg, solParam[3]);
    }
    else
    {
      G4String ErrMessage1 = "Solid type " + stype;
      G4String ErrMessage2 = " should have 11 or 4 parameters,\n";
      G4String ErrMessage3 =
        "and it has " + G4UIcommand::ConvertToString(G4int(solParam.size()));
      G4String ErrMessage = ErrMessage1 + ErrMessage2 + ErrMessage3 + " !";
      G4Exception("G4tgbVolume::FindOrConstructG4Solid()", "InvalidSetup",
                  FatalException, ErrMessage);
      return 0;
    }
  }
  else if(stype == "TRD")
  {
    CheckNoSolidParams(stype, 5, (G4int)solParam.size());
    solid = new G4Trd(sname, solParam[0], solParam[1], solParam[2], solParam[3],
                      solParam[4]);
  }
  else if(stype == "PARA")
  {
    CheckNoSolidParams(stype, 6, (G4int)solParam.size());
    solid = new G4Para(sname, solParam[0], solParam[1], solParam[2],
                       solParam[3], solParam[4], solParam[5]);
  }
  else if(stype == "CONE")
  {
    CheckNoSolidParams(stype, 5, (G4int)solParam.size());
    solid = new G4Cons(sname, solParam[0], solParam[1], solParam[2],
                       solParam[3], solParam[4], 0., 360. * deg);
  }
  else if(stype == "CONS")
  {
    CheckNoSolidParams(stype, 7, (G4int)solParam.size());
    G4double phiDelta = solParam[6];
    if(std::fabs(phiDelta - twopi) < angularTolerance)
    {
      phiDelta = twopi;
    }
    solid = new G4Cons(sname, solParam[0], solParam[1], solParam[2],
                       solParam[3], solParam[4], solParam[5], phiDelta);
  }
  else if(stype == "SPHERE")
  {
    CheckNoSolidParams(stype, 6, (G4int)solParam.size());
    G4double phiDelta = solParam[3];
    if(std::fabs(phiDelta - twopi) < angularTolerance)
    {
      phiDelta = twopi;
    }
    G4double thetaDelta = solParam[5];
    if(std::fabs(thetaDelta - pi) < angularTolerance)
    {
      thetaDelta = pi;
    }
    solid = new G4Sphere(sname, solParam[0], solParam[1], solParam[2], phiDelta,
                         solParam[4], thetaDelta);
  }
  else if(stype == "ORB")
  {
    CheckNoSolidParams(stype, 1, (G4int)solParam.size());
    solid = new G4Orb(sname, solParam[0]);
  }
  else if(stype == "TORUS")
  {
    CheckNoSolidParams(stype, 5, (G4int)solParam.size());
    G4double phiDelta = solParam[4];
    if(std::fabs(phiDelta - twopi) < angularTolerance)
    {
      phiDelta = twopi;
    }
    solid = new G4Torus(sname, solParam[0], solParam[1], solParam[2],
                        solParam[3], phiDelta);
  }
  else if(stype == "POLYCONE" 
	  || stype == "GENERICPOLYCONE")
  {
    std::size_t nplanes = std::size_t(solParam[2]);
    G4bool genericPoly = false;
    if(solParam.size() == 3 + nplanes * 3)
    {
      genericPoly = false;
    }
    else if(solParam.size() == 3 + nplanes * 2)
    {
      genericPoly = true;
    }
    else
    {
      G4String Err1 = "Solid type " + stype + " should have ";
      G4String Err2 = G4UIcommand::ConvertToString(G4int(3 + nplanes * 3)) +
                      " (Z,Rmin,Rmax)\n";
      G4String Err3 =
        "or " + G4UIcommand::ConvertToString(G4int(3 + nplanes * 2));
      G4String Err4 = " (RZ corners) parameters,\n";
      G4String Err5 =
        "and it has " + G4UIcommand::ConvertToString(G4int(solParam.size()));
      G4String ErrMessage = Err1 + Err2 + Err3 + Err4 + Err5 + " !";
      G4Exception("G4tgbVolume::FindOrConstructG4Solid()", "InvalidSetup",
                  FatalException, ErrMessage);
      return nullptr;
    }

    if(!genericPoly)
    {
      std::vector<G4double>* z_p    = new std::vector<G4double>;
      std::vector<G4double>* rmin_p = new std::vector<G4double>;
      std::vector<G4double>* rmax_p = new std::vector<G4double>;
      for(std::size_t ii = 0; ii < nplanes; ++ii)
      {
        (*z_p).push_back(solParam[3 + 3 * ii]);
        (*rmin_p).push_back(solParam[3 + 3 * ii + 1]);
        (*rmax_p).push_back(solParam[3 + 3 * ii + 2]);
      }
      G4double phiTotal = solParam[1];
      if(std::fabs(phiTotal - twopi) < angularTolerance)
      {
        phiTotal = twopi;
      }
      solid = new G4Polycone(sname, solParam[0], phiTotal,  // start,delta-phi
                             (G4int)nplanes,                // sections
                             &((*z_p)[0]), &((*rmin_p)[0]), &((*rmax_p)[0]));
    }
    else
    {
      std::vector<G4double>* R_c = new std::vector<G4double>;
      std::vector<G4double>* Z_c = new std::vector<G4double>;
      for(size_t ii = 0; ii < nplanes; ii++)
      {
        (*R_c).push_back(solParam[3 + 2 * ii]);
        (*Z_c).push_back(solParam[3 + 2 * ii + 1]);
      }
      G4double phiTotal = solParam[1];
      if(std::fabs(phiTotal - twopi) < angularTolerance)
      {
        phiTotal = twopi;
      }
      solid =
        new G4GenericPolycone(sname, solParam[0], phiTotal,  // start,delta-phi
                              (G4int)nplanes,                // sections
                              &((*R_c)[0]), &((*Z_c)[0]));
    }
  }
  else if(stype == "POLYHEDRA")
  {
    std::size_t nplanes = std::size_t(solParam[3]);
    G4bool genericPoly = false;
    if(solParam.size() == 4 + nplanes * 3)
    {
      genericPoly = false;
    }
    else if(solParam.size() == 4 + nplanes * 2)
    {
      genericPoly = true;
    }
    else
    {
      G4String Err1 = "Solid type " + stype + " should have ";
      G4String Err2 = G4UIcommand::ConvertToString(G4int(4 + nplanes * 3)) +
                      " (Z,Rmin,Rmax)\n";
      G4String Err3 =
        "or " + G4UIcommand::ConvertToString(G4int(4 + nplanes * 2));
      G4String Err4 = " (RZ corners) parameters,\n";
      G4String Err5 =
        "and it has " + G4UIcommand::ConvertToString(G4int(solParam.size()));
      G4String ErrMessage = Err1 + Err2 + Err3 + Err4 + Err5 + " !";
      G4Exception("G4tgbVolume::FindOrConstructG4Solid()", "InvalidSetup",
                  FatalException, ErrMessage);
      return nullptr;
    }

    if(!genericPoly)
    {
      std::vector<G4double>* z_p    = new std::vector<G4double>;
      std::vector<G4double>* rmin_p = new std::vector<G4double>;
      std::vector<G4double>* rmax_p = new std::vector<G4double>;
      for(std::size_t ii = 0; ii < nplanes; ++ii)
      {
        (*z_p).push_back(solParam[4 + 3 * ii]);
        (*rmin_p).push_back(solParam[4 + 3 * ii + 1]);
        (*rmax_p).push_back(solParam[4 + 3 * ii + 2]);
      }
      G4double phiTotal = solParam[1];
      if(std::fabs(phiTotal - twopi) < angularTolerance)
      {
        phiTotal = twopi;
      }
      solid = new G4Polyhedra(sname, solParam[0], phiTotal, G4int(solParam[2]),
                              (G4int)nplanes, &((*z_p)[0]), &((*rmin_p)[0]),
                              &((*rmax_p)[0]));
    }
    else
    {
      std::vector<G4double>* R_c = new std::vector<G4double>;
      std::vector<G4double>* Z_c = new std::vector<G4double>;
      for(std::size_t ii = 0; ii < nplanes; ++ii)
      {
        (*R_c).push_back(solParam[4 + 2 * ii]);
        (*Z_c).push_back(solParam[4 + 2 * ii + 1]);
      }
      G4double phiTotal = solParam[1];
      if(std::fabs(phiTotal - twopi) < angularTolerance)
      {
        phiTotal = twopi;
      }
      solid = new G4Polyhedra(sname, solParam[0], phiTotal, G4int(solParam[2]),
                              (G4int)nplanes, &((*R_c)[0]), &((*Z_c)[0]));
    }
  }
  else if(stype == "ELLIPTICALTUBE")
  {
    CheckNoSolidParams(stype, 3, (G4int)solParam.size());
    solid = new G4EllipticalTube(sname, solParam[0], solParam[1], solParam[2]);
  }
  else if(stype == "ELLIPSOID")
  {
    CheckNoSolidParams(stype, 5, (G4int)solParam.size());
    solid = new G4Ellipsoid(sname, solParam[0], solParam[1], solParam[2],
                            solParam[3], solParam[4]);
  }
  else if(stype == "ELLIPTICALCONE")
  {
    CheckNoSolidParams(stype, 4, (G4int)solParam.size());
    solid = new G4EllipticalCone(sname, solParam[0], solParam[1], solParam[2],
                                 solParam[3]);
  }
  else if(stype == "HYPE")
  {
    CheckNoSolidParams(stype, 5, (G4int)solParam.size());
    solid = new G4Hype(sname, solParam[0], solParam[1], solParam[2],
                       solParam[3], solParam[4]);
  }
  else if(stype == "TET")
  {
    CheckNoSolidParams(stype, 12, (G4int)solParam.size());
    G4ThreeVector anchor(solParam[0], solParam[1], solParam[2]);
    G4ThreeVector p2(solParam[3], solParam[4], solParam[5]);
    G4ThreeVector p3(solParam[6], solParam[7], solParam[8]);
    G4ThreeVector p4(solParam[9], solParam[10], solParam[11]);
    solid = new G4Tet(sname, anchor, p2, p3, p4);
  }
  else if(stype == "TWISTEDBOX")
  {
    CheckNoSolidParams(stype, 4, (G4int)solParam.size());
    solid = new G4TwistedBox(sname, solParam[0], solParam[1], solParam[2],
                             solParam[3]);
  }
  else if(stype == "TWISTEDTRAP")
  {
    CheckNoSolidParams(stype, 11, (G4int)solParam.size());
    solid =
      new G4TwistedTrap(sname, solParam[0], solParam[1], solParam[2],
                        solParam[3], solParam[4], solParam[5], solParam[6],
                        solParam[7], solParam[8], solParam[9], solParam[10]);
  }
  else if(stype == "TWISTEDTRD")
  {
    CheckNoSolidParams(stype, 6, (G4int)solParam.size());
    solid = new G4TwistedTrd(sname, solParam[0], solParam[1], solParam[2],
                             solParam[3], solParam[4], solParam[5]);
  }
    else if(stype == "SCALED")
  {
    const G4tgrSolidScaled* tgrSol = dynamic_cast<const G4tgrSolidScaled*>(sol);
    if(tgrSol == nullptr)
    {
      G4Exception("G4tgbVolume::FindOrConstructG4Solid()", "InvalidSetup",
                  FatalException, "Invalid Solid pointer");
      return nullptr;
    }
    G4VSolid* sol0   = FindOrConstructG4Solid(tgrSol->GetOrigSolid());    
    G4Scale3D scale  = tgrSol->GetScale3d();
    solid  = new G4ScaledSolid(sname, sol0, scale);
  }
  else if(stype == "TWISTEDTUBS")
  {
    CheckNoSolidParams(stype, 5, (G4int)solParam.size());
    G4double phiTotal = solParam[4];
    if(std::fabs(phiTotal - twopi) < angularTolerance)
    {
      phiTotal = twopi;
    }
    solid = new G4TwistedTubs(sname, solParam[0], solParam[1], solParam[2],
                              solParam[3], phiTotal);
  }
  else if(stype == "TESSELLATED")
  {
    G4int nFacets               = G4int(solParam[0]);
    G4int jj                    = 0;
    solid                       = new G4TessellatedSolid(sname);
    G4TessellatedSolid* solidTS = (G4TessellatedSolid*) (solid);
    G4VFacet* facet             = nullptr;

    for(G4int ii = 0; ii < nFacets; ++ii)
    {
      G4int nPoints = G4int(solParam[jj + 1]);
      if(G4int(solParam.size()) < jj + nPoints * 3 + 2)
      {
        G4String Err1 = "Too small number of parameters in tesselated solid, "
                        "it should be at least " +
                        G4UIcommand::ConvertToString(jj + nPoints * 3 + 2);
        G4String Err2 = " facet number " + G4UIcommand::ConvertToString(ii);
        G4String Err3 = " number of parameters is " +
                        G4UIcommand::ConvertToString(G4int(solParam.size()));
        G4String ErrMessage = Err1 + Err2 + Err3 + " !";
        G4Exception("G4tgbVolume::FindOrConstructG4Solid()", "InvalidSetup",
                    FatalException, ErrMessage);
        return nullptr;
      }

      if(nPoints == 3)
      {
        G4ThreeVector pt0(solParam[jj + 2], solParam[jj + 3], solParam[jj + 4]);
        G4ThreeVector vt1(solParam[jj + 5], solParam[jj + 6], solParam[jj + 7]);
        G4ThreeVector vt2(solParam[jj + 8], solParam[jj + 9],
                          solParam[jj + 10]);
        G4FacetVertexType vertexType = ABSOLUTE;
        if(solParam[jj + 11] == 0)
        {
          vertexType = ABSOLUTE;
        }
        else if(solParam[jj + 11] == 1)
        {
          vertexType = RELATIVE;
        }
        else
        {
          G4String Err1 = "Wrong number of vertex type in tesselated solid, it "
                          "should be 0 =ABSOLUTE) or 1 (=RELATIVE)";
          G4String Err2 =
            " facet number " + G4UIcommand::ConvertToString(G4int(ii));
          G4String Err3 = " vertex type is " +
                          G4UIcommand::ConvertToString(solParam[jj + 11]);
          G4String ErrMessage = Err1 + Err2 + Err3 + " !";
          G4Exception("G4tgbVolume::FindOrConstructG4Solid()", "InvalidSetup",
                      FatalException, ErrMessage);
          return nullptr;
        }
        facet = new G4TriangularFacet(pt0, vt1, vt2, vertexType);
      }
      else if(nPoints == 4)
      {
        G4ThreeVector pt0(solParam[jj + 2], solParam[jj + 3], solParam[jj + 4]);
        G4ThreeVector vt1(solParam[jj + 5], solParam[jj + 6], solParam[jj + 7]);
        G4ThreeVector vt2(solParam[jj + 8], solParam[jj + 9],
                          solParam[jj + 10]);
        G4ThreeVector vt3(solParam[jj + 11], solParam[jj + 12],
                          solParam[jj + 13]);
        G4FacetVertexType vertexType = ABSOLUTE;
        if(solParam[jj + 14] == 0)
        {
          vertexType = ABSOLUTE;
        }
        else if(solParam[jj + 14] == 1)
        {
          vertexType = RELATIVE;
        }
        else
        {
          G4String Err1 = "Wrong number of vertex type in tesselated solid, it "
                          "should be 0 =ABSOLUTE) or 1 (=RELATIVE)";
          G4String Err2 =
            " facet number " + G4UIcommand::ConvertToString(G4int(ii));
          G4String Err3 = " vertex type is " +
                          G4UIcommand::ConvertToString(solParam[jj + 14]);
          G4String ErrMessage = Err1 + Err2 + Err3 + " !";
          G4Exception("G4tgbVolume::FindOrConstructG4Solid()", "InvalidSetup",
                      FatalException, ErrMessage);
          return nullptr;
        }
        facet = new G4QuadrangularFacet(pt0, vt1, vt2, vt3, vertexType);
      }
      else
      {
        G4String Err1 =
          "Wrong number of points in tesselated solid, it should be 3 or 4";
        G4String Err2 =
          " facet number " + G4UIcommand::ConvertToString(G4int(ii));
        G4String Err3 = " number of points is " +
                        G4UIcommand::ConvertToString(G4int(nPoints));
        G4String ErrMessage = Err1 + Err2 + Err3 + " !";
        G4Exception("G4tgbVolume::FindOrConstructG4Solid()", "InvalidSetup",
                    FatalException, ErrMessage);
        return nullptr;
      }

      solidTS->AddFacet(facet);
      jj += nPoints * 3 + 2;
    }
  }
  else if(stype == "EXTRUDED")
  {
    std::vector<G4TwoVector> polygonList;
    std::vector<G4ExtrudedSolid::ZSection> zsectionList;
    G4int nPolygons = G4int(solParam[0]);
    G4int ii        = 1;
    G4int nMax      = nPolygons * 2 + 1;
    for(; ii < nMax; ii += 2)
    {
      polygonList.push_back(G4TwoVector(solParam[ii], solParam[ii + 1]));
    }
    G4int nZSections = G4int(solParam[ii]);
    nMax             = nPolygons * 2 + nZSections * 4 + 2;
    ++ii;
    for(; ii < nMax; ii += 4)
    {
      G4TwoVector offset(solParam[ii + 1], solParam[ii + 2]);
      zsectionList.push_back(
        G4ExtrudedSolid::ZSection(solParam[ii], offset, solParam[ii + 3]));
    }
    solid = new G4ExtrudedSolid(sname, polygonList, zsectionList);
  }
  else if(stype.substr(0, 7) == "Boolean")
  {
    const G4tgrSolidBoolean* solb = dynamic_cast<const G4tgrSolidBoolean*>(sol);
    if(solb == nullptr)
    {
      G4Exception("G4tgbVolume::FindOrConstructG4Solid()", "InvalidSetup",
                  FatalException, "Invalid Solid pointer");
      return nullptr;
    }
    G4VSolid* sol1 = FindOrConstructG4Solid(solb->GetSolid(0));
    G4VSolid* sol2 = FindOrConstructG4Solid(solb->GetSolid(1));
    G4RotationMatrix* relRotMat =
      G4tgbRotationMatrixMgr::GetInstance()->FindOrBuildG4RotMatrix(
        sol->GetRelativeRotMatName());
    G4ThreeVector relPlace = solb->GetRelativePlace();

    if(stype == "Boolean_UNION")
    {
      solid = new G4UnionSolid(sname, sol1, sol2, relRotMat, relPlace);
    }
    else if(stype == "Boolean_SUBTRACTION")
    {
      solid = new G4SubtractionSolid(sname, sol1, sol2, relRotMat, relPlace);
    }
    else if(stype == "Boolean_INTERSECTION")
    {
      solid = new G4IntersectionSolid(sname, sol1, sol2, relRotMat, relPlace);
    }
    else
    {
      G4String ErrMessage = "Unknown Boolean type " + stype;
      G4Exception("G4tgbVolume::FindOrConstructG4Solid()", "InvalidSetup",
                  FatalException, ErrMessage);
      return nullptr;
    }
  }
  else if(stype == "MULTIUNION")
  {
    const G4tgrSolidMultiUnion* tgrSol = dynamic_cast<const G4tgrSolidMultiUnion*>(sol);
    if(tgrSol == nullptr)
    {
      G4Exception("G4tgbVolume::FindOrConstructG4Solid()", "InvalidSetup",
                  FatalException, "Invalid Solid pointer");
      return nullptr;
    }

    G4int nsol =  tgrSol->GetNSolid();
    G4VSolid*     soli;
    G4Transform3D tri;
    G4MultiUnion* solidu = new G4MultiUnion(sname);

    for (G4int i=0; i<nsol; ++i)
    {
      soli = FindOrConstructG4Solid(tgrSol->GetSolid(i));
      tri  = tgrSol->GetTransformation(i);
      solidu->AddNode(*soli, tri);
    }
    solidu->Voxelize();
    solid = dynamic_cast<G4VSolid*>(solidu);
  }
  else
  {
    G4String ErrMessage =
      "Solids of type " + stype + " not implemented yet, sorry...";
    G4Exception("G4tgbVolume::FindOrConstructG4Solid()", "NotImplemented",
                FatalException, ErrMessage);
    return nullptr;
  }

#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 2)
  {
    G4cout << " G4tgbVolume::FindOrConstructG4Solid()" << G4endl
           << "   Created solid " << sname << " of type "
           << solid->GetEntityType() << G4endl;
  }
#endif

#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 1)
  {
    G4cout << " Constructing new G4Solid: " << *solid << G4endl;
  }
#endif

  return solid;
}

// --------------------------------------------------------------------
void G4tgbVolume::CheckNoSolidParams(const G4String& solidType,
                                     const unsigned int NoParamExpected,
                                     const unsigned int NoParam)
{
  if(NoParamExpected != NoParam)
  {
    G4String Err1 = "Solid type " + solidType + " should have ";
    G4String Err2 =
      G4UIcommand::ConvertToString(G4int(NoParamExpected)) + " parameters,\n";
    G4String Err3 =
      "and it has " + G4UIcommand::ConvertToString(G4int(NoParam));
    G4String ErrMessage = Err1 + Err2 + Err3 + " !";
    G4Exception("G4tgbVolume::CheckNoSolidParams()", "InvalidSetup",
                FatalException, ErrMessage);
  }
}

// --------------------------------------------------------------------
G4LogicalVolume* G4tgbVolume::ConstructG4LogVol(const G4VSolid* solid)
{
  G4LogicalVolume* logvol;

#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 2)
  {
    G4cout << " G4tgbVolume::ConstructG4LogVol() - " << GetName() << G4endl;
  }
#endif

  //----------- Get the material first
  G4Material* mate = G4tgbMaterialMgr::GetInstance()->FindOrBuildG4Material(
    theTgrVolume->GetMaterialName());
  if(mate == nullptr)
  {
    G4String ErrMessage = "Material not found " +
                          theTgrVolume->GetMaterialName() + " for volume " +
                          GetName() + ".";
    G4Exception("G4tgbVolume::ConstructG4LogVol()", "InvalidSetup",
                FatalException, ErrMessage);
  }
#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 2)
  {
    G4cout << " G4tgbVolume::ConstructG4LogVol() -"
           << " Material constructed: " << mate->GetName() << G4endl;
  }
#endif

  //---------- Construct the LV
  logvol = new G4LogicalVolume(const_cast<G4VSolid*>(solid),
                               const_cast<G4Material*>(mate), GetName());

#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 1)
  {
    G4cout << " Constructing new G4LogicalVolume: " << logvol->GetName()
           << " mate " << mate->GetName() << G4endl;
  }
#endif

  //---------- Set visibility and colour
  if(!GetVisibility() || GetColour()[0] != -1)
  {
    G4VisAttributes* visAtt = new G4VisAttributes();
#ifdef G4VERBOSE
    if(G4tgrMessenger::GetVerboseLevel() >= 1)
    {
      G4cout << " Constructing new G4VisAttributes: " << *visAtt << G4endl;
    }
#endif

    if(!GetVisibility())
    {
      visAtt->SetVisibility(GetVisibility());
    }
    else if(GetColour()[0] != -1)
    {
      // this else should not be necessary, because if the visibility
      // is set to off, colour should have no effect. But it does not
      // work: if you set colour and vis off, it is visualized!?!?!?

      const G4double* col = GetColour();
      if(col[3] == -1.)
      {
        visAtt->SetColour(G4Colour(col[0], col[1], col[2]));
      }
      else
      {
        visAtt->SetColour(G4Colour(col[0], col[1], col[2], col[3]));
      }
    }
    logvol->SetVisAttributes(visAtt);
  }

#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 2)
  {
    G4cout << " G4tgbVolume::ConstructG4LogVol() -"
           << " Created logical volume: " << GetName() << G4endl;
  }
#endif

  return logvol;
}

// --------------------------------------------------------------------
G4VPhysicalVolume*
G4tgbVolume::ConstructG4PhysVol(const G4tgrPlace* place,
                                const G4LogicalVolume* currentLV,
                                const G4LogicalVolume* parentLV)
{
  G4VPhysicalVolume* physvol = nullptr;
  G4int copyNo;

  //----- Case of placement of top volume
  if(place == nullptr)
  {
#ifdef G4VERBOSE
    if(G4tgrMessenger::GetVerboseLevel() >= 2)
    {
      G4cout << " G4tgbVolume::ConstructG4PhysVol() - World: " << GetName()
             << G4endl;
    }
#endif
    physvol = new G4PVPlacement(
      nullptr, G4ThreeVector(), const_cast<G4LogicalVolume*>(currentLV),
      GetName(), 0, false, 0, theTgrVolume->GetCheckOverlaps());
#ifdef G4VERBOSE
    if(G4tgrMessenger::GetVerboseLevel() >= 1)
    {
      G4cout << " Constructing new : G4PVPlacement " << physvol->GetName()
             << G4endl;
    }
#endif
  }
  else
  {
    copyNo = place->GetCopyNo();

#ifdef G4VERBOSE
    if(G4tgrMessenger::GetVerboseLevel() >= 2)
    {
      G4cout << " G4tgbVolume::ConstructG4PhysVol() - " << GetName() << G4endl
             << "   inside " << parentLV->GetName() << " copy No: " << copyNo
             << " of type= " << theTgrVolume->GetType() << G4endl
             << "   placement type= " << place->GetType() << G4endl;
    }
#endif

    if(theTgrVolume->GetType() == "VOLSimple")
    {
      //----- Get placement
#ifdef G4VERBOSE
      if(G4tgrMessenger::GetVerboseLevel() >= 2)
      {
        G4cout << " G4tgbVolume::ConstructG4PhysVol() - Placement type = "
               << place->GetType() << G4endl;
      }
#endif

      //--------------- If it is  G4tgrPlaceSimple
      if(place->GetType() == "PlaceSimple")
      {
        //----- Get rotation matrix
        G4tgrPlaceSimple* placeSimple = (G4tgrPlaceSimple*) place;
        G4String rmName               = placeSimple->GetRotMatName();

        G4RotationMatrix* rotmat =
          G4tgbRotationMatrixMgr::GetInstance()->FindOrBuildG4RotMatrix(rmName);
        //----- Place volume in mother
        G4double check =
          (rotmat->colX().cross(rotmat->colY())) * rotmat->colZ();
        G4double tol = 1.0e-3;
        //---- Check that matrix is ortogonal
        if(1 - std::abs(check) > tol)
        {
          G4cerr << " Matrix : " << rmName << " " << rotmat->colX() << " "
                 << rotmat->colY() << " " << rotmat->colZ() << G4endl
                 << " product x X y * z = " << check << " x X y "
                 << rotmat->colX().cross(rotmat->colY()) << G4endl;
          G4String ErrMessage = "Rotation is not ortogonal " + rmName + " !";
          G4Exception("G4tgbVolume::ConstructG4PhysVol()", "InvalidSetup",
                      FatalException, ErrMessage);
          //---- Check if it is reflection
        }
        else if(1 + check <= tol)
        {
          G4Translate3D transl = place->GetPlacement();
          G4Transform3D trfrm  = transl * G4Rotate3D(*rotmat);
          physvol =
            (G4ReflectionFactory::Instance()->Place(
               trfrm, GetName(), const_cast<G4LogicalVolume*>(currentLV),
               const_cast<G4LogicalVolume*>(parentLV), false, copyNo, false))
              .first;
        }
        else
        {
#ifdef G4VERBOSE
          if(G4tgrMessenger::GetVerboseLevel() >= 1)
          {
            G4cout << "Construction new G4VPhysicalVolume"
                   << " through G4ReflectionFactory " << GetName()
                   << " in volume " << parentLV->GetName() << " copyNo "
                   << copyNo << " position " << place->GetPlacement() << " ROT "
                   << rotmat->colX() << " " << rotmat->colY() << " "
                   << rotmat->colZ() << G4endl;
          }
#endif
          physvol =
            new G4PVPlacement(rotmat, place->GetPlacement(),
                              const_cast<G4LogicalVolume*>(currentLV),
                              GetName(), const_cast<G4LogicalVolume*>(parentLV),
                              false, copyNo, theTgrVolume->GetCheckOverlaps());
        }

        //--------------- If it is G4tgrPlaceParam
      }
      else if(place->GetType() == "PlaceParam")
      {
        G4tgrPlaceParameterisation* dp = (G4tgrPlaceParameterisation*) (place);

        //----- See what parameterisation type
#ifdef G4VERBOSE
        if(G4tgrMessenger::GetVerboseLevel() >= 2)
        {
          G4cout << " G4tgbVolume::ConstructG4PhysVol() -" << G4endl
                 << "   param: " << GetName() << " in " << parentLV->GetName()
                 << " param type= " << dp->GetParamType() << G4endl;
        }
#endif

        G4tgbPlaceParameterisation* param = nullptr;

        if((dp->GetParamType() == "CIRCLE") ||
           (dp->GetParamType() == "CIRCLE_XY") ||
           (dp->GetParamType() == "CIRCLE_XZ") ||
           (dp->GetParamType() == "CIRCLE_YZ"))
        {
          param = new G4tgbPlaceParamCircle(dp);
        }
        else if((dp->GetParamType() == "LINEAR") ||
                (dp->GetParamType() == "LINEAR_X") ||
                (dp->GetParamType() == "LINEAR_Y") ||
                (dp->GetParamType() == "LINEAR_Z"))
        {
          param = new G4tgbPlaceParamLinear(dp);
        }
        else if((dp->GetParamType() == "SQUARE") ||
                (dp->GetParamType() == "SQUARE_XY") ||
                (dp->GetParamType() == "SQUARE_XZ") ||
                (dp->GetParamType() == "SQUARE_YZ"))
        {
          param = new G4tgbPlaceParamSquare(dp);
        }
        else
        {
          G4String ErrMessage = "Parameterisation has wrong type, TYPE: " +
                                G4String(dp->GetParamType()) + " !";
          G4Exception("G4tgbVolume::ConstructG4PhysVol", "WrongArgument",
                      FatalException, ErrMessage);
          return nullptr;
        }
#ifdef G4VERBOSE
        if(G4tgrMessenger::GetVerboseLevel() >= 1)
        {
          G4cout << " G4tgbVolume::ConstructG4PhysVol() -" << G4endl
                 << "   New G4PVParameterised: " << GetName() << " vol "
                 << currentLV->GetName() << " in vol " << parentLV->GetName()
                 << " axis " << param->GetAxis() << " nCopies "
                 << param->GetNCopies() << G4endl;
        }
#endif
        physvol = new G4PVParameterised(
          GetName(), const_cast<G4LogicalVolume*>(currentLV),
          const_cast<G4LogicalVolume*>(parentLV), EAxis(param->GetAxis()),
          param->GetNCopies(), param);
#ifdef G4VERBOSE
        if(G4tgrMessenger::GetVerboseLevel() >= 1)
        {
          G4cout << " Constructing new G4PVParameterised: "
                 << physvol->GetName() << " in volume " << parentLV->GetName()
                 << " N copies " << param->GetNCopies() << " axis "
                 << param->GetAxis() << G4endl;
        }
#endif
      }
      else if(place->GetType() == "PlaceReplica")
      {
        //--------------- If it is  PlaceReplica
        G4tgrPlaceDivRep* dpr = (G4tgrPlaceDivRep*) place;

#ifdef G4VERBOSE
        if(G4tgrMessenger::GetVerboseLevel() >= 2)
        {
          G4cout << " G4tgbVolume::ConstructG4PhysVol() -" << G4endl
                 << "   replica"
                 << " " << currentLV->GetName() << " in " << parentLV->GetName()
                 << " NDiv " << dpr->GetNDiv() << " Width " << dpr->GetWidth()
                 << " offset " << dpr->GetOffset() << G4endl;
        }
#endif
        physvol = new G4PVReplica(
          GetName(), const_cast<G4LogicalVolume*>(currentLV),
          const_cast<G4LogicalVolume*>(parentLV), EAxis(dpr->GetAxis()),
          dpr->GetNDiv(), dpr->GetWidth(), dpr->GetOffset());
#ifdef G4VERBOSE
        if(G4tgrMessenger::GetVerboseLevel() >= 1)
        {
          G4cout << " Constructing new G4PVReplica: " << currentLV->GetName()
                 << " in " << parentLV->GetName() << " NDiv " << dpr->GetNDiv()
                 << " Width " << dpr->GetWidth() << " offset "
                 << dpr->GetOffset() << G4endl;
        }
#endif
      }
    }
    else if(theTgrVolume->GetType() == "VOLDivision")
    {
      G4tgrVolumeDivision* volr  = (G4tgrVolumeDivision*) theTgrVolume;
      G4tgrPlaceDivRep* placeDiv = volr->GetPlaceDivision();
      G4VSolid* solid =
        BuildSolidForDivision(parentLV->GetSolid(), placeDiv->GetAxis());
      G4Material* mate = G4tgbMaterialMgr::GetInstance()->FindOrBuildG4Material(
        theTgrVolume->GetMaterialName());
      G4LogicalVolume* divLV =
        new G4LogicalVolume(solid, const_cast<G4Material*>(mate), GetName());
#ifdef G4VERBOSE
      if(G4tgrMessenger::GetVerboseLevel() >= 1)
      {
        G4cout << " Constructed new G4LogicalVolume for division: "
               << divLV->GetName() << " mate " << mate->GetName() << G4endl;
      }
#endif

      G4DivType divType = placeDiv->GetDivType();
      switch(divType)
      {
        case DivByNdiv:
          physvol = new G4PVDivision(GetName(), (G4LogicalVolume*) divLV,
                                     const_cast<G4LogicalVolume*>(parentLV),
                                     placeDiv->GetAxis(), placeDiv->GetNDiv(),
                                     placeDiv->GetOffset());
#ifdef G4VERBOSE
          if(G4tgrMessenger::GetVerboseLevel() >= 1)
          {
            G4cout << " Constructing new G4PVDivision by number of divisions: "
                   << GetName() << " in " << parentLV->GetName() << " axis "
                   << placeDiv->GetAxis() << " Ndiv " << placeDiv->GetNDiv()
                   << " offset " << placeDiv->GetOffset() << G4endl;
          }
#endif
          break;
        case DivByWidth:
          physvol = new G4PVDivision(GetName(), (G4LogicalVolume*) divLV,
                                     const_cast<G4LogicalVolume*>(parentLV),
                                     placeDiv->GetAxis(), placeDiv->GetWidth(),
                                     placeDiv->GetOffset());
#ifdef G4VERBOSE
          if(G4tgrMessenger::GetVerboseLevel() >= 1)
          {
            G4cout << " Constructing new G4PVDivision by width: " << GetName()
                   << " in " << parentLV->GetName() << " axis "
                   << placeDiv->GetAxis() << " width " << placeDiv->GetWidth()
                   << " offset " << placeDiv->GetOffset() << G4endl;
          }
#endif
          break;
        case DivByNdivAndWidth:
          physvol = new G4PVDivision(
            GetName(), (G4LogicalVolume*) divLV,
            const_cast<G4LogicalVolume*>(parentLV), placeDiv->GetAxis(),
            placeDiv->GetNDiv(), placeDiv->GetWidth(), placeDiv->GetOffset());
#ifdef G4VERBOSE
          if(G4tgrMessenger::GetVerboseLevel() >= 1)
          {
            G4cout << " Constructing new G4PVDivision by width"
                   << " and number of divisions: " << GetName() << " in "
                   << parentLV->GetName() << " axis " << placeDiv->GetAxis()
                   << " Ndiv " << placeDiv->GetNDiv() << " width "
                   << placeDiv->GetWidth() << " offset "
                   << placeDiv->GetOffset() << G4endl;
          }
#endif
          break;
      }
    }
    else if(theTgrVolume->GetType() == "VOLAssembly")
    {
      // Define one layer as one assembly volume
      G4tgrVolumeAssembly* tgrAssembly = (G4tgrVolumeAssembly*) theTgrVolume;

      if(!theG4AssemblyVolume)
      {
        theG4AssemblyVolume = new G4AssemblyVolume;
#ifdef G4VERBOSE
        if(G4tgrMessenger::GetVerboseLevel() >= 1)
        {
          G4cout << " Constructing new G4AssemblyVolume: "
                 << " number of assembly components "
                 << tgrAssembly->GetNoComponents() << G4endl;
        }
#endif
        G4tgbVolumeMgr* g4vmgr = G4tgbVolumeMgr::GetInstance();
        for(G4int ii = 0; ii < tgrAssembly->GetNoComponents(); ++ii)
        {
          // Rotation and translation of a plate inside the assembly

          G4ThreeVector transl = tgrAssembly->GetComponentPos(ii);
          G4String rmName      = tgrAssembly->GetComponentRM(ii);
          G4RotationMatrix* rotmat =
            G4tgbRotationMatrixMgr::GetInstance()->FindOrBuildG4RotMatrix(
              rmName);

          //----- Get G4LogicalVolume of component
          G4String lvname         = tgrAssembly->GetComponentName(ii);
          G4LogicalVolume* logvol = g4vmgr->FindG4LogVol(lvname);
          if(logvol == nullptr)
          {
            g4vmgr->FindVolume(lvname)->ConstructG4Volumes(0, 0);
            logvol = g4vmgr->FindG4LogVol(lvname, true);
          }
          // Fill the assembly by the plates
          theG4AssemblyVolume->AddPlacedVolume(logvol, transl, rotmat);
#ifdef G4VERBOSE
          if(G4tgrMessenger::GetVerboseLevel() >= 1)
          {
            G4cout << " G4AssemblyVolume->AddPlacedVolume " << ii << " "
                   << logvol->GetName() << " translation " << transl
                   << " rotmat " << rotmat->colX() << " " << rotmat->colY()
                   << " " << rotmat->colZ() << G4endl;
          }
#endif
        }
      }

      // Rotation and Translation of the assembly inside the world

      G4tgrPlaceSimple* placeSimple = (G4tgrPlaceSimple*) place;
      G4String rmName               = placeSimple->GetRotMatName();
      G4RotationMatrix* rotmat =
        G4tgbRotationMatrixMgr::GetInstance()->FindOrBuildG4RotMatrix(rmName);
      G4ThreeVector transl = place->GetPlacement();

      G4LogicalVolume* parentLV_nonconst =
        const_cast<G4LogicalVolume*>(parentLV);
      theG4AssemblyVolume->MakeImprint(parentLV_nonconst, transl, rotmat);
    }
    else  // If it is G4tgrVolumeAssembly
    {
      G4String ErrMessage =
        "Volume type not supported: " + theTgrVolume->GetType() + ", sorry...";
      G4Exception("G4tgbVolume::ConstructG4PhysVol()", "NotImplemented",
                  FatalException, ErrMessage);
    }
  }

  return physvol;
}

// --------------------------------------------------------------------
G4VSolid* G4tgbVolume::BuildSolidForDivision(G4VSolid* parentSolid, EAxis axis)
{
  G4VSolid* solid = nullptr;
  G4double redf =
    (parentSolid->GetExtent().GetXmax() - parentSolid->GetExtent().GetXmin());
  redf = std::min(redf, parentSolid->GetExtent().GetYmax() -
                          parentSolid->GetExtent().GetYmin());
  redf = std::min(redf, parentSolid->GetExtent().GetZmax() -
                          parentSolid->GetExtent().GetZmin());
  redf *= 0.001;  // make daugther much smaller, to fit in parent

  if(parentSolid->GetEntityType() == "G4Box")
  {
    G4Box* psolid = (G4Box*) (parentSolid);
    solid         = new G4Box(GetName(), psolid->GetXHalfLength() * redf,
                      psolid->GetZHalfLength() * redf,
                      psolid->GetZHalfLength() * redf);
  }
  else if(parentSolid->GetEntityType() == "G4Tubs")
  {
    G4Tubs* psolid = (G4Tubs*) (parentSolid);
    solid          = new G4Tubs(GetName(), psolid->GetInnerRadius() * redf,
                       psolid->GetOuterRadius() * redf,
                       psolid->GetZHalfLength() * redf,
                       psolid->GetStartPhiAngle(), psolid->GetDeltaPhiAngle());
  }
  else if(parentSolid->GetEntityType() == "G4Cons")
  {
    G4Cons* psolid = (G4Cons*) (parentSolid);
    solid = new G4Cons(GetName(), psolid->GetInnerRadiusMinusZ() * redf,
                       psolid->GetOuterRadiusMinusZ() * redf,
                       psolid->GetInnerRadiusPlusZ() * redf,
                       psolid->GetOuterRadiusPlusZ() * redf,
                       psolid->GetZHalfLength() * redf,
                       psolid->GetStartPhiAngle(), psolid->GetDeltaPhiAngle());
  }
  else if(parentSolid->GetEntityType() == "G4Trd")
  {
    G4Trd* psolid  = (G4Trd*) (parentSolid);
    G4double mpDx1 = psolid->GetXHalfLength1();
    G4double mpDx2 = psolid->GetXHalfLength2();

    if(axis == kXAxis &&
       std::fabs(mpDx1 - mpDx2) >
         G4GeometryTolerance::GetInstance()->GetSurfaceTolerance())
    {
      solid = new G4Trap(GetName(), psolid->GetZHalfLength() * redf,
                         psolid->GetYHalfLength1() * redf,
                         psolid->GetXHalfLength2() * redf,
                         psolid->GetXHalfLength1() * redf);
    }
    else
    {
      solid = new G4Trd(
        GetName(), psolid->GetXHalfLength1() * redf,
        psolid->GetXHalfLength2() * redf, psolid->GetYHalfLength1() * redf,
        psolid->GetYHalfLength2() * redf, psolid->GetZHalfLength() * redf);
    }
  }
  else if(parentSolid->GetEntityType() == "G4Para")
  {
    G4Para* psolid = (G4Para*) (parentSolid);
    solid          = new G4Para(
      GetName(), psolid->GetXHalfLength() * redf,
      psolid->GetYHalfLength() * redf, psolid->GetZHalfLength() * redf,
      std::atan(psolid->GetTanAlpha()), psolid->GetSymAxis().theta(),
      psolid->GetSymAxis().phi());
  }
  else if(parentSolid->GetEntityType() == "G4Polycone")
  {
    G4Polycone* psolid             = (G4Polycone*) (parentSolid);
    G4PolyconeHistorical origParam = *(psolid->GetOriginalParameters());
    for(G4int ii = 0; ii < origParam.Num_z_planes; ++ii)
    {
      origParam.Rmin[ii] = origParam.Rmin[ii] * redf;
      origParam.Rmax[ii] = origParam.Rmax[ii] * redf;
    }
    solid = new G4Polycone(GetName(), psolid->GetStartPhi(),
                           psolid->GetEndPhi(), origParam.Num_z_planes,
                           origParam.Z_values, origParam.Rmin, origParam.Rmax);
  }
  else if(parentSolid->GetEntityType() == "G4GenericPolycone")
  {
    G4GenericPolycone* psolid = (G4GenericPolycone*) (parentSolid);
    const G4int numRZ         = psolid->GetNumRZCorner();
    G4double* r               = new G4double[numRZ];
    G4double* z               = new G4double[numRZ];
    for(G4int ii = 0; ii < numRZ; ++ii)
    {
      r[ii] = psolid->GetCorner(ii).r;
      z[ii] = psolid->GetCorner(ii).z;
    }
    solid = new G4GenericPolycone(GetName(), psolid->GetStartPhi(),
                                  psolid->GetEndPhi() - psolid->GetStartPhi(),
                                  numRZ, r, z);
    delete[] r;
    delete[] z;
  }
  else if(parentSolid->GetEntityType() == "G4Polyhedra")
  {
    G4Polyhedra* psolid             = (G4Polyhedra*) (parentSolid);
    G4PolyhedraHistorical origParam = *(psolid->GetOriginalParameters());
    for(G4int ii = 0; ii < origParam.Num_z_planes; ++ii)
    {
      origParam.Rmin[ii] = origParam.Rmin[ii] * redf;
      origParam.Rmax[ii] = origParam.Rmax[ii] * redf;
    }
    solid =
      new G4Polyhedra(GetName(), psolid->GetStartPhi(), psolid->GetEndPhi(),
                      psolid->GetNumSide(), origParam.Num_z_planes,
                      origParam.Z_values, origParam.Rmin, origParam.Rmax);
  }
  else
  {
    G4String ErrMessage = "Solid type not supported. VOLUME= " + GetName() +
                          " Solid type= " + parentSolid->GetEntityType() +
                          "\n" +
                          "Only supported types are: G4Box, G4Tubs, G4Cons," +
                          " G4Trd, G4Para, G4Polycone, G4Polyhedra.";
    G4Exception("G4tgbVolume::BuildSolidForDivision()", "NotImplemented",
                FatalException, ErrMessage);
    return nullptr;
  }

#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 1)
  {
    G4cout << " Constructing new G4Solid for division: " << *solid << G4endl;
  }
#endif
  return solid;
}
