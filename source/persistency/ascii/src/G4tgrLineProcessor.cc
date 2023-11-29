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
// G4tgrLineProcessor
//
// Author: P.Arce, CIEMAT (November 2007)
// --------------------------------------------------------------------

#include "G4tgrLineProcessor.hh"

#include "G4SystemOfUnits.hh"
#include "G4tgrParameterMgr.hh"
#include "G4tgrFileIn.hh"
#include "G4tgrElementSimple.hh"
#include "G4tgrElementFromIsotopes.hh"
#include "G4tgrVolume.hh"
#include "G4tgrVolumeDivision.hh"
#include "G4tgrVolumeAssembly.hh"
#include "G4tgrPlaceDivRep.hh"
#include "G4tgrPlaceParameterisation.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMaterialFactory.hh"
#include "G4tgrRotationMatrixFactory.hh"
#include "G4tgrMessenger.hh"

// --------------------------------------------------------------------
G4tgrLineProcessor::G4tgrLineProcessor()
{
  volmgr = G4tgrVolumeMgr::GetInstance();
}

// --------------------------------------------------------------------
G4tgrLineProcessor::~G4tgrLineProcessor()
{
}

// --------------------------------------------------------------------
G4bool G4tgrLineProcessor::ProcessLine(const std::vector<G4String>& wl)
{
#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 1)
  {
    G4tgrUtils::DumpVS(wl, "@@@ Processing input line");
  }
#endif

  G4String wl0 = wl[0];
  for(G4int ii = 0; ii < (G4int)wl0.length(); ++ii)
  {
    wl0[ii] = (char)std::toupper(wl0[ii]);
  }

  //------------------------------- parameter number
  if(wl0 == ":P")
  {
    G4tgrParameterMgr::GetInstance()->AddParameterNumber(wl);

    //------------------------------- parameter string
  }
  else if(wl0 == ":PS")
  {
    G4tgrParameterMgr::GetInstance()->AddParameterString(wl);

    //------------------------------- isotope
  }
  else if(wl0 == ":ISOT")
  {
    G4tgrIsotope* isot = G4tgrMaterialFactory::GetInstance()->AddIsotope(wl);
    volmgr->RegisterMe(isot);

    //------------------------------- element
  }
  else if(wl0 == ":ELEM")
  {
    G4tgrElementSimple* elem =
      G4tgrMaterialFactory::GetInstance()->AddElementSimple(wl);
    volmgr->RegisterMe(elem);

    //------------------------------- element from isotopes
  }
  else if(wl0 == ":ELEM_FROM_ISOT")
  {
    //:ELEM_FROM_ISOT NAME SYMBOL N_ISOT (ISOT_NAME ISOT_ABUNDANCE)
    G4tgrElementFromIsotopes* elem =
      G4tgrMaterialFactory::GetInstance()->AddElementFromIsotopes(wl);
    volmgr->RegisterMe(elem);

    //------------------------------- material
  }
  else if(wl0 == ":MATE")
  {
    G4tgrMaterialSimple* mate =
      G4tgrMaterialFactory::GetInstance()->AddMaterialSimple(wl);
    volmgr->RegisterMe(mate);

    //------------------------------- material mixtures & by weight
  }
  else if((wl0 == ":MIXT") || (wl0 == ":MIXT_BY_WEIGHT"))
  {
    G4tgrMaterialMixture* mate =
      G4tgrMaterialFactory::GetInstance()->AddMaterialMixture(
        wl, "MaterialMixtureByWeight");
    volmgr->RegisterMe(mate);

    //------------------------------- material mixture by number of atoms
  }
  else if(wl0 == ":MIXT_BY_NATOMS")
  {
    G4tgrMaterialMixture* mate =
      G4tgrMaterialFactory::GetInstance()->AddMaterialMixture(
        wl, "MaterialMixtureByNoAtoms");
    volmgr->RegisterMe(mate);

    //------------------------------- material mixture by volume
  }
  else if(wl0 == ":MIXT_BY_VOLUME")
  {
    G4tgrMaterialMixture* mate =
      G4tgrMaterialFactory::GetInstance()->AddMaterialMixture(
        wl, "MaterialMixtureByVolume");
    volmgr->RegisterMe(mate);

    //------------------------------- material Mean Excitation Energy of
    //                                Ionisation Potential
  }
  else if(wl0 == ":MATE_MEE")
  {
    G4tgrMaterial* mate = G4tgrMaterialFactory::GetInstance()->FindMaterial(
      G4tgrUtils::GetString(wl[1]));
    if(mate == 0)
    {
      G4Exception("G4tgrLineProcessor::ProcessLine()", "Material not found",
                  FatalException, G4tgrUtils::GetString(wl[1]));
      return false;
    }
    mate->SetIonisationMeanExcitationEnergy(G4tgrUtils::GetDouble(wl[2]));

    //------------------------------- material
  }
  else if(wl0 == ":MATE_STATE")
  {
    G4tgrMaterial* mate = G4tgrMaterialFactory::GetInstance()->FindMaterial(
      G4tgrUtils::GetString(wl[1]));
    if(mate == 0)
    {
      G4Exception("G4tgrLineProcessor::ProcessLine()", "Material not found",
                  FatalException, wl[1]);
    }
    mate->SetState(wl[2]);

    //------------------------------- material
  }
  else if(wl0 == ":MATE_TEMPERATURE")
  {
    G4tgrMaterial* mate = G4tgrMaterialFactory::GetInstance()->FindMaterial(
      G4tgrUtils::GetString(wl[1]));
    if(mate == 0)
    {
      G4Exception("G4tgrLineProcessor::ProcessLine()", "Material not found",
                  FatalException, wl[1]);
    }
    mate->SetTemperature(G4tgrUtils::GetDouble(wl[2], kelvin));

    //------------------------------- material
  }
  else if(wl0 == ":MATE_PRESSURE")
  {
    G4tgrMaterial* mate = G4tgrMaterialFactory::GetInstance()->FindMaterial(
      G4tgrUtils::GetString(wl[1]));
    if(mate == 0)
    {
      G4Exception("G4tgrLineProcessor::ProcessLine()", "Material not found",
                  FatalException, wl[1]);
    }
    mate->SetPressure(G4tgrUtils::GetDouble(wl[2], atmosphere));

    //------------------------------- solid
  }
  else if(wl0 == ":SOLID")
  {  // called from here or from G4tgrVolume::G4tgrVolume
    volmgr->CreateSolid(wl, 0);

    //------------------------------- volume
  }
  else if(wl0 == ":VOLU")
  {
    G4tgrVolume* vol = new G4tgrVolume(wl);
    volmgr->RegisterMe(vol);

    //--------------------------------- single placement
  }
  else if(wl0 == ":PLACE")
  {
    G4tgrVolume* vol = FindVolume(G4tgrUtils::GetString(wl[1]));
    G4tgrPlace* vpl  = vol->AddPlace(wl);
    volmgr->RegisterMe(vpl);

    //--------------------------------- parameterisation
  }
  else if(wl0 == ":PLACE_PARAM")
  {
    G4tgrVolume* vol                = FindVolume(G4tgrUtils::GetString(wl[1]));
    G4tgrPlaceParameterisation* vpl = vol->AddPlaceParam(wl);
    volmgr->RegisterMe(vpl);

    //--------------------------------- division
  }
  else if((wl0 == ":DIV_NDIV") || (wl0 == ":DIV_WIDTH") ||
          (wl0 == ":DIV_NDIV_WIDTH"))
  {
    //---------- Create G4tgrVolumeDivision and fill the volume params
    G4tgrVolumeDivision* vol = new G4tgrVolumeDivision(wl);
    volmgr->RegisterMe(vol);

    //--------------------------------- replica
  }
  else if(wl0 == ":REPL")
  {
    G4tgrVolume* vol      = FindVolume(G4tgrUtils::GetString(wl[1]));
    G4tgrPlaceDivRep* vpl = vol->AddPlaceReplica(wl);
    volmgr->RegisterMe(vpl);

    //----------------------------- assembly volume: definition of components
  }
  else if(wl0 == ":VOLU_ASSEMBLY")
  {
    G4tgrVolumeAssembly* vol = new G4tgrVolumeAssembly(wl);
    volmgr->RegisterMe(vol);

    //----------------------------- assembly volume: definition of components
  }
  else if(wl0 == ":PLACE_ASSEMBLY")
  {
    G4tgrVolume* vol = FindVolume(G4tgrUtils::GetString(wl[1]));
    G4tgrPlace* vpl  = vol->AddPlace(wl);
    volmgr->RegisterMe(vpl);

    //---------------------------------  rotation matrix
  }
  else if(wl0 == ":ROTM")
  {
    //---------- When second word is ':NEXT/:MNXT' it is used for defining a
    //           rotation matrix that will be used for the next placement/s
    G4tgrRotationMatrix* rm =
      G4tgrRotationMatrixFactory::GetInstance()->AddRotMatrix(wl);
    volmgr->RegisterMe(rm);

    //------------------------------- visualisation
  }
  else if(wl0 == ":VIS")
  {
    std::vector<G4tgrVolume*> vols =
      volmgr->FindVolumes(G4tgrUtils::GetString(wl[1]), 1);
    for(std::size_t ii = 0; ii < vols.size(); ++ii)
    {
      vols[ii]->AddVisibility(wl);
    }

    //--------------------------------- colour
  }
  else if((wl0 == ":COLOUR") || (wl0 == ":COLOR"))
  {
    std::vector<G4tgrVolume*> vols =
      volmgr->FindVolumes(G4tgrUtils::GetString(wl[1]), 1);
    for(std::size_t ii = 0; ii < vols.size(); ++ii)
    {
      vols[ii]->AddRGBColour(wl);
    }

    //--------------------------------- check overlaps
  }
  else if(wl0 == ":CHECK_OVERLAPS")
  {
    std::vector<G4tgrVolume*> vols =
      volmgr->FindVolumes(G4tgrUtils::GetString(wl[1]), 1);
    for(std::size_t ii = 0; ii < vols.size(); ++ii)
    {
      vols[ii]->AddCheckOverlaps(wl);
    }
    //--------------------------------- ERROR
  }
  else
  {
    return false;
  }

  return true;
}

// --------------------------------------------------------------------
G4tgrVolume* G4tgrLineProcessor::FindVolume(const G4String& volname)
{
  G4tgrVolume* vol = volmgr->FindVolume(volname, 1);

  if(vol->GetType() == "VOLDivision")
  {
    G4Exception("G4tgrLineProcessor::FindVolume()", "InvalidSetup",
                FatalException,
                "Using 'PLACE' for a volume created by a division !");
  }

  return vol;
}
