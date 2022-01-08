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
// G4tgrBorderSurface
//
#include "G4tgrBorderSurface.hh"
#include "G4tgrUtils.hh"

//---------------------------------------------------------------------
G4tgrBorderSurface::G4tgrBorderSurface(const std::vector<G4String>& wl)
{
  // :SURF name physVolu1:CopyNo            // :SURF name skin logVolu
  //   physVolu2:CopyNo                     //   const_key1 str1
  //   const_key1 str1                      //   array_key1 array_value1
  //   array_key1 array_value1

  theName = wl[1];
  V1Name  = wl[2].substr(0,wl[2].find(":"));
  V2Name  = wl[3].substr(0,wl[3].find(":"));
  copyNo1 = atoi(wl[2].substr(wl[2].find(":")+1).data()); // if not return 0
  copyNo2 = atoi(wl[3].substr(wl[3].find(":")+1).data());

  theOptic = new G4OpticalSurface(wl[1]);
  
  size_t i=4;
  while (i<wl.size()) {
    G4String setting = wl[i]; 
    G4String value = wl[i+1];
    G4StrUtil::to_lower(setting);
    G4StrUtil::to_lower(value);

    if (setting=="property")    {i++; break;} 
    else if (setting=="type")   theOptic->SetType(GetType(value));
    else if (setting=="model")  theOptic->SetModel(glisur);
    else if (setting=="finish") theOptic->SetFinish(GetFinish(value));
    else if (setting=="sigma_alpha")
      theOptic->SetSigmaAlpha(G4tgrUtils::GetInt(value));
    else
      G4cout<<"TextGeom: Ignore unknown surface setting "<<value<<G4endl;
    i+=2;
  }
  if (i<wl.size()) { // break while loop because of "property"
    theMPT = new G4tgrMaterialPropertiesTable(wl);
    theMPT->G4tgrPopMaterialPropertiesTable(i);
  }
}

//---------------------------------------------------------------------
G4tgrBorderSurface::~G4tgrBorderSurface()
{
}

//---------------------------------------------------------------------
G4SurfaceType G4tgrBorderSurface::GetType(G4String typeStr)
{
       if (typeStr == "dielectric_metal")      return dielectric_metal;
  else if (typeStr == "dielectric_dielectric") return dielectric_dielectric;
  else if (typeStr == "dielectric_lut")        return dielectric_LUT;
  else if (typeStr == "dielectric_lutdavis")   return dielectric_LUTDAVIS;
  else if (typeStr == "dielectric_dichroic")   return dielectric_dichroic;
  else if (typeStr == "firsov")                return firsov;
  else if (typeStr == "x_ray")                 return x_ray;
  else {
    G4cout << "TextGeom: Unknown surface type: " << typeStr 
          << " Dielectric_dielectric will be used." << G4endl;
    return dielectric_dielectric;
  }
}

//---------------------------------------------------------------------
G4OpticalSurfaceFinish G4tgrBorderSurface::GetFinish(G4String finishStr)
{
       if (finishStr == "polished")              return polished;
  else if (finishStr == "polishedfrontpainted")  return polishedfrontpainted;
  else if (finishStr == "polishedbackpainted")   return polishedbackpainted;
  else if (finishStr == "ground")                return ground;
  else if (finishStr == "groundfrontpainted")    return groundfrontpainted;
  else if (finishStr == "groundbackpainted")     return groundbackpainted;
  else if (finishStr == "polishedlumirrorair")   return polishedlumirrorair;
  else if (finishStr == "polishedlumirrorglue")  return polishedlumirrorglue;
  else if (finishStr == "polishedair")           return polishedair;
  else if (finishStr == "polishedteflonair")     return polishedteflonair;
  else if (finishStr == "polishedtioair")        return polishedtioair;
  else if (finishStr == "polishedtyvekair")      return polishedtyvekair;
  else if (finishStr == "polishedvm2000air")     return polishedvm2000air;
  else if (finishStr == "polishedvm2000glue")    return polishedvm2000glue;
  else if (finishStr == "etchedlumirrorair")     return etchedlumirrorair;
  else if (finishStr == "etchedlumirrorglue")    return etchedlumirrorglue;
  else if (finishStr == "etchedair")             return etchedair;
  else if (finishStr == "etchedteflonair")       return etchedteflonair;
  else if (finishStr == "etchedtioair")          return etchedtioair;
  else if (finishStr == "etchedtyvekair")        return etchedtyvekair;
  else if (finishStr == "etchedvm2000air")       return etchedvm2000air;
  else if (finishStr == "etchedvm2000glue")      return etchedvm2000glue;
  else if (finishStr == "groundlumirrorair")     return groundlumirrorair;
  else if (finishStr == "groundlumirrorglue")    return groundlumirrorglue;
  else if (finishStr == "groundair")             return groundair;
  else if (finishStr == "groundteflonair")       return groundteflonair;
  else if (finishStr == "groundtioair")          return groundtioair;
  else if (finishStr == "groundtyvekair")        return groundtyvekair;
  else if (finishStr == "groundvm2000air")       return groundvm2000air;
  else if (finishStr == "groundvm2000glue")      return groundvm2000glue;
  else if (finishStr == "rough_lut")             return Rough_LUT;
  else if (finishStr == "roughteflon_lut")       return RoughTeflon_LUT;
  else if (finishStr == "roughesr_lut")          return RoughESR_LUT;
  else if (finishStr == "roughesrgrease_lut")    return RoughESRGrease_LUT;
  else if (finishStr == "polished_lut")          return Polished_LUT;
  else if (finishStr == "polishedTeflon_lut")    return PolishedTeflon_LUT;
  else if (finishStr == "polishedESR_lut")       return PolishedESR_LUT;
  else if (finishStr == "polishedESRGrease_lut") return PolishedESRGrease_LUT;
  else if (finishStr == "detector_lut")          return Detector_LUT;
  else {
    G4cout << "TextGeom: Unknown surface finish: " << finishStr 
      << " Polished will be used." << G4endl;
    return polished;
  }
}

//---------------------------------------------------------------------
G4OpticalSurfaceModel G4tgrBorderSurface::GetModel(G4String modelStr)
{
       if (modelStr == "glisur")   return glisur;
  else if (modelStr == "unified")  return unified;
  else if (modelStr == "lut")      return LUT;
  else if (modelStr == "davis")    return DAVIS;
  else if (modelStr == "dichroic") return dichroic;
  else {
    G4cout << "TextGeom: Unknown surface model: " << modelStr 
      << " Glisur will be used." << G4endl;
    return glisur;
  }
}