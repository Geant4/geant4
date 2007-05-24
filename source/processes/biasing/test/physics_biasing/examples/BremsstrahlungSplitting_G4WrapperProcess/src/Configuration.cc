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
// $Id: Configuration.cc,v 1.1 2007-05-24 21:57:03 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation - built job configurations
//
#include "Configuration.hh"

#include "ConfigData.hh"
#include "G4NistManager.hh"

namespace Be_15pt18MeV_1_10_degrees
{
  void Initialise()
  {
    ConfigData::SetPrimaryEnergy(15.8*MeV);
    ConfigData::SetTargetDistance(6.31*cm);
    ConfigData::SetChamberWindowDistance(-0.9*cm);
    ConfigData::SetAirGap1Distance(-0.9051*cm);
    ConfigData::SetMonitorDistance(-2.2*cm);
    ConfigData::SetAirGap2Distance(-2.21*cm);
    ConfigData::SetBeamWindowDistance(-2.6*cm);
    ConfigData::SetBeamPipeDistance(-2.613*cm);
    ConfigData::SetBeamPipeEndDistance(-13*cm);

    ConfigData::SetTargetMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_Be"));
    ConfigData::SetChamberWindowMaterial(G4Material::GetMaterial("User_Steel"));
  }
}

namespace Be_15pt18MeV_30_90_degrees
{
  void Initialise()
  {
    Be_15pt18MeV_1_10_degrees::Initialise();
    ConfigData::SetChamberWindowMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"));
  }
}

namespace Al_15pt18MeV_1_10_degrees
{
  void Initialise()
  {
    ConfigData::SetTargetDistance(3.61*cm);
    ConfigData::SetChamberWindowDistance(-0.9*cm);
    ConfigData::SetAirGap1Distance(-0.9051*cm);
    ConfigData::SetMonitorDistance(-2.2*cm);
    ConfigData::SetAirGap2Distance(-2.21*cm);
    ConfigData::SetBeamWindowDistance(-2.6*cm);
    ConfigData::SetBeamPipeDistance(-2.613*cm);
    ConfigData::SetBeamPipeEndDistance(-13*cm);

    ConfigData::SetTargetMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_Al"));
    ConfigData::SetChamberWindowMaterial(G4Material::GetMaterial("User_Steel"));

  }
}

namespace Al_15pt18MeV_30_90_degrees
{
  void Initialise()
  {
    Al_15pt18MeV_1_10_degrees::Initialise();
    ConfigData::SetChamberWindowMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"));
  }
}

namespace Pb_15pt18MeV_1_10_degrees
{
  void Initialise()
  {
    ConfigData::SetTargetDistance(0.805*cm);
    ConfigData::SetTargetRadius(1.583*cm);
    ConfigData::SetChamberWindowDistance(-1.6*cm);
    ConfigData::SetAirGap1Distance(-1.6051*cm);
    ConfigData::SetMonitorDistance(-2.9*cm);
    ConfigData::SetAirGap2Distance(-2.91*cm);
    ConfigData::SetBeamWindowDistance(-3.3*cm);
    ConfigData::SetBeamPipeDistance(-3.313*cm);
    ConfigData::SetBeamPipeEndDistance(-13*cm);

    ConfigData::SetTargetMaterial(G4Material::GetMaterial("User_Lead"));
    ConfigData::SetChamberWindowMaterial(G4Material::GetMaterial("User_Steel"));

  }
}

namespace Pb_15pt18MeV_30_90_degrees
{
  void Initialise()
  {
    Pb_15pt18MeV_1_10_degrees::Initialise();
    ConfigData::SetChamberWindowMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"));
  }
}

/////////////////////////////////////////////

namespace Al_10pt09MeV_0_degrees
{
  void Initialise() {
    ConfigData::SetPrimaryEnergy(10.09*MeV);
    ConfigData::SetChamberWindowDistance(-0.9*cm);
    ConfigData::SetAirGap1Distance(-0.9051*cm);
    ConfigData::SetMonitorDistance(-2.2*cm);
    ConfigData::SetAirGap2Distance(-2.215*cm);
    ConfigData::SetBeamWindowDistance(-2.6*cm);
    ConfigData::SetBeamPipeDistance(-2.613*cm);
    ConfigData::SetBeamPipeEndDistance(-13*cm);
    ConfigData::SetChamberWindowMaterial(G4Material::GetMaterial("User_Steel"));

    G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    ConfigData::SetTargetMaterial(material);

    ConfigData::SetTargetDistance((6.48*g/(cm*cm))/material->GetDensity());
  }
}

namespace Al_15pt18MeV_0_degrees
{
  void Initialise() 
  {
    Al_10pt09MeV_0_degrees::Initialise();
    ConfigData::SetPrimaryEnergy(15.18*MeV);
    ConfigData::SetTargetDistance((9.73*g/(cm*cm))/G4NistManager::Instance()->FindOrBuildMaterial("G4_Al")->GetDensity());
  }
}

namespace Al_20pt28MeV_0_degrees
{
  void Initialise()
  {
    Al_10pt09MeV_0_degrees::Initialise();
    ConfigData::SetPrimaryEnergy(20.28*MeV);
    ConfigData::SetTargetDistance((11.63*g/(cm*cm))/G4NistManager::Instance()->FindOrBuildMaterial("G4_Al")->GetDensity());
  }
}

namespace Al_25pt38MeV_0_degrees
{
  void Initialise()
  {
    Al_10pt09MeV_0_degrees::Initialise();
    ConfigData::SetPrimaryEnergy(25.38*MeV);
    ConfigData::SetTargetDistance((15.14*g/(cm*cm))/G4NistManager::Instance()->FindOrBuildMaterial("G4_Al")->GetDensity());
  }
}

namespace Al_30pt45MeV_0_degrees
{
  void Initialise()
  {
    Al_10pt09MeV_0_degrees::Initialise();
    ConfigData::SetPrimaryEnergy(30.45*MeV);
    ConfigData::SetTargetDistance((16.21*g/(cm*cm))/G4NistManager::Instance()->FindOrBuildMaterial("G4_Al")->GetDensity());
  }
}

namespace Pb_10pt09MeV_0_degrees
{
  void Initialise()
  {
    ConfigData::SetTargetRadius(1.583*cm);
    ConfigData::SetChamberWindowDistance(-1.6*cm);
    ConfigData::SetAirGap1Distance(-1.6051*cm);
    ConfigData::SetMonitorDistance(-2.9*cm);
    ConfigData::SetAirGap2Distance(-2.915*cm);
    ConfigData::SetBeamWindowDistance(-3.3*cm);
    ConfigData::SetBeamPipeDistance(-3.313*cm);
    ConfigData::SetBeamPipeEndDistance(-13.3*cm);
    ConfigData::SetChamberWindowMaterial(G4Material::GetMaterial("User_Steel"));

    G4Material* material = G4Material::GetMaterial("User_Lead");
    ConfigData::SetTargetMaterial(material);

    ConfigData::SetTargetDistance((6.85*g/(cm*cm))/material->GetDensity());
  }
}

namespace Pb_15pt18MeV_0_degrees
{
  void Initialise()
  {
    Pb_10pt09MeV_0_degrees::Initialise();

    ConfigData::SetTargetDistance((9.13*g/(cm*cm))/G4Material::GetMaterial("User_Lead")->GetDensity());
  }
}

namespace Pb_20pt28MeV_0_degrees
{
  void Initialise()
  {
    Pb_10pt09MeV_0_degrees::Initialise();
    
    ConfigData::SetTargetDistance((11.43*g/(cm*cm))/G4Material::GetMaterial("User_Lead")->GetDensity());

  }
}

namespace Pb_25pt38MeV_0_degrees
{
  void Initialise()
  {
    Pb_10pt09MeV_0_degrees::Initialise();
    
    ConfigData::SetTargetDistance((11.43*g/(cm*cm))/G4Material::GetMaterial("User_Lead")->GetDensity());

  }
}

namespace Pb_30pt45MeV_0_degrees
{
  void Initialise()
  {
    Pb_10pt09MeV_0_degrees::Initialise();

    ConfigData::SetTargetDistance((13.68*g/(cm*cm))/G4Material::GetMaterial("User_Lead")->GetDensity());

  }
}
