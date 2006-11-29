#include "GmGeometryFromText.hh"
#include <vector>
#include <map>
//#include "GamosCore/GamosUtils/include/GmFileIn.hh"
//#include "GamosCore/GamosUtils/include/GmGenUtils.hh"
//#include "GamosCore/GamosBase/include/GmParameterMgr.hh"
#include "G4tgbVolumeMgr.hh"
#include "G4tgrMessenger.hh"


GmGeometryFromText::GmGeometryFromText()
{
  new G4tgrMessenger;
}


G4VPhysicalVolume* GmGeometryFromText::Construct()
{
  //------------------- construct g4 geometry
  std::string filename = "cmsgeom.txt";
    
  G4tgbVolumeMgr* volmgr = G4tgbVolumeMgr::GetInstance();
  volmgr->AddTextFile(filename);
  G4VPhysicalVolume* physiWorld = volmgr->ReadAndConstructDetector();

  return physiWorld;
}


