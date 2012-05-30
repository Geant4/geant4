#include "G4ITTransportationManager.hh"
#include "G4TransportationManager.hh"
#include "G4ITNavigator.hh"

G4ITTransportationManager* G4ITTransportationManager::fpInstance (0);

G4ITTransportationManager::G4ITTransportationManager()
{
    Initialize();
}

G4ITTransportationManager::~G4ITTransportationManager()
{
    if(fpNavigator) delete fpNavigator;
}

void G4ITTransportationManager::Initialize()
{
    fpNavigator = new G4ITNavigator();
    fpNavigator->Activate(true);
    G4Navigator* navForTracking = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    G4VPhysicalVolume* world = navForTracking->GetWorldVolume();
    fpNavigator->SetWorldVolume(world);
}

G4ITTransportationManager* G4ITTransportationManager::GetTransportationManager()
{
    if(fpInstance == 0) fpInstance = new G4ITTransportationManager;
    return fpInstance;
}

G4ITNavigator* G4ITTransportationManager::GetNavigatorForTracking()
{
    return fpNavigator;
}
