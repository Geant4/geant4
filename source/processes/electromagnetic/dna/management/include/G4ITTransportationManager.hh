#ifndef G4ITTRANSPORTATIONMANAGER_HH
#define G4ITTRANSPORTATIONMANAGER_HH


#include "globals.hh"

class G4ITNavigator;

class G4ITTransportationManager
{
public:
    G4ITTransportationManager();
    static void DeleteInstance();

    static G4ITTransportationManager* GetTransportationManager();

    G4ITNavigator* GetNavigatorForTracking();
private:
    ~G4ITTransportationManager();
    static G4ITTransportationManager* fpInstance;
    void Initialize();
    G4ITNavigator* fpNavigator;
};

#endif // G4ITTRANSPORTATIONMANAGER_HH
