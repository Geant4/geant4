#include"G24PhysicsList.hpp"

#include "G4PhysicsListHelper.hh"
#include "G4VUserPhysicsList.hh"
#include"G4Electron.hh"
#include"G4Geantino.hh"
#include"G4ChargedGeantino.hh"
#include"G4Transportation.hh"

#include "G4StepLimiterPhysics.hh"
#include "G4UserSpecialCuts.hh"
class G4Electron;
class G4Geantino;

G24PhysicsList:: G24PhysicsList()
{
    G4StepLimiterPhysics* stepLimitPhys = new G4StepLimiterPhysics();
    RegisterPhysics(stepLimitPhys);
    
    

    
}

G24PhysicsList:: ~G24PhysicsList()
{}

void G24PhysicsList::ConstructParticle()
{
    G4VModularPhysicsList::ConstructParticle();
    G4Electron::ElectronDefinition();
   
    //G4Geantino::GeantinoDefinition();
    //G4ChargedGeantino::ChargedGeantinoDefinition();
}

void G24PhysicsList::ConstructProcess()
{
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    G4VModularPhysicsList::ConstructProcess();
   AddTransportation();  
   ph->RegisterProcess(new G4UserSpecialCuts(),G4Electron::ElectronDefinition());
}