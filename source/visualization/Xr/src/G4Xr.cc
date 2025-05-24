#include "G4Xr.hh"
#include "G4XrSceneHandler.hh"

#include "G4UIQt.hh"
#include "G4UIbatch.hh"
#include "G4UImanager.hh"

G4Xr::G4Xr()
  : G4VGraphicsSystem("Xr", "Xr", "Web delivery of XR file", G4VGraphicsSystem::noFunctionality)
{
  std::cout << "G4Xr::G4Xr()" << std::endl;
}

G4VSceneHandler* G4Xr::CreateSceneHandler(const G4String& name)
{
  std::cout << "G4Xr::CreateSceneHandler(name=" << name << ")" << std::endl;
  return new G4XrSceneHandler(*this, name);
}

G4VViewer* G4Xr::CreateViewer(G4VSceneHandler& scene, const G4String& name)
{
  std::cout << "G4Xr::CreateViewer(scene=" << &scene << ")" << std::endl;
  return new G4XrViewer(scene, name);
}

G4bool G4Xr::IsUISessionCompatible() const
{
  // Qt windows require a Qt session.
  G4UIsession* baseSession = G4UImanager::GetUIpointer()->GetBaseSession();
  return dynamic_cast<G4UIQt*>(baseSession) != nullptr;
}