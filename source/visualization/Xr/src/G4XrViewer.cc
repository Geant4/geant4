#include "G4XrViewer.hh"
#include "G4VSceneHandler.hh"

// tinygltf
// Define these only in *one* .cc file.
#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
// #define TINYGLTF_NOEXCEPTION // optional. disable exception handling.

G4XrViewer::G4XrViewer(G4VSceneHandler& sceneHandler, const G4String& name)
  : G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name)
{
  // Set default and current view parameters
  fVP.SetAutoRefresh(true);
  fDefaultVP.SetAutoRefresh(true);
}

void G4XrViewer::Initialise()
{
  std::cout << "G4XrViewer::Initialise()" << std::endl;

  httplib::Server svr;
  svr.Get("/hi", [](const httplib::Request &, httplib::Response &res) {
    res.set_content("Hello World!", "text/plain");
  });

  svr.listen("0.0.0.0", 8080);
}

G4XrViewer::~G4XrViewer()
{
}

void G4XrViewer::SetView()
{
}

void G4XrViewer::DrawView()
{
  NeedKernelVisit();

  ProcessView();

  FinishView();
}

void G4XrViewer::ShowView() {
}

void G4XrViewer::ClearView()
{
}

void G4XrViewer::FinishView()
{
}
