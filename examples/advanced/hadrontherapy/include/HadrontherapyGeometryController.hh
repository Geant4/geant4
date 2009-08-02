#ifndef HadrontherapyGeometryController_hh
#define HadrontherapyGeometryController_hh 1

#include "globals.hh"
#include "G4String.hh"
#include "G4VUserDetectorConstruction.hh"

/**
 * Controller for geometry selection
 *
 * This controller is called by the geometry messenger and used to
 * select the geometry. Each available geometry must have unique name
 * and it must be known by the geometry controller.
 */
class HadrontherapyGeometryController
{
public:
  HadrontherapyGeometryController();
  ~HadrontherapyGeometryController();

  /**
   * Select a geometry by name.
   */
  void SetGeometry(G4String);

private:
  void registerGeometry(G4VUserDetectorConstruction *detector);
};

#endif
