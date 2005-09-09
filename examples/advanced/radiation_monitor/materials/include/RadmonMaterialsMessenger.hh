//
// File name:     RadmonMaterialsMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonMaterialsMessenger.hh,v 1.1 2005-09-09 08:27:13 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   UI commands for managing materials
//

#ifndef   RADMONMATERIALMESSENGER_HH
 #define  RADMONMATERIALMESSENGER_HH

 // Include files
 #include "G4UImessenger.h"

 class RadmonMaterialMessenger : public G4UImessenger
 {
  public:
                                                RadmonMaterialMessenger();
                                               ~RadmonMaterialMessenger();

    virtual G4String                            GetCurrentValue(G4UIcommand * command);
    virtual void                                SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonMaterialMessenger(const RadmonMaterialMessenger & copy);
    RadmonMaterialMessenger &                   operator=(const RadmonMaterialMessenger & copy);
 };
#endif /* RADMONMATERIALMESSENGER_HH */
