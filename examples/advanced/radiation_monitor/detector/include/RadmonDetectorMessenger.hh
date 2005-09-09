//
// File name:     RadmonDetectorMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMessenger.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   UI commands for managing a layout
//

#ifndef   RADMONDETECTORMESSENGER_HH
 #define  RADMONDETECTORMESSENGER_HH

 // Include files
 #include "G4UImessenger.hh"
 #include "G4String.hh"
 #include <map>

 // Forward declarations
 class RadmonVDetectorLayout;

 class RadmonDetectorMessenger : public G4UImessenger
 {
  public:
                                                RadmonDetectorMessenger(RadmonVDetectorLayout * layout);
                                               ~RadmonDetectorMessenger();

    virtual G4String                            GetCurrentValue(G4UIcommand * command);
    virtual void                                SetNewValue(G4UIcommand * command, G4String newValue);

  private:
                                                RadmonDetectorMessenger();
                                                RadmonDetectorMessenger(const RadmonDetectorMessenger & copy);
    RadmonDetectorMessenger &                   operator=(const RadmonDetectorMessenger & copy);

    RadmonVDetectorLayout *                     detectorLayout;
    std::map<G4String, G4String>                environmentVariables;
 };
#endif /* RADMONDETECTORMESSENGER_HH */
