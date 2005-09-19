//
// File name:     RadmonMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonMessenger.hh,v 1.1 2005-09-19 19:39:43 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   UI commands for managing a layout
//

#ifndef   RADMONMESSENGER_HH
 #define  RADMONMESSENGER_HH

 // Include files
 #include "G4UImessenger.hh"
 #include "G4UIdirectory.hh"
 #include "RadmonMessengersSupport.hh"

 class RadmonMessenger : public G4UImessenger
 {
  public:
   static G4bool                                ProcessArguments(const G4String & rawArguments, G4int nArgs, G4String * arguments);
   static G4double                              GetUnit(const G4String & unitStr, const char * cathegory);
   
  protected:
                                                RadmonMessenger(const char * path, const char * guidance=0);
   inline                                      ~RadmonMessenger();

   std::istream *                               OpenForInput(const G4String & fileName) const;
   std::ostream *                               OpenForOutput(const G4String & fileName) const;

  private:
  // Hidden constructors and operators
                                                RadmonMessenger();
                                                RadmonMessenger(const RadmonMessenger & copy);
   RadmonMessenger &                            operator=(const RadmonMessenger & copy);

  // Private Attributes
   G4UIdirectory *                              directory;
 };

 // Inline implementation
 #include "RadmonMessenger.icc"
#endif /* RADMONMESSENGER_HH */
