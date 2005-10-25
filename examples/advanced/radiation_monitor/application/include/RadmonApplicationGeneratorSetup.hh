//
// File name:     RadmonApplicationGeneratorSetup.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationGeneratorSetup.hh,v 1.1 2005-10-25 16:39:12 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon application generators algoritms setup
//

#ifndef   RADMONAPPLICATIONGENERATORSETUP_HH
 #define  RADMONAPPLICATIONGENERATORSETUP_HH

 // Inglude files
 #include "globals.hh"

 // Forward declarations
 class RadmonApplicationOptions;
 class RadmonGeneratorsWithLabelFactory;

 class RadmonApplicationGeneratorSetup
 {
  public:
   inline                                       RadmonApplicationGeneratorSetup(const RadmonApplicationOptions & options);
   inline                                      ~RadmonApplicationGeneratorSetup();

  protected:
   G4bool                                       CreateGenerators(RadmonGeneratorsWithLabelFactory * factory);

  // Hidden constructors and operators
                                                RadmonApplicationGeneratorSetup();
                                                RadmonApplicationGeneratorSetup(const RadmonApplicationGeneratorSetup & copy);
   RadmonApplicationGeneratorSetup &            operator=(const RadmonApplicationGeneratorSetup & copy);
   
  // Private attributes
   const RadmonApplicationOptions &             currentOptions;
 };
 
 #include "RadmonApplicationGeneratorSetup.icc"
#endif /* RADMONAPPLICATIONGENERATORSETUP_HH */
