//
// File name:     RadmonApplicationDetectorSetup.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationDetectorSetup.hh,v 1.1 2005-10-25 16:39:12 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon application detector entity constructors setup
//

#ifndef   RADMONAPPLICATIONDETECTORSETUP_HH
 #define  RADMONAPPLICATIONDETECTORSETUP_HH

 // Inglude files
 #include "globals.hh"

 // Forward declarations
 class RadmonApplicationOptions;
 class RadmonDetectorLabelledEntitiesConstructorsFactory;

 class RadmonApplicationDetectorSetup
 {
  public:
   inline                                       RadmonApplicationDetectorSetup(const RadmonApplicationOptions & options);
   inline                                      ~RadmonApplicationDetectorSetup();

  protected:
   G4bool                                       CreateDetectorEntityConstructors(RadmonDetectorLabelledEntitiesConstructorsFactory * factory);

  // Hidden constructors and operators
                                                RadmonApplicationDetectorSetup();
                                                RadmonApplicationDetectorSetup(const RadmonApplicationDetectorSetup & copy);
   RadmonApplicationDetectorSetup &             operator=(const RadmonApplicationDetectorSetup & copy);
   
  // Private attributes
   const RadmonApplicationOptions &             currentOptions;
 };
 
 #include "RadmonApplicationDetectorSetup.icc"
#endif /* RADMONAPPLICATIONDETECTORSETUP_HH */
