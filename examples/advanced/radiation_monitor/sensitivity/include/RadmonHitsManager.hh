//
// File name:     RadmonHitsManager.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonHitsManager.hh,v 1.1 2005-11-24 02:31:47 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Global features for hits
//

#ifndef   RADMONHITSMANAGER_HH
 #define  RADMONHITSMANAGER_HH
 
 // Include files
 #include "G4Allocator.hh"
 #include "globals.hh"
 #include <map>
 
 class RadmonHit;
 
 class RadmonHitsManager
 {
  public:
   typedef G4Allocator<RadmonHit>                Allocator;

   inline static RadmonHitsManager *             Instance();

   inline Allocator &                            HitsAllocator();
   G4int                                         ReserveIdByLabel(const G4String & label);
   
  private:
                                                RadmonHitsManager();
   inline                                      ~RadmonHitsManager();

  // Hidden constructors and operators
                                                RadmonHitsManager(const RadmonHitsManager & copy);
   RadmonHitsManager &                          operator=(const RadmonHitsManager & copy);
   
   typedef std::map<G4String, G4int>            IdLabels;
   
   static RadmonHitsManager *                   instance;

   IdLabels                                     idLabels;
   G4int                                        idCounter;
   Allocator *                                  allocator;
 };
 
 // Inline implementations
 #include "RadmonHitsManager.icc"
#endif /* RADMONHITSMANAGER_HH */
