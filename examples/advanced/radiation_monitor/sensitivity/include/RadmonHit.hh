//
// File name:     RadmonHit.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonHit.hh,v 1.1 2005-11-24 02:31:47 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Radmon hit
//

#ifndef   RADMONHIT_HH
 #define  RADMONHIT_HH
 
 // Include files
 #include "G4VHit.hh"
 #include "G4THitsCollection.hh"
 #include "RadmonHitsManager.hh"
 #include "globals.hh"
 #include <map>
 
 class RadmonHit : public G4VHit
 {
  public:
   inline                                       RadmonHit();
   inline virtual                              ~RadmonHit();

   void *                                       operator new(size_t size);
   void                                         operator delete(void * hit);
   
   G4double                                     RetrieveById(G4int id) const;
   void                                         StoreById(G4int id, G4double value);

   inline G4double                              RetrieveByLabel(const G4String & label) const;
   inline void                                  StoreByLabel(const G4String & label, G4double value);
   
  private:
  // Hidden constructors and operators
                                                RadmonHit(const RadmonHit & copy);
   RadmonHit &                                  operator=(const RadmonHit & copy);
   
   typedef std::map<G4int, G4double>            HitData;
   
   HitData                                      hitData;
 };
 
 typedef G4THitsCollection<RadmonHit>           RadmonHitsCollection;

 // Inline implementations
 #include "RadmonHit.icc"
#endif /* RADMONHIT_HH */
