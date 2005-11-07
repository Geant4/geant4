//
// File name:     RadmonSubPhysicsListLayout.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSubPhysicsListLayout.hh,v 1.1 2005-11-07 17:52:36 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Internal class that describes layer informations
//

#ifndef   RADMONSUBPHYSICSLISTLAYOUT_HH
 #define  RADMONSUBPHYSICSLISTLAYOUT_HH

 // Include files
 #include "RadmonLayoutEntityWithAttributes.hh"

 class RadmonSubPhysicsListLayout : public RadmonLayoutEntityWithAttributes
 {
  public:
   inline                                       RadmonSubPhysicsListLayout();
   inline                                       RadmonSubPhysicsListLayout(const RadmonSubPhysicsListLayout & copy);
   inline                                      ~RadmonSubPhysicsListLayout();

   inline RadmonSubPhysicsListLayout &          operator=(const RadmonSubPhysicsListLayout & copy);

   inline const G4String &                      GetLabel(void) const;
   inline void                                  SetLabel(const G4String & label);

   void                                         DumpLayout(std::ostream & out, const G4String &indent=G4String()) const;

  private:
  // Private attributes
   G4String                                     layerLabel;
 };

 // Inline implementations
 #include "RadmonSubPhysicsListLayout.icc"
#endif /* RADMONSUBPHYSICSLISTLAYOUT_HH */
