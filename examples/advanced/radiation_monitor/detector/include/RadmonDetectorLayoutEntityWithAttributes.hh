//
// File name:     RadmonDetectorLayoutEntityWithAttributes.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayoutEntityWithAttributes.hh,v 1.2 2005-09-12 17:14:17 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Provides attributes to other detector classes
//

#ifndef   RADMONDETECTORLAYOUTENTITYWITHATTRIBUTES_HH
 #define  RADMONDETECTORLAYOUTENTITYWITHATTRIBUTES_HH
 
 // Include files
 #include "G4String.hh"
 #include "globals.hh"

 #include <map>
 
 class RadmonDetectorLayoutEntityWithAttributes
 {
  public:
   G4String                                     GetAttribute(const G4String & attributeName, const G4String & defaultValue = G4String("")) const;
   G4bool                                       ExistsAttribute(const G4String & attributeName) const;
   void                                         SetAttribute(const G4String & attributeName, const G4String & value);
   void                                         ClearAttribute(const G4String & attributeName);
   void                                         ClearAllAttributes(void);

  protected:
   inline                                       RadmonDetectorLayoutEntityWithAttributes();
   inline                                       RadmonDetectorLayoutEntityWithAttributes(const RadmonDetectorLayoutEntityWithAttributes & copy);
   inline                                      ~RadmonDetectorLayoutEntityWithAttributes();

   void                                         DumpAttributesLayout(std::ostream & out, const G4String & indent=G4String()) const;

   inline RadmonDetectorLayoutEntityWithAttributes & operator=(const RadmonDetectorLayoutEntityWithAttributes & copy);

  private:
  // Private methods
   void                                         CopyFrom(const RadmonDetectorLayoutEntityWithAttributes & copy);
   
  // Private attributes
   typedef std::map<G4String, G4String>         AttributesMap;
   AttributesMap                                attributesMap;
 };
 
 // Inline implementations
 #include "RadmonDetectorLayoutEntityWithAttributes.icc"
#endif /* RADMONDETECTORLAYOUTENTITYWITHATTRIBUTES_HH */
