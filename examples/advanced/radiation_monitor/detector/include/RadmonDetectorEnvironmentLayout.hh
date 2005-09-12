//
// File name:     RadmonDetectorEnvironmentLayout.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorEnvironmentLayout.hh,v 1.2 2005-09-12 17:14:17 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Internal class to describe the environment layout
//

#ifndef   RADMONDETECTORENVIRONMENTLAYOUT_HH
 #define  RADMONDETECTORENVIRONMENTLAYOUT_HH
 
 // Include files
 #include "RadmonDetectorLayoutEntityWithAttributes.hh"
 
 class RadmonDetectorEnvironmentLayout : public RadmonDetectorLayoutEntityWithAttributes
 {
  public:
   inline                                     RadmonDetectorEnvironmentLayout();
   inline                                    ~RadmonDetectorEnvironmentLayout();

   inline G4bool                              IsEnabled(void) const;
   inline void                                Enable(void);
   inline void                                Disable(void);

   inline const G4String &                    GetType(void) const;
   inline void                                SetType(const G4String & type);

   void                                       DumpLayout(std::ostream & out, const G4String & indent=G4String("")) const;

  private:
                                              RadmonDetectorEnvironmentLayout(const RadmonDetectorEnvironmentLayout & copy);
   RadmonDetectorEnvironmentLayout &          operator=(const RadmonDetectorEnvironmentLayout & copy);

   G4bool                                     enabled;
   G4String                                   environmentType;
 };

 // Inline implementations
 #include "RadmonDetectorEnvironmentLayout.icc"
#endif /* RADMONDETECTORENVIRONMENTLAYOUT_HH */

