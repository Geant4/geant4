//
// File name:     RadmonDetectorEnvironmentLayout.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorEnvironmentLayout.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
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
                                              RadmonDetectorEnvironmentLayout();
                                             ~RadmonDetectorEnvironmentLayout();

   G4bool                                     IsEnabled(void) const;
   void                                       Enable(void);
   void                                       Disable(void);

   const G4String &                           GetType(void) const;
   void                                       SetType(const G4String & type);

   void                                       DumpLayout(std::ostream & out) const;

  private:
                                              RadmonDetectorEnvironmentLayout(const RadmonDetectorEnvironmentLayout & copy);
   RadmonDetectorEnvironmentLayout &          operator=(const RadmonDetectorEnvironmentLayout & copy);

   G4bool                                     enabled;
   G4String                                   environmentType;
 };
#endif /* RADMONDETECTORENVIRONMENTLAYOUT_HH */

