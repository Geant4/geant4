//
// File name:     RadmonApplicationOptions.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationOptions.hh,v 1.1 2005-09-09 08:26:54 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Supporting class for application startup
//

#ifndef   RADMONAPPLICATIONOPTIONS_HH
 #define  RADMONAPPLICATIONOPTIONS_HH

 // Include files
 #include "globals.hh"
 #include "G4String.hh"

 class RadmonApplicationOptions
 {
  public:
                                                RadmonApplicationOptions(int argc, char * * argv);
   inline                                      ~RadmonApplicationOptions();

   inline G4bool                                Valid(void) const;

   inline G4bool                                Interactive(void) const;
   inline G4bool                                Help(void) const;
   inline G4bool                                Verbose(void) const;
   inline const char *                          ApplicationName(void) const;
   inline const char *                          FileName(void) const;
   inline const char *                          StartupFileName(void) const;
   
   void                                         DumpHelp(void) const;

  private:
  // Hidden constructors and operators
                                                RadmonApplicationOptions(const RadmonApplicationOptions & copy);
   RadmonApplicationOptions &                   operator=(const RadmonApplicationOptions & copy);

  // Private attributes
   G4bool                                       valid;

   G4bool                                       interactive;
   G4bool                                       verbose;
   G4bool                                       help;

   const char *                                 applicationName;
   const char *                                 fileName;
   const char *                                 startupFileName;
 };
 
 #include "RadmonApplicationOptions.icc"
#endif /* RADMONAPPLICATIONOPTIONS_HH */
