//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// File name:     RadmonApplicationOptions.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationOptions.hh,v 1.2 2006-06-28 13:45:20 gunter Exp $
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
