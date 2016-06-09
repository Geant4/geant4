//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// File name:     RadmonApplicationOptions.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationOptions.hh,v 1.3 2006/06/29 16:08:17 gunter Exp $
// Tag:           $Name: geant4-08-02 $
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
