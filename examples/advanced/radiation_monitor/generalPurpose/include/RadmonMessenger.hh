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
// File name:     RadmonMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonMessenger.hh,v 1.3 2006/06/29 16:14:23 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//
// Description:   UI commands for managing a layout
//

#ifndef   RADMONMESSENGER_HH
 #define  RADMONMESSENGER_HH

 // Include files
 #include "G4UImessenger.hh"
 #include "G4UIdirectory.hh"
 #include "RadmonMessengersSupport.hh"

 class RadmonMessenger : public G4UImessenger
 {
  public:
   static G4bool                                ProcessArguments(const G4String & rawArguments, G4int nArgs, G4String * arguments);
   static G4double                              GetUnit(const G4String & unitStr, const char * cathegory);
   
  protected:
                                                RadmonMessenger(const char * path, const char * guidance=0);
   inline                                      ~RadmonMessenger();

   std::istream *                               OpenForInput(const G4String & fileName) const;
   std::ostream *                               OpenForOutput(const G4String & fileName) const;

  private:
  // Hidden constructors and operators
                                                RadmonMessenger();
                                                RadmonMessenger(const RadmonMessenger & copy);
   RadmonMessenger &                            operator=(const RadmonMessenger & copy);

  // Private Attributes
   G4UIdirectory *                              directory;
 };

 // Inline implementation
 #include "RadmonMessenger.icc"
#endif /* RADMONMESSENGER_HH */
