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
// File name:     RadmonMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonMessenger.hh,v 1.2 2006-06-28 13:52:02 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
