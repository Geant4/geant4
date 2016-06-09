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
// File name:     RadmonSubPhysicsListLayout.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSubPhysicsListLayout.hh,v 1.5 2006/06/29 16:18:13 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//
// Description:   Internal class that describes sub physics list informations
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
   G4String                                     subPhysicsListLabel;
 };

 // Inline implementations
 #include "RadmonSubPhysicsListLayout.icc"
#endif /* RADMONSUBPHYSICSLISTLAYOUT_HH */
