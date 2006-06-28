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
// File name:     RadmonSubPhysicsListLayout.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSubPhysicsListLayout.hh,v 1.4 2006-06-28 13:55:39 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
