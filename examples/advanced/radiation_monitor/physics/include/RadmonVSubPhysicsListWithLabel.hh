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
// File name:     RadmonVSubPhysicsListWithLabel.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVSubPhysicsListWithLabel.hh,v 1.3 2006-06-28 13:56:01 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class of a physics list piece with label
//                and attributes
//

#ifndef   RADMONVSUBPHYSICSLISTWITHLABEL_HH
 #define  RADMONVSUBPHYSICSLISTWITHLABEL_HH
 
 // Include files
 #include "RadmonLayoutEntityWithAttributes.hh"
 #include "RadmonVSubPhysicsList.hh"
 
 class RadmonVSubPhysicsListWithLabel : public RadmonVSubPhysicsList, public RadmonLayoutEntityWithAttributes
 {
  public:
   inline virtual                              ~RadmonVSubPhysicsListWithLabel();

   inline const G4String &                      GetLabel(void) const;
   inline virtual void                          SetPhysicsListAttribute(const G4String & attributeName, const G4String & value);

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const = 0;

  protected:
   inline                                       RadmonVSubPhysicsListWithLabel(const G4String & label);
   
  private:
  // Hidden constructors and operators
                                                RadmonVSubPhysicsListWithLabel();
                                                RadmonVSubPhysicsListWithLabel(const RadmonVSubPhysicsListWithLabel & copy);
   RadmonVSubPhysicsListWithLabel &             operator=(const RadmonVSubPhysicsListWithLabel & copy);

  // Private attributes
   G4String                                     physiscListLabel;
 };
 
 #include "RadmonVSubPhysicsListWithLabel.icc"
#endif /* RADMONVSUBPHYSICSLISTWITHLABEL_HH */
