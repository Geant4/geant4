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
// File name:     RadmonVSubPhysicsListWithLabel.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVSubPhysicsListWithLabel.hh,v 1.4 2006/06/29 16:18:35 gunter Exp $
// Tag:           $Name: geant4-08-02 $
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
