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
// File name:     RadmonVDetectorEntityConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorEntityConstructor.hh,v 1.3 2006-06-28 13:49:55 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class of a detector-entity constructor
//

#ifndef   RADMONVDETECTORENTITYCONSTRUCTOR_HH
 #define  RADMONVDETECTORENTITYCONSTRUCTOR_HH
 
 // Forward declarations
 class G4LogicalVolume;
 class G4String;
 
 class RadmonVDetectorEntityConstructor
 {
  public:
   inline virtual                              ~RadmonVDetectorEntityConstructor();
    
   virtual void                                 SetEntityAttribute(const G4String & attributeName, const G4String &value) = 0;
   virtual G4LogicalVolume *                    ConstructLogicalVolume() = 0;

  protected:
   inline                                       RadmonVDetectorEntityConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonVDetectorEntityConstructor(const RadmonVDetectorEntityConstructor & copy);
   RadmonVDetectorEntityConstructor &           operator=(const RadmonVDetectorEntityConstructor & copy);

 };
 
 // Inline implementations
 #include "RadmonVDetectorEntityConstructor.icc"
#endif /* RADMONVDETECTORENTITYCONSTRUCTOR_HH */
