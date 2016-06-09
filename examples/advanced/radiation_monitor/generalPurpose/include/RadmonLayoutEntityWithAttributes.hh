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
// File name:     RadmonLayoutEntityWithAttributes.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonLayoutEntityWithAttributes.hh,v 1.3.2.2 2006/06/29 16:14:19 gunter Exp $
// Tag:           $Name: geant4-09-01 $
//
// Description:   Provides attributes to other detector classes
//

#ifndef   RADMONLAYOUTENTITYWITHATTRIBUTES_HH
 #define  RADMONLAYOUTENTITYWITHATTRIBUTES_HH
 
 // Include files
 #include "G4String.hh"
 #include "globals.hh"
 #include "G4ThreeVector.hh"
 #include "G4RotationMatrix.hh"

 #include <vector>
 #include <utility>
 
 class RadmonLayoutEntityWithAttributes
 {
  public:
   G4int                                        GetNAttributes(void) const;
   const G4String &                             GetAttributeName(G4int index) const;

   G4String                                     GetAttribute(const G4String & attributeName, const G4String & defaultValue = G4String("")) const;
   G4bool                                       ExistsAttribute(const G4String & attributeName) const;
   void                                         SetAttribute(const G4String & attributeName, const G4String & value);
   void                                         ClearAttribute(const G4String & attributeName);
   void                                         ClearAllAttributes(void);

   G4double                                     GetAttributeAsDouble(const G4String & attributeName, double defaultValue) const;
   G4double                                     GetAttributeAsMeasure(const G4String & attributeName, const char * category, double defaultValue) const;
   G4int                                        GetAttributeAsInteger(const G4String & attributeName, G4int defaultValue) const;
   G4ThreeVector                                GetAttributeAsThreeVector(const G4String & attributeName, const G4ThreeVector & defaultValue) const;
   G4ThreeVector                                GetAttributeAsThreeVectorWithMeasure(const G4String & attributeName, const char * category, const G4ThreeVector & defaultValue) const;
   G4ThreeVector                                GetAttributeAsDirection(const G4String & attributeName, const G4ThreeVector & defaultValue) const;
   G4RotationMatrix                             GetAttributeAsRotationMatrix(const G4String & attributeName, const G4RotationMatrix & defaultValue) const;
   
  protected:
   inline                                       RadmonLayoutEntityWithAttributes();
   inline                                       RadmonLayoutEntityWithAttributes(const RadmonLayoutEntityWithAttributes & copy);
   inline                                      ~RadmonLayoutEntityWithAttributes();

   void                                         DumpAttributesLayout(std::ostream & out, const G4String & indent=G4String()) const;

   inline RadmonLayoutEntityWithAttributes &    operator=(const RadmonLayoutEntityWithAttributes & copy);

  private:
  // Private methods
   void                                         CopyFrom(const RadmonLayoutEntityWithAttributes & copy);
   
  // Private attributes
   typedef std::pair<G4String, G4String>        AttributeItem;
   typedef std::vector<AttributeItem>           AttributesVector;
   AttributesVector                             attributesVector;
 };
 
 // Inline implementations
 #include "RadmonLayoutEntityWithAttributes.icc"
#endif /* RADMONLAYOUTENTITYWITHATTRIBUTES_HH */
