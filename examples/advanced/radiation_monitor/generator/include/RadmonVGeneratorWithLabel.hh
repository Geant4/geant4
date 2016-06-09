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
// File name:     RadmonVGeneratorWithLabel.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVGeneratorWithLabel.hh,v 1.3 2006/06/29 16:15:58 gunter Exp $
// Tag:           $Name: geant4-08-02 $
//
// Description:   Primary generators with a label abstract class
//

#ifndef   RADMONVGENERATORWITHLABEL_HH
 #define  RADMONVGENERATORWITHLABEL_HH

 // Include files
 #include "RadmonVGenerator.hh"
 #include "RadmonLayoutEntityWithAttributes.hh"

 class RadmonVGeneratorWithLabel : public RadmonVGenerator, public RadmonLayoutEntityWithAttributes
 {
  public:
   inline virtual                              ~RadmonVGeneratorWithLabel();

   inline const G4String &                      GetLabel(void) const;

   inline virtual void                          SetGeneratorAttribute(const G4String & attribute, const G4String & value);

   virtual RadmonVGeneratorWithLabel *          New(void) const = 0;

  protected:
   inline                                       RadmonVGeneratorWithLabel(const G4String & label);

  // Hidden constructors and operators
  private:
                                                RadmonVGeneratorWithLabel(const RadmonVGeneratorWithLabel & copy);
   RadmonVGeneratorWithLabel &                   operator=(const RadmonVGeneratorWithLabel & copy);

  // Private attributes
   G4String                                     generatorLabel;
 };
 
 #include "RadmonVGeneratorWithLabel.icc"   
#endif /* RADMONVGENERATORWITHLABEL_HH */
