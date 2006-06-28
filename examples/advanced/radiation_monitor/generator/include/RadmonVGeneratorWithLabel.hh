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
// File name:     RadmonVGeneratorWithLabel.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVGeneratorWithLabel.hh,v 1.2 2006-06-28 13:53:35 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
