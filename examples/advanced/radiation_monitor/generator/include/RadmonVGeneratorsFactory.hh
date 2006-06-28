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
// File name:     RadmonVGeneratorsFactory.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVGeneratorsFactory.hh,v 1.2 2006-06-28 13:53:39 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract factory for a primary generators
//

#ifndef   RADMONVGENERATORSFACTORY_HH
 #define  RADMONVGENERATORSFACTORY_HH

 // Forward declarations
 class RadmonVGenerator;
 class G4String;
 
 class RadmonVGeneratorsFactory
 {
  public:
   inline                                     RadmonVGeneratorsFactory();
   inline virtual                            ~RadmonVGeneratorsFactory();

   virtual RadmonVGenerator *                 GetGenerator(const G4String & generatorType) = 0;

  // Hidden constructors and operators
  private:
                                              RadmonVGeneratorsFactory(const RadmonVGeneratorsFactory & copy);
   RadmonVGeneratorsFactory &                 operator=(const RadmonVGeneratorsFactory & copy);
 };
 
 // Inline implementations
 #include "RadmonVGeneratorsFactory.icc"
#endif /* RADMONVGENERATORSFACTORY_HH */
