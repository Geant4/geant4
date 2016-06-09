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
// File name:     RadmonGeneratorsWithLabelFactory.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorsWithLabelFactory.hh,v 1.3 2006/06/29 16:15:44 gunter Exp $
// Tag:           $Name: geant4-08-02 $
//
// Description:   Factory of primary generators with a label
//

#ifndef   RADMONGENERATORSWITHLABELFACTORY_HH
 #define  RADMONGENERATORSWITHLABELFACTORY_HH
 
 // Include files
 #include "RadmonVGeneratorsFactory.hh"
 #include <list>
 
 // Forward declarations
 class RadmonVGeneratorWithLabel;

 class RadmonGeneratorsWithLabelFactory : public RadmonVGeneratorsFactory
 {
  public:
                                                RadmonGeneratorsWithLabelFactory();
   virtual                                     ~RadmonGeneratorsWithLabelFactory();

   virtual RadmonVGenerator *                   GetGenerator(const G4String & generatorType);

   void                                         AppendGenerator(RadmonVGeneratorWithLabel * generator);

  // Hidden constructors and operators
  private:
                                                RadmonGeneratorsWithLabelFactory(const RadmonGeneratorsWithLabelFactory & copy);
   RadmonGeneratorsWithLabelFactory &           operator=(const RadmonGeneratorsWithLabelFactory & copy);

  // Private data types
   typedef std::list<RadmonVGeneratorWithLabel *> GeneratorsList;

  // Private attributes
   GeneratorsList                                generatorsList;
 };
#endif /* RADMONGENERATORSWITHLABELFACTORY_HH */
