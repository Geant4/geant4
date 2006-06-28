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
// File name:     RadmonGeneratorsWithLabelFactory.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorsWithLabelFactory.hh,v 1.2 2006-06-28 13:53:21 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
