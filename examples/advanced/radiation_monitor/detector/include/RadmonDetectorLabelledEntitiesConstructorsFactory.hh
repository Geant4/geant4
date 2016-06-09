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
// File name:     RadmonDetectorLabelledEntitiesConstructorsFactory.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLabelledEntitiesConstructorsFactory.hh,v 1.4 2006/06/29 16:10:07 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//
// Description:   Concrete factory that constructs 
//                RadmonVDetectorLabelledEntityConstructor objects
//

#ifndef   RADMONDETECTORLABELLEDENTITIESCONSTRUCTORSFACTORY_HH
 #define  RADMONDETECTORLABELLEDENTITIESCONSTRUCTORSFACTORY_HH

 // Include files
 #include "RadmonVDetectorEntitiesConstructorsFactory.hh"
 #include <list>

 // Forward declarations
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonVDetectorEntityConstructor;

 class RadmonDetectorLabelledEntitiesConstructorsFactory : public RadmonVDetectorEntitiesConstructorsFactory
 {
  public:
   inline                                       RadmonDetectorLabelledEntitiesConstructorsFactory();
   virtual                                     ~RadmonDetectorLabelledEntitiesConstructorsFactory();

   virtual RadmonVDetectorEntityConstructor *   CreateEntityConstructor(const G4String & entityName);

   void                                         AppendLabelledEntityConstructor(RadmonVDetectorLabelledEntityConstructor * constructor);

  private:
  // Hidden constructors and operators
                                                RadmonDetectorLabelledEntitiesConstructorsFactory(const RadmonDetectorLabelledEntitiesConstructorsFactory & copy);
   RadmonDetectorLabelledEntitiesConstructorsFactory & operator=(const RadmonDetectorLabelledEntitiesConstructorsFactory & copy);

  // Private attributes
   typedef std::list<RadmonVDetectorLabelledEntityConstructor *> EntityConstructorsList;
   EntityConstructorsList                       entityConstructorsList;
 };
 
 // Inline implementations
 #include "RadmonDetectorLabelledEntitiesConstructorsFactory.icc"
#endif /* RADMONDETECTORLABELLEDENTITIESCONSTRUCTORSFACTORY_HH */
