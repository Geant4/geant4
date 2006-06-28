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
// File name:     RadmonDetectorLabelledEntitiesConstructorsFactory.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLabelledEntitiesConstructorsFactory.hh,v 1.3 2006-06-28 13:47:26 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
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
