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
// File name:     RadmonVDetectorEntitiesConstructorsFactory.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorEntitiesConstructorsFactory.hh,v 1.3 2006-06-28 13:49:51 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class of a factory of detector-entity constructor
//

#ifndef   RADMONVDETECTORENTITIESCONSTRUCTORSFACTORY_HH
 #define  RADMONVDETECTORENTITIESCONSTRUCTORSFACTORY_HH

 // Forward declaration
 class RadmonVDetectorEntityConstructor;
 class G4String;

 class RadmonVDetectorEntitiesConstructorsFactory
 {
  public:
   inline virtual                              ~RadmonVDetectorEntitiesConstructorsFactory();

   virtual RadmonVDetectorEntityConstructor *   CreateEntityConstructor(const G4String & entityName) = 0;

  protected:
   inline                                       RadmonVDetectorEntitiesConstructorsFactory();

  private:
  // Hidden constructors and operators
                                                RadmonVDetectorEntitiesConstructorsFactory(const RadmonVDetectorEntitiesConstructorsFactory & copy);
   RadmonVDetectorEntitiesConstructorsFactory & operator=(const RadmonVDetectorEntitiesConstructorsFactory & copy);
 };
 
 // Inline implementations
 #include "RadmonVDetectorEntitiesConstructorsFactory.icc"
#endif /* RADMONVDETECTORENTITIESCONSTRUCTORSFACTORY_HH */
