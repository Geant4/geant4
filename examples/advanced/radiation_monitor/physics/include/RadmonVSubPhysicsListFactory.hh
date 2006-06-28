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
// File name:     RadmonVSubPhysicsListFactory.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVSubPhysicsListFactory.hh,v 1.3 2006-06-28 13:55:57 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class of a factory of sub physics lists
//

#ifndef   RADMONVSUBPHYSICSLISTFACTORY_HH
 #define  RADMONVSUBPHYSICSLISTFACTORY_HH

 // Forward declaration
 class RadmonVSubPhysicsList;
 class G4String;

 class RadmonVSubPhysicsListFactory
 {
  public:
   inline virtual                              ~RadmonVSubPhysicsListFactory();

   virtual RadmonVSubPhysicsList *              CreateSubPhysicsList(const G4String & subPhysicsListName) = 0;

  protected:
   inline                                       RadmonVSubPhysicsListFactory();

  private:
  // Hidden constructors and operators
                                                RadmonVSubPhysicsListFactory(const RadmonVSubPhysicsListFactory & copy);
   RadmonVSubPhysicsListFactory &               operator=(const RadmonVSubPhysicsListFactory & copy);
 };
 
 // Inline implementations
 #include "RadmonVSubPhysicsListFactory.icc"
#endif /* RADMONVSUBPHYSICSLISTFACTORY_HH */
