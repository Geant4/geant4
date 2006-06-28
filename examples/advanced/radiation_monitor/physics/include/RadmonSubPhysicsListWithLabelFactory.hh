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
// File name:     RadmonSubPhysicsListWithLabelFactory.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSubPhysicsListWithLabelFactory.hh,v 1.2 2006-06-28 13:55:43 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Concrete factory that constructs 
//                RadmonSubPhysiscsListWithLabel objects
//

#ifndef   RADMONSUBPHYSICSLISTWITHLABELFACTORY_HH
 #define  RADMONSUBPHYSICSLISTWITHLABELFACTORY_HH

 // Include files
 #include "RadmonVSubPhysicsListFactory.hh"
 #include <list>

 // Forward declarations
 class RadmonVSubPhysicsList;
 class RadmonVSubPhysicsListWithLabel;

 class RadmonSubPhysicsListWithLabelFactory : public RadmonVSubPhysicsListFactory
 {
  public:
   inline                                       RadmonSubPhysicsListWithLabelFactory();
   virtual                                     ~RadmonSubPhysicsListWithLabelFactory();

   virtual RadmonVSubPhysicsList *              CreateSubPhysicsList(const G4String & subPhysicsListName);

   void                                         AppendSubPhysicsListWithLabel(RadmonVSubPhysicsListWithLabel * physicsList);

  private:
  // Hidden constructors and operators
                                                RadmonSubPhysicsListWithLabelFactory(const RadmonSubPhysicsListWithLabelFactory & copy);
   RadmonSubPhysicsListWithLabelFactory &       operator=(const RadmonSubPhysicsListWithLabelFactory & copy);

  // Private attributes
   typedef std::list<RadmonVSubPhysicsListWithLabel *> SubPhysicsLists;
   SubPhysicsLists                              subPhysicsLists;
 };
 
 // Inline implementations
 #include "RadmonSubPhysicsListWithLabelFactory.icc"
#endif /* RADMONSUBPHYSICSLISTWITHLABELFACTORY_HH */
