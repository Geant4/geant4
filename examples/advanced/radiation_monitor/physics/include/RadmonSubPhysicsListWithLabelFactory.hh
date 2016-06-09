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
// File name:     RadmonSubPhysicsListWithLabelFactory.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSubPhysicsListWithLabelFactory.hh,v 1.3 2006/06/29 16:18:17 gunter Exp $
// Tag:           $Name: geant4-09-02 $
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
