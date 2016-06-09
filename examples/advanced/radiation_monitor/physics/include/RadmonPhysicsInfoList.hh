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
// File name:     RadmonPhysicsInfoList.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsInfoList.hh,v 1.3 2006/06/29 16:17:23 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//
// Description:   Provides informations about a process
//

#ifndef   RADMONPHYSICSINFOLIST_HH
 #define  RADMONPHYSICSINFOLIST_HH
 
 // Include files
 #include "RadmonPhysicsInfo.hh"
 #include "globals.hh"
 #include <vector>
 
 class RadmonPhysicsInfoList
 {
  public:
   inline                                       RadmonPhysicsInfoList();
                                                RadmonPhysicsInfoList(const RadmonPhysicsInfoList & copy);
   inline                                      ~RadmonPhysicsInfoList();
   
   RadmonPhysicsInfoList &                      operator=(const RadmonPhysicsInfoList & copy);
   
   G4bool                                       CollidesWith(const RadmonPhysicsInfoList & other) const;
   
   void                                         InsertPhysicsInfo(const RadmonPhysicsInfo & info);
   
   inline G4int                                 GetNPhysicsInfos(void) const;
   inline const RadmonPhysicsInfo &             GetPhysicsInfo(G4int index) const;
   
  private:
  // Private data types
   typedef std::vector<RadmonPhysicsInfo>       InfoVector;

  // Private attributes
   InfoVector                                   infoVector;
 };
 
 std::ostream &                                 operator<<(std::ostream & out, const RadmonPhysicsInfoList & infoList);
 
 // Inline implementations
 #include "RadmonPhysicsInfoList.icc"
#endif /* RADMONPHYSICSINFOLIST_HH */
