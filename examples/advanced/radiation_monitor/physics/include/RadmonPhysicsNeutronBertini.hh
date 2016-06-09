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
// File name:     RadmonPhysicsNeutronBertini.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsNeutronBertini.hh,v 1.3 2006/06/29 16:17:37 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//
// Description:   Neutron processes based on bertini cascade model
//

#ifndef   RADMONPHYSICSNEUTRONBERTINI_HH
 #define  RADMONPHYSICSNEUTRONBERTINI_HH
 
 // Include files
 #include "RadmonVSubPhysicsListWithLabel.hh"
 #include "RadmonPhysicsInfoList.hh"
 
 class RadmonPhysicsNeutronBertini : public RadmonVSubPhysicsListWithLabel
 {
  public:
   inline                                       RadmonPhysicsNeutronBertini();
   inline virtual                              ~RadmonPhysicsNeutronBertini();

   virtual RadmonVSubPhysicsListWithLabel *     New(void) const;

   virtual void                                 ConstructParticle(void);
   virtual void                                 ConstructProcess(void);
   virtual void                                 SetCuts(void);
   
   virtual const RadmonPhysicsInfoList &        Provides(void) const;

  protected:
   
  private:
  // Hidden constructors and operators
                                                RadmonPhysicsNeutronBertini(const RadmonPhysicsNeutronBertini & copy);
   RadmonPhysicsNeutronBertini &                operator=(const RadmonPhysicsNeutronBertini & copy);
   
   mutable RadmonPhysicsInfoList                infoList;
 };
 
 #include "RadmonPhysicsNeutronBertini.icc"
#endif /* RADMONPHYSICSNEUTRONBERTINI_HH */
