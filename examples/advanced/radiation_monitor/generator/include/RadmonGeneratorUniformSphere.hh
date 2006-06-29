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
// File name:     RadmonGeneratorUniformSphere.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorUniformSphere.hh,v 1.3 2006-06-29 16:15:40 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class for a primary generator 
//

#ifndef   RADMONGENERATORUNIFORMSPHERE_HH
 #define  RADMONGENERATORUNIFORMSPHERE_HH

 // Include files
 #include "RadmonVGeneratorWithLabel.hh"
 
 // Forward declaration
 class G4String;
 class G4ParticleGun;
 
 class RadmonGeneratorUniformSphere : public RadmonVGeneratorWithLabel
 {
  public:
   inline                                       RadmonGeneratorUniformSphere();
   inline virtual                              ~RadmonGeneratorUniformSphere();

   virtual void                                 ConvolveParticleGun(G4ParticleGun & gun);
   virtual RadmonVGeneratorWithLabel *          New(void) const;

  // Hidden constructors and operators
  private:
                                                RadmonGeneratorUniformSphere(const RadmonGeneratorUniformSphere & copy);
   RadmonGeneratorUniformSphere &               operator=(const RadmonGeneratorUniformSphere & copy);
 };
 
 // Inline implementations
 #include "RadmonGeneratorUniformSphere.icc"
#endif /* RADMONGENERATORUNIFORMSPHERE_HH */
