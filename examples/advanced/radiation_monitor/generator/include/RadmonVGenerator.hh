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
// File name:     RadmonVGenerator.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVGenerator.hh,v 1.3 2006/06/29 16:15:50 gunter Exp $
// Tag:           $Name: geant4-09-01 $
//
// Description:   Abstract class for a primary generator 
//

#ifndef   RADMONVGENERATOR_HH
 #define  RADMONVGENERATOR_HH

 // Forward declaration
 class G4String;
 class G4ParticleGun;
 
 class RadmonVGenerator
 {
  public:
   inline virtual                              ~RadmonVGenerator();

   virtual void                                 SetGeneratorAttribute(const G4String & attribute, const G4String & value) = 0;
   virtual void                                 ConvolveParticleGun(G4ParticleGun & gun) = 0;

  protected:
   inline                                       RadmonVGenerator();

  // Hidden constructors and operators
  private:
                                                RadmonVGenerator(const RadmonVGenerator & copy);
   RadmonVGenerator &                           operator=(const RadmonVGenerator & copy);
 };
 
 // Inline implementations
 #include "RadmonVGenerator.icc"
#endif /* RADMONVGENERATOR_HH */
