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
// File name:     RadmonGeneratorFixedParticle.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorFixedParticle.hh,v 1.2 2006-06-28 13:52:51 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class for a primary generator 
//

#ifndef   RADMONGENERATORFIXEDPARTICLE_HH
 #define  RADMONGENERATORFIXEDPARTICLE_HH

 // Include files
 #include "RadmonVGeneratorWithLabel.hh"
 
 // Forward declaration
 class G4String;
 class G4ParticleGun;
 
 class RadmonGeneratorFixedParticle : public RadmonVGeneratorWithLabel
 {
  public:
   inline                                       RadmonGeneratorFixedParticle();
   inline virtual                              ~RadmonGeneratorFixedParticle();

   virtual void                                 ConvolveParticleGun(G4ParticleGun & gun);
   virtual RadmonVGeneratorWithLabel *          New(void) const;

  // Hidden constructors and operators
  private:
                                                RadmonGeneratorFixedParticle(const RadmonGeneratorFixedParticle & copy);
   RadmonGeneratorFixedParticle &               operator=(const RadmonGeneratorFixedParticle & copy);
 };
 
 // Inline implementations
 #include "RadmonGeneratorFixedParticle.icc"
#endif /* RADMONGENERATORFIXEDPARTICLE_HH */
