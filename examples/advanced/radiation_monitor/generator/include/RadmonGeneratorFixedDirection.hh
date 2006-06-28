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
// File name:     RadmonGeneratorFixedDirection.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorFixedDirection.hh,v 1.2 2006-06-28 13:52:43 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class for a primary generator 
//

#ifndef   RADMONGENERATORFIXEDDIRECTION_HH
 #define  RADMONGENERATORFIXEDDIRECTION_HH

 // Include files
 #include "RadmonVGeneratorWithLabel.hh"
 
 // Forward declaration
 class G4String;
 class G4ParticleGun;
 
 class RadmonGeneratorFixedDirection : public RadmonVGeneratorWithLabel
 {
  public:
   inline                                       RadmonGeneratorFixedDirection();
   inline virtual                              ~RadmonGeneratorFixedDirection();

   virtual void                                 ConvolveParticleGun(G4ParticleGun & gun);
   virtual RadmonVGeneratorWithLabel *          New(void) const;

  // Hidden constructors and operators
  private:
                                                RadmonGeneratorFixedDirection(const RadmonGeneratorFixedDirection & copy);
   RadmonGeneratorFixedDirection &              operator=(const RadmonGeneratorFixedDirection & copy);
 };
 
 // Inline implementations
 #include "RadmonGeneratorFixedDirection.icc"
#endif /* RADMONGENERATORFIXEDDIRECTION_HH */
