//
// File name:     RadmonGeneratorFixedPosition.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorFixedPosition.hh,v 1.1 2005-10-25 16:36:41 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class for a primary generator 
//

#ifndef   RADMONGENERATORFIXEDPOSITION_HH
 #define  RADMONGENERATORFIXEDPOSITION_HH

 // Include files
 #include "RadmonVGeneratorWithLabel.hh"
 
 // Forward declaration
 class G4String;
 class G4ParticleGun;
 
 class RadmonGeneratorFixedPosition : public RadmonVGeneratorWithLabel
 {
  public:
   inline                                       RadmonGeneratorFixedPosition();
   inline virtual                              ~RadmonGeneratorFixedPosition();

   virtual void                                 ConvolveParticleGun(G4ParticleGun & gun);
   virtual RadmonVGeneratorWithLabel *          New(void) const;

  // Hidden constructors and operators
  private:
                                                RadmonGeneratorFixedPosition(const RadmonGeneratorFixedPosition & copy);
   RadmonGeneratorFixedPosition &               operator=(const RadmonGeneratorFixedPosition & copy);
 };
 
 // Inline implementations
 #include "RadmonGeneratorFixedPosition.icc"
#endif /* RADMONGENERATORFIXEDPOSITION_HH */
