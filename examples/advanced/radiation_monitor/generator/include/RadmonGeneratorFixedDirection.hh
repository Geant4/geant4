//
// File name:     RadmonGeneratorFixedDirection.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorFixedDirection.hh,v 1.1 2005-10-25 16:36:41 capra Exp $
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
