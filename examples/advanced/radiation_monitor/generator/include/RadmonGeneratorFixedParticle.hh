//
// File name:     RadmonGeneratorFixedParticle.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorFixedParticle.hh,v 1.1 2005-10-25 16:36:41 capra Exp $
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
