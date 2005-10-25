//
// File name:     RadmonGeneratorUniformSphere.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorUniformSphere.hh,v 1.1 2005-10-25 16:36:41 capra Exp $
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
