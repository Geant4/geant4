//
// File name:     RadmonGeneratorUniformPlane.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorUniformPlane.hh,v 1.1 2005-10-25 16:36:41 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class for a primary generator 
//

#ifndef   RADMONGENERATORUNIFORMPLANE_HH
 #define  RADMONGENERATORUNIFORMPLANE_HH

 // Include files
 #include "RadmonVGeneratorWithLabel.hh"
 
 // Forward declaration
 class G4String;
 class G4ParticleGun;
 
 class RadmonGeneratorUniformPlane : public RadmonVGeneratorWithLabel
 {
  public:
   inline                                       RadmonGeneratorUniformPlane();
   inline virtual                              ~RadmonGeneratorUniformPlane();

   virtual void                                 ConvolveParticleGun(G4ParticleGun & gun);
   virtual RadmonVGeneratorWithLabel *          New(void) const;

  // Hidden constructors and operators
  private:
                                                RadmonGeneratorUniformPlane(const RadmonGeneratorUniformPlane & copy);
   RadmonGeneratorUniformPlane &                operator=(const RadmonGeneratorUniformPlane & copy);
 };
 
 // Inline implementations
 #include "RadmonGeneratorUniformPlane.icc"
#endif /* RADMONGENERATORUNIFORMPLANE_HH */
