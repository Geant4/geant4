//
// File name:     RadmonVGenerator.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVGenerator.hh,v 1.1 2005-10-25 16:36:41 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
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
