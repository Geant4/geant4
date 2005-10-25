//
// File name:     RadmonGeneratorFixedEnergy.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorFixedEnergy.hh,v 1.1 2005-10-25 16:36:41 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Abstract class for a primary generator 
//

#ifndef   RADMONGENERATORFIXEDENERGY_HH
 #define  RADMONGENERATORFIXEDENERGY_HH

 // Include files
 #include "RadmonVGeneratorWithLabel.hh"
 
 // Forward declaration
 class G4String;
 class G4ParticleGun;
 
 class RadmonGeneratorFixedEnergy : public RadmonVGeneratorWithLabel
 {
  public:
   inline                                       RadmonGeneratorFixedEnergy();
   inline virtual                              ~RadmonGeneratorFixedEnergy();

   virtual void                                 ConvolveParticleGun(G4ParticleGun & gun);
   virtual RadmonVGeneratorWithLabel *          New(void) const;

  // Hidden constructors and operators
  private:
                                                RadmonGeneratorFixedEnergy(const RadmonGeneratorFixedEnergy & copy);
   RadmonGeneratorFixedEnergy &                 operator=(const RadmonGeneratorFixedEnergy & copy);
 };
 
 // Inline implementations
 #include "RadmonGeneratorFixedEnergy.icc"
#endif /* RADMONGENERATORFIXEDENERGY_HH */
