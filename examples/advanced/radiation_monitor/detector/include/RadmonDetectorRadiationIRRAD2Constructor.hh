//
// File name:     RadmonDetectorRadiationIRRAD2Constructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorRadiationIRRAD2Constructor.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds the IRRAD2 environment
//

#ifndef   RADMONDETECTORRADIATIONIRRAD2CONSTRUCTOR_HH
 #define  RADMONDETECTORRADIATIONIRRAD2CONSTRUCTOR_HH

 // Include files
 #include "RadmonVDetectorLabelledEntityConstructor.hh"

 class RadmonDetectorRadiationIRRAD2Constructor : public RadmonVDetectorLabelledEntityConstructor
 {
  public:
                                                RadmonDetectorRadiationIRRAD2Constructor();
   virtual                                     ~RadmonDetectorRadiationIRRAD2Constructor();
   
   virtual G4LogicalVolume *                    ConstructLogicalVolume(void);
   virtual RadmonVDetectorLabelledEntityConstructor * New(void);

  private:
  // Hidden constructors and operators
                                                RadmonDetectorRadiationIRRAD2Constructor(const RadmonDetectorRadiationIRRAD2Constructor & copy);
   RadmonDetectorRadiationIRRAD2Constructor &   operator=(const RadmonDetectorRadiationIRRAD2Constructor & copy);
 };
#endif /* RADMONDETECTORRADIATIONIRRAD2CONSTRUCTOR_HH */
