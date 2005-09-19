//
// File name:     RadmonDetectorSimpleBoxConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorSimpleBoxConstructor.hh,v 1.1 2005-09-19 19:39:52 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds the IRRAD2 environment
//

#ifndef   RADMONDETECTORSIMPLEBOXCONSTRUCTOR_HH
 #define  RADMONDETECTORSIMPLEBOXCONSTRUCTOR_HH

 // Include files
 #include "RadmonVDetectorLabelledEntityConstructor.hh"
 
 // Forward declarations
 class G4Box;
 class G4LogicalVolume;

 class RadmonDetectorSimpleBoxConstructor : public RadmonVDetectorLabelledEntityConstructor
 {
  public:
   inline                                       RadmonDetectorSimpleBoxConstructor();
   virtual                                     ~RadmonDetectorSimpleBoxConstructor();
   
   virtual G4LogicalVolume *                    ConstructLogicalVolume(void);
   virtual RadmonVDetectorLabelledEntityConstructor * New(void);

  private:
  // Hidden constructors and operators
                                                RadmonDetectorSimpleBoxConstructor(const RadmonDetectorSimpleBoxConstructor & copy);
   RadmonDetectorSimpleBoxConstructor &         operator=(const RadmonDetectorSimpleBoxConstructor & copy);
   
  // Private attributes
   G4Box *                                      box;
   G4VisAttributes *                            visAttributes;
   G4LogicalVolume *                            logicalVolume;
 };
 
 // Inline implementations
 #include "RadmonDetectorSimpleBoxConstructor.icc"
#endif /* RADMONDETECTORSIMPLEBOXCONSTRUCTOR_HH */
