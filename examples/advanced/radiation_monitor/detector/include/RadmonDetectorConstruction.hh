//
// File name:     RadmonDetectorConstruction.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorConstruction.hh,v 1.2 2005-09-19 19:42:13 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Implementation of the G4VUserDetectorConstruction
//

#ifndef   RADMONDETECTORCONSTRUCTION_HH
 #define  RADMONDETECTORCONSTRUCTION_HH

 // Include files
 #include "globals.hh"
 #include "G4VUserDetectorConstruction.hh"
 #include "RadmonVDetectorLayoutObserver.hh"
 #include <stack>
 #include <utility>
 
 // Forward declaration
 class G4VPhysicalVolume;
 class G4VSolid;
 class G4LogicalVolume;
 class RadmonVDetectorLayout;
 class RadmonVDetectorEntityConstructor;
 class RadmonVDetectorEntitiesConstructorsFactory;
 class G4String;

 class RadmonDetectorConstruction : public G4VUserDetectorConstruction, public RadmonVDetectorLayoutObserver
 {
  public:
                                                RadmonDetectorConstruction(RadmonVDetectorLayout * layout, RadmonVDetectorEntitiesConstructorsFactory * factory);
   virtual                                     ~RadmonDetectorConstruction();

   virtual G4VPhysicalVolume *                  Construct(void);

   virtual void                                 OnLayoutChange(void);

  private:
  // Private methods
   void                                         Destruct(void);
   
   void                                         BuildEnvironmentFromType(const G4String & type);
   void                                         BuildEnvironmentSphere(void);
   void                                         BuildMultilayer(G4int index);

  // Hidden constructors and operators
                                                RadmonDetectorConstruction();
                                                RadmonDetectorConstruction(const RadmonDetectorConstruction & copy);
   RadmonDetectorConstruction &                 operator=(const RadmonDetectorConstruction & copy);

  // Private data types
   typedef std::pair<RadmonVDetectorEntityConstructor *, G4VPhysicalVolume *> LayerItem;
   typedef std::stack<LayerItem>                LayersStack;

  // Private attributes
   RadmonVDetectorLayout *                      detectorLayout;
   RadmonVDetectorEntitiesConstructorsFactory * constructorsFactory;
    
   LayersStack                                  layersStack;

   RadmonVDetectorEntityConstructor *           environmentConstructor;
   G4VPhysicalVolume *                          environmentPhysicalVolume;
   G4LogicalVolume *                            environmentLogicalVolume;
   G4VSolid *                                   environmentSolid;
 };
#endif /* RADMONDETECTORCONSTRUCTION_HH */
