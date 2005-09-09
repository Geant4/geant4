//
// File name:     RadmonDetectorConstruction.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorConstruction.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Implementation of the G4VUserDetectorConstruction
//

#ifndef   RADMONDETECTORCONSTRUCTION_HH
 #define  RADMONDETECTORCONSTRUCTION_HH

 // Include files
 #include "G4VUserDetectorConstruction.hh"
 #include "RadmonVDetectorLayoutObserver.hh"
 #include <list>
 
 // Forward declaration
 class G4VPhysicalVolume;
 class G4LogicalVolume;
 class RadmonVDetectorLayout;
 class RadmonVDetectorEntityConstructor;
 class RadmonVDetectorEntitiesConstructorsFactory;

 class RadmonDetectorConstruction : public G4VUserDetectorConstruction, public RadmonVDetectorLayoutObserver
 {
  public:
                                                RadmonDetectorConstruction(RadmonVDetectorLayout * layout, RadmonVDetectorEntitiesConstructorsFactory * factory);
    virtual                                    ~RadmonDetectorConstruction();

    virtual G4VPhysicalVolume *                 Construct(void);

    virtual void                                OnLayoutChange(void);

  private:
                                                RadmonDetectorConstruction();
                                                RadmonDetectorConstruction(const RadmonDetectorConstruction & copy);
    RadmonDetectorConstruction &                operator=(const RadmonDetectorConstruction & copy);

    void                                        Destruct(void);

    RadmonVDetectorLayout *                     detectorLayout;
    RadmonVDetectorEntitiesConstructorsFactory * constructorsFactory;
    
    std::list<RadmonVDetectorEntityConstructor *> constructorsList;
    G4VPhysicalVolume *                         motherPhysicalVolume;
    G4LogicalVolume *                           motherLogicalVolume;
    std::list<G4VPhysicalVolume *>              daughtersPhysicalVolumesList;
 };
#endif /* RADMONDETECTORCONSTRUCTION_HH */
