
#ifndef TstVADetectorConstruction_h
#define TstVADetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class TstVADetectorMessenger;
class G4AssemblyVolume;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <vector>


class TstVADetectorConstruction : public G4VUserDetectorConstruction
{
//  public:
    struct sClassic
    {
      G4LogicalVolume*                    caloLV;
      G4std::vector<G4VPhysicalVolume*>   PVs;
    };
    
  public:
    TstVADetectorConstruction();
    ~TstVADetectorConstruction();

  public:
    G4VPhysicalVolume*      Construct();
    void                    SwitchDetector();
    void                    SelectDetector(G4String val);
    void                    SelectMaterial(G4String val);

//  private:
    G4VPhysicalVolume*      SelectDetector();
    void                    SelectMaterialPointer();
    void                    ConstructClassic();
    void                    ConstructAssembly();
    void                    CleanClassic();
    void                    CleanAssembly();

  private:
    G4VPhysicalVolume*      worldVol;
    G4Material*             Air;
    G4Material*             Al;
    G4Material*             Pb;
    G4Material*             selectedMaterial;
    G4int                   detectorChoice;
    G4String                materialChoice;
    TstVADetectorMessenger* detectorMessenger;
     
  private:
    // Very private data
    G4LogicalVolume*        plateLV;
    sClassic                classicDetector;
    G4AssemblyVolume*       assemblyDetector;
};

#endif

