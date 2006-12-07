#include "RandomCaloDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"

#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
      
#include "RandomCaloDetectorMessenger.hh"

#include "G4SDManager.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4NistManager.hh"
#include "G4StableIsotopes.hh"
           

RandomCaloDetectorConstruction::RandomCaloDetectorConstruction() :  
  experimentalHall_log( 0 ), experimentalHall_phys( 0 ),
  fieldMgr( 0 ), uniformMagField( 0 ), 
  detectorMessenger( 0 ), nistManager( 0 ),

  thicknessLayer( 10.0*cm ),                                //***LOOKHERE***
  radiusCalo( 50.0*cm )                                     //***LOOKHERE***
  // These two constant parameters specify, respectively,
  // the radius of the cylindrical layers, and thickness 
  // of each of them, which is assumed to be the same
  // for all layers, independently of the material.

{
  fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  DefineMaterials();
  detectorMessenger = new RandomCaloDetectorMessenger( this );
}


RandomCaloDetectorConstruction::~RandomCaloDetectorConstruction() {
  delete detectorMessenger;
  delete uniformMagField;
}


G4VPhysicalVolume* RandomCaloDetectorConstruction::Construct() {
  return ConstructCalorimeter();
}


void RandomCaloDetectorConstruction::DefineMaterials() { 
  // This method fills the (empty) vector  vecNistMaterials  with
  // the pointers to all NIST materials which satisfy two conditions
  // (necessary to be able to simulate Geant4 interactions a
  // material): first, the atomic number (Z) of each element
  // of the material must not exceed 92; second, for each 
  // element of the material, there must be at least one
  // stable isotope. The class G4StableIsotopes is used to
  // check the number of stable isotopes for a given element.

  G4cout << G4endl
         << " -------------------------------------------------------- " << G4endl
         << " BEGIN  RandomCaloDetectorConstruction::DefineMaterials()" << G4endl
         << G4endl;

  nistManager = G4NistManager::Instance();
  G4StableIsotopes stableIsotopesObject;
  if ( nistManager ) {
    std::vector< G4String > vecNistMaterialNames = nistManager->GetNistMaterialNames();
    G4cout << " Number of NIST materials = " << vecNistMaterialNames.size() << G4endl; 
    for ( std::vector<G4String>::const_iterator cit = vecNistMaterialNames.begin();
	  cit != vecNistMaterialNames.end(); ++cit ) {
      G4cout << " ---> materialString = " << *cit << G4endl;
      G4Material* nistMaterial = nistManager->FindOrBuildMaterial( *cit ); 
      if ( nistMaterial ) {
	G4cout << "\t materialName = " << nistMaterial->GetName() << G4endl
	       << "\t density = " << nistMaterial->GetDensity() / (gram/cm3) 
               << "  g/cm^3 " << G4endl
	       << "\t X0 = " << nistMaterial->GetRadlen() / cm << " cm " << G4endl
	       << "\t Lambda = " << nistMaterial->GetNuclearInterLength() / cm << " cm "
	       << G4endl
	       << "\t numberOfElements = " << nistMaterial->GetNumberOfElements() 
	       << G4endl;
	bool isMaterialOK = true;
	for ( int iElement = 0; 
	      iElement < static_cast< int >( nistMaterial->GetNumberOfElements() ); 
	      iElement++ ) {
	  const G4Element* pElement = nistMaterial->GetElement( iElement );
	  if ( pElement ) {
	    G4int numStableIsotopes =
	      stableIsotopesObject.GetNumberOfIsotopes( static_cast< int >( pElement->GetZ() ) ); 
	    G4cout << "\t \t element=" << pElement->GetName() 
	           << "  Z=" << pElement->GetZ()
                   << "  #Nucleons=" << pElement->GetN() 
	           << "  #StableIsotopes=" << numStableIsotopes
		   << G4endl;
	    if ( pElement->GetZ() > 92  ||  numStableIsotopes == 0 ) {
	      isMaterialOK = false;
            }
	  } else {
	    G4cout << "\t \t NO ELEMENT!" << G4endl;
	  }
	}
	if ( isMaterialOK ) {
	  vecNistMaterials.push_back( nistMaterial ); 
	} else {
	  G4cout << "\t NOT CONSIDERED! " << G4endl;
        }
      } else {
	G4cout << "\t Something wrong with the NIST Material!" << G4endl;
      }
    }
  } else {
    G4cout << "\t Cannot instantiate the NIST material manager!" << G4endl;
  }

  G4cout << G4endl
	 << " Number of considered NIST materials = " << vecNistMaterials.size() 
	 << G4endl << G4endl
         << " END  RandomCaloDetectorConstruction::DefineMaterials()" << G4endl
         << " ------------------------------------------------------" << G4endl
	 << G4endl;
}


G4VPhysicalVolume* RandomCaloDetectorConstruction::ConstructCalorimeter() {

  //G4cout << " BEGIN  RandomCaloDetectorConstruction::ConstructCalorimeter()" 
  //       << G4endl; //***DEBUG***

  // Clean old geometry, if any.
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  G4Material* material = 0;

  // --- Experimental hall (world volume)    ***LOOKHERE***
  // Beam line along the Z-axis.
  // The total length of the calorimeter is given by the product
  // of the number of layers, which is assumed to be equal to
  // the number of selected NIST materials, and the thickness
  // of each layer, assumed to be the same for all layers
  // independently of the material.
  // The experimental hall is a box slightly larger, by 1 cm,
  // with respect to the cylindrical calorimeter.

  G4double zMaxCalorimeter = ( thicknessLayer / 2.0 ) * vecNistMaterials.size();
  G4double expHall_x = radiusCalo + 1.0*cm;         // half dimension along x 
  G4double expHall_y = radiusCalo + 1.0*cm;         // half dimension along y
  G4double expHall_z = zMaxCalorimeter + 1.0*cm;    // half dimension along z

  G4cout << G4endl
         << " SIZE Experimental Hall: " << G4endl
         << "\t X-axis : " << -expHall_x/m << " : " << expHall_x/m << " m " << G4endl
         << "\t Y-axis : " << -expHall_y/m << " : " << expHall_y/m << " m " << G4endl
         << "\t X-axis : " << -expHall_z/m << " : " << expHall_z/m << " m " << G4endl
         << G4endl;

  G4Box* experimentalHall_box
    = new G4Box( "expHall_box", expHall_x, expHall_y, expHall_z );

  // Air is assigned as material of the experimental hall.
  if ( nistManager ) {
    material = nistManager->FindOrBuildMaterial( "G4_AIR" );
  }

  experimentalHall_log = new G4LogicalVolume( experimentalHall_box, // solid 
                                              material,             // material
                                              "expHall_log",        // name
                                              0,                    // field manager
                                              0,                    // sensitive detector
                                              0 );                  // user limits

  experimentalHall_phys = new G4PVPlacement( 0,                     // rotation
                                             G4ThreeVector(),       // translation
                                             "expHall",             // name
                                             experimentalHall_log,  // logical volume
                                             0,                     // mother phy. vol
                                             false,                 // boolean operation
                                             0 );                   // copy number

  // All layers have the same thickness, thicknessLayer, and
  // the same radius, radiusCalo.
  // Each layer have a different material. 
  // The order of the materials follows the order in which they
  // appear in the NIST table, from which we extract the subset
  // vecNistMaterials .

  G4Tubs* solidLayer = 
    new G4Tubs( "solidLayer", 
		0.0,                 // inner radius
		radiusCalo,          // outer radius
		thicknessLayer/2.0,  // half cylinder length in z
		0.0,                 // starting phi angle in rad
		2.0*pi );            // final phi angle in rad
  
  G4LogicalVolume* logicLayer = 0;
  G4VPhysicalVolume* physicLayer = 0;
  G4double zPosition = -zMaxCalorimeter + thicknessLayer/2.0;

  G4cout << " POSITIONS and MATERIALS of the Calorimeter Layers:" << G4endl;
  for ( std::vector< G4Material* >::const_iterator cit = vecNistMaterials.begin();
	cit != vecNistMaterials.end(); ++cit ) {
      material = *cit; 
      if ( material ) {
	G4cout << "\t Layer of: " << material->GetName() 
               << "  zMin=" << ( zPosition - thicknessLayer/2.0 ) / cm 
               << "  zMax=" << ( zPosition + thicknessLayer/2.0 ) / cm 
	       << " cm " << G4endl;
	logicLayer = new G4LogicalVolume( solidLayer,     // solid 
					  material     ,  // material
					  "logicLayer",   // name
					  0,              // field manager
					  0,              // sensitive det.
					  0 );            // user limits

	physicLayer = new G4PVPlacement( 0,                     // rotation
                                                                // translation
					 G4ThreeVector( 0.0, 0.0, zPosition ), 
					 "physiLayer",          // its name
					 logicLayer,            // logical volume
					 experimentalHall_phys, // mother physical volume
					 false,                 // boolean operation
					 1 );                   // copy number

	zPosition += thicknessLayer;
      }
  }
  G4cout << G4endl;

  // --- Visualization attributes

  experimentalHall_log->SetVisAttributes( G4VisAttributes::Invisible );
  // The World is not visualized.

  G4VisAttributes* theVisAttLayer = new G4VisAttributes( G4Colour( 1.0, 1.0, 1.0 ) );
  theVisAttLayer->SetVisibility( true );
  theVisAttLayer->SetForceWireframe( true );

  logicLayer->SetVisAttributes( theVisAttLayer );
  // The layer will appear in white colour.
  // (the order of colours is: (red, green, blue) )

  return experimentalHall_phys;
}


void RandomCaloDetectorConstruction::SetMagField( const G4double fieldValue ) {
  if ( uniformMagField ) {
    delete uniformMagField;
  }
  if ( std::abs( fieldValue ) > 0.0 ) {
    // Apply a global uniform magnetic field along the Y axis.
    // Notice that only if the magnetic field is not zero, the Geant4
    // transportion in field gets activated.

    uniformMagField = new G4UniformMagField( G4ThreeVector( 0.0, fieldValue, 0.0 ) );

    fieldMgr->SetDetectorField( uniformMagField );
    fieldMgr->CreateChordFinder( uniformMagField );

  } 
}


void RandomCaloDetectorConstruction::UpdateGeometry() {
  G4RunManager::GetRunManager()->DefineWorldVolume( ConstructCalorimeter() );
  PrintParameters();
}


void RandomCaloDetectorConstruction::PrintParameters() {
  G4cout << G4endl << G4endl
         << " ------  RandomCaloDetectorConstruction::PrintParameters() ------ " 
	 << G4endl;

  G4cout << G4endl << " Magnetic field [T]    = ";
  if ( uniformMagField ) {
    G4cout << uniformMagField->GetConstantFieldValue() / tesla;
  } else {
    G4cout << "(0,0,0)";
  }

  G4cout << G4endl << " -------------------------------------------------------- "
         << G4endl << G4endl;
}
