///////////////////////////////////////////////////////////////////////////////
// File: CCalG4Hall.cc
// Description: CCalG4Hall Geometry factory class to construct G4 geometry 
//              of the experimental hall
///////////////////////////////////////////////////////////////////////////////
#include "CCalG4Hall.hh"

#include "CCalMaterialFactory.hh"

#include "CCalG4Hcal.hh"
#include "CCalG4Ecal.hh"

#include "CCalutils.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"

//#define debug
//#define pdebug

////////////////////////////////////////////////////////////////////
// CCalG4Hall constructor & destructor...
////////////////////////////////////////////////////////////////////

CCalG4Hall::CCalG4Hall(const G4String &name):
  CCalHall(name),CCalG4Able(name) {}

CCalG4Hall::~CCalG4Hall() {}

////////////////////////////////////////////////////////////////////
// CCalG4Hall methods...
////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* CCalG4Hall::constructIn(G4VPhysicalVolume* mother) {

  cout << "==>> Constructing CCalG4Hall..." << endl;

  ///////////////////////////////////////////////////////////////
  //Pointers to the Rotation Matrices and to the Materials
  CCalMaterialFactory* matfact       = CCalMaterialFactory::getInstance();

  //Experimental Hall. The mother volume....
#ifdef debug
  cout << tab << "Building experimental Hall geometry...." << endl;
#endif

  G4Box* solid = new G4Box(Name(), getDx_2Hall()*mm, getDy_2Hall()*mm,
			   getDy_2Hall()*mm);
#ifdef debug
  cout << solid->GetName() << " made of " << getMaterial() << " of dimension " 
       << getDx_2Hall()*mm << " " << getDy_2Hall()*mm << " " 
       << getDy_2Hall()*mm << " (all in mm)" << endl;
#endif

  G4Material* matter = matfact->findMaterial(getMaterial());
  G4LogicalVolume* glog = new G4LogicalVolume (solid, matter, Name(), 0, 0, 0);
  setVisType(CCalVisualisable::PseudoVolumes,glog);

  G4VPhysicalVolume* volume = new G4PVPlacement(0,G4ThreeVector(),Name(),
                                                glog,mother,false,0);
#ifdef pdebug
  G4String name("Null");
  if (mother != 0) name = mother->GetName();
  cout << glog->GetName() << " number 1 positioned in " << name
       << " at (0,0,0) with no rotation" << endl;
#endif  

  cout << "<<== End of CCalG4Hall construction ..." << endl;

  return volume;
}

void CCalG4Hall::constructDaughters(){
  //Hadron Calorimeter
  CCalG4Hcal* hcal = new CCalG4Hcal("HadronCalorimeter");
  addDetector(hcal);
  AddCCalG4Able(hcal);

  //Crystal matrix
  CCalG4Ecal* ecal = new CCalG4Ecal("CrystalMatrixModule");
  ecal->setType(CCalG4Ecal::module1);
  addDetector(ecal);
  AddCCalG4Able(ecal);
}
