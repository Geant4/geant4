///////////////////////////////////////////////////////////////////////////////
// File: G4HcalTB96.cc
// Date: 08/00 S.Banerjee
// Modifications: 
///////////////////////////////////////////////////////////////////////////////
#include "G4HcalTB96.hh"

#include "CCalMaterialFactory.hh"

#include "G4HcalTB96HCal.hh"
#include "G4CrystalMatrix.hh"

#include "utils.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"

//#define debug
//#define pdebug

////////////////////////////////////////////////////////////////////
// G4HcalTB96 constructor & destructor...
////////////////////////////////////////////////////////////////////

G4HcalTB96::G4HcalTB96(const G4String &name):
  HcalTB96(name),G4Able(name) {}

G4HcalTB96::~G4HcalTB96() {}

////////////////////////////////////////////////////////////////////
// G4HcalTB96 methods...
////////////////////////////////////////////////////////////////////

G4VPhysicalVolume* G4HcalTB96::constructIn(G4VPhysicalVolume* mother) {

  cout << "==>> Constructing G4HcalTB96..." << endl;

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
  setVisType(Visualisable::PseudoVolumes,glog);

  G4VPhysicalVolume* volume = new G4PVPlacement(0,G4ThreeVector(),Name(),
                                                glog,mother,false,0);
#ifdef pdebug
  G4String name("Null");
  if (mother != 0) name = mother->GetName();
  cout << glog->GetName() << " number 1 positioned in " << name
       << " at (0,0,0) with no rotation" << endl;
#endif  

  cout << "<<== End of G4HcalTB96 construction ..." << endl;

  return volume;
}

void G4HcalTB96::constructDaughters(){
  //Hadron Calorimeter
  G4HcalTB96HCal* hcal = new G4HcalTB96HCal("HadronCalorimeter");
  addDetector(hcal);
  AddG4Able(hcal);

  //Crystal matrix
  G4CrystalMatrix* xtalmod = new G4CrystalMatrix("CrystalMatrixModule");
  xtalmod->setType(G4CrystalMatrix::module1);
  addDetector(xtalmod);
  AddG4Able(xtalmod);
}
