#include "Tst05DetectorConstruction.hh"

#include "Tst05DetectorMessenger.hh"

Tst05DetectorConstruction::Tst05DetectorConstruction():detectorChoice(0)
{
  expHall_x = 1111600.*cm;
  expHall_y = 1111600.*cm;
  expHall_z = 1111600.*cm;

  calBox_x = 50.*cm;
  calBox_y = 50.*cm;
  calBox_z = 50.*cm;
  rotAngle = 30.*deg;
  calPos = 200.*cm;

  trackerRadius = 50.*cm;
  trackerHight = 100.*cm;
  trackerPos = -200.*cm;

  detectorMessenger = new Tst05DetectorMessenger(this);
}

Tst05DetectorConstruction::~Tst05DetectorConstruction()
{;}

G4VPhysicalVolume* Tst05DetectorConstruction::Construct()
{

  //------------------------------------------------------ materials

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", iz=8., a);

  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  a = 207.19*g/mole;
  density = 11.35*g/cm3;
  G4Material* Lead = new G4Material(name="Lead", z=82., a, density);

  a = 39.95*g/mole;
  density = 1.782e-03*g/cm3;
  G4Material* Ar = new G4Material(name="ArgonGas", z=18., a, density);

  //------------------------------------------------------ volumes

  //------------------------------ experimental hall

  G4Box * experimentalHall_box
    = new G4Box("expHall_b",expHall_x,expHall_y,expHall_z);
  G4LogicalVolume * experimentalHall_log
    = new G4LogicalVolume(experimentalHall_box,Air,"expHall_L",0,0,0);
  G4VPhysicalVolume * experimentalHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),"expHall_P",
                        experimentalHall_log,0,false,0);


  // STEP assembly creation
  const char * stepfile;

  switch(detectorChoice)
  { 
    
    default:
      stepfile = "stepfiles/G4rod_place_asm.stp";
      break;

    case 1:
      stepfile = "stepfiles/CMSTracker.step";
      break;  
  } 

  G4AssemblyCreator MyAC(stepfile);
  MyAC.ReadStepFile();
  STEPentity* ent=0;
  MyAC.CreateG4Geometry(*ent);

  void *tmp =  MyAC.GetCreatedObject();
  G4Assembly* assembly = new G4Assembly();
  assembly->SetPlacedVector(*(G4PlacedVector*)tmp);

  G4PlacedSolid* ps=0;
  G4int solids = assembly->GetNumberOfSolids();
  for(G4int c=0;c<solids;c++)
    {
      ps = assembly->GetPlacedSolid(c);
      G4LogicalVolume* lv = new G4LogicalVolume(ps->GetSolid(),
						Lead,
						"STEPlog");
      HepRotation* hr = ps->GetRotation();
      G4ThreeVector* tr = ps->GetTranslation();
      new G4PVPlacement(hr,
			*tr,
		        ps->GetSolid()->GetName(),
			lv,
			experimentalHall_phys,
			false,
			c);
    }

  return experimentalHall_phys;
}


void Tst05DetectorConstruction::SelectDetector(G4String val)
{
  if(val=="CMSTracker") 
  { detectorChoice = 1; }
  else
  { detectorChoice = 0; }
  G4cout << "Now Detector is " << val << G4endl;
}














