//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include <iomanip>

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "NTSTDetectorConstruction.hh"
#include "NTSTDetectorMessenger.hh"
#include "NTSTRotationMatrix.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4Transform3D.hh"
#include "G4Point3D.hh"
#include "NTSTFileRead.hh"

#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4SimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4HelixMixedStepper.hh"
#include "G4NystromRK4.hh"

#include "G4DELPHIMagField.hh"
#include "G4PropagatorInField.hh"

NTSTDetectorConstruction::NTSTDetectorConstruction() 
  : _FileRead(0), debug(false), radius(19*cm), NSubLayer(0),
    disableSVT(false), disableDCH(false),
    field( 1.5*tesla, 0, 0 ),
    fpChordFinder( 0 ), 
    fMinChordStep( 0.1 )  // was 0.001 *mm )
{
  _FileRead = new NTSTFileRead("SVT.dat");

  // create commands necessary for the definition of the SVT
    
  DetectorMessenger = new NTSTDetectorMessenger(this);
}

NTSTDetectorConstruction::~NTSTDetectorConstruction()
{
  delete _FileRead;
  delete fpChordFinder;
  delete DetectorMessenger;
}

void NTSTDetectorConstruction::SetInputFileName(G4String FileName)
{
  delete _FileRead;
  _FileRead = new NTSTFileRead(FileName);
}


void NTSTDetectorConstruction::SetDebugCmd(G4int NewDebug)
{
  debug = NewDebug;
  if (debug) {
    G4cout << "Reset debug flag to true" << G4endl;
  }
  else {
    G4cout << "Reset debug flag to false" << G4endl;
  }
}


void
NTSTDetectorConstruction::SetNSubLayer(G4int NewNSubLayer)
{
  NSubLayer = NewNSubLayer;
  G4cout << "Reset number of sublayers to " << NSubLayer << G4endl;
}

void
NTSTDetectorConstruction::SetOuterRadius(G4double NewRadius)
{
  radius = NewRadius;
  G4cout << "Reset SVT mother volume outer radius parameter to " << radius
	 << G4endl;
}

void
NTSTDetectorConstruction::DisableDetector(G4String theDetector)
{
  if (theDetector == "SVT") {
    G4cout << "Disable " << theDetector << " detector" << G4endl;
    disableSVT=true;
  }
  else if (theDetector == "DCH") {
    G4cout << "Disable " << theDetector << " detector" << G4endl;
    disableDCH=true;
  }
  else if (theDetector == "all") {
    G4cout << "Disable SVT and DCH" << G4endl;
    disableSVT=true;
    disableDCH=true;
  }
  else if (theDetector == "none") {
    G4cout << "Enable SVT and DCH" << G4endl;
    disableSVT=false;
    disableDCH=false;
  }
}


void
NTSTDetectorConstruction::PrintCorners(const G4Transform3D& theT,
                                       G4LogicalVolume* theLV)
{
  G4VSolid* theSolid=theLV->GetSolid();
  G4Trd* theTRD=(G4Trd*)theSolid;
        
  G4double x1=theTRD->GetXHalfLength1();
  G4double x2=theTRD->GetXHalfLength2();
  G4double y1=theTRD->GetYHalfLength1();
  G4double y2=theTRD->GetYHalfLength2();
  G4double z =theTRD->GetZHalfLength();

  G4Point3D t1(-x1, -y1, -z);
  G4Point3D t2(+x1, -y1, -z);
  G4Point3D t3(-x1, +y1, -z);
  G4Point3D t4(+x1, +y1, -z);
  G4Point3D t5(-x2, -y2,  z);
  G4Point3D t6(+x2, -y2,  z);
  G4Point3D t7(-x2, +y2,  z);
  G4Point3D t8(+x2, +y2,  z);

  G4Point3D u1 = theT*t1;
  G4Point3D u2 = theT*t2;
  G4Point3D u3 = theT*t3;
  G4Point3D u4 = theT*t4;
  G4Point3D u5 = theT*t5;
  G4Point3D u6 = theT*t6;
  G4Point3D u7 = theT*t7;
  G4Point3D u8 = theT*t8;

  G4cout << std::setw(9) << u1.z() << std::setw(9) << u2.z() << std::setw(9) << u3.z()
	 << std::setw(9) << u4.z() << std::setw(9) << u5.z() << std::setw(9) << u6.z()
	 << std::setw(9) << u7.z() << std::setw(9) << u8.z() << G4endl;
}


G4VPhysicalVolume*
NTSTDetectorConstruction::Construct()
{
  //------------------------------------------------------ field
  G4Mag_UsualEqRhs *pEquation;
  G4MagIntegratorStepper *pStepper;
  G4FieldManager *globalFieldManager; 

  globalFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  G4PropagatorInField *
  globalPropagatorInField= G4TransportationManager::GetTransportationManager()->GetPropagatorInField();

  globalPropagatorInField->SetMaxLoopCount( 10000 ); 
  G4cout 
    << "PropagatorInField parameter(s) are: " << G4endl
    << " SetMaxLoopCount=" << globalPropagatorInField->GetMaxLoopCount()
    << " minEpsilonStep= " << globalPropagatorInField->GetMinimumEpsilonStep() << " "
    << " maxEpsilonStep= " << globalPropagatorInField->GetMaximumEpsilonStep() << " " 
    << G4endl;

  globalFieldManager->SetDetectorField( (G4MagneticField *)&field );

  // globalFieldManager->SetMinimumEpsilonStep( 5.0e-7 );    // Old value
  // globalFieldManager->SetMaximumEpsilonStep( 0.05 );      // FIX - old value
  // globalFieldManager->SetDeltaOneStep( 0.25 * mm );       // original value
  // globalFieldManager->SetDeltaIntersection( 0.10 * mm );  // original value

  G4cout << "Field Manager's parameters are " 
	 << " minEpsilonStep= " << globalFieldManager->GetMinimumEpsilonStep() << " "
	 << " maxEpsilonStep= " << globalFieldManager->GetMaximumEpsilonStep() << " " 
	 << " deltaOneStep=   " << globalFieldManager->GetDeltaOneStep() << " "
	 << " deltaIntersection= " << globalFieldManager->GetDeltaIntersection() 
	 << G4endl;

  pEquation = new G4Mag_UsualEqRhs( &field); 
 
  // pStepper =   
  //           new G4ClassicalRK4( pEquation ); G4cout << "Stepper is " << "ClassicalRK4" << G4endl;
  //           new G4RKG3_Stepper( pEquation );  // Nystrom, like Geant3 
  // pStepper= new G4SimpleRunge( pEquation ); G4cout << "Stepper is " << "CashKarpRKF45" << G4endl;
  // pStepper= new G4CashKarpRKF45( pEquation ); G4cout << "Stepper is " << "CashKarpRKF45" << G4endl;
  // pStepper= new G4HelixMixedStepper( pEquation ); G4cout << "Stepper is " << "HelixMixed" << G4endl;
  // pStepper=  StepperFactory::CreateStepper( order ); 

  pStepper= new G4NystromRK4( pEquation ); G4cout << "Stepper is " << "NystromRK4" << G4endl;

    // G4cout << "Stepper is " << "CashKarpRKF45" << G4endl;
    //	 << "ClassicalRK4" << G4endl;
    //   << " G4HelixMixedStepper " << G4endl;

  // globalFieldManager->CreateChordFinder( (G4MagneticField *)&field );
  fpChordFinder= new G4ChordFinder( (G4MagneticField *)&field, 
				    fMinChordStep,
                                    pStepper );
  fpChordFinder->SetVerbose(1); 
  globalFieldManager->SetChordFinder( fpChordFinder );

  //------------------------------------------------------ materials

  G4double a;  // atomic mass
  G4double zn; // atomic number
  G4double density;
  G4String name;

  a = 39.95*g/mole;
  density = 1.782e-03*g/cm3;
  G4Material* Ar = new G4Material(name="ArgonGas", zn=18., a, density);

  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  // G4Material* Al = 
  new G4Material(name="Aluminum", zn=13., a, density);

  a = 28.0855*g/mole;
  density = 2.33*g/cm3;
  G4Material* Si = new G4Material(name="Silicon", zn=14., a, density);

  //------------------------------------------------------ volumes

  //------------------------------ experimental hall (world volume)

  G4double expHall_x = 1000*cm;
  G4double expHall_y = 1000*cm;
  G4double expHall_z = 2000*cm;
  G4Box* experimentalHall_box
    = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  G4LogicalVolume* experimentalHall_log
    = new G4LogicalVolume(experimentalHall_box,Ar,"expHall_log",0,0,0);

  experimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);
  
  G4VPhysicalVolume* experimentalHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),"expHall",
                        experimentalHall_log,0,false,0);

  G4double innerRadiusOfTheSvt = 2.9*cm;
  G4double outerRadiusOfTheSvt = radius;
  G4double lengthOfTheSvt = 40.*cm;
  G4double startAngleOfTheSvt = 0*deg;
  G4double spanningAngleOfTheSvt = 360.*deg;

  G4double SvtPos_x = 0.*m;
  G4double SvtPos_y = 0.*m;
  G4double SvtPos_z = 0.*m;

  G4double innerRadiusOfTheDch = 24*cm;
  G4double outerRadiusOfTheDch = 81*cm;
  G4double lengthOfTheDch = 250*cm;
  G4double startAngleOfTheDch = 0*deg;
  G4double spanningAngleOfTheDch = 360.*deg;
  
  G4double DchPos_x = 0.*m;
  G4double DchPos_y = 0.*m;
  G4double DchPos_z = 0.*m;

  disableSVT=false;
  
  if (disableSVT == false){
  
    //------------------------------ SVT tracker volume
      
    G4Tubs* Svt_tube
      = new G4Tubs("Svt_tube",innerRadiusOfTheSvt,
		   outerRadiusOfTheSvt,lengthOfTheSvt,
		   startAngleOfTheSvt,spanningAngleOfTheSvt);
    G4LogicalVolume* Svt_log
      = new G4LogicalVolume(Svt_tube,Ar,"Svt_log",0,0,0);

    Svt_log -> SetVisAttributes(G4VisAttributes::Invisible);
      
    // G4VPhysicalVolume* Svt_phys =
      new G4PVPlacement(0,
			G4ThreeVector(SvtPos_x,
				      SvtPos_y,
				      SvtPos_z),
			Svt_log,"Svt",experimentalHall_log,false,0);
      
    if (debug)
      G4cout << "Placed SVT mother of length: "
             << std::setw(7) << lengthOfTheSvt/cm
	     << " and radii (cm): "
	     << std::setw(7) << innerRadiusOfTheSvt/cm
	     << std::setw(7) << outerRadiusOfTheSvt/cm << G4endl;
      
    //------------------------------ SVT guts
      
    // read in parameters of the wafers
      
    int NwafType=0;
    _FileRead->StreamLine() >> NwafType;
    if (debug) G4cout << "Number of wafer types: " << NwafType << G4endl;
      
    G4LogicalVolume** theWafer_log = new G4LogicalVolume*[NwafType];
      
    // define wafer vis attributes
      
    G4Color red(1,0,0);
    // G4Color green(0,1,0);
    // G4Color blue(0,0,1);
    G4VisAttributes* vAttr = new G4VisAttributes(red);
    // make solid
    vAttr->SetForceSolid(true);
      
    // define wafer shapes and create logical volumes indexed by Wafer type
      
    for (int ind=0; ind<NwafType; ind++){
          
      G4double Hmin, Hmax, Hzlen, Hthick;
      G4int IwafType;
      _FileRead->StreamLine() >> IwafType >> Hmin >> Hmax >> Hzlen >> Hthick;
      if (debug) G4cout << "Wafer type " << std::setw(3) << IwafType
			<< " Hmin " << std::setw(10) << Hmin/cm
			<< " Hmax " << std::setw(10) << Hmax/cm
			<< " Hzlen " << std::setw(10) << Hzlen/cm
			<< " Hthick " << std::setw(6) << Hthick/cm << G4endl;
          
      G4Trd* aWafer = new G4Trd("aWafer", Hthick*mm, Hthick*mm,
				Hmin*mm, Hmax*mm, Hzlen*mm);
      theWafer_log[IwafType-1]
	= new G4LogicalVolume(aWafer, Si, "aWafer_log", 0,0,0);
          
      theWafer_log[IwafType-1] -> SetVisAttributes(vAttr);
    }
      
    // get number of layers
      
    G4int Nsublayer=NSubLayer;
      
    _FileRead->StreamLine() >> Nsublayer;
    if (debug) G4cout << "Number of layers " << Nsublayer << G4endl;
      
    if (NSubLayer>0 && NSubLayer<=7){
      Nsublayer = NSubLayer;
    }
      
    // loop over the number of layers
      
    for (G4int Isublay=0; Isublay<Nsublayer;Isublay++){
      G4int Ilayer, Isublayer, Nmodule;
      _FileRead->StreamLine() >> Ilayer >> Isublayer >> Nmodule;
      if (debug) G4cout << "Number of modules for layer "
                        << std::setw(3) << Ilayer
			<< " sublayer " << std::setw(3) << Isublayer << " = "
			<< std::setw(3) << Nmodule << G4endl;
          
      // loop over the number of modules
          
      for (G4int Imod=0; Imod<Nmodule; Imod++){
	G4int Imodule, Nwafer;
	_FileRead->StreamLine() >> Imodule >> Nwafer;
	if (debug) G4cout << "Number of wafers in module "
                          << std::setw(3) << Imodule
			  << " = " << std::setw(3) << Nwafer << G4endl;
              
	// loop over the number of wafers in a module
              
	for (G4int Iwaf=0; Iwaf < Nwafer; Iwaf++){
	  G4int Iwafer, IwaferType;
	  _FileRead->StreamLine() >> Iwafer >> IwaferType;
	  if (debug) G4cout << "Wafer " << std::setw(3) << Iwafer
                            << " type " << std::setw(3)
			    << IwaferType << G4endl;
	  G4double x,y,z;
	  _FileRead->StreamLine() >> x >> y >> z;
	  G4ThreeVector WafPos(x*mm, y*mm, z*mm);
	  if (debug) G4cout << " position " << std::setw(9) << x << " "
                            << std::setw(9) << y
			    << " " << std::setw(9) << z << G4endl;
                  
	  _FileRead->StreamLine() >> x >> y >> z;
	  if (debug) G4cout << "Rotation Matrix:" << G4endl;
                  
	  G4ThreeVector row1(x,y,z);
	  if (debug) G4cout << row1 << G4endl;
                  
	  _FileRead->StreamLine() >> x >> y >> z;
                  
	  G4ThreeVector row2(x,y,z);
	  if (debug) G4cout << row2 << G4endl;
                  
	  _FileRead->StreamLine() >> x >> y >> z;
                  
	  G4ThreeVector row3(x,y,z);
	  if (debug) G4cout << row3 << G4endl;
                  
	  NTSTRotationMatrix WafMat;
                  
	  WafMat.SetRotationMatrixByRow(row1,row2,row3);
                  
	  G4Transform3D theTransform(WafMat, WafPos);              
                  
	  // G4VPhysicalVolume* wafer_phys =
	    new G4PVPlacement(theTransform, theWafer_log[IwaferType-1],
			      "WaferPos",Svt_log,false,0);
	  if (Imod==0 && debug) {
	    G4cout << "lay " << std::setw(3) << Ilayer << " Waf "
		   << Iwafer;
	    PrintCorners(theTransform, theWafer_log[IwaferType-1]);
	  }
	}
      }
    }
  } // end SVT block
  if (disableDCH == false) {
    G4Tubs* Dch_tube
      = new G4Tubs("Dch_tube",innerRadiusOfTheDch,
		   outerRadiusOfTheDch,lengthOfTheDch,
		   startAngleOfTheDch,spanningAngleOfTheDch);
    G4LogicalVolume* Dch_log
      = new G4LogicalVolume(Dch_tube,Ar,"Dch_log",0,0,0);
      
    Dch_log -> SetVisAttributes(G4VisAttributes::Invisible);
      
    // G4VPhysicalVolume* Dch_phys =
      new G4PVPlacement(0,
			G4ThreeVector(DchPos_x,
				      DchPos_y,
				      DchPos_z),
			Dch_log,"Dch",experimentalHall_log,false,0);
      
    if (debug)
      G4cout << "Placed DCH mother of length: "
             << std::setw(7) << lengthOfTheDch/cm
	     << " and radii (cm): "
	     << std::setw(7) << innerRadiusOfTheDch/cm
	     << std::setw(7) << outerRadiusOfTheDch/cm << G4endl;

    G4double r[41] = {25, 26, 27, 28, 30, 32, 33, 34, 35, 37, 38, 39, 41, 42,
		      43, 45, 46, 48, 49, 50, 52, 53, 54, 56, 57, 59, 60, 61,
		      62, 64, 66, 67, 68, 70, 71, 72, 73, 75, 76, 77, 78
    };

    for (int lay=0; lay < 40; lay++){
      G4double innerRadiusOfTheLayer=r[lay]*cm;
      G4double outerRadiusOfTheLayer=r[lay+1]*cm;
      G4double lengthOfTheLayer=lengthOfTheDch;
      G4double startAngleOfTheLayer=0*deg;
      G4double spanningAngleOfTheLayer=360*deg;

      G4Tubs* LayTub
	= new G4Tubs("Lay_tube",innerRadiusOfTheLayer,
		     outerRadiusOfTheLayer,lengthOfTheLayer,
		     startAngleOfTheLayer,spanningAngleOfTheLayer);
      G4LogicalVolume* Layer_log
	= new G4LogicalVolume(LayTub,Ar,"Layer_log",0,0,0);
          
      // G4VPhysicalVolume* Layer_phys =
	new G4PVPlacement(0,
			  G4ThreeVector(0),
			  Layer_log,"Layer", Dch_log,false,0);
      if (debug)
	G4cout << "Placed LAYER mother of length: "
               << std::setw(7) << lengthOfTheLayer/cm
	       << " and radii (cm): "
	       << std::setw(7) << innerRadiusOfTheLayer/cm
	       << std::setw(7) << outerRadiusOfTheLayer/cm << G4endl;
    }

  } // end DCH block
  

  //------------------------------------------------------------------

  return experimentalHall_phys;
}


