#include "G4tgbGeometryDumper.hh"

#include "G4tgrMessenger.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4BooleanSolid.hh"
#include "G4ReflectionFactory.hh"
#include "G4ReflectedSolid.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include <iomanip>

using namespace CLHEP;

//------------------------------------------------------------------------
G4tgbGeometryDumper* G4tgbGeometryDumper::theInstance = 0;

//------------------------------------------------------------------------
G4tgbGeometryDumper::G4tgbGeometryDumper()
{

  theRotationNumber = 0;
}

//------------------------------------------------------------------------
G4tgbGeometryDumper* G4tgbGeometryDumper::GetInstance()
{
  if( theInstance == 0 ){
    theInstance = new G4tgbGeometryDumper;
  }

  return theInstance;

}

//------------------------------------------------------------------------
void G4tgbGeometryDumper::DumpGeometry( const G4String& fname )
{
  
  theFile = new std::ofstream(fname);

  G4VPhysicalVolume* pv = GetTopPhysVol();
  DumpPhysVol( pv ); // dump volume and recursively it will dump all hierarchy
}

//---------------------------------------------------------------------
G4VPhysicalVolume* G4tgbGeometryDumper::GetTopPhysVol()
{
  G4PhysicalVolumeStore* pvstore = G4PhysicalVolumeStore::GetInstance();
  G4PhysicalVolumeStore::const_iterator ite;
  G4VPhysicalVolume* pv = *(pvstore->begin());
  for( ;; ){
    G4LogicalVolume* lv = pv->GetMotherLogical();
    if( lv == 0 ) break;
    //----- look for one PV of this LV
    for( ite = pvstore->begin(); ite != pvstore->end(); ite++ ){
      pv = (*ite);
      if( pv->GetLogicalVolume() == lv ) {
	break;
      }
    }
  }

  return pv;
}


//------------------------------------------------------------------------
G4tgbGeometryDumper::~G4tgbGeometryDumper()
{
}

//------------------------------------------------------------------------
void G4tgbGeometryDumper::DumpPhysVol( G4VPhysicalVolume* pv )
{

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " G4tgbGeometryDumper::DumpPhysVol " << pv->GetName() << G4endl;
#endif

  //---- Dump logical volume first
  G4LogicalVolume* lv = pv->GetLogicalVolume();

  G4ReflectionFactory* reffact = G4ReflectionFactory::Instance();

  //---- It is not needed to dump _refl volumes created when parent is reflected 
  // !!WARNING : it must be avoided to reflect a volume hierarchy if children has also been reflected, as both will have same name
  if( reffact->IsReflected( lv ) && reffact->IsReflected( pv->GetMotherLogical() ) ) return;

  G4bool bVolExists =  !CheckIfLogVolExists( lv->GetName(), lv );
  if( bVolExists ){
    DumpLogVol( lv );
  }

  //---- Construct this PV
  if( pv->GetMotherLogical() != 0 ) { // WORLD volume
  
    if( pv->IsReplicated() || pv->IsParameterised() ){
      G4Exception("G4tgbGeometryDumper::DumpPlacements   only G4PVPlacement is supported yet " + pv->GetName() );
    }
    
    G4String pvName = pv->GetName();

    /*
      if( pv->GetMotherLogical() != 0 ) {
      (*theFile) << " GPV: " << pv->GetName() << " # " << pv->GetCopyNo() << " in " << pv->GetMotherLogical()->GetName();
      if( pv->GetRotation() ){
      (*theFile) << " ROT: " <<  pv->GetRotation()->colX() << " " <<  pv->GetRotation()->colY() << " " << pv->GetRotation()->colZ();
      }
      (*theFile)  << G4endl;
      }
    */
    G4RotationMatrix* rotMat = pv->GetRotation();
    if( !rotMat ) rotMat = new G4RotationMatrix();
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
      G4cout << " PV RotationMatrix " << rotMat << G4endl;
#endif

    //  G4PVPlacement* pvpl = static_cast<G4PVPlacement*>(pv);
    //---- Check if it is reflected
    if( reffact->IsReflected( lv ) ){
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
	G4cout << " DumpPhysVol: reflected volume " << pv->GetName() << G4endl;
#endif
      size_t irefl = lv->GetName().find("_refl");
      Hep3Vector colx = rotMat->colX();
      Hep3Vector coly = rotMat->colY();
      Hep3Vector colz = rotMat->colZ();
      // apply a Z reflection (reflection matrix is decomposed in new reflection-free rotation + z-reflection)
      colz *= -1.;
      HepRep3x3 rottemp(colx.x(),coly.x(),colz.x(),
                  colx.y(),coly.y(),colz.y(),
		  colx.z(),coly.z(),colz.z()); //matrix representation (inverted)
      *rotMat = G4RotationMatrix(rottemp);
      pvName += "_refl";
    }
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
      G4cout << " calling DumpRotationMatrix " << rotMat << G4endl;
#endif
    G4String rotName = DumpRotationMatrix( rotMat );
    G4ThreeVector pos = pv->GetTranslation();
  
    G4String fullname = pvName
      +"#"+itoa(pv->GetCopyNo())
      +"/"+pv->GetMotherLogical()->GetName();
    if( !CheckIfPhysVolExists(fullname, pv )) {
      (*theFile) << ":PLACE " << SubstituteRefl(AddQuotes(pvName)) << " "
		 << pv->GetCopyNo() << " "
		 << SubstituteRefl(AddQuotes(pv->GetMotherLogical()->GetName())) << " " 
		 << AddQuotes(rotName) << " " 
		 << pos.x() << " " 
		 << pos.y() << " " 
		 << pos.z() << G4endl;
      
      thePhysVols[fullname] = pv;
    }
  }

  if( bVolExists ){
    //---- Construct PV's who has this LV as mother
    std::vector<G4VPhysicalVolume*> pvChildren = GetPVChildren( lv );
    std::vector<G4VPhysicalVolume*>::const_iterator ite;
    for( ite = pvChildren.begin(); ite != pvChildren.end(); ite++ ){
      DumpPhysVol( *ite );
    }
  }

}

//------------------------------------------------------------------------
void G4tgbGeometryDumper::DumpLogVol( G4LogicalVolume* lv )
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << "Log vol: " << lv->GetName() << lv->GetSolid()->GetName() << G4endl;
#endif
  G4VSolid* solid;
  //--- take out the '_refl' in the name
  G4String lvName = lv->GetName();

  /* G4ReflectionFactory* reffact = G4ReflectionFactory::Instance();
  if( reffact->IsReflected( lv ) ){
    solid = reffact->GetReflectedLV( lv )->GetSolid();
  } else {
    solid = lv->GetSolid();
    }*/
  
  solid = lv->GetSolid();
  //---- Dump solid 
  DumpSolid( solid );

  //---- Dump material
  G4Material* mate = lv->GetMaterial();
  DumpMaterial( mate );

  //---- Dump logical volume (solid + material)
  (*theFile) << ":VOLU " << SubstituteRefl(AddQuotes(lvName)) << " "
	     << SupressRefl(AddQuotes(solid->GetName())) << " " << AddQuotes(mate->GetName()) << G4endl;

  theLogVols[lvName] = lv;

}


//------------------------------------------------------------------------
void G4tgbGeometryDumper::DumpMaterial( G4Material* mat )
{
  if( CheckIfMaterialExists( mat->GetName(), mat ) ) return;

  //  if( theElements.find( ele ) != theElements.end() ) return; //do not repeat material names

  size_t numElements           = mat->GetNumberOfElements();
  G4double density             = mat->GetDensity()/g*cm3;
  const G4ElementVector* elems = mat->GetElementVector();
  const G4double* fractions    = mat->GetFractionVector();
  
  // start tag
  if (numElements == 1) {
    (*theFile) << ":MATE " << AddQuotes(mat->GetName()) << " " << mat->GetZ() << " " << mat->GetA()/(g/mole) << " " << density << G4endl;
  } else {
    (*theFile) << ":MIXT "<< AddQuotes(mat->GetName()) << " " << density << " " << numElements << G4endl;
    // close start element tag and get ready to do composit "parts"
    for (size_t ii = 0; ii < numElements; ii++) {
      (*theFile) << "   " << AddQuotes((*elems)[ii]->GetName()) << " " << fractions[ii] << G4endl;
    }

    for (size_t ii = 0; ii < numElements; ii++) {
      DumpElement( (*elems)[ii] );
    }
  }

  theMaterials[mat->GetName()] = mat;

}

//------------------------------------------------------------------------
void G4tgbGeometryDumper::DumpElement( G4Element* ele)
{
  if( CheckIfElementExists( ele->GetName(),ele ) ) return;

  if( ele->GetNumberOfIsotopes() != 0 ){
    G4Exception("G4tgbGeometryDumper::DumpElement Elements from isotopes not supported yet ");
  }

  //----- Material mixtures store the components as elements (even if the input are materials), but without symbol
  G4String symbol = ele->GetSymbol();
  if( symbol == "" || symbol == " " ) {
    symbol = ele->GetName();
  }
  (*theFile) << ":ELEM " << AddQuotes(ele->GetName()) << " "  << AddQuotes(symbol) << " " << ele->GetZ() << " " << ele->GetA()/(g/mole) << " " << G4endl;
  
  theElements[ele->GetName()] = ele;

}


//------------------------------------------------------------------------
void G4tgbGeometryDumper::DumpSolid( G4VSolid* solid )
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " G4tgbGeometryDumper::DumpSolid ? " << solid->GetName() << G4endl;
#endif
  if( CheckIfSolidExists( solid->GetName(), solid) ) return;

  G4String solidType = solid->GetEntityType();
  solidType = GetTGSolidType( solidType );

  if (solidType == "UNIONSOLID") {
    DumpBooleanVolume( "UNION", solid );

  } else if (solidType == "SUBTRACTIONSOLID") {
    DumpBooleanVolume( "SUBS", solid );

  } else if (solidType == "INTERSECTIONSOLID") {
    DumpBooleanVolume( "INTERS", solid );

  } else if (solidType == "REFLECTEDSOLID") {
    G4ReflectedSolid* solidrefl = dynamic_cast<G4ReflectedSolid*>(solid);
    G4VSolid* solidori = solidrefl->GetConstituentMovedSolid();
    DumpSolid( solidori );
  } else {
    (*theFile) << ":SOLID " << AddQuotes(solid->GetName()) << " ";
    (*theFile) << AddQuotes(solidType) << " ";
    DumpSolidParams( solid );

    theSolids[solid->GetName()] = solid;
  }

}

//------------------------------------------------------------------------
void G4tgbGeometryDumper::DumpBooleanVolume( const G4String& solidType, G4VSolid* so )
{
  G4BooleanSolid * bso = dynamic_cast < G4BooleanSolid * > (so);
  G4VSolid* solid0 = bso->GetConstituentSolid( 0 );
  G4VSolid* solid1 = bso->GetConstituentSolid( 1 );
  G4DisplacedSolid* solid1Disp;
  G4bool displaced = dynamic_cast<G4DisplacedSolid*>(solid1);
  if( displaced ) {
    solid1Disp = dynamic_cast<G4DisplacedSolid*>(solid1);
    solid1 = solid1Disp->GetConstituentMovedSolid();
  }
  DumpSolid( solid0 );
  DumpSolid( solid1 );

  G4String rotName;
  G4ThreeVector pos;
  if( displaced ) {
    rotName = DumpRotationMatrix( &(solid1Disp->GetTransform().NetRotation()) );// rotation is of mother frame, no need of invert
    pos = -solid1Disp->GetTransform().NetTranslation(); // translation is of mother frame
  } else { // no displacement
    rotName = DumpRotationMatrix( new G4RotationMatrix );
    pos = G4ThreeVector();
  }

  if( CheckIfSolidExists( bso->GetName(), bso) ) return;

  (*theFile) << ":SOLID " 
	     << AddQuotes(bso->GetName()) << " " 
	     << AddQuotes(solidType) << " " 
	     << AddQuotes(solid0->GetName()) << " " 
	     << AddQuotes(solid1->GetName()) << " " 
	     << AddQuotes(rotName) << " " 
	     << approxTo0(pos.x()) << " " 
	     << approxTo0(pos.y()) << " " 
	     << approxTo0(pos.z()) << " " << G4endl;

  theSolids[bso->GetName()] = bso;

}


//------------------------------------------------------------------------
void G4tgbGeometryDumper::DumpSolidParams( G4VSolid * so) 
{
  G4String solidType = so->GetEntityType();
  solidType = GetTGSolidType( solidType );

  if (solidType == "BOX")  {
    G4Box * sb = dynamic_cast<G4Box*>(so);
    (*theFile) << sb->GetXHalfLength() << " " 
	    << sb->GetYHalfLength() << " " 
	    << sb->GetZHalfLength() << " " << G4endl;

  } else if (solidType == "TUBS") {
    G4Tubs * tu = dynamic_cast < G4Tubs * > (so);
    (*theFile) << tu->GetInnerRadius()   << " "
	    << tu->GetOuterRadius()   << " "
	    << tu->GetZHalfLength()   << " "
	    << tu->GetStartPhiAngle()/deg << " "
	    << tu->GetDeltaPhiAngle()/deg << G4endl;

  } else if (solidType == "CONS") {
    G4Cons * cn = dynamic_cast < G4Cons * > (so);
    (*theFile) << cn->GetInnerRadiusMinusZ() << " " 
	       << cn->GetOuterRadiusMinusZ() << " "
	       << cn->GetInnerRadiusPlusZ()  << " "	    
	       << cn->GetOuterRadiusPlusZ()  << " "
	       << cn->GetZHalfLength() << " "
	       << cn->GetStartPhiAngle()/deg  << " "
	       << cn->GetDeltaPhiAngle()/deg  << G4endl;
  } else if (solidType == "POLYCONE") {
    //--- Dump RZ corners, as original parameters will not be present if it was build from RZ corners
    G4Polycone * pc = dynamic_cast < G4Polycone * > (so);
    
    double angphi = pc->GetStartPhi()/deg;
    if( angphi > 180*deg ) angphi -= 360*deg;
    G4int ncor = pc->GetNumRZCorner();
    (*theFile) << angphi << " "
	       << pc->GetOriginalParameters()->Opening_angle/deg << " " << ncor << G4endl;
  
    for( size_t ii = 0; ii < ncor; ii++ ){
      (*theFile) << pc->GetCorner(ii).r << " " << pc->GetCorner(ii).z << G4endl;
    }

    //    G4cout << " POLYCONE " << pc->GetName() << *pc << G4endl;

  }else if (solidType == "POLYHEDRA") {
    //--- Dump RZ corners, as original parameters will not be present if it was build from RZ corners
    //      bool isOpen = (dynamic_cast < G4Polyhedra * > (so))->IsOpen();
    G4Polyhedra * ph = (dynamic_cast < G4Polyhedra * > (so));
    
    double angphi = ph->GetStartPhi()/deg;
    if( angphi > 180*deg ) angphi -= 360*deg;

    G4int ncor = ph->GetNumRZCorner();

    (*theFile) << angphi << " "
	    << ph->GetOriginalParameters()->Opening_angle/deg << " " 
	       << ph->GetNumSide() << " " << ncor << G4endl;

    for( size_t ii = 0; ii < ncor; ii++ ){
      (*theFile) << ph->GetCorner(ii).r << " " << ph->GetCorner(ii).z <<G4endl;
    }

  } else if (solidType == "TRAP") {
    G4Trap * trp = dynamic_cast < G4Trap * > (so);
    G4ThreeVector symAxis(trp->GetSymAxis());
    double theta, phi;
    theta = symAxis.theta()/deg;
    phi = symAxis.phi()/deg;
    (*theFile) << trp->GetZHalfLength()     << " "
	       << theta << " " 
	       << phi<< " "
	       << trp->GetYHalfLength1() << " "
	       << trp->GetXHalfLength1() << " "
	       << trp->GetXHalfLength2() << " "	    
	       << atan(trp->GetTanAlpha1())/deg << " " 
	       << trp->GetYHalfLength2()    << " "
	       << trp->GetXHalfLength3()    << " "
	       << trp->GetXHalfLength4()    << " "	    
	       << atan(trp->GetTanAlpha2())/deg << " "
	       << G4endl;
  } else if (solidType == "TRD") {
    G4Trd * tr = dynamic_cast < G4Trd * > (so);
    (*theFile) << tr->GetZHalfLength()<< " "
	    << tr->GetYHalfLength1() << " "
            << tr->GetYHalfLength2() << " " 
	    << tr->GetXHalfLength1() << " "
	    << tr->GetXHalfLength2() << G4endl;
    
  } else if (solidType == "ORB") {
    G4Orb * orb = dynamic_cast < G4Orb * > (so);
    (*theFile) << orb->GetRadius()  << G4endl;

  }else {
    G4Exception("G4tgbGeometryDumpe::DumpSolidParams   solid type not supported " + solidType );
  }
    
}   


//------------------------------------------------------------------------
G4String G4tgbGeometryDumper::DumpRotationMatrix( G4RotationMatrix* rotm )
{
  G4String rotName = LookForExistingRotation( rotm );
  if( rotName != "" ) return rotName;

  if (!rotm) {
    rotm = new G4RotationMatrix();
  } 

  G4ThreeVector v(1.,1.,1.);
  double de = MatDeterminant(rotm);
  if (de < -0.9 ) { // a reflection ....
    (*theFile) << ":ROTM ";
    rotName = "RRM";
    rotName += itoa(theRotationNumber++);
 
    (*theFile) << AddQuotes(rotName) << std::setprecision(9) << " " 
	       << approxTo0(rotm->xx())  << " "
	       << approxTo0(rotm->xy())  << " "
	       << approxTo0(rotm->xz())  << " "
	       << approxTo0(rotm->yx())  << " "
	       << approxTo0(rotm->yy())  << " "
	       << approxTo0(rotm->yz())  << " "
	       << approxTo0(rotm->zx())  << " "
	       << approxTo0(rotm->zy())  << " "
	       << approxTo0(rotm->zz())  << G4endl;
  } else if(de > 0.9 ) { // a rotation
    (*theFile) << ":ROTM ";
    rotName = "RM";
    rotName += itoa(theRotationNumber++);
    
    (*theFile) << AddQuotes(rotName) << std::setprecision(9) << " " 
	       << approxTo0(rotm->xx())  << " "
	       << approxTo0(rotm->xy())  << " "
	       << approxTo0(rotm->xz())  << " "
	       << approxTo0(rotm->yx())  << " "
	       << approxTo0(rotm->yy())  << " "
	       << approxTo0(rotm->yz())  << " "
	       << approxTo0(rotm->zx())  << " "
	       << approxTo0(rotm->zy())  << " "
	       << approxTo0(rotm->zz())  << G4endl;
    /*t    (*theFile) << AddQuotes(rotName) << " " 
	       << approxTo0(rotm->thetaX()/deg)  << " "
	       << approxTo0(rotm->phiX()/deg)    << " "
	       << approxTo0(rotm->thetaY()/deg)  << " "
	       << approxTo0(rotm->phiY()/deg)    << " "
	       << approxTo0(rotm->thetaZ()/deg)  << " "
	       << approxTo0(rotm->phiZ()/deg)    << G4endl;
    */
  }
  
  
  theRotMats[rotName] = rotm;

  return rotName;

}

//------------------------------------------------------------------------
std::vector<G4VPhysicalVolume*> G4tgbGeometryDumper::GetPVChildren( G4LogicalVolume* lv )
{
  G4PhysicalVolumeStore* pvstore = G4PhysicalVolumeStore::GetInstance();
  G4PhysicalVolumeStore::const_iterator ite;
  std::vector<G4VPhysicalVolume*> children;
  for( ite = pvstore->begin(); ite != pvstore->end(); ite++ ){
    if( (*ite)->GetMotherLogical() == lv ) {
      children.push_back( *ite );
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " adding children " << (*ite)->GetName() << " of " << lv->GetName() <<  G4endl;
#endif
    }
  }

  return children;
}

//------------------------------------------------------------------------
G4String G4tgbGeometryDumper::GetTGSolidType( G4String& solidType )
{
  G4String newsolidType = solidType.substr(2,solidType.length() );
  for( size_t ii = 0; ii < newsolidType.length(); ii++ ){
    newsolidType[ii] = toupper(newsolidType[ii] );
  }
  return newsolidType;
}

//------------------------------------------------------------------------
double G4tgbGeometryDumper::MatDeterminant(G4RotationMatrix * ro) 
{
   HepRep3x3 r = ro->rep3x3();
   return       r.xx_*(r.yy_*r.zz_ - r.zy_*r.yz_)
              - r.yx_*(r.xy_*r.zz_ - r.zy_*r.xz_)
	      + r.zx_*(r.xy_*r.yz_ - r.yy_*r.xz_);
}

//------------------------------------------------------------------------
G4String G4tgbGeometryDumper::itoa(int current)
{
  const char theDigits[11] = "0123456789";
  G4String result;
  int digit;
  do
    {
      digit = current-10*(current/10);
      result=theDigits[digit]+result;
      current/=10;
    }
  while(current!=0);
  return result;
}


//-----------------------------------------------------------------------
G4double G4tgbGeometryDumper::approxTo0( double val )
{
  double precision = kCarTolerance;
  if( fabs(val) < precision ) val = 0;
  return val;
}


//-----------------------------------------------------------------------
G4String G4tgbGeometryDumper::AddQuotes( const G4String& str )
{
  //--- look if there is a separating blank
  G4bool bBlank = FALSE;
  size_t siz = str.length();
  for( size_t ii = 0; ii < siz; ii++ ){
    if( str.substr(ii,1) == " " ) {
      bBlank = TRUE;
      break;
    }
  }
  G4String str2 = str;
  if( bBlank ) {
    str2 = G4String("\"") + str2 + G4String("\"");
  }
  return str2;
}


//------------------------------------------------------------------------
G4String G4tgbGeometryDumper::SupressRefl( G4String name )
{
  size_t irefl = name.rfind("_refl");
  if( irefl != -1 ) {
    name = name.substr( 0, irefl );
  }
  return name;
}

//------------------------------------------------------------------------
G4String G4tgbGeometryDumper::SubstituteRefl( G4String name )
{
  size_t irefl = name.rfind("_refl");
  if( irefl != -1 ) {
    name = name.substr( 0, irefl ) + "_REFL";
  }
  return name;
}

//------------------------------------------------------------------------
G4bool G4tgbGeometryDumper::CheckIfElementExists( const G4String& name, G4Element* pt )
{
  if( theElements.find( name ) != theElements.end() ){
    if( (*(theElements.find(name))).second != pt ){
      G4Exception("G4tgbGeometryDumper::CheckIfElementExists. Element found but not same as before : " + name );
    }
    return 1;
  } else {
    return 0;
  }
}

//------------------------------------------------------------------------
G4bool G4tgbGeometryDumper::CheckIfMaterialExists( const G4String& name, G4Material* pt )
{
  if( theMaterials.find( name ) != theMaterials.end() ){
    if( (*(theMaterials.find(name))).second != pt ){
      //t TEMPORARY: DUMMY material created for boolean solids, until :SOLID exists
      if( name != "DUMMY" ) {
	G4Exception("G4tgbGeometryDumper::CheckIfMaterialExists. Material found but not same as before : " + name );
      }
    }
    return 1;
  } else {
    return 0;
  } 
}

//------------------------------------------------------------------------
G4bool G4tgbGeometryDumper::CheckIfLogVolExists( const G4String& name, G4LogicalVolume* pt )
{
  if( theLogVols.find( name ) != theLogVols.end() ){
    G4LogicalVolume* lvnew = (*(theLogVols.find(name))).second;
    if( lvnew != pt ){
      /*      //---- Reflected volumes are repeated
      G4ReflectionFactory* reffact = G4ReflectionFactory::Instance();
      if( !reffact->IsReflected( pt ) && !reffact->IsReflected( lvnew ) ){
	G4Exception("G4tgbGeometryDumper::CheckIfLogVolExists. LogVol found but not same as before : " + name );
	} */
    }
    return 1;
  } else {
    return 0;
  }
}

//------------------------------------------------------------------------
G4bool G4tgbGeometryDumper::CheckIfSolidExists( const G4String& name, G4VSolid* pt )
{
  if( theSolids.find( name ) != theSolids.end() ){
    if( (*(theSolids.find(name))).second != pt ){
      uint nn = 2;
      G4String solidName = pt->GetName();
      G4String addName = "_ADDED_";
      for( ;; ){
	if(theSolids.find( solidName+ addName + itoa(nn) ) == theSolids.end() ) {
	  addName += itoa(nn);
	  break;
	}
	nn++;
      }
      G4cerr << "G4tgbGeometryDumper::CheckIfSolidExists. Solid found but not same as before : " << name << " Changed to " << name + addName << G4endl;
      pt->SetName(solidName+addName);
      //---- Change name of log vol that contains this solid
      G4LogicalVolumeStore* lvstore = G4LogicalVolumeStore::GetInstance();
      G4LogicalVolumeStore::const_iterator ite;
      for( ite = lvstore->begin(); ite != lvstore->end(); ite++ ){
	if( (*ite)->GetSolid() == pt ) (*ite)->SetName( (*ite)->GetName() +addName );
      }
      return 0;
    }
    return 1;
  } else {
    return 0;
  }
}

//-----------------------------------------------------------------------
G4bool G4tgbGeometryDumper::CheckIfPhysVolExists( const G4String& name, G4VPhysicalVolume* pt )
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    G4cout << " G4tgbGeometryDumper::CheckIfPhysVolExists " << name << G4endl;
#endif
  if( thePhysVols.find( name ) != thePhysVols.end() ){
    if( (*(thePhysVols.find(name))).second != pt ){
      //      G4Exception("G4tgbGeometryDumper::CheckIfPhysVolExists. Placement found but not same as before : " + name );
      G4cerr << " G4tgbGeometryDumper::CheckIfPhysVolExists  Placement found but not same as before : " << name << G4endl;
    }
    return 1;
  } else {
    return 0;
  }
}

//-----------------------------------------------------------------------
G4String G4tgbGeometryDumper::LookForExistingRotation( const G4RotationMatrix* rotm )
{
  G4String rmName = "";

  std::map<G4String,G4RotationMatrix*>::const_iterator ite;
  for( ite = theRotMats.begin(); ite != theRotMats.end(); ite++ ){
    if( (*ite).second->isNear( *rotm ) ){
      rmName = (*ite).first;
      /*      (*theFile) << " near matrices " << rmName << " " 
	     << " new " << rotm->colX() << " " << rotm->colY() << " " << rotm->colZ() << " " << G4endl
	     << " old " << (*ite).second->colX() << " " << (*ite).second->colY() << " " << (*ite).second->colZ() << " " << G4endl;
      */
      break;
    }
  }
  return rmName;
}

