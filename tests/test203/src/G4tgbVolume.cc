#include "G4tgbVolume.hh"
#include "G4tgbVolumeMgr.hh"
#include "G4tgbMaterialMgr.hh"
#include "G4tgbRotationMatrixMgr.hh"
#include "G4tgbPlaceParamLinear.hh"
#include "G4tgbPlaceParamCircle.hh"

#include "G4tgrSolid.hh"
#include "G4tgrSolidBoolean.hh"
#include "G4tgrVolume.hh"
#include "G4tgrVolumeDivision.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgrPlace.hh"
#include "G4tgrPlaceDivRep.hh"
#include "G4tgrPlaceParameterisation.hh"
#include "G4tgrUtils.hh"

#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4Trd.hh"
#include "G4Para.hh"
#include "G4Torus.hh"
#include "G4Hype.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4EllipticalTube.hh"
#include "G4Ellipsoid.hh"
#include "G4EllipticalCone.hh"
#include "G4Hype.hh"
#include "G4Tet.hh"
#include "G4TwistedBox.hh"
#include "G4TwistedTrap.hh"
#include "G4TwistedTrd.hh"
#include "G4TwistedTubs.hh"

#include "G4Material.hh"
#include "G4RotationMatrix.hh"
#include "G4ReflectionFactory.hh"

#include "G4VisAttributes.hh"
#include "G4RegionStore.hh"
#include "G4tgrMessenger.hh"

using namespace CLHEP;

//-------------------------------------------------------------------
G4tgbVolume::G4tgbVolume( G4tgrVolume* vol)
{
  theTgrVolume = vol;
}


//-------------------------------------------------------------------
void G4tgbVolume::ConstructG4Volumes( const G4tgrPlace* place, const G4LogicalVolume* parentLV )
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << theTgrVolume->GetName() << " G4tgbVolume::ConstructG4Volumes " << parentLV << endl;
#endif
  G4tgbVolumeMgr* g4vmgr = G4tgbVolumeMgr::GetInstance();
  G4LogicalVolume* logvol = g4vmgr->FindG4LogVol( GetName() );
  if( logvol == 0 ) {
    //--- If first time build solid and LogVol
    //    G4cout << " calling ConstrucG4Solid " << GetName() << G4endl;
    G4VSolid* solid = FindOrConstructG4Solid( theTgrVolume->GetSolid() ); 
    g4vmgr->RegisterMe( solid );
    
    logvol = ConstructG4LogVol( solid );
    g4vmgr->RegisterMe( logvol );
    g4vmgr->RegisterChildParentLVs( logvol, parentLV ); 
    
    //--- If first copy build children placements in this LogVol
    pair<mmapspl::iterator, mmapspl::iterator> children = G4tgrVolumeMgr::GetInstance()->GetChildren( GetName() );
    mmapspl::iterator cite; 
    //-    cout << "children " << &(*(children.second))  << " - " << &(*(children.first)) << endl;
    for( cite = children.first; cite != children.second; cite++ ) {
      //----- Call G4tgrPlace ->constructG4Volumes 
      //---- find G4tgbVolume corresponding to the G4tgrVolume pointed by G4tgrPlace
      G4tgrPlace* pl = const_cast<G4tgrPlace*>((*cite).second);
      //-      cout << " pl " << pl << (pl->vol() ) << endl;
      //-  cout << " pl " << pl->vol()->GetName() << endl;
      G4tgbVolume* svol = g4vmgr->FindVolume( pl->GetVolume()->GetName() );
      //--- find copyNo
      svol->ConstructG4Volumes( pl, logvol );
    }
  }

  //--- Construct PhysVol
  G4VPhysicalVolume* physvol = ConstructG4PhysVol( place, logvol, parentLV );
  //  cout << "onstructG4PV( called " <<physvol <<  endl;
  //   cout << "onstructG4PV( called " << physvol->GetName() << endl;
  g4vmgr->RegisterMe( physvol );

}


//-------------------------------------------------------------------
G4VSolid* G4tgbVolume::FindOrConstructG4Solid( const G4tgrSolid* sol ) 
{

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) { 
    G4cout << this << "G4tgbVolume::FindOrConstructG4Solid  SOLID= " << sol;
    G4cout << " " << sol->GetName() << " of type " << sol->GetType() << G4endl;
  }
#endif 

  //----- Check if solid exists already
  G4VSolid* solid = G4tgbVolumeMgr::GetInstance()->FindG4Solid( sol->GetName() );
  if( solid ) return solid;

  // Give sol as boolean solids needs to call this method twice

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbVolume::FindOrConstructG4Solid " << " " << sol->GetSolidParams().size() << endl;
#endif
  
  vector<double> solParam;
  // In case of BOOLEAN solids, solidParams are taken from components
  if( sol->GetSolidParams().size() == 1) { 
    //-   cout << "1  solidpara " << endl;
    solParam = * sol->GetSolidParams()[ 0 ];
    //-    cout << " paramt hs "  <<  sol->solidParams()[ 0 ] << " " <<  (sol->solidParams()[ 0 ])->size() << endl;
  }

  //----------- instantiate the appropiate G4VSolid type
  G4String stype = sol->GetType();
  G4String sname = sol->GetName();
  if( stype == "BOX" ) {
    //    G4cout << " box " <<  stype << " " << solParam.size() << G4endl;
    CheckNoSolidParams( stype, 3, solParam.size() );
    solid = new G4Box( sname, solParam[0], solParam[1], solParam[2] ); 
  } else if( stype == "TUBE" ) {
    CheckNoSolidParams( stype, 3, solParam.size() );
    solid = new G4Tubs( sname, solParam[0], solParam[1], solParam[2], 0.*deg, 360.*deg ); 
  } else if( stype == "TUBS" ) {
    CheckNoSolidParams( stype, 5, solParam.size() );
    solid = new G4Tubs( sname, solParam[0], solParam[1], solParam[2], solParam[3], solParam[4] ); 
  } else if( stype == "TRAP" ) {
    CheckNoSolidParams( stype, 11, solParam.size() );
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 3 ) 
      G4cout << " new G4Trap " << sname <<" "<< solParam[0]<<" "<< solParam[1]<<" "<< solParam[2]
	     <<" "<< solParam[3]<<" "<< solParam[4]<<" "<< solParam[5]<<" "<< solParam[6]
	     <<" "<< solParam[7]<<" "<< solParam[8]<<" "<< solParam[9]<<" "<< solParam[10] << G4endl;
#endif
    /*    for( size_t ii = 0; ii < 11; ii++ ){
      if( solParam[ii] == 0. ) solParam[ii] = 1.e-6*mm;
      } */
    solid = new G4Trap( sname, solParam[0], solParam[1], solParam[2]
			, solParam[3], solParam[4], solParam[5], solParam[6]
			, solParam[7], solParam[8], solParam[9], solParam[10] );
  } else if( stype == "TRD" ) {
    CheckNoSolidParams( stype, 5, solParam.size() );
    solid = new G4Trd( sname, solParam[0], solParam[1], solParam[2]
			, solParam[3], solParam[4] );
  } else if( stype == "PARA" ) {
    CheckNoSolidParams( stype, 6, solParam.size() );
    solid = new G4Para( sname, solParam[0], solParam[1], solParam[2]
			, solParam[3], solParam[4], solParam[5] );
  } else if( stype == "CONE" ) {
    CheckNoSolidParams( stype, 5, solParam.size() );
    solid = new G4Cons( sname, solParam[0], solParam[1], solParam[2]
			, solParam[3], solParam[4], 0., 360.*deg);
  } else if( stype == "CONS" ) {
    CheckNoSolidParams( stype, 7, solParam.size() );
    solid = new G4Cons( sname, solParam[0], solParam[1], solParam[2]
			, solParam[3], solParam[4], solParam[5], solParam[6]);
  } else if( stype == "SPHERE" ) {
    CheckNoSolidParams( stype, 6, solParam.size() );
    solid = new G4Sphere( sname, solParam[0], solParam[1], solParam[2]
			, solParam[3], solParam[4], solParam[5]);
  } else if( stype == "ORB" ) {
    CheckNoSolidParams( stype, 1, solParam.size() );
    solid = new G4Orb( sname, solParam[0] );
  } else if( stype == "TORUS" ) {
    CheckNoSolidParams( stype, 5, solParam.size() );
    solid = new G4Torus( sname, solParam[0], solParam[1], solParam[2]
			, solParam[3], solParam[4] );
  } else if( stype == "POLYCONE" ) {
    int nplanes = int(solParam[2]);
    //    G4cout << " POLYCONE nplanes " << nplanes << " phimin " << solParam[0] << " phimax " << solParam[1] << " size " << solParam.size() << G4endl;
    G4bool genericPoly = false;
    if( solParam.size() == 3+nplanes*3 ){ 
      genericPoly = true;
    }else if( solParam.size() == 3+nplanes*2 ){ 
      genericPoly = false;
    } else {
      G4cerr << "G4tgbVolume::CheckNoSolidParams, solid type " << stype << " should have " << 3+nplanes*3 << " (Z,Rmin,Rmax) or " << 3+nplanes*2 << " (RZ corners) parameters, and it has " <<  solParam.size() << G4endl;
      G4Exception("");
    }

    if( genericPoly ) {
      vector<double>* z_p = new vector<double>;
      vector<double>* rmin_p = new vector<double>;
      vector<double>* rmax_p = new vector<double>;
      //    vector<double>::const_iterator i = solParam.begin()+2;
      for( uint ii = 0; ii < nplanes; ii++ ){
	(*z_p).push_back( solParam[3+3*ii] );
	(*rmin_p).push_back( solParam[3+3*ii+1] );
	(*rmax_p).push_back(  solParam[3+3*ii+2] );
      }
      solid = new G4Polycone( sname, solParam[0], solParam[1], // start,delta-phi
			      nplanes, // sections
			      &((*z_p)[0]),
			      &((*rmin_p)[0]),
			      &((*rmax_p)[0]));
    }else {
      vector<double>* R_c = new vector<double>;
      vector<double>* Z_c = new vector<double>;
      for( uint ii = 0; ii < nplanes; ii++ ){
	(*R_c).push_back( solParam[3+2*ii] );
	(*Z_c).push_back( solParam[3+2*ii+1] );
      }
      solid = new G4Polycone( sname, solParam[0], solParam[1], // start,delta-phi
			      nplanes, // sections
			      &((*R_c)[0]),
			      &((*Z_c)[0]));
    }
  } else if( stype == "POLYHEDRA" ) {
    int nplanes = int(solParam[3]);
    G4bool genericPoly = false;
    if( solParam.size() == 4+nplanes*3 ){ 
      genericPoly = true;
    }else if( solParam.size() == 4+nplanes*2 ){ 
      genericPoly = false;
    } else {
      G4cerr << "G4tgbVolume::CheckNoSolidParams, solid type " << stype << " should have " << 4+nplanes*3 << " (Z,Rmin,Rmax) or " << 4+nplanes*2 << " (RZ corners) parameters, and it has " <<  solParam.size() << G4endl;
      G4Exception("");
    }
    
    if( genericPoly ) {
      vector<double>* z_p = new vector<double>;
      vector<double>* rmin_p = new vector<double>;
      vector<double>* rmax_p = new vector<double>;
      //    vector<double>::const_iterator i = solParam.begin()+2;
      for( uint ii = 0; ii < nplanes; ii++ ){
	(*z_p).push_back( solParam[4+3*ii] );
	(*rmin_p).push_back( solParam[4+3*ii+1] );
	(*rmax_p).push_back(  solParam[4+3*ii+2] );
      }
      solid = new G4Polyhedra( sname, solParam[0], solParam[1], int(solParam[2]), nplanes,
			       &((*z_p)[0]),
			       &((*rmin_p)[0]),
			       &((*rmax_p)[0]));
    }else {
      vector<double>* R_c = new vector<double>;
      vector<double>* Z_c = new vector<double>;
      for( uint ii = 0; ii < nplanes; ii++ ){
	(*R_c).push_back( solParam[4+2*ii] );
	(*Z_c).push_back( solParam[4+2*ii+1] );
      }
      solid = new G4Polyhedra( sname, solParam[0], solParam[1], int(solParam[2]), nplanes,
			       &((*R_c)[0]),
			       &((*Z_c)[0]));
    }

  } else if( stype == "ELLIPTICAL_TUBE" ) {
    CheckNoSolidParams( stype, 3, solParam.size() );
    solid = new G4EllipticalTube( sname, solParam[0], solParam[1], solParam[2]);
  } else if( stype == "ELLIPSOID" ) {
    CheckNoSolidParams( stype, 5, solParam.size() );
    solid = new G4Ellipsoid( sname, solParam[0], solParam[1], solParam[2]
			, solParam[3], solParam[4] );
  } else if( stype == "ELLIPTICAL_CONE" ) {
    CheckNoSolidParams( stype, 4, solParam.size() );
    solid = new G4EllipticalCone( sname, solParam[0], solParam[1], solParam[2]
			, solParam[3] );
  } else if( stype == "HYPE" ) {
    CheckNoSolidParams( stype, 5, solParam.size() );
    solid = new G4Hype( sname, solParam[0], solParam[1], solParam[2] 
			, solParam[3], solParam[4] );
  } else if( stype == "TET" ) {
    CheckNoSolidParams( stype, 12, solParam.size() );
    G4ThreeVector anchor(solParam[0], solParam[1], solParam[2]);
    G4ThreeVector p2(solParam[3], solParam[4], solParam[5]);
    G4ThreeVector p3(solParam[6], solParam[7], solParam[8]);
    G4ThreeVector p4(solParam[9], solParam[10], solParam[11]);
    solid = new G4Tet( sname, anchor, p2, p3, p4 );
  } else if( stype == "TWISTED_BOX" ) {
    CheckNoSolidParams( stype, 4, solParam.size() );
    solid = new G4TwistedBox( sname, solParam[0], solParam[1], solParam[2]
			      , solParam[3]);
  } else if( stype == "TWISTED_TRAP" ) {
    CheckNoSolidParams( stype, 11, solParam.size() );
    solid = new G4TwistedTrap( sname, solParam[0], solParam[1], solParam[2]
			, solParam[3], solParam[4], solParam[5], solParam[6]
			, solParam[7], solParam[8], solParam[9], solParam[10] );
  } else if( stype == "TWISTED_TRD" ) {
    CheckNoSolidParams( stype, 6, solParam.size() );
    solid = new G4TwistedTrd( sname, solParam[0], solParam[1], solParam[2]
			, solParam[3], solParam[4], solParam[5]);
  } else if( stype == "TWISTED_TUBS" ) {
    CheckNoSolidParams( stype, 5, solParam.size() );
    solid = new G4TwistedTubs( sname, solParam[0], solParam[1], solParam[2]
			, solParam[3], solParam[4]);
  } else if( stype.substr(0,7) == "Boolean" ) {
    const G4tgrSolidBoolean* solb = dynamic_cast<const G4tgrSolidBoolean*>(sol);
    //-    cout << " first solid " << sol->getSOL(0)->GetName() << endl;
    G4VSolid* sol1 = FindOrConstructG4Solid( solb->GetSolid(0));
    //-    cout << " second solid " << sol->getSOL(1)->GetName() << endl;
    G4VSolid* sol2 = FindOrConstructG4Solid( solb->GetSolid(1));
    G4RotationMatrix* relRotMat = G4tgbRotationMatrixMgr::GetInstance()
      ->FindOrBuildG4RotMatrix( sol->GetRelativeRotMatName() );
    //-    cout << " relrotmat " << relRotMat << endl;
    Hep3Vector relPlace = solb->GetRelativePlace();
    //-    cout << " relplace " << relPlace << endl;

    if( stype == "Boolean_UNION" ){
      solid = new G4UnionSolid( sname, sol1, sol2, relRotMat, relPlace );
    } else if( stype == "Boolean_SUBS" ){
      solid = new G4SubtractionSolid( sname, sol1, sol2, relRotMat, relPlace );
    } else if( stype == "Boolean_INTERS" ){
      solid = new G4IntersectionSolid( sname, sol1, sol2, relRotMat, relPlace );
    } else { 
      G4Exception("G4tgbVolume::FindOrConstructG4Solid:BOOLEAN  unknown boolean type " + stype);
    }

  } else {
    G4Exception("G4tgbVolume::FindOrConstructG4Solid. type of solid " 
	 + stype + " not implemented yet \n"
	 + " Please, if you think it is necessary, implement it yourself or ask the responsible person to do it ");
  } 
  
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbVolume::FindOrConstructG4Solid: created solid " << sname << " of type " << solid->GetEntityType() << endl; 
#endif
 
  return solid;
  
}


//-------------------------------------------------------------------
void G4tgbVolume::CheckNoSolidParams( const G4String& solidType, const int NoParamExpected, const uint NoParam )
{
  if( NoParamExpected != NoParam ) {
    G4cerr << "G4tgbVolume::CheckNoSolidParams, solid type " << solidType << " should have " << NoParamExpected << " parameters, and it has " << NoParam << G4endl;
    G4Exception("");
  }
}


//-------------------------------------------------------------------
G4LogicalVolume* G4tgbVolume::ConstructG4LogVol( const G4VSolid* solid )
{
  G4LogicalVolume* logvol;

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbVolume::ConstructG4LogVol " << theTgrVolume->GetName() << endl;
#endif

  //----------- get the material first
  G4Material* mate = G4tgbMaterialMgr::GetInstance()->FindOrBuildG4Material( theTgrVolume->GetMaterialName() );
  if( mate == 0 ) {
    G4Exception("G4tgbVolume::ConstructG4LogVol. material not found " + theTgrVolume->GetMaterialName() + " for volume " + theTgrVolume->GetName());
  }
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << "G4tgbVolume::ConstructG4LogVol. mate constructed " << mate->GetName() << endl; 
#endif
 
  //---------- construct the LV
  logvol = new G4LogicalVolume( const_cast<G4VSolid*>(solid), const_cast<G4Material*>(mate), theTgrVolume->GetName() );

  //---------- Set visibility and colour
  if( !GetVisibility() || GetColour()[0] != -1 ) {
    G4VisAttributes* visAtt = new G4VisAttributes();
    if( !GetVisibility() ) {
      visAtt->SetVisibility( GetVisibility() );
    } else if( GetColour()[0] != -1 ) { //this else should not be necessary, because if the visibility is set to off, colour should have no efect. But it does not work for G4: if you set colour and vis off, it is visualized!?!?!?
      const double* col = GetColour();
      visAtt->SetColour( G4Colour(col[0],col[1],col[2]));
    }
    logvol->SetVisAttributes(visAtt);
  }


#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbVolume::ConstructG4LogVol: created logical volume " << theTgrVolume->GetName() << endl;
#endif

  return logvol;
}


//-------------------------------------------------------------------
G4VPhysicalVolume* G4tgbVolume::ConstructG4PhysVol( const G4tgrPlace* place, const G4LogicalVolume* currentLV, const G4LogicalVolume* parentLV )
{
  G4VPhysicalVolume* physvol = 0;
  int copyNo;
  
  //----- Case of placement of top volume
  if( place == 0 ) {
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 ) cout << " G4tgbVolume::ConstructG4PhysVol world: " << GetName() << endl;
#endif
    physvol = new G4PVPlacement(0, G4ThreeVector(), const_cast<G4LogicalVolume*>(currentLV), GetName(), 0, false, 0);
    
  } else { 
    copyNo = place->GetCopyNo();

#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
      cout << " G4tgbVolume::ConstructG4PhysVol " << theTgrVolume->GetName() << " inside " << parentLV->GetName() << " copy No " << copyNo << endl;
#endif
    
    if( theTgrVolume->GetType() == "VOLSimple" ) {
      //----- get placement
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
	cout << " G4tgbVolume::ConstructG4PhysVol Placement type " <<  place->GetType() << G4endl;
#endif
      
      //--------------- If it is  G4tgrPlaceSimple
      if( place->GetType() == "PlaceSimple" ) {
	//----- get rotation matrix
	G4String rmName = place->GetRotMatName();

	G4RotationMatrix* rotmat = G4tgbRotationMatrixMgr::GetInstance()
	  ->FindOrBuildG4RotMatrix( rmName );
	//-	*rotmat = (*rotmat).inverse();
	//----- place volume in mother
	double check = (rotmat->colX().cross(rotmat->colY()))*rotmat->colZ();
	double tol = 1.0e-3;
	//---- Check that matrix is ortogonal
	if (1-abs(check)>tol) {
	  G4cerr << " Matrix : " << rmName << " " << rotmat->colX() << " " << rotmat->colY() << " " << rotmat->colZ() << G4endl
		 << " product x X y * z = " << check << " x X y " << rotmat->colX().cross(rotmat->colY()) << G4endl;
	  G4Exception("G4tgbVolume::ConstructG4PhysVol  Rotation is not ortogonal " + rmName );
	  //---- Check if it is reflection
	} else if (1+check<=tol) {
	  G4Translate3D transl = place->GetPlacement();
	  // G3 convention of defining rot-matrices ...
	  G4Transform3D trfrm  = transl * G4Rotate3D(*rotmat);
	  physvol = (G4ReflectionFactory::Instance()
	    ->Place(trfrm, GetName(), const_cast<G4LogicalVolume*>(currentLV), const_cast<G4LogicalVolume*>(parentLV), false, copyNo, false )).first;
	  //	  G4cout << " building reflected pv " << physvol << G4endl;
	} else {
#ifdef G4VERBOSE
	  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) {
	    G4cout << "PLACE " << GetName() << " # " << copyNo 
		   << " ROT " << rotmat->colX() 
		   << " " << rotmat->colY() 
		   << " " << rotmat->colZ() << G4endl;
	  }
#endif
	  physvol = new G4PVPlacement( rotmat, place->GetPlacement(),
				     const_cast<G4LogicalVolume*>(currentLV), GetName(), const_cast<G4LogicalVolume*>(parentLV), false, copyNo);
	}
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << physvol << " G4tgbVolume::ConstructG4PhysVol  place= " << place->GetPlacement() << " " << GetName() << " in " <<  parentLV->GetName() << endl;
#endif
  
  //--------------- If it is G4tgrPlaceParam
      } else if( place->GetType() == "PlaceParam" ) {
	G4tgrPlaceParameterisation* dp = (G4tgrPlaceParameterisation*)(place);
	//	G4tgrPlaceParam* dp = (G4tgrPlaceParam)(dp);
	//----- see what parameterisation type
#ifdef G4VERBOSE
	if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
	  cout << physvol << " G4tgbVolume::ConstructG4PhysVol param " << GetName() << " in " <<  parentLV->GetName() << endl;
#endif
	G4VPVParameterisation * param;
	
	if(dp->GetParamType() == "CIRCLE" ) { 
#ifdef G4VERBOSE
	  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
	    cout << " noCopies " << dp->GetNumberOfCopies() 
		 << " offset " << dp->GetOffset() << " step " << dp->GetStep() 
		 << " radius "<< dp->GetExtraData()[0] << endl;
#endif
	  param = new G4tgbPlaceParamCircle(dp->GetNumberOfCopies(), dp->GetOffset()
					    , dp->GetStep(), EAxis(dp->GetAxis())
					    , dp->GetExtraData()[0]); 
	  //	  physvol = new G4PVPlacement(0, G4ThreeVector(0,100,0.), currentLV, GetName(),  parentLV, false, 0);
	} else if(dp->GetParamType() == "LINEAR_X" || dp->GetParamType() == "LINEAR_Y" || dp->GetParamType() == "LINEAR_Z" ) { 
#ifdef G4VERBOSE
	  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
	    cout << " noCopies " << dp->GetNumberOfCopies() 
		 << " offset " << dp->GetOffset() << " step " << dp->GetStep()  << endl;
#endif
	  param = new G4tgbPlaceParamLinear(dp->GetNumberOfCopies(), dp->GetOffset()
					    , dp->GetStep(), EAxis(dp->GetAxis())); 
	}
	physvol = new G4PVParameterised(GetName(), const_cast<G4LogicalVolume*>(currentLV), const_cast<G4LogicalVolume*>(parentLV), EAxis(dp->GetAxis()), dp->GetNumberOfCopies(), param);
      }
      
      //--------------- If it is  G4tgrVolumeDivision
    } else if( theTgrVolume->GetType() == "PlaceReplica" ) {
      G4tgrVolumeDivision* volr = (G4tgrVolumeDivision*)theTgrVolume;
      G4tgrPlaceDivRep* dpr = volr->GetPlaceDivision() ;
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
	cout << physvol << " G4tgbVolume::ConstructG4PhysVol replica" << " " << GetName() << " in " <<  parentLV->GetName() 
	     << " " << dpr->GetNDiv() << " " << dpr->GetWidth()
	     << " " << dpr->GetOffset() << endl;
#endif
      physvol = new G4PVReplica(GetName(),const_cast<G4LogicalVolume*>(currentLV),const_cast<G4LogicalVolume*>(parentLV),EAxis(dpr->GetAxis())
				,dpr->GetNDiv(),dpr->GetWidth(),dpr->GetOffset());
    } else if( theTgrVolume->GetType() == "PlaceDivision" ) {
      G4Exception("G4tgbVolume::ConstructG4PhysVol  placement type not supported: " +  theTgrVolume->GetType() + "  . Please contact GEANT4 expert");
    } else {
      G4Exception("G4tgbVolume::ConstructG4PhysVol  placement type not supported: " +  theTgrVolume->GetType() + "  . Please contact GEANT4 expert");
   }
    
  } 

  //  G4cout << " physvol constructed " << physvol << G4endl;
  return physvol;
}

