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
//
// $Id: G4tgbGeometryDumper.cc,v 1.3 2008-10-23 16:20:01 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4tgbGeometryDumper

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

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
#include "G4Sphere.hh"
#include "G4Ellipsoid.hh"
#include "G4PVPlacement.hh"
#include "G4BooleanSolid.hh"
#include "G4ReflectionFactory.hh"
#include "G4ReflectedSolid.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GeometryTolerance.hh"

#include <iomanip>

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
  for( ;; )
  {
    G4LogicalVolume* lv = pv->GetMotherLogical();
    if( lv == 0 ) { break; }

    //----- look for one PV of this LV
    for( ite = pvstore->begin(); ite != pvstore->end(); ite++ )
    {
      pv = (*ite);
      if( pv->GetLogicalVolume() == lv )
      {
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
  {
    G4cout << " G4tgbGeometryDumper::DumpPhysVol() - Volume: "
           << pv->GetName() << G4endl;
  }
#endif

  //---- Dump logical volume first
  G4LogicalVolume* lv = pv->GetLogicalVolume();

  G4ReflectionFactory* reffact = G4ReflectionFactory::Instance();

  //---- It is not needed to dump _refl volumes created when parent is reflected 
  // !!WARNING : it must be avoided to reflect a volume hierarchy if children
  //             has also been reflected, as both will have same name
  if( reffact->IsReflected( lv )
   && reffact->IsReflected( pv->GetMotherLogical() ) )  { return; }

  G4bool bVolExists = !CheckIfLogVolExists( lv->GetName(), lv );
  if( bVolExists )
  {
    DumpLogVol( lv );
  }

  //---- Construct this PV
  if( pv->GetMotherLogical() != 0 )   // WORLD volume
  {
    if( pv->IsReplicated() || pv->IsParameterised() )
    {
      G4String ErrMessage = "Only G4PVPlacement is supported yet, sorry... "
                          + pv->GetName();
      G4Exception("G4tgbGeometryDumper::DumpPlacements()", "NotImplemented",
                  FatalException, ErrMessage);
    }
    
    G4String pvName = pv->GetName();

    G4RotationMatrix* rotMat = pv->GetRotation();
    if( !rotMat ) rotMat = new G4RotationMatrix();
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    {
      G4cout << " G4tgbGeometryDumper::DumpPhysVol() - PV RotationMatrix: "
             << rotMat << G4endl;
    }
#endif

    //---- Check if it is reflected
    if( reffact->IsReflected( lv ) )
    {
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 2 )
      {
        G4cout << " G4tgbGeometryDumper::DumpPhysVol() - Reflected volume: "
               << pv->GetName() << G4endl;
      }
#endif
      CLHEP::Hep3Vector colx = rotMat->colX();
      CLHEP::Hep3Vector coly = rotMat->colY();
      CLHEP::Hep3Vector colz = rotMat->colZ();
      // apply a Z reflection (reflection matrix is decomposed in new
      // reflection-free rotation + z-reflection)
      colz *= -1.;
      CLHEP::HepRep3x3 rottemp(colx.x(),coly.x(),colz.x(),
                               colx.y(),coly.y(),colz.y(),
                               colx.z(),coly.z(),colz.z());
        // matrix representation (inverted)
      *rotMat = G4RotationMatrix(rottemp);
      *rotMat = (*rotMat).inverse();
      pvName += "_refl";
    }
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    {
      G4cout << " G4tgbGeometryDumper::DumpPhysVol() -"
             << " Calling DumpRotationMatrix :" << rotMat << G4endl;
    }
#endif
    G4String rotName = DumpRotationMatrix( rotMat );
    G4ThreeVector pos = pv->GetTranslation();
  
    G4String fullname = pvName
      +"#"+itoa(pv->GetCopyNo())
      +"/"+pv->GetMotherLogical()->GetName();

    if( !CheckIfPhysVolExists(fullname, pv ))
    {
      (*theFile)
           << ":PLACE "
           << SubstituteRefl(AddQuotes(pv->GetLogicalVolume()->GetName()))
           << " " << pv->GetCopyNo() << " "
           << SubstituteRefl(AddQuotes(pv->GetMotherLogical()->GetName()))
           << " " << AddQuotes(rotName) << " " 
           << pos.x() << " " << pos.y() << " " << pos.z() << G4endl;

      thePhysVols[fullname] = pv;
    }
  }

  if( bVolExists )
  {
    //---- Construct PV's who has this LV as mother
    std::vector<G4VPhysicalVolume*> pvChildren = GetPVChildren( lv );
    std::vector<G4VPhysicalVolume*>::const_iterator ite;
    for( ite = pvChildren.begin(); ite != pvChildren.end(); ite++ )
    {
      DumpPhysVol( *ite );
    }
  }
}

//------------------------------------------------------------------------
void G4tgbGeometryDumper::DumpLogVol( G4LogicalVolume* lv )
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgbGeometryDumper::DumpLogVol() - Log vol: "
           << lv->GetName() << lv->GetSolid()->GetName() << G4endl;
  }
#endif
  G4VSolid* solid;
  //--- take out the '_refl' in the name
  G4String lvName = lv->GetName();

  solid = lv->GetSolid();
  //---- Dump solid 
  DumpSolid( solid );

  //---- Dump material
  G4Material* mate = lv->GetMaterial();
  DumpMaterial( mate );

  //---- Dump logical volume (solid + material)
  (*theFile) << ":VOLU " << SubstituteRefl(AddQuotes(lvName)) << " "
             << SupressRefl(AddQuotes(solid->GetName()))
             << " " << AddQuotes(mate->GetName()) << G4endl;

  theLogVols[lvName] = lv;
}


//------------------------------------------------------------------------
void G4tgbGeometryDumper::DumpMaterial( G4Material* mat )
{
  if( CheckIfMaterialExists( mat->GetName(), mat ) )  { return; }

  size_t numElements           = mat->GetNumberOfElements();
  G4double density             = mat->GetDensity()/g*cm3;
  const G4ElementVector* elems = mat->GetElementVector();
  const G4double* fractions    = mat->GetFractionVector();
  
  // start tag
  if (numElements == 1)
  {
    (*theFile) << ":MATE " << AddQuotes(mat->GetName()) << " "
               << mat->GetZ() << " " << mat->GetA()/(g/mole) << " "
               << density << G4endl;
  }
  else
  {
    (*theFile) << ":MIXT "<< AddQuotes(mat->GetName()) << " "
               << density << " " << numElements << G4endl;
    // close start element tag and get ready to do composit "parts"
    for (size_t ii = 0; ii < numElements; ii++)
    {
      (*theFile) << "   " << AddQuotes((*elems)[ii]->GetName()) << " "
                 << fractions[ii] << G4endl;
    }

    for (size_t ii = 0; ii < numElements; ii++)
    {
      DumpElement( (*elems)[ii] );
    }
  }

  theMaterials[mat->GetName()] = mat;

}


//------------------------------------------------------------------------
void G4tgbGeometryDumper::DumpElement( G4Element* ele)
{
  if( CheckIfElementExists( ele->GetName(),ele ) )  { return; }

  if( ele->GetNumberOfIsotopes() != 0 )
  {
    G4Exception("G4tgbGeometryDumper::DumpElement()",
                "NotImplemented", FatalException,
                "Elements from isotopes not supported yet, sorry...");
  }

  //----- Material mixtures store the components as elements
  //      (even if the input are materials), but without symbol
  G4String symbol = ele->GetSymbol();
  if( symbol == "" || symbol == " " )
  {
    symbol = ele->GetName();
  }
  (*theFile) << ":ELEM " << AddQuotes(ele->GetName()) << " "
             << AddQuotes(symbol) << " " << ele->GetZ() << " "
             << ele->GetA()/(g/mole) << " " << G4endl;
  
  theElements[ele->GetName()] = ele;
}


//------------------------------------------------------------------------
void G4tgbGeometryDumper::DumpSolid( G4VSolid* solid )
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgbGeometryDumper::DumpSolid() - Solid: "
           << solid->GetName() << G4endl;
  }
#endif
  if( CheckIfSolidExists( solid->GetName(), solid) )  { return; }

  G4String solidType = solid->GetEntityType();
  solidType = GetTGSolidType( solidType );

  if (solidType == "UNIONSOLID")
  {
    DumpBooleanVolume( "UNION", solid );

  } else if (solidType == "SUBTRACTIONSOLID")  {
    DumpBooleanVolume( "SUBS", solid );

  } else if (solidType == "INTERSECTIONSOLID") {
    DumpBooleanVolume( "INTERS", solid );

  } else if (solidType == "REFLECTEDSOLID") {
    G4ReflectedSolid* solidrefl = dynamic_cast<G4ReflectedSolid*>(solid);
    G4VSolid* solidori = solidrefl->GetConstituentMovedSolid();
    DumpSolid( solidori );
  }
  else
  {
    (*theFile) << ":SOLID " << AddQuotes(solid->GetName()) << " ";
    (*theFile) << AddQuotes(solidType) << " ";
    DumpSolidParams( solid );

    theSolids[solid->GetName()] = solid;
  }
}


//------------------------------------------------------------------------
void G4tgbGeometryDumper::DumpBooleanVolume( const G4String& solidType,
                                                   G4VSolid* so )
{
  G4BooleanSolid * bso = dynamic_cast < G4BooleanSolid * > (so);
  G4VSolid* solid0 = bso->GetConstituentSolid( 0 );
  G4VSolid* solid1 = bso->GetConstituentSolid( 1 );
  G4DisplacedSolid* solid1Disp = 0;
  G4bool displaced = dynamic_cast<G4DisplacedSolid*>(solid1);
  if( displaced )
  {
    solid1Disp = dynamic_cast<G4DisplacedSolid*>(solid1);
    solid1 = solid1Disp->GetConstituentMovedSolid();
  }
  DumpSolid( solid0 );
  DumpSolid( solid1 );

  G4String rotName;
  G4ThreeVector pos;
  if( displaced )
  {
    pos = solid1Disp->GetObjectTranslation(); // translation is of mother frame
    rotName = DumpRotationMatrix( new G4RotationMatrix( (solid1Disp->
                                  GetTransform().NetRotation()).inverse() ) );
  }
  else  // no displacement
  {
    rotName = DumpRotationMatrix( new G4RotationMatrix );
    pos = G4ThreeVector();
  }

  if( CheckIfSolidExists( bso->GetName(), bso) )  { return; }

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
    //--- Dump RZ corners, as original parameters will not be present
    //    if it was build from RZ corners
    G4Polycone * pc = dynamic_cast < G4Polycone * > (so);
    
    G4double angphi = pc->GetStartPhi()/deg;
    if( angphi > 180*deg ) angphi -= 360*deg;
    G4int ncor = pc->GetNumRZCorner();
    (*theFile) << angphi << " "
               << pc->GetOriginalParameters()->Opening_angle/deg << " "
               << ncor << G4endl;
  
    for( G4int ii = 0; ii < ncor; ii++ )
    {
      (*theFile) << pc->GetCorner(ii).r << " " << pc->GetCorner(ii).z << G4endl;
    }

    // transform rz points into z,rmin,rmax
    //
    std::map<G4double,std::vector<G4double> > rzs;
    std::map<G4double,std::vector<G4double> >::iterator mite,mite2,mite3; 
    for( G4int ii = 0; ii < ncor; ii++ )
    {
      mite = rzs.find( pc->GetCorner(ii).z );
      if( mite == rzs.end() )
      {
        std::vector<G4double> val; 
        val.push_back( pc->GetCorner(ii).r );
        rzs[pc->GetCorner(ii).z] = val;
      }
      else
      {
        std::vector<G4double> val = (*mite).second;
        val.push_back( pc->GetCorner(ii).r );
        rzs[pc->GetCorner(ii).z] = val;
      }
    }

    // transform 3 R's into 2 (get min and max)
    //
    for( mite = rzs.begin(); mite != rzs.end(); mite++ )
    {
      std::vector<G4double> val = (*mite).second; 
      G4cout << " PG z " << (*mite).first << " N R's " << val.size() << G4endl;
      if( val.size() > 2 )
      {
        //take min and max values only
        G4double vmin = DBL_MAX;
        G4double vmax = -DBL_MAX; 
        std::vector<G4double> valnew;
        for( size_t ii = 0; ii < val.size(); ii++ )
        {
          if( val[ii] < vmin )  { vmin = val[ii]; }
          if( val[ii] > vmax )  { vmax = val[ii]; }
          G4cout << " PG min R " << vmin << " max R " << vmax << G4endl;
        }
        valnew.push_back( vmin );
        valnew.push_back( vmax );
        (*mite).second = valnew;
      }
    }

    
    mite2 = rzs.end(); mite2--;
    // transform 1 R into 2. First do end points (easier logic later): min=max
    //
    for( mite = rzs.begin(); mite != rzs.end(); mite++)
    {
      if( mite == rzs.begin() || mite == mite2 )
      {
        std::vector<G4double> val = (*mite).second; 
        G4cout << " PG 2nd: z " << (*mite).first
               << " N R's " << val.size() << G4endl;
        if( val.size() == 1 )
        {
          std::vector<G4double> valnew;
          valnew.push_back( val[0] );
          valnew.push_back( val[0] );
          (*mite).second = valnew;
        }
      }
    }

    // transform 1 R into 2. No end points: interpolate between two neighbours
    //
    for( mite = rzs.begin(); mite != mite2; mite++ )
    {
      if( mite == rzs.begin() )  { continue; }
      std::vector<G4double> val = (*mite).second; 
      G4cout << " PG 2nd: z " << (*mite).first
             << " N R's " << val.size() << G4endl;
      if( val.size() == 1 )
      {
        // Check that neighbours do not also 1
        mite3 = mite; mite3--;
        std::vector<G4double> valleft = (*(mite3)).second;
        G4double zleft = (*(mite3)).first;
        mite3++; mite3++;
        std::vector<G4double> valright = (*(mite3)).second;
        G4double zright = (*(mite3)).first;
        if( valleft.size() == 1 || valright.size() == 1 )
        {
          G4String ErrMessage = "RZ polygone with only two neighbour Z's"
                              + G4String(" having only 1 R !");
          G4Exception("G4tgbGeometryDumper::DumpSolidParams()",
                      "InvalidSetup", FatalException, ErrMessage);
        }
        // Interpolate min and max and then determine if val[0] is Rmin or Rmax
        G4double rmin = valleft[0]+(valright[0]- valleft[0])/(zright-zleft)
                                  *((*mite).first-zleft);
        G4double rmax = valleft[1]+(valright[1]- valleft[1])/(zright-zleft)
                                  *((*mite).first-zleft);
        if( val[0] < rmin ) {
          val.push_back(rmax);
        } else if ( val[0] > rmax ) {
          val.push_back(rmin);
        }
        else
        {
          G4String ErrMessage = "RZ polygone with only 1 R cannot be "
                 + G4String("determined\n if it is Rmin or Rmax !");
          G4Exception("G4tgbGeometryDumper::DumpSolidParams()",
                      "InvalidSetup", FatalException, ErrMessage);
        }
        (*mite).second = val;
      }
    }

    // Print result

    for( mite = rzs.begin(); mite != rzs.end(); mite++ )
    {
      std::vector<G4double> val = (*mite).second; 
      G4cout << " POLYCONE z " << (*mite).first
             << " rmin " << val[0] << " rmax  " << val[1] << G4endl;
    }

  } else if (solidType == "POLYHEDRA") {
    //--- Dump RZ corners, as original parameters will not be present
    //    if it was build from RZ corners
    G4Polyhedra * ph = (dynamic_cast < G4Polyhedra * > (so));
    
    G4double angphi = ph->GetStartPhi()/deg;
    if( angphi > 180*deg ) angphi -= 360*deg;

    G4int ncor = ph->GetNumRZCorner();

    (*theFile) << angphi << " "
               << ph->GetOriginalParameters()->Opening_angle/deg << " " 
               << ph->GetNumSide() << " " << ncor << G4endl;

    for( G4int ii = 0; ii < ncor; ii++ )
    {
      (*theFile) << ph->GetCorner(ii).r << " "
                 << ph->GetCorner(ii).z << G4endl;
    }

  } else if (solidType == "TRAP") {
    G4Trap * trp = dynamic_cast < G4Trap * > (so);
    G4ThreeVector symAxis(trp->GetSymAxis());
    G4double theta, phi;
    theta = symAxis.theta()/deg;
    phi = symAxis.phi()/deg;
    (*theFile) << trp->GetZHalfLength()  << " "
               << theta << " " 
               << phi<< " "
               << trp->GetYHalfLength1() << " "
               << trp->GetXHalfLength1() << " "
               << trp->GetXHalfLength2() << " "    
               << std::atan(trp->GetTanAlpha1())/deg << " " 
               << trp->GetYHalfLength2()    << " "
               << trp->GetXHalfLength3()    << " "
               << trp->GetXHalfLength4()    << " "    
               << std::atan(trp->GetTanAlpha2())/deg << " "
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

  } else if (solidType == "SPHERE") {
    G4Sphere * sphe = dynamic_cast < G4Sphere * > (so);
    (*theFile) << sphe->GetInsideRadius() << G4endl;
    (*theFile) << sphe->GetOuterRadius() << G4endl;
    (*theFile) << sphe->GetStartPhiAngle() << G4endl;
    (*theFile) << sphe->GetDeltaPhiAngle() << G4endl;
    (*theFile) << sphe->GetStartThetaAngle() << G4endl;
    (*theFile) << sphe->GetDeltaThetaAngle() << G4endl;

  } else if (solidType == "ELLIPSOID" ){
    G4Ellipsoid* dso = dynamic_cast < G4Ellipsoid * > (so);
    (*theFile) << dso->GetSemiAxisMax(0)  << " "
               << dso->GetSemiAxisMax(1)  << " "
               << dso->GetSemiAxisMax(2)  << " "
               << dso->GetZBottomCut()   << " "
               << dso->GetZTopCut() << G4endl;
  }
  else
  {
    G4String ErrMessage = "Solid type not supported, sorry... " + solidType;
    G4Exception("G4tgbGeometryDumpe::DumpSolidParams()",
                "NotImplemented", FatalException, ErrMessage);
  }
}


//------------------------------------------------------------------------
std::vector<G4double> G4tgbGeometryDumper::GetSolidParams( const G4VSolid * so) 
{
  std::vector<G4double> params;

  G4String solidType = so->GetEntityType();
  solidType = GetTGSolidType( solidType );

  if (solidType == "BOX")  {
    const G4Box * sb = dynamic_cast < const G4Box*>(so);
    params.push_back( sb->GetXHalfLength() ); 
    params.push_back( sb->GetYHalfLength() ); 
    params.push_back( sb->GetZHalfLength() ); 

  } else if (solidType == "TUBS") {
    const G4Tubs * tu = dynamic_cast < const G4Tubs * > (so);
    params.push_back( tu->GetInnerRadius()   );
    params.push_back( tu->GetOuterRadius()   );
    params.push_back( tu->GetZHalfLength()   );
    params.push_back( tu->GetStartPhiAngle()/deg );
    params.push_back( tu->GetDeltaPhiAngle()/deg );
    
  } else if (solidType == "CONS") {
    const G4Cons * cn = dynamic_cast < const G4Cons * > (so);
    params.push_back( cn->GetInnerRadiusMinusZ() ); 
    params.push_back( cn->GetOuterRadiusMinusZ() );
    params.push_back( cn->GetInnerRadiusPlusZ()  );    
    params.push_back( cn->GetOuterRadiusPlusZ()  );
    params.push_back( cn->GetZHalfLength() );
    params.push_back( cn->GetStartPhiAngle()/deg  );
    params.push_back( cn->GetDeltaPhiAngle()/deg  );
  } else if (solidType == "POLYCONE") {
    //--- Dump RZ corners, as original parameters will not be present
    //    if it was build from RZ corners
    const G4Polycone * pc = dynamic_cast < const G4Polycone * > (so);
    
    G4double angphi = pc->GetStartPhi()/deg;
    if( angphi > 180*deg )  { angphi -= 360*deg; }
    G4int ncor = pc->GetNumRZCorner();
    params.push_back( angphi );
    params.push_back( pc->GetOriginalParameters()->Opening_angle/deg ); 
    params.push_back( ncor );
    
    for( G4int ii = 0; ii < ncor; ii++ )
    {
      params.push_back( pc->GetCorner(ii).r ); 
      params.push_back( pc->GetCorner(ii).z );
    }
    
    // transform rz points into z,rmin,rmax
    //
    std::map<G4double,std::vector<G4double> > rzs;
    std::map<G4double,std::vector<G4double> >::iterator mite,mite2,mite3; 
    for( G4int ii = 0; ii < ncor; ii++ )
    {
      mite = rzs.find( pc->GetCorner(ii).z );
      if( mite == rzs.end() )
      {
        std::vector<G4double> val; 
        val.push_back( pc->GetCorner(ii).r );
        rzs[pc->GetCorner(ii).z] = val;
      }
      else
      {
        std::vector<G4double> val = (*mite).second;
        val.push_back( pc->GetCorner(ii).r );
        rzs[pc->GetCorner(ii).z] = val;
      }
    }
    
    // transform 3 R's into 2 (get min and max)
    //
    for( mite = rzs.begin(); mite != rzs.end(); mite++ )
    {
      std::vector<G4double> val = (*mite).second; 
      G4cout << " PG z " << (*mite).first << " N R's " << val.size() << G4endl;
      if( val.size() > 2 ) {
        // take min and max values only
        G4double vmin = DBL_MAX;
        G4double vmax = -DBL_MAX; 
        std::vector<G4double> valnew;
        for( size_t ii = 0; ii < val.size(); ii++ )
        {
          if( val[ii] < vmin )  { vmin = val[ii]; }
          if( val[ii] > vmax )  { vmax = val[ii]; }
          G4cout << " PG min R " << vmin << " max R " << vmax << G4endl;
        }
        valnew.push_back( vmin );
        valnew.push_back( vmax );
        (*mite).second = valnew;
      }
    }
    mite2 = rzs.end(); mite2--;

    // transform 1 R into 2. First do end points (easier logic later): min=max
    //
    for( mite = rzs.begin(); mite != rzs.end(); mite++)
    {
      if( mite == rzs.begin() || mite == mite2 )
      {
        std::vector<G4double> val = (*mite).second; 
        G4cout << " PG 2nd: z " << (*mite).first
               << " N R's " << val.size() << G4endl;
        if( val.size() == 1 )
        {
          std::vector<G4double> valnew;
          valnew.push_back( val[0] );
          valnew.push_back( val[0] );
          (*mite).second = valnew;
        }
      }
    }
    
    // transform 1 R into 2. No end points: interpolate between two neighbours
    //
    for( mite = rzs.begin(); mite != mite2; mite++ )
    {
      if( mite == rzs.begin() )  { continue; }
      std::vector<G4double> val = (*mite).second; 
      G4cout << " PG 2nd: z " << (*mite).first
             << " N R's " << val.size() << G4endl;
      if( val.size() == 1 )
      {
        // Check that neighbours do not also 1
        mite3 = mite; mite3--;
        std::vector<G4double> valleft = (*(mite3)).second;
        G4double zleft = (*(mite3)).first;
        mite3++; mite3++;
        std::vector<G4double> valright = (*(mite3)).second;
        G4double zright = (*(mite3)).first;
        if( valleft.size() == 1 || valright.size() == 1 )
        {
          G4String ErrMessage = "RZ polygone with only two neighbours Z's \n"
                              + G4String("that have only 1 R !");
          G4Exception("G4tgbGeometryDumper::GetSolidParams()",
                      "InvalidSetup", FatalException, ErrMessage);
        }
        // Interpolate min and max and then determine if val[0] ir Rmin or Rmax
        //
        G4double rmin = valleft[0]+(valright[0]- valleft[0])/(zright-zleft)
                                  *((*mite).first-zleft);
        G4double rmax = valleft[1]+(valright[1]- valleft[1])/(zright-zleft)
                                  *((*mite).first-zleft);
        if( val[0] < rmin ) {
          val.push_back(rmax);
        } else if ( val[0] > rmax ) {
          val.push_back(rmin);
        }
        else
        {
          G4String ErrMessage = "RZ polygone with only 1 R cannot be "
                 + G4String("determined \n if it is Rmin or Rmax !");
          G4Exception("G4tgbGeometryDumper::GetSolidParams()",
                      "InvalidSetup", FatalException, ErrMessage);
        }
        (*mite).second = val;
      }
    }

    // Print result

    for( mite = rzs.begin(); mite != rzs.end(); mite++ )
    {
      std::vector<G4double> val = (*mite).second; 
      G4cout << " POLYCONE z " << (*mite).first
             << " rmin " << val[0] << " rmax  " << val[1] << G4endl;
    }
    //    G4cout << " POLYCONE " << pc->GetName() << *pc );
    
  } else if (solidType == "POLYHEDRA") {
    //--- Dump RZ corners, as original parameters will not be present
    //    if it was build from RZ corners
    const G4Polyhedra * ph = (dynamic_cast < const G4Polyhedra * > (so));
    
    G4double angphi = ph->GetStartPhi()/deg;
    if( angphi > 180*deg ) angphi -= 360*deg;

    G4int ncor = ph->GetNumRZCorner();
    
    params.push_back( angphi );
    params.push_back( ph->GetOriginalParameters()->Opening_angle/deg ); 
    params.push_back( ph->GetNumSide() ); 
    params.push_back( ncor );

    for( G4int ii = 0; ii < ncor; ii++ )
    {
       params.push_back( ph->GetCorner(ii).r ); 
       params.push_back( ph->GetCorner(ii).z );
    }

  } else if (solidType == "TRAP") {
    const G4Trap * trp = dynamic_cast < const G4Trap * > (so);
    G4ThreeVector symAxis(trp->GetSymAxis());
    G4double theta, phi;
    theta = symAxis.theta()/deg;
    phi = symAxis.phi()/deg;
    params.push_back( trp->GetZHalfLength() );
    params.push_back( theta ); 
    params.push_back( phi);
    params.push_back( trp->GetYHalfLength1() );
    params.push_back( trp->GetXHalfLength1() );
    params.push_back( trp->GetXHalfLength2() );    
    params.push_back( std::atan(trp->GetTanAlpha1())/deg ); 
    params.push_back( trp->GetYHalfLength2()    );
    params.push_back( trp->GetXHalfLength3()    );
    params.push_back( trp->GetXHalfLength4()    );    
    params.push_back( std::atan(trp->GetTanAlpha2())/deg );
  } else if (solidType == "TRD") {
    const G4Trd * tr = dynamic_cast < const G4Trd * > (so);
    params.push_back( tr->GetZHalfLength());
    params.push_back( tr->GetYHalfLength1() );
    params.push_back( tr->GetYHalfLength2() ); 
    params.push_back( tr->GetXHalfLength1() );
    params.push_back( tr->GetXHalfLength2() );
    
  } else if (solidType == "ORB") {
    const G4Orb * orb = dynamic_cast < const G4Orb * > (so);
    params.push_back( orb->GetRadius()  );
     
  } else if (solidType == "SPHERE") {
    const G4Sphere * sphe = dynamic_cast < const G4Sphere * > (so);
    params.push_back( sphe->GetInsideRadius() );
    params.push_back( sphe->GetOuterRadius() );
    params.push_back( sphe->GetStartPhiAngle() );
    params.push_back( sphe->GetDeltaPhiAngle() );
    params.push_back( sphe->GetStartThetaAngle() );
    params.push_back( sphe->GetDeltaThetaAngle() );

  } else if (solidType == "ELLIPSOID" ){
    const G4Ellipsoid* dso = dynamic_cast < const G4Ellipsoid * > (so);
     params.push_back( dso->GetSemiAxisMax(0)  );
     params.push_back( dso->GetSemiAxisMax(1)  );
     params.push_back( dso->GetSemiAxisMax(2)  );
     params.push_back( dso->GetZBottomCut()   );
     params.push_back( dso->GetZTopCut() );
  }
  else
  {
    G4String ErrMessage = "Solid type not supported, sorry... " + solidType;
    G4Exception("G4tgbGeometryDumpe::DumpSolidParams()",
                "NotImplemented", FatalException, ErrMessage);
  }
   
  return params;
}   


//------------------------------------------------------------------------
G4String G4tgbGeometryDumper::DumpRotationMatrix( G4RotationMatrix* rotm )
{
  G4double de = MatDeterminant(rotm);
 
  G4String rotName = LookForExistingRotation( rotm );
  if( rotName != "" )  { return rotName; }

  if (!rotm)  { rotm = new G4RotationMatrix(); } 

  G4ThreeVector v(1.,1.,1.);
  if (de < -0.9 )  // a reflection ....
  {
    (*theFile) << ":ROTM ";
    rotName = "RRM";
    rotName += itoa(theRotationNumber++);
 
    (*theFile) << AddQuotes(rotName) << std::setprecision(9) << " " 
               << approxTo0(rotm->xx())  << " "
               << approxTo0(rotm->yx())  << " "
               << approxTo0(rotm->zx())  << " "
               << approxTo0(rotm->xy())  << " "
               << approxTo0(rotm->yy())  << " "
               << approxTo0(rotm->zy())  << " "
               << approxTo0(rotm->xz())  << " "
               << approxTo0(rotm->yz())  << " "
               << approxTo0(rotm->zz())  << G4endl;
  }
  else if(de > 0.9 )  // a rotation ....
  {
    (*theFile) << ":ROTM ";
    rotName = "RM";
    rotName += itoa(theRotationNumber++);
    
    (*theFile) << AddQuotes(rotName) << " " 
               << approxTo0(rotm->thetaX()/deg)  << " "
               << approxTo0(rotm->phiX()/deg)    << " "
               << approxTo0(rotm->thetaY()/deg)  << " "
               << approxTo0(rotm->phiY()/deg)    << " "
               << approxTo0(rotm->thetaZ()/deg)  << " "
               << approxTo0(rotm->phiZ()/deg)    << G4endl;
  }
  
  theRotMats[rotName] = rotm;

  return rotName;
}

//------------------------------------------------------------------------
std::vector<G4VPhysicalVolume*>
G4tgbGeometryDumper::GetPVChildren( G4LogicalVolume* lv )
{
  G4PhysicalVolumeStore* pvstore = G4PhysicalVolumeStore::GetInstance();
  G4PhysicalVolumeStore::const_iterator ite;
  std::vector<G4VPhysicalVolume*> children;
  for( ite = pvstore->begin(); ite != pvstore->end(); ite++ )
  {
    if( (*ite)->GetMotherLogical() == lv )
    {
      children.push_back( *ite );
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 2 )
      {
        G4cout << " G4tgbGeometryDumper::GetPVChildren() - adding children: "
               << (*ite)->GetName() << " of " << lv->GetName() <<  G4endl;
      }
#endif
    }
  }

  return children;
}

//------------------------------------------------------------------------
G4String G4tgbGeometryDumper::GetTGSolidType( G4String& solidType )
{
  G4String newsolidType = solidType.substr(2,solidType.length() );
  for( size_t ii = 0; ii < newsolidType.length(); ii++ )
  {
    newsolidType[ii] = toupper(newsolidType[ii] );
  }
  return newsolidType;
}

//------------------------------------------------------------------------
G4double G4tgbGeometryDumper::MatDeterminant(G4RotationMatrix * ro) 
{
   CLHEP::HepRep3x3 r = ro->rep3x3();
   return       r.xx_*(r.yy_*r.zz_ - r.zy_*r.yz_)
              - r.yx_*(r.xy_*r.zz_ - r.zy_*r.xz_)
              + r.zx_*(r.xy_*r.yz_ - r.yy_*r.xz_);
}

//------------------------------------------------------------------------
G4String G4tgbGeometryDumper::itoa(G4int current)
{
  const char theDigits[11] = "0123456789";
  G4String result;
  G4int digit;
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
G4double G4tgbGeometryDumper::approxTo0( G4double val )
{
  G4double precision = G4GeometryTolerance::GetInstance()
                       ->GetSurfaceTolerance();

  if( std::fabs(val) < precision )  { val = 0; }
  return val;
}


//-----------------------------------------------------------------------
G4String G4tgbGeometryDumper::AddQuotes( const G4String& str )
{
  //--- look if there is a separating blank

  G4bool bBlank = FALSE;
  size_t siz = str.length();
  for( size_t ii = 0; ii < siz; ii++ )
  {
    if( str.substr(ii,1) == " " )
    {
      bBlank = TRUE;
      break;
    }
  }
  G4String str2 = str;
  if( bBlank )
  {
    str2 = G4String("\"") + str2 + G4String("\"");
  }
  return str2;
}


//------------------------------------------------------------------------
G4String G4tgbGeometryDumper::SupressRefl( G4String name )
{
  G4int irefl = name.rfind("_refl");
  if( irefl != -1 )
  {
    name = name.substr( 0, irefl );
  }
  return name;
}

//------------------------------------------------------------------------
G4String G4tgbGeometryDumper::SubstituteRefl( G4String name )
{
  G4int irefl = name.rfind("_refl");
  if( irefl != -1 )
  {
    name = name.substr( 0, irefl ) + "_REFL";
  }
  return name;
}


//------------------------------------------------------------------------
G4bool G4tgbGeometryDumper::CheckIfElementExists( const G4String& name,
                                                        G4Element* pt )
{
  if( theElements.find( name ) != theElements.end() )
  {
    if( pt != (*(theElements.find(name))).second )
    {
      G4cerr << "G4tgbGeometryDumper::CheckIfElementExists() - "
             << "Element found but not same as before: " << name << G4endl;
      if( !Same2G4Elements(pt, (*(theElements.find(name))).second))
      {
        G4String ErrMessage = "Element found but with different A or Z \n"
                            + G4String("as as before: ") + name;
        G4Exception("G4tgbGeometryDumper::CheckIfElementExists()",
                    "InvalidSetup", FatalException, ErrMessage);
      }
    }
    return 1;
  }
  else
  {
    return 0;
  }
}


//------------------------------------------------------------------------
G4bool G4tgbGeometryDumper::CheckIfMaterialExists( const G4String& name,
                                                         G4Material* pt )
{
  if( theMaterials.find( name ) != theMaterials.end() )
  {
    if( (*(theMaterials.find(name))).second != pt )
    {
      G4cerr << "G4tgbGeometryDumper::CheckIfMaterialExists() -"
             << " Material found but not same as before : " << name << G4endl;
      if(Same2G4Materials(pt, (*(theMaterials.find(name))).second ))
      {
        G4String ErrMessage = "Material found but not same as before: " + name;
        G4Exception("G4tgbGeometryDumper::CheckIfMaterialExists()",
                    "InvalidSetup", FatalException, ErrMessage);
      }
    }
    return 1;
  }
  else
  {
    return 0;
  } 
}


//------------------------------------------------------------------------
G4bool G4tgbGeometryDumper::CheckIfLogVolExists( const G4String& name,
                                                       G4LogicalVolume* pt )
{
  if( theLogVols.find( name ) != theLogVols.end() )
  {
    G4LogicalVolume* lvnew = (*(theLogVols.find(name))).second;
    if( lvnew != pt )
    {
      /*
      //---- Reflected volumes are repeated

      G4ReflectionFactory* reffact = G4ReflectionFactory::Instance();
      if( !reffact->IsReflected( pt ) && !reffact->IsReflected( lvnew ) )
      {
        G4String ErrMessage = "LogVol found but not same as before: " + name;
        G4Exception("G4tgbGeometryDumper::CheckIfLogVolExists()",
                    "InvalidSetup", FatalException, ErrMessage);
      }
      */
    }
    return 1;
  }
  else
  {
    return 0;
  }
}


//------------------------------------------------------------------------
G4bool G4tgbGeometryDumper::CheckIfSolidExists( const G4String& name,
                                                      G4VSolid* pt )
{
  if( theSolids.find( name ) != theSolids.end() )
  {
    if( (*(theSolids.find(name))).second != pt )
    {
      unsigned int nn = 2;
      G4String solidName = pt->GetName();
      G4String addName = "_ADDED_";
      for( ;; )
      {
        if(theSolids.find( solidName+ addName + itoa(nn) ) == theSolids.end() )
        {
          addName += itoa(nn);
          break;
        }
        nn++;
      }
      G4cerr << "G4tgbGeometryDumper::CheckIfSolidExists() -"
             << " Solid found but not same as before : " << name
             << " Changed to " << name + addName << G4endl;
      pt->SetName(solidName+addName);
      //---- Change name of log vol that contains this solid
      G4LogicalVolumeStore* lvstore = G4LogicalVolumeStore::GetInstance();
      G4LogicalVolumeStore::const_iterator ite;
      for( ite = lvstore->begin(); ite != lvstore->end(); ite++ )
      {
        if( (*ite)->GetSolid() == pt )
        {
          (*ite)->SetName( (*ite)->GetName() +addName );
        }
      }
      return 0;
    }
    return 1;
  }
  else
  {
    return 0;
  }
}


//-----------------------------------------------------------------------
G4bool G4tgbGeometryDumper::CheckIfPhysVolExists( const G4String& name,
                                                        G4VPhysicalVolume* pt )
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgbGeometryDumper::CheckIfPhysVolExists() - "
           << name << G4endl;
  }
#endif
  if( thePhysVols.find( name ) != thePhysVols.end() )
  {
    if( (*(thePhysVols.find(name))).second != pt )
    {
      // G4String ErrMessage = "Placement found but not same as before: "
      //                     + name;
      // G4Exception("G4tgbGeometryDumper::CheckIfPhysVolExists()",
      //             "InvalidSetup", FatalException, ErrMessage);
      G4cerr << " G4tgbGeometryDumper::CheckIfPhysVolExists () -"
             << " Placement found but not same as before : " << name << G4endl;
    }
    return 1;
  }
  else
  {
    return 0;
  }
}


//-----------------------------------------------------------------------
G4String
G4tgbGeometryDumper::LookForExistingRotation( const G4RotationMatrix* rotm )
{
  G4String rmName = "";

  std::map<G4String,G4RotationMatrix*>::const_iterator ite;
  for( ite = theRotMats.begin(); ite != theRotMats.end(); ite++ )
  {
    if( (*ite).second->isNear( *rotm ) )
    {
      rmName = (*ite).first;
      break;
    }
  }
  return rmName;
}


//------------------------------------------------------------------------
G4bool G4tgbGeometryDumper::Same2G4Elements( G4Element* ele1, G4Element* ele2 )
{
  if( ele1->GetZ() != ele2->GetZ() || ele1->GetA() != ele2->GetA() )
  {
    return 0;
  }
  else
  {
    return 1;
  }
}


//------------------------------------------------------------------------
G4bool G4tgbGeometryDumper::Same2G4Materials( G4Material* mat1,
                                              G4Material* mat2 )
{
  G4bool bSame = 1;
  if( mat1->GetDensity() != mat2->GetDensity() ) { bSame = 0;}
  if( mat1->GetNumberOfElements() != mat2->GetNumberOfElements() ) { bSame = 0;}
  G4int nele = mat1->GetNumberOfElements();
  for( G4int ii=0; ii<nele; ii++)
  {
    if( mat1->GetFractionVector()[ii] != mat2->GetFractionVector()[ii] )
    {
      bSame = 0;
    }
    if( mat1->GetElement(ii) != mat2->GetElement(ii) )
    {
      bSame = 0;
    }
  }

  return bSame;
}
