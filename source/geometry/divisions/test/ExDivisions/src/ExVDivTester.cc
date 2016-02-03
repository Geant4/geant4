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
// $Id$
//
// class ExVDivTester Implementation file
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************

#include "ExVDivTester.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4Box.hh"
#include "G4PVPlacement.hh"

#include "Randomize.hh"

#include "G4PVDivisionFactory.hh"
#include "G4ReflectionFactory.hh"
#include "G4PVReplica.hh"
#include "G4PVDivision.hh"

#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include <sstream>
#include <fstream>

G4ThreadLocal G4bool ExVDivTester::bDivCylindrical = 0;

//--------------------------------------------------------------------------
ExVDivTester::
ExVDivTester( PVType& pvtype, PlaceType& postype,
              std::vector<G4String>& extraPars )
  : thePVType( pvtype ), thePosType( postype ), theExtraPars( extraPars )
{
  verbose = 2;
  //--- set the number of replicas
  theWorldLengthXY = 1*m;
  theWorldLengthZ = 8*theWorldLengthXY;
  theWorldGap = 10*cm;
  //--- default material
  theMate = new G4Material("Pb", 82., 207.19*g/mole, 11.35*g/cm3);

  theDivType = DivNDIV;
  theOffsetFactor = 0.;
  theWidthFactor = 1.;
  theNDiv = 3;

  theStartPhi = 0.*deg;
  theDeltaPhi = 360.*deg;

  // Create the division factory singleton
  //
  G4PVDivisionFactory::GetInstance();

  //---- extract PVType from extraPars (NREPLICAS, WIDTH, NREPLICAS&WIDTH)
  if( theExtraPars.size() != 0 )
  {
    readExtraPars();
  }
}


//--------------------------------------------------------------------------
void ExVDivTester::readExtraPars()
{
  G4int ii, siz = theExtraPars.size();
  G4bool ifound;
  for( ii = 0; ii < siz; ii++ ){
    ifound = 0;
    if( theExtraPars[ii].substr(0,6) == "offset" ) {
      theOffsetFactor = getValueFromExtraPar( theExtraPars[ii] );
      ifound = 1;
    }else if( theExtraPars[ii].substr(0,5) == "width" ) {
      theWidthFactor = getValueFromExtraPar( theExtraPars[ii] );
      ifound = 1;
    } else if( theExtraPars[ii].substr(0,7) == "divtype" ) {
      G4int ival = G4int(getValueFromExtraPar( theExtraPars[ii] ));
      if( ival == 0 ){
        theDivType = DivNDIVandWIDTH;
        ifound = 1;
      }else if( ival == 1 ){
        theDivType = DivNDIV;
        ifound = 1;
      }else if( ival == 2 ){
        theDivType = DivWIDTH;
        ifound = 1;
      }
    }else if( theExtraPars[ii].substr(0,8) == "startphi" ) {
      theStartPhi = getValueFromExtraPar( theExtraPars[ii] )*deg;
      ifound = 1;
    }else if( theExtraPars[ii].substr(0,8) == "deltaphi" ) {
      theDeltaPhi = getValueFromExtraPar( theExtraPars[ii] )*deg;
      ifound = 1;
    }

    if( !ifound ) {
      G4String msg = "Extra Parameter not in list! : "+theExtraPars[ii];
      G4Exception("ExVDivTester::readExtraPars()",
                  "IllegalConstruct", FatalException, msg.c_str() );
    }
  }
}


//---------------------------------------------------------------------------
G4double ExVDivTester::getValueFromExtraPar( G4String& epar )
{
  G4double val = 0.;

  G4int ieq = epar.find("=");
  if( ieq == -1 ) {
    G4String msg ="Extra Parameter does not have an '='! :"+epar;
    G4Exception("ExVDivTester::getValueFromExtraPar()",
                "IllegalConstruct", FatalException, msg.c_str() );

  } else {
    G4String eparValue = epar.substr(ieq+1,epar.size() );
    val = atof( eparValue.c_str() );
  }

  return val;
}


//--------------------------------------------------------------------------
ExVDivTester::~ExVDivTester()
{
  delete theMate;
}

//--------------------------------------------------------------------------
void ExVDivTester::GenerateRandomPoints()
{
  std::ofstream fout("points.lis");
  G4double posX, posY, posZ;

  numberOfPoints = 1000;
  for( G4int ii = 0; ii < numberOfPoints; ii++ ) {
    posX = G4RandFlat::shoot(-theWorldLengthXY, theWorldLengthXY);
    posY = G4RandFlat::shoot(-theWorldLengthXY, theWorldLengthXY);
    posZ = G4RandFlat::shoot(-theWorldLengthZ, theWorldLengthZ);
    fout << posX << " " << posY << " " << posZ << G4endl;
  }
}

//--------------------------------------------------------------------------
G4VPhysicalVolume* ExVDivTester::BuildGeometry()
{
  // Build the world volume
  //
  G4Box *worldSolid =
    new G4Box ("Big Cube", theWorldLengthXY+theWorldGap,
                           theWorldLengthXY+theWorldGap,
                           theWorldLengthZ+theWorldGap);

  // Air
  G4double a, density;
  a = 14.01*g/mole;
  G4Element* elN = new G4Element("Nitrogen", "N", 7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element("Oxigen", "O", 8., a);
  density = 1.29*mg/cm3;
  G4Material* Air = new G4Material("Air", density, 2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  G4LogicalVolume *worldLog=new G4LogicalVolume(worldSolid,Air,
                                                "World",0,0,0);
  
  G4PVPlacement *worldPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
                                             "World",worldLog,
                                             0,false,0);
  G4VisAttributes* worldVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  worldLog->SetVisAttributes(worldVisAtt);

  if( verbose >= 1 )
    G4cout << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@ BuildParentSolids "<< G4endl;
  BuildParentSolids();
  if( verbose >= 1 )
    G4cout << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@ BuildParentVolumes "<< G4endl;
  BuildParentVolumes( worldLog );

  if( verbose >= 1 )
    G4cout << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@ BuildChildrenSolids "<< G4endl;
  BuildChildrenSolids();
  if( verbose >= 1 )
    G4cout << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@ SetOffsets "<< G4endl;
  SetOffsets();
  if( verbose >= 1 )
    G4cout << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@ BuildChildrenVolumes "<< G4endl;
  BuildChildrenVolumes(); //LV and divisions as PV


  return worldPhys;
}

//--------------------------------------------------------------------------
void ExVDivTester::BuildParentVolumes( G4LogicalVolume* worldLog )
{
  // build three volumes
  G4LogicalVolume* parent1Log =
    new G4LogicalVolume(theParentSolids[0],theMate,"parentLog-1",0,0,0);
  G4LogicalVolume* parent2Log =
    new G4LogicalVolume(theParentSolids[1],theMate,"parentLog-2",0,0,0);
  G4LogicalVolume* parent3Log =
    new G4LogicalVolume(theParentSolids[2],theMate,"parentLog-3",0,0,0);

  G4VisAttributes* parentVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  parent1Log->SetVisAttributes(parentVisAtt);
  parent2Log->SetVisAttributes(parentVisAtt);
  parent3Log->SetVisAttributes(parentVisAtt);
  if ( thePosType == pvReflected )
  {
    G4LogicalVolume* parent1LogRefl
      = G4ReflectionFactory::Instance()->GetReflectedLV(parent1Log);
    G4LogicalVolume* parent2LogRefl
      = G4ReflectionFactory::Instance()->GetReflectedLV(parent2Log);
    G4LogicalVolume* parent3LogRefl
      = G4ReflectionFactory::Instance()->GetReflectedLV(parent2Log);
    if(parent1LogRefl) parent1LogRefl->SetVisAttributes(parentVisAtt);
    if(parent2LogRefl) parent2LogRefl->SetVisAttributes(parentVisAtt);
    if(parent3LogRefl) parent3LogRefl->SetVisAttributes(parentVisAtt);
  }

  theParentLogs.push_back(parent1Log);
  theParentLogs.push_back(parent2Log);
  theParentLogs.push_back(parent3Log);
  
  //parent physical volumes positions do not depend on SolidType, PVType
  G4int ii;
  G4int nParents = theAxis.size();
  for( ii = 0; ii < nParents; ii++ )
  {
    G4String parentstr = "parent-";
    std::ostringstream os;
    os << ii;
    G4String buf = os.str();
    parentstr += buf;
    if ( thePosType == pvReflected )
    {
      // With reflection
      G4Translate3D translate(0.,0.,(ii-1)*7*theWorldLengthXY);
      G4ReflectZ3D  reflect;
      G4Transform3D transform = translate * reflect; 

      G4PhysicalVolumesPair pvPair
        = G4ReflectionFactory::Instance()
          ->Place(transform, parentstr, theParentLogs[ii], worldLog, FALSE, 0);  
      theParentPhyss.push_back(pvPair.first);
    }
    else
    {
      // Without reflection
      theParentPhyss.push_back( new G4PVPlacement( 0,
                              G4ThreeVector(0.,0.,(ii-1)*7*theWorldLengthXY),
                              theParentLogs[ii] , parentstr,
                              worldLog, FALSE, 0 ) );
    }
  }
}

//--------------------------------------------------------------------------
void ExVDivTester::SetOffsets()
{
  G4int ii;
  G4int nParents = theAxis.size();
  for( ii = 0; ii < nParents; ii++ )
  {
    theOffsets.push_back( theWidths[ii]*theOffsetFactor );
    if( verbose >= 2 )
      G4cout << ii << " offsets " << theOffsets[ii] << " "
             << theWidths[ii] << G4endl;
  }
}


//--------------------------------------------------------------------------
void ExVDivTester::BuildChildrenVolumes()
{
  G4int ii;
  G4VisAttributes* childVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  G4VisAttributes* childVisAtt2 = new G4VisAttributes(G4Colour(0.0,1.0,0.0));

  G4int nParents = theAxis.size();
  for( ii = 0; ii < nParents; ii++ )
  {
    G4String childlogstr = "childLog-";
    std::ostringstream os;
    os << ii;
    G4String buf = os.str();
    childlogstr += buf;
    G4LogicalVolume* childLog = new G4LogicalVolume(theChildSolids[ii],
                                theMate, childlogstr,0,0,0);
    childLog->SetVisAttributes(childVisAtt);
    theChildLogs.push_back(childLog);
  }

  //----- Build divisions
  if( thePVType == pvDivision )
  {
    for( ii = 0; ii < nParents; ii++ )
    {
      if( verbose >= 1 )
        G4cout << " @@@@ Building Child volume " << ii << G4endl;
      G4String childstr = "child-";
      std::ostringstream os;
      os << ii;
      G4String buf = os.str();
      childstr += buf;

      if ( thePosType == pvReflected )
      {
        // With reflection
        //
        G4ReflectionFactory::Instance()->SetVerboseLevel(2);
        if( theDivType == DivNDIVandWIDTH )
        {
          G4ReflectionFactory::Instance()
            ->Divide(childstr, theChildLogs[ii], theParentLogs[ii],
                               theAxis[ii],
                               theNDiv, 
                               theWidths[ii]*theWidthFactor,
                               theOffsets[ii] );

          if( verbose >= 1 ) G4cout << "division NDIVandWIDTH " << theNDiv
                 << " " << theWidths[ii]*theWidthFactor
                 << " " << theOffsets[ii]*theOffsetFactor << G4endl;
        }
        else if( theDivType == DivNDIV )
        {
          G4ReflectionFactory::Instance()
            ->Divide(childstr, theChildLogs[ii], theParentLogs[ii],
                               theAxis[ii],
                               theNDiv, 
                               theOffsets[ii] );
          if( verbose >= 1 )
            G4cout << "division NDIV " << theNDiv << " "
                   << theOffsets[ii] << G4endl;
        }
        else if( theDivType == DivWIDTH )
        {
          G4ReflectionFactory::Instance()
            ->Divide(childstr, theChildLogs[ii], theParentLogs[ii],
                               theAxis[ii],
                               theWidths[ii]*theWidthFactor,
                               theOffsets[ii] );
          if( verbose >= 1 )
            G4cout << "division WIDTH " << theWidths[ii]*theWidthFactor
                   << " " << theOffsets[ii]*theOffsetFactor << G4endl;
        }

        // set vis attributes to reflected volumes
        G4LogicalVolume* childLogRefl
          = G4ReflectionFactory::Instance()->GetReflectedLV(theChildLogs[ii]);
        if (childLogRefl) childLogRefl->SetVisAttributes(childVisAtt2);
      }
      else
      {
        // Without reflection
        //
        if( theDivType == DivNDIVandWIDTH )
        {
          new G4PVDivision(childstr, theChildLogs[ii], theParentLogs[ii], 
                           theAxis[ii],
                           theNDiv, 
                           theWidths[ii]*theWidthFactor,
                           theOffsets[ii] );
          if( verbose >= 1 ) G4cout << "division NDIVandWIDTH " << theNDiv
                 << " " << theWidths[ii]*theWidthFactor
                 << " " << theOffsets[ii]*theOffsetFactor << G4endl;
        }
        else if( theDivType == DivNDIV )
        {
          new G4PVDivision(childstr, theChildLogs[ii], theParentLogs[ii], 
                           theAxis[ii],
                           theNDiv, 
                           theOffsets[ii] );
          if( verbose >= 1 )
            G4cout << "division NDIV " << theNDiv << " "
                   << theOffsets[ii] << G4endl;
        }
        else if( theDivType == DivWIDTH )
        {
          new G4PVDivision(childstr, theChildLogs[ii], theParentLogs[ii], 
                           theAxis[ii],
                           theWidths[ii]*theWidthFactor, 
                           theOffsets[ii] );
          if( verbose >= 1 )
            G4cout << "division WIDTH " << theWidths[ii]*theWidthFactor
                   << " " << theOffsets[ii]*theOffsetFactor << G4endl;
        }
      }
    }
  }
  else if( thePVType == pvReplica )
  {
    for( ii = 0; ii < nParents; ii++ )
    {
      G4String childstr = "child-";
      std::ostringstream os;
      os << ii;
      G4String buf = os.str();
      childstr += buf;
      if ( thePosType == pvReflected )
      {
        G4ReflectionFactory::Instance()
          ->Replicate(childstr, theChildLogs[ii], theParentLogs[ii],
                                theAxis[ii],
                                theNDiv, 
                                theWidths[ii]*theWidthFactor,
                                theOffsets[ii] );

        // set vis attributes to reflected volumes
        G4LogicalVolume* childLogRefl
          = G4ReflectionFactory::Instance()->GetReflectedLV(theChildLogs[ii]);
        if (childLogRefl) childLogRefl->SetVisAttributes(childVisAtt2);
      }
      else
      {
        new G4PVReplica(childstr, theChildLogs[ii], theParentLogs[ii], 
                                  theAxis[ii],
                                  theNDiv, 
                                  theWidths[ii]*theWidthFactor, 
                                  theOffsets[ii] );
      }
      if( verbose >= 1 )
        G4cout << "replica " <<  theNDiv << " "
               << theWidths[ii]*theWidthFactor << " "
               << theOffsets[ii]*theOffsetFactor << G4endl;
    }
  }
}

//------------------------------------------------------------------------
void ExVDivTester::PrintChildrenSolids( std::ostream& out )
{
  for( size_t ii = 0; ii < theChildSolids.size(); ii++)
  {
    out << theChildSolids[0] << "child solid built  0 "
        << *theChildSolids[0] << G4endl;
    out << " theChildSolid built " <<  ii << G4endl;
    out <<  *theChildSolids[ii] << G4endl;
  }
}

//------------------------------------------------------------------------
void ExVDivTester::PrintParentSolid( std::ostream& out )
{
  out << " PARENT SOLID " << G4endl;

  G4int nParents = theAxis.size();
  for( G4int ii = 0; ii < nParents; ii++ ){
    out << *(theParentSolids[ii]) << G4endl;
  }
}


//--------------------------------------------------------------------------
void ExVDivTester::GenerateScanPointsAsBox()
{
  std::ofstream fout("points.lis");
  G4int ii;

  G4int nPointsPerDiv = 2;
  G4double disppoint = 0*mm;
  numberOfPoints = theNDiv * nPointsPerDiv;
  // For division along X
  G4ThreeVector centre(0.,0.,-theWorldLengthZ+theWorldLengthXY);
  for( ii = 0; ii < numberOfPoints; ii++ )
  {
    // any Z, any Y
    G4ThreeVector pR( 0., disppoint, disppoint );
    G4double X = -theWorldLengthXY
               + (ii+0.5) * 2*theWorldLengthXY/numberOfPoints;
    pR.setX( X );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }

  // For division along Y
  centre = G4ThreeVector(0.,0.,0.);
  for( ii = 0; ii < numberOfPoints; ii++ )
  {
    // any X, any Z
    G4ThreeVector pR(disppoint, 0., disppoint );
    G4double Y = -theWorldLengthXY
               + (ii+0.5) * 2*theWorldLengthXY/numberOfPoints;
    pR.setY( Y );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }


  // For division along Z
  centre = G4ThreeVector(0.,0.,theWorldLengthZ-theWorldLengthXY);
  for( ii = 0; ii < numberOfPoints; ii++ )
  {
    // any X, any Y
    G4ThreeVector pR( disppoint, disppoint, 0. );
    G4double Z = -theWorldLengthXY
               + (ii+0.5) * 2*theWorldLengthXY/numberOfPoints;
    pR.setZ( Z );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }
}


//--------------------------------------------------------------------------
void ExVDivTester::GenerateScanPointsAsTubs()
{
  std::ofstream fout("points.lis");
  G4int ii;

  G4int nPointsPerDiv = 5;
  numberOfPoints = theNDiv*nPointsPerDiv;
  if( verbose >= 1 ) G4cout << " numberOfPoints " << numberOfPoints << G4endl;
  // For division along R
  G4ThreeVector centre(0.,0.,-theWorldLengthZ+theWorldLengthXY);
  for( ii = 0; ii < numberOfPoints; ii++ )
  {
    // any Z, phi close initial phi
    G4double phi = 0.*deg;
    //  if( theExtraPars.size() > 0 ) phi = (theExtraPars[0]+1.)*deg;
    G4ThreeVector pR( std::cos(phi), std::sin(phi), 0. );
    G4double rho = (ii+0.5) * theWorldLengthXY/numberOfPoints;
    pR.setRho( rho );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }
  
  // For division along phi
  centre = G4ThreeVector(theWorldLengthXY*0.7,0.,0.);
  for( ii = 0; ii < numberOfPoints; ii++ )
  {
    G4double phi = (ii+0.5) * 360*deg/numberOfPoints;
    // any Z
    //    G4ThreeVector pR( std::cos(phi), std::sin(phi), theWorldLengthXY/100. );
    //    pR += centre;
    G4ThreeVector pR( centre );
    pR.setPhi( phi );
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }
  
  // For division along Z
  centre = G4ThreeVector(theWorldLengthXY*0.7,
                         0., theWorldLengthZ-theWorldLengthXY);
  for( ii = 0; ii < numberOfPoints; ii++ )
  {
    //any rho (=1), phi close initial phi
    G4ThreeVector pR( centre);
    G4double Z = centre.z() - theWorldLengthXY
               + (ii+0.5) * 2*theWorldLengthXY / numberOfPoints;
    pR.setZ( Z );
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }
}
