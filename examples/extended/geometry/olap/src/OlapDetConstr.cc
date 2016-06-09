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
/// \file geometry/olap/src/OlapDetConstr.cc
/// \brief Implementation of the OlapDetConstr class
//
//
// $Id$
//
// 
// --------------------------------------------------------------
// OlapDetConstr
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#include "OlapDetConstr.hh"
#include "OlapGenerator.hh"

#include "SolidAnalyser.hh"
#include "OlapManager.hh"

#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4VisExtent.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Box.hh"
#include "G4Polyline.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4ThreeVector.hh"

OlapDetConstr::OlapDetConstr(G4VUserDetectorConstruction * aGeometry,
                             G4VPhysicalVolume *aWorld) :
  theTheta(0),
  thePhi(0),
  theAlpha(0),
  theNewWorldRot(0),
  theFullGeometry(aGeometry), 
  theWorld(aWorld),
  theNewWorld(0),
  theNewLV(0),
  theNewWorldBox(new G4Box("world",1.,1.,1.)),
  syncVis(true),
  nrLV(0)
{  
   //theNewWorldLV = new G4LogicalVolume 
   visMother =   new G4VisAttributes(G4Colour(7.0, 0.0, 0.0, 0.4));   // red
   visMother->SetForceWireframe(true);
   visDaughterA = new G4VisAttributes(G4Colour(.0, 0.0, 1.0, 0.3)); // blue
   visDaughterA->SetForceWireframe(false);
   visDaughterB = new G4VisAttributes(G4Colour(0., 0., 1.0, 0.3)); // light-blue
   visDaughterB->SetForceWireframe(false);
   visWorld =    new G4VisAttributes(G4Colour(1., 0, 0, 1.)); // red
   visWorld->SetForceWireframe(false);
   visWorld->SetVisibility(true);
   
   visFullWorld = new G4VisAttributes(G4Colour(1.,1.,1.,0.9));
   visFullWorld->SetForceWireframe(true);
   visWorld->SetVisibility(false);
   
   visInvisible = new G4VisAttributes(G4Colour(1.,1.,1.,0.9));
   visInvisible->SetForceWireframe(true);
   visInvisible->SetVisibility(false);
    
   visFirstLevel = new G4VisAttributes(G4Colour(1.,1.,1,1));
   visFirstLevel->SetForceWireframe(true);
   visFirstLevel->SetVisibility(true);
}


OlapDetConstr::~OlapDetConstr()
{
}


void OlapDetConstr::SetRotation(G4double t, G4double p, G4double a)
{
   theTheta = t;
   thePhi = p;
   theAlpha = a;
}


G4VPhysicalVolume * OlapDetConstr::Construct()
{
   if (!theWorld)
   {
     theWorld = theFullGeometry->Construct();
   }
   if (!theWorld)
   {
     G4Exception("OlapDetConstr::Construct()", "InvalidSetup", FatalException,
                 "Cannot create full world. Exiting...");
   }  
   
   // logical volume which will have a box with its dimensions derived from
   // the extent of the 'new-mother' volume positioned inside
   // the size of this NewWorld will be set in 'SetNewWorld(..)'
   theNewWorldLV = 
      new G4LogicalVolume(theNewWorldBox,
                          theWorld->GetLogicalVolume()->GetMaterial(),
                          "NewWorld"
                          );

   // physical volume which will serve as the 'small' NewWorld to be
   // used for overlap checking
   theNewWorld =
      new G4PVPlacement( 0,               // rotation
                         G4ThreeVector(), // translation
                         "NewWorld",      // name of phys = name of logical!!!
                         theNewWorldLV,   // own logical vol
                         0, false, 0      // no mother, not MANY, cpNr=0
                       );  
   
   theNewWorldLV->SetVisAttributes(visWorld);                       

   //ResetColors();
   
   
   theNewLV = theWorld->GetLogicalVolume();
   nrLV = G4LogicalVolumeStore::GetInstance()->size();

   return theWorld;
}


G4VPhysicalVolume * OlapDetConstr::SetNewWorld(G4LogicalVolume * aMotherLV,
                                               G4bool debugFlag)
{
   if (debugFlag)
   {
      G4cout << "Mother: " << aMotherLV->GetName() << " : " 
             << aMotherLV->GetSolid()->GetName() << G4endl;
   }

   //ML: ResetColors(theNewLV);
   theNewLV = aMotherLV;
   
   ConstructNewWorld();
   
   // set new world to Geant4-Kernel
   G4RunManager::GetRunManager()->DefineWorldVolume(theNewWorld);
   
   // try to set the OlapGenerator   
   const OlapGenerator * aGen = dynamic_cast<const OlapGenerator *>
         (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
   if ( aGen )
   {
     OlapGenerator * aNonConstGen = const_cast<OlapGenerator*>(aGen);
        
     // extent which allows arbitrary rotations ...
     G4Box * bx = dynamic_cast<G4Box*>
                (theNewWorld->GetLogicalVolume()->GetSolid());
     aNonConstGen->SetExtent(2.*bx->GetXHalfLength());
   }
   else
   {
      G4cout << "Warning: Primary generator is not OlapGenerator!" << G4endl
             << "Overlap Detection will not work!" << G4endl;
   }
   return theNewWorld;
}


G4VPhysicalVolume * OlapDetConstr::GetNewWorld()
{
   return theNewWorld; 
}


G4VPhysicalVolume * OlapDetConstr::GetFullWorld()
{
   return theWorld;
}


void OlapDetConstr::DrawPolyOutline()
{
   SolidAnalyser::GetSolidAnalyser();
   G4VSolid * solid = theNewLV->GetSolid();
   G4Polyline line, endline;
   G4double ex(0), first_r(0), first_z(0);
   G4double maxz(0), maxr(0);
   G4int nr(0);  
   G4Polycone  * s  = dynamic_cast<G4Polycone*>(solid);
   G4Polyhedra * sh = dynamic_cast<G4Polyhedra*>(solid);
   if (s) 
   { 
      nr = s->GetNumRZCorner();
      
      first_r = s->GetCorner(0).r;
      first_z = s->GetCorner(0).z;
      for (G4int i=0; i<nr; i++)
      {
          G4double r = s->GetCorner(i).r;
          G4double z = s->GetCorner(i).z;
          //line.push_back(G4ThreeVector(r,0.,z));  
          line.push_back(G4ThreeVector(r,z,0));  
          if (r>maxr)
            maxr=r;
          if (z>maxz)
            maxz=z;  
          G4cout << G4ThreeVector(r,0,z) << G4endl;
          //endline.push_back(G4ThreeVector(s->GetCorner(nr-1).r, 0.,
          //                                s->GetCorner(nr-1).z));
          endline.push_back(G4ThreeVector(s->GetCorner(nr-1).r,
                                          s->GetCorner(nr-1).z,0));
       }
    }
    else if (sh) 
    { 
      nr = sh->GetNumRZCorner();

      first_r = sh->GetCorner(0).r;
      first_z = sh->GetCorner(0).z;
      for (G4int i=0; i<nr; i++)
      {
          G4double r = sh->GetCorner(i).r;
          G4double z = sh->GetCorner(i).z;
          line.push_back(G4ThreeVector(r,0.,z));  
          if (r>maxr)
            maxr=r;
          if (z>maxz)
            maxz=z;  
          G4cout << G4ThreeVector(r,0,z) << G4endl;
          endline.push_back(G4ThreeVector(sh->GetCorner(nr-1).r, 0.,
                                          sh->GetCorner(nr-1).z));
      }
    }
    else
      return;

     //endline.push_back(G4ThreeVector(first_r,0.,first_z));
     endline.push_back(G4ThreeVector(first_r,first_z,0.));
     if (maxz>=maxr)
       ex = maxz;
     else
       ex = maxr;
     
     G4cout << "ex=" << ex << G4endl;  
     theNewWorldBox->SetXHalfLength(ex + ex/30.);
     theNewWorldBox->SetYHalfLength(ex + ex/30.);
     theNewWorldBox->SetZHalfLength(ex + ex/30.);
   
   G4Colour colour;
   colour = G4Colour(1.,0.,0.);
   line.SetVisAttributes(theNewLV->GetVisAttributes()->GetColor());
   G4VVisManager::GetConcreteInstance()->Draw(line);
   colour = G4Colour(0.,1.,0.);
   endline.SetVisAttributes(theNewLV->GetVisAttributes()->GetColor());
   //endline.SetVisAttributes(attribs3);
   G4VVisManager::GetConcreteInstance()->Draw(endline);
   colour = G4Colour(0.,0.,1.);
   G4Polyline ax;
   //ax.push_back(G4ThreeVector(0.,0.,-ex));
   //ax.push_back(G4ThreeVector(0.,0.,ex));
   ax.push_back(G4ThreeVector(0.,-ex,0.));
   ax.push_back(G4ThreeVector(0.,ex,0.));
   G4VisAttributes attribs2(colour);
   ax.SetVisAttributes(attribs2);
   G4VVisManager::GetConcreteInstance()->Draw(ax);
}

void OlapDetConstr::ConstructNewWorld()
{ 
   // delete the current NewWorld
   DeleteNewWorld();

   G4int nr = theNewLV->GetNoDaughters();
   
   // create a NewWorld with a 2% larger extent than the
   // solid of the NewMother
   
   const G4VisExtent & ext = theNewLV->GetSolid()->GetExtent();
   G4double extX = (ext.GetXmax()-ext.GetXmin())/2.;
   G4double extY = (ext.GetYmax()-ext.GetYmin())/2.;
   G4double extZ = (ext.GetZmax()-ext.GetZmin())/2.;
   G4ThreeVector pos( (ext.GetXmax()+ext.GetXmin())/2. ,
                      (ext.GetYmax()+ext.GetYmin())/2. ,
                      (ext.GetZmax()+ext.GetZmin())/2.
                     ); 
   pos=-pos;                     
  
   /*
   extX += extX/1000.;
   extY += extY/1000.;
   extZ += extZ/1000.;
   */
   G4double worldDim = std::sqrt(extX*extX + extY*extY + extZ*extZ);
   worldDim += worldDim/100.;
   // this automatically sets the dimensions of the 
   // solid in theNewWorldLV, because it's a ptr to the same solid
   theNewWorldBox->SetXHalfLength(worldDim);
   theNewWorldBox->SetYHalfLength(worldDim);
   theNewWorldBox->SetZHalfLength(worldDim);
      
   // create the mother with one layer of daughters
   G4LogicalVolume * aNewMotherLV =
     new G4LogicalVolume(theNewLV->GetSolid(),
                         theNewLV->GetMaterial(),
                         theNewLV->GetName()
                         );
   
   // below the mother is fixed colored instead as taken from vis-atts of the
   // full detector
   aNewMotherLV->SetVisAttributes(theNewLV->GetVisAttributes());

   delete theNewWorldRot;
   G4ThreeVector rotAxis(std::cos(thePhi)*std::sin(theTheta),
                         std::sin(thePhi)*std::sin(theTheta),
                         std::cos(theTheta));
   
   theNewWorldRot = new G4RotationMatrix(rotAxis,theAlpha);
   
   // G4VPhysicalVolume * aNewMotherPV =
   new G4PVPlacement( theNewWorldRot,   // rotation
                      pos,              // translation
                      aNewMotherLV,     // own logical vol
                      theNewLV->GetName(), // name of phys = name of logical!!!
                      
                      theNewWorldLV,    // mother logical vol
                      //aInbetweenLV,
                      
                      false, 0          // not MANY, cpNr=0
                      );  

  
  // loop over daughters an copy them into the new world ...
  //      - but don't iterate recursively over deeper layers of daughters!
  for (G4int i=0; i<nr; i++)
  {
    G4VPhysicalVolume * pv = theNewLV->GetDaughter(i);

    // ---------------------------------------------------------------
    // only do the stuff, if OlapManager.NoOlapMap[lv] is set to false!
    // otherwise the volume should not be considered for olap detection
    
    OlapManager * m = OlapManager::GetOlapManager();
    if (m->NoOlapMap[pv->GetLogicalVolume()])
      continue;  
    // ---------------------------------------------------------------                
    
    G4LogicalVolume * myLogical = 
      new G4LogicalVolume(pv->GetLogicalVolume()->GetSolid(),
                          pv->GetLogicalVolume()->GetMaterial(),
                          pv->GetLogicalVolume()->GetName()
                          );
                          
    // take vis atts from whole detector instead of fixed values!
    myLogical->SetVisAttributes(pv->GetLogicalVolume()->GetVisAttributes());
      
    // select between placements, replicas and parameterizations
    G4PVPlacement * pvPlace; 
    if (dynamic_cast<G4PVPlacement*>(pv)) 
    {
       pvPlace = dynamic_cast<G4PVPlacement*>(pv);
       new G4PVPlacement( pvPlace->GetRotation(),
                          pvPlace->GetTranslation(),
                          myLogical,     // my own logical vol
                          pvPlace->GetName(),
                          aNewMotherLV,  // my mother logical vol
                          false,         // not many
                          pvPlace->GetCopyNo()
                        );  
    } 
    else if ( pv->IsReplicated() && ! pv->GetParameterisation() ) // Replicas
    {
       // retrieve replication data
       EAxis aAxis;
       G4int nReplicas;
       G4double aWidth, anOffset;
       G4bool aDummy;
       pv->GetReplicationData(aAxis, nReplicas, aWidth, anOffset, aDummy);
       
       new G4PVReplica( pv->GetName(),
                        myLogical,
                        aNewMotherLV,
                        aAxis,
                        nReplicas,
                        aWidth,
                        anOffset
                      );        
    }
    else if ( pv->IsReplicated() && pv->GetParameterisation() )
    {
       // retrieve replication/parameterisation data
       EAxis aAxis;
       G4int nReplicas;
       G4double aWidth, anOffset;
       G4bool aDummy;
       pv->GetReplicationData(aAxis, nReplicas, aWidth, anOffset, aDummy);
       G4VPVParameterisation *pParam = pv->GetParameterisation();
       
       new G4PVParameterised( pv->GetName(),
                              myLogical,
                              aNewMotherLV,
                              aAxis,
                              nReplicas,
                              pParam
                            );
    }// end:placement?replica?parameterisation
    
  } // end:Loop over Daughters                                    
}


void OlapDetConstr::DeleteNewWorld()
{
   G4int nr = theNewWorldLV->GetNoDaughters();
   
   if (nr==0)
   {
     G4cout << "OlapDetConstr::DeleteNewWorld(): no daughter in NewWorld!"
            << G4endl << G4endl;
     return;
   }  
     
   if (nr>1 || nr<0)  // only one daughter is allowed!!!!
   {  
     G4Exception("OlapDetConstr::DeleteNewWorld()",
                 "InvalidSetup", FatalException,
                 "Too many daughters in NewWorldLV! Exiting...");
   }
      
   G4LogicalVolume * aMother;
   G4VPhysicalVolume * aDaughter;
   
   aMother = theNewWorldLV->GetDaughter(0)->GetLogicalVolume();
   
   G4int nrd = aMother->GetNoDaughters();
   for (G4int i=0; i<nrd; i++)
   {
     aDaughter = aMother->GetDaughter(i);
     G4LogicalVolume * tmp = aDaughter->GetLogicalVolume();
     delete tmp;
     delete aDaughter;
   }  
   
   aDaughter = theNewWorldLV->GetDaughter(0);
   theNewWorldLV->RemoveDaughter(aDaughter);
   delete aDaughter;
   delete aMother;
}


// set everything in the full geometry to invisible
void OlapDetConstr::ResetColors()
{
/*
   G4ColorMap * aColMap = G4ColorMap::GetColorMap();
   aColMap->SetAllInvisible();   
*/   
   //ColorFirstLevel();   
}


// colors the first layer of daughters of the full geometry
void OlapDetConstr::ColorFirstLevel()
{
/*  
   G4ColorMap * aColMap = G4ColorMap::GetColorMap();

   G4LogicalVolume * lvWorld = theWorld->GetLogicalVolume();
   aColMap->SetVisible(lvWorld,true);
   G4int nr = lvWorld->GetNoDaughters();
   for (G4int i = 0; i<nr; i++) {
     aColMap->SetVisible(lvWorld->GetDaughter(i)->GetLogicalVolume(),true);
   }  
*/   
}


// sets NewWorld invisible
void OlapDetConstr::ResetColors(G4LogicalVolume*)
{
/*
   G4ColorMap * aColMap = G4ColorMap::GetColorMap();
   G4int nr = lv->GetNoDaughters();
   for (G4int i =0; i<nr; i++) {
     aColMap->SetVisible(lv->GetDaughter(i)->GetLogicalVolume(),false);
   }
   aColMap->SetVisible(lv,false);
*/   
   //ColorFirstLevel();

}


G4VPhysicalVolume * OlapDetConstr::SetFullWorld()
{
   G4RunManager::GetRunManager()->DefineWorldVolume(theWorld);
   
   return theWorld;
}


G4LogicalVolume * OlapDetConstr::GetOriginalWorld()
{
   return theNewLV;
}   
