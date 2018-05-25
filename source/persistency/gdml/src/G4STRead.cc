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
// $Id: G4STRead.cc 110108 2018-05-15 11:46:54Z gcosmo $
//
// class G4STRead Implementation
//
// History:
// - Created.                                  Zoltan Torzsok, November 2007
// -------------------------------------------------------------------------

#include <fstream>

#include "G4STRead.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4QuadrangularFacet.hh"
#include "G4TriangularFacet.hh"
#include "G4TessellatedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

void G4STRead::TessellatedRead(const std::string& line)
{
   if (tessellatedList.size()>0)
   {
     tessellatedList.back()->SetSolidClosed(true);
       // Finish the previous solid at first!
   }

   std::istringstream stream(line.substr(2));
   
   G4String name;
   stream >> name;
   
   G4TessellatedSolid* tessellated = new G4TessellatedSolid(name);
   volumeMap[tessellated] =
     new G4LogicalVolume(tessellated, solid_material, name+"_LV" , 0, 0, 0);
   tessellatedList.push_back(tessellated);

#ifdef G4VERBOSE
   G4cout << "G4STRead: Reading solid: " << name << G4endl;
#endif
}

void G4STRead::FacetRead(const std::string& line)
{
   if (tessellatedList.size()==0)
   {
     G4Exception("G4STRead::FacetRead()", "ReadError", FatalException,
                 "A solid must be defined before defining a facet!");
   }

   if (line[2]=='3')  // Triangular facet
   {
      G4double x1,y1,z1;
      G4double x2,y2,z2;
      G4double x3,y3,z3;
      
      std::istringstream stream(line.substr(4));
      stream >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
      tessellatedList.back()->
        AddFacet(new G4TriangularFacet(G4ThreeVector(x1,y1,z1),
                                       G4ThreeVector(x2,y2,z2),
                                       G4ThreeVector(x3,y3,z3), ABSOLUTE));
   }
   else if (line[2]=='4')  // Quadrangular facet
   {
      G4double x1,y1,z1;
      G4double x2,y2,z2;
      G4double x3,y3,z3;
      G4double x4,y4,z4;

      std::istringstream stream(line.substr(4));
      stream >> x1 >> y1 >> z1 >> x2 >> y2 >> z2
             >> x3 >> y3 >> z3 >> x4 >> y4 >> z4;
      tessellatedList.back()->
        AddFacet(new G4QuadrangularFacet(G4ThreeVector(x1,y1,z1),
                                         G4ThreeVector(x2,y2,z2),
                                         G4ThreeVector(x3,y3,z3),
                                         G4ThreeVector(x4,y4,z4), ABSOLUTE));
   }
   else
   {
     G4Exception("G4STRead::FacetRead()", "ReadError", FatalException,
                 "Number of vertices per facet should be either 3 or 4!");
   }
}

void G4STRead::PhysvolRead(const std::string& line)
{
   G4int level;
   G4String name;
   G4double r1,r2,r3;
   G4double r4,r5,r6;
   G4double r7,r8,r9;
   G4double pX,pY,pZ;
   G4double n1,n2,n3,n4,n5;

   std::istringstream stream(line.substr(2));
   stream >> level >> name >> r1 >> r2 >> r3 >> n1 >> r4 >> r5 >> r6
          >> n2 >> r7 >> r8 >> r9 >> n3 >> pX >> pY >> pZ >> n4 >> n5;
   std::string::size_type idx = name.rfind("_");
   if (idx!=std::string::npos)
   {
     name.resize(idx);
   }
   else
   {
     G4Exception("G4STRead::PhysvolRead()", "ReadError",
                 FatalException, "Invalid input stream!");
     return;
   }

   G4cout << "G4STRead: Placing tessellated solid: " << name << G4endl;

   G4TessellatedSolid* tessellated = 0;

   for (size_t i=0; i<tessellatedList.size(); i++)
   {                                    // Find the volume for this physvol!
      if (tessellatedList[i]->GetName() == G4String(name))
      {
         tessellated = tessellatedList[i];
	 break;
      }
   }

   if (tessellated == 0)
   {
     G4String error_msg = "Referenced solid '" + name + "' not found!";
     G4Exception("G4STRead::PhysvolRead()", "ReadError",
                 FatalException, error_msg);
   }
   if (volumeMap.find(tessellated) == volumeMap.end())
   {
     G4String error_msg = "Referenced solid '" + name
                        + "' is not associated with a logical volume!";
     G4Exception("G4STRead::PhysvolRead()", "InvalidSetup",
                 FatalException, error_msg);
   }
   const G4RotationMatrix rot(G4ThreeVector(r1,r2,r3),
                              G4ThreeVector(r4,r5,r6),
                              G4ThreeVector(r7,r8,r9));
   const G4ThreeVector pos(pX,pY,pZ);

   new G4PVPlacement(G4Transform3D(rot.inverse(),pos),
                     volumeMap[tessellated], name+"_PV", world_volume, 0, 0);
     // Note: INVERSE of rotation is needed!!!

   G4double minx,miny,minz;
   G4double maxx,maxy,maxz;
   const G4VoxelLimits limits;

   tessellated->CalculateExtent(kXAxis,limits,
                G4AffineTransform(rot,pos),minx,maxx);
   tessellated->CalculateExtent(kYAxis,limits,
                G4AffineTransform(rot,pos),miny,maxy);
   tessellated->CalculateExtent(kZAxis,limits,
                G4AffineTransform(rot,pos),minz,maxz);

   if (world_extent.x() < std::fabs(minx))
     { world_extent.setX(std::fabs(minx)); }
   if (world_extent.y() < std::fabs(miny))
     { world_extent.setY(std::fabs(miny)); }
   if (world_extent.z() < std::fabs(minz))
     { world_extent.setZ(std::fabs(minz)); }
   if (world_extent.x() < std::fabs(maxx))
     { world_extent.setX(std::fabs(maxx)); }
   if (world_extent.y() < std::fabs(maxy))
     { world_extent.setY(std::fabs(maxy)); }
   if (world_extent.z() < std::fabs(maxz))
     { world_extent.setZ(std::fabs(maxz)); }
}

void G4STRead::ReadGeom(const G4String& name)
{
#ifdef G4VERBOSE
   G4cout << "G4STRead: Reading '" << name << "'..." << G4endl;
#endif
   std::ifstream GeomFile(name);
   
   if (!GeomFile)
   {
      G4String error_msg = "Cannot open file: " + name;
      G4Exception("G4STRead::ReadGeom()", "ReadError",
                  FatalException, error_msg);
   }

   tessellatedList.clear();
   volumeMap.clear();
   std::string line;
   
   while (getline(GeomFile,line))
   {
      if (line[0] == 'f') { TessellatedRead(line); } else
      if (line[0] == 'p') { FacetRead(line); }
   }

   if (tessellatedList.size()>0)   // Finish the last solid!
   {
      tessellatedList.back()->SetSolidClosed(true);
   }

   G4cout << "G4STRead: Reading '" << name << "' done." << G4endl;
}

void G4STRead::ReadTree(const G4String& name)
{
#ifdef G4VERBOSE
   G4cout << "G4STRead: Reading '" << name << "'..." << G4endl;
#endif
   std::ifstream TreeFile(name);

   if (!TreeFile)
   {
      G4String error_msg = "Cannot open file: " + name;
      G4Exception("G4STRead::ReadTree()", "ReadError",
                  FatalException, error_msg);
   }

   std::string line;
   
   while (getline(TreeFile,line))
   {
      if (line[0] == 'g')  { PhysvolRead(line); }
   }

   G4cout << "G4STRead: Reading '" << name << "' done." << G4endl;
}

G4LogicalVolume*
G4STRead::Read(const G4String& name, G4Material* mediumMaterial,
                                     G4Material* solidMaterial)
{
   if (mediumMaterial == 0)
   {
     G4Exception("G4STRead::Read()", "InvalidSetup", FatalException,
                 "Pointer to medium material is not valid!");
   }
   if (solidMaterial == 0)
   {
     G4Exception("G4STRead::Read()", "InvalidSetup", FatalException,
                 "Pointer to solid material is not valid!");
   }

   solid_material = solidMaterial;

   world_box = new G4Box("TessellatedWorldBox",kInfinity,kInfinity,kInfinity);
     // We don't know the extent of the world yet!
   world_volume = new G4LogicalVolume(world_box, mediumMaterial,
                                      "TessellatedWorldLV", 0, 0, 0);
   world_extent = G4ThreeVector(0,0,0);

   ReadGeom(name+".geom");
   ReadTree(name+".tree");

   // Now setting the world extent ...
   //
   if (world_box->GetXHalfLength() > world_extent.x())
     { world_box->SetXHalfLength(world_extent.x()); }
   if (world_box->GetYHalfLength() > world_extent.y())
     { world_box->SetYHalfLength(world_extent.y()); }
   if (world_box->GetZHalfLength() > world_extent.z())
     { world_box->SetZHalfLength(world_extent.z()); }

   return world_volume;
}
