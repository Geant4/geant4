#include "G4STEPRead.hh"

void G4STEPRead::tessellatedRead(const std::string& line) {

   if (tessellatedList.size()>0) tessellatedList.back()->SetSolidClosed(true); // Finish the previous solid at first!

   std::istringstream stream(line.substr(2));
   
   G4String name;
   stream >> name;
   
   G4TessellatedSolid* tessellated = new G4TessellatedSolid(name);
   volumeMap[tessellated] = new G4LogicalVolume(tessellated,solid_material,"volume"+name,0,0,0);
   tessellatedList.push_back(tessellated);

   G4cout << "G4STEPRead: Reading solid: " << name << G4endl;
}

void G4STEPRead::facetRead(const std::string& line) {

   if (tessellatedList.size()==0) G4Exception("G4STEPRead: ERROR! A solid must be defined before defining a facet!");

   if (line[2]=='3') { // Triangular facet

      G4double x1,y1,z1;
      G4double x2,y2,z2;
      G4double x3,y3,z3;
      
      std::istringstream stream(line.substr(4));
      stream >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
      tessellatedList.back()->AddFacet(new G4TriangularFacet(G4ThreeVector(x1,y1,z1),G4ThreeVector(x2,y2,z2),G4ThreeVector(x3,y3,z3),ABSOLUTE));
   } else
   if (line[2]=='4') { // Quadrangular facet

      G4double x1,y1,z1;
      G4double x2,y2,z2;
      G4double x3,y3,z3;
      G4double x4,y4,z4;

      std::istringstream stream(line.substr(4));
      stream >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3 >> x4 >> y4 >> z4;
      tessellatedList.back()->AddFacet(new G4QuadrangularFacet(G4ThreeVector(x1,y1,z1),G4ThreeVector(x2,y2,z2),G4ThreeVector(x3,y3,z3),G4ThreeVector(x4,y4,z4),ABSOLUTE));
   } else
   G4Exception("G4STEPRead: ERROR! Number of vertices per facet should be either 3 or 4!");
}

void G4STEPRead::physvolRead(const std::string& line) {

   G4int level;
   G4String name;
   G4double r1,r2,r3;
   G4double r4,r5,r6;
   G4double r7,r8,r9;
   G4double pX,pY,pZ;
   G4double n1,n2,n3,n4,n5;

   std::istringstream stream(line.substr(2));
   stream >> level >> name >> r1 >> r2 >> r3 >> n1 >> r4 >> r5 >> r6 >> n2 >> r7 >> r8 >> r9 >> n3 >> pX >> pY >> pZ >> n4 >> n5;

   name.resize(name.rfind("_"));

   G4cout << "G4STEPRead: Placing solid: " << name << G4endl;

   G4TessellatedSolid* tessellated = 0;

   for (size_t i=0;i<tessellatedList.size();i++) { // Find the volume for this physvol!
   
      if (tessellatedList[i]->GetName() == G4String(name)) {
      
         tessellated = tessellatedList[i];
	 break;
      }
   }

   if (tessellated == 0) G4Exception("G4STEPRead: ERROR! Referenced solid '"+name+"' not found!");
   if (volumeMap.find(tessellated) == volumeMap.end()) G4Exception("G4STEPRead: ERROR! Referenced solid '"+name+"' is not associated with a logical volume!");

   const G4RotationMatrix rot(G4ThreeVector(r1,r2,r3),G4ThreeVector(r4,r5,r6),G4ThreeVector(r7,r8,r9));
   const G4ThreeVector pos(pX,pY,pZ);

   new G4PVPlacement(G4Transform3D(rot.inverse(),pos),volumeMap[tessellated],"physvol"+name,world_volume,0,0); // Note: INVERSE of rotation is needed!!!

   G4double minx,miny,minz;
   G4double maxx,maxy,maxz;
   const G4VoxelLimits limits;

   tessellated->CalculateExtent(kXAxis,limits,G4AffineTransform(rot,pos),minx,maxx);
   tessellated->CalculateExtent(kYAxis,limits,G4AffineTransform(rot,pos),miny,maxy);
   tessellated->CalculateExtent(kZAxis,limits,G4AffineTransform(rot,pos),minz,maxz);

   if (world_extent.x() < fabs(minx)) world_extent.setX(fabs(minx));
   if (world_extent.y() < fabs(miny)) world_extent.setY(fabs(miny));
   if (world_extent.z() < fabs(minz)) world_extent.setZ(fabs(minz));
   if (world_extent.x() < fabs(maxx)) world_extent.setX(fabs(maxx));
   if (world_extent.y() < fabs(maxy)) world_extent.setY(fabs(maxy));
   if (world_extent.z() < fabs(maxz)) world_extent.setZ(fabs(maxz));
}

void G4STEPRead::ReadGeom(const G4String& name) {

   G4cout << "G4STEPRead: Reading '" << name << "'..." << G4endl;

   std::ifstream GeomFile(name);
   
   if (!GeomFile) G4Exception("G4STEPRead: ERROR! Can not open file: "+name);

   tessellatedList.clear();
   volumeMap.clear();
   std::string line;
   
   while (getline(GeomFile,line)) {
   
      if (line[0] == 'f') tessellatedRead(line); else
      if (line[0] == 'p') facetRead(line);
   }

   if (tessellatedList.size()>0) tessellatedList.back()->SetSolidClosed(true); // Finish the last solid!

   G4cout << "G4STEPRead: Reading '" << name << "' done." << G4endl;
}

void G4STEPRead::ReadTree(const G4String& name) {

   G4cout << "G4STEPRead: Reading '" << name << "'..." << G4endl;

   std::ifstream TreeFile(name);

   if (!TreeFile) G4Exception("G4STEPRead: ERROR! Can not open file: "+name);

   std::string line;
   
   while (getline(TreeFile,line)) {
   
      if (line[0] == 'g') physvolRead(line);
   }

   G4cout << "G4STEPRead: Reading '" << name << "' done." << G4endl;
}

G4VPhysicalVolume* G4STEPRead::Read(const G4String& name,G4Material* mediumMaterial,G4Material* solidMaterial) {

   if (mediumMaterial == 0) G4Exception("G4STEPRead: ERROR! Pointer to medium material is nod valid!");
   if (solidMaterial == 0) G4Exception("G4STEPRead: ERROR! Pointer to solid material is nod valid!");

   solid_material = solidMaterial;

   world_box = new G4Box("WorldBox",kInfinity,kInfinity,kInfinity); // We don't know the extent of the world yet!
   world_volume = new G4LogicalVolume(world_box,mediumMaterial,"WorldLogicalVolume",0,0,0);
   world_extent = G4ThreeVector(0,0,0);

   ReadGeom(name+".geom");
   ReadTree(name+".tree");

   if (world_box->GetXHalfLength() > world_extent.x()) { world_box->SetXHalfLength(world_extent.x()); } // Now we know the world extent!
   if (world_box->GetYHalfLength() > world_extent.y()) { world_box->SetYHalfLength(world_extent.y()); }
   if (world_box->GetZHalfLength() > world_extent.z()) { world_box->SetZHalfLength(world_extent.z()); }

   return new G4PVPlacement(0,G4ThreeVector(0,0,0),world_volume,"WorldPhysicalVolume",0,0,0);
}
