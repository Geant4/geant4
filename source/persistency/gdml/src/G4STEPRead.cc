#include "G4STEPRead.hh"

G4TessellatedSolid* G4STEPRead::ReadTessellatedSolid(const std::string& line) {

   std::istringstream stream(line.substr(2));
   
   G4String name;
   stream >> name;
   
   return new G4TessellatedSolid(name);
}

G4VFacet* G4STEPRead::ReadFacet(const std::string& line) {

   if (line[2]=='3') { // Triangular facet

      G4double x1,y1,z1;
      G4double x2,y2,z2;
      G4double x3,y3,z3;
      
      std::istringstream stream(line.substr(4));
      stream >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
      return new G4TriangularFacet(G4ThreeVector(x1,y1,z1),G4ThreeVector(x2,y2,z2),G4ThreeVector(x3,y3,z3),ABSOLUTE);
   }

   if (line[2]=='4') { // Quadrangular facet

      G4double x1,y1,z1;
      G4double x2,y2,z2;
      G4double x3,y3,z3;
      G4double x4,y4,z4;

      std::istringstream stream(line.substr(4));
      stream >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3 >> x4 >> y4 >> z4;
      return new G4QuadrangularFacet(G4ThreeVector(x1,y1,z1),G4ThreeVector(x2,y2,z2),G4ThreeVector(x3,y3,z3),G4ThreeVector(x4,y4,z4),ABSOLUTE);
   }

   G4Exception("G4STEPRead: ERROR! Number of vertices per facet should be either 3 or 4!");
   return 0;
}

void G4STEPRead::ReadGeom(const G4String& name) {

   G4cout << "G4STEPRead: Reading '" << name << "'..." << G4endl;

   std::ifstream GeomFile(name);
   
   if (!GeomFile) G4Exception("G4STEPRead: ERROR! Can not open file: "+name);

   std::string line;
   
   G4TessellatedSolid* tessellated = 0;
   
   while (getline(GeomFile,line)) {
   
      if (line[0] == 'f') { tessellated = ReadTessellatedSolid(line); } else
      if (line[0] == 'p') { tessellated->AddFacet(ReadFacet(line)); }
   }

   G4cout << "G4STEPRead: Reading '" << name << "' done." << G4endl;
}
