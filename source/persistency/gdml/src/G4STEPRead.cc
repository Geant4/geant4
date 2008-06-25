#include "G4STEPRead.hh"

void G4STEPRead::tessellatedRead(const std::string& line) {

   std::istringstream stream(line.substr(2));
   
   G4String name;
   stream >> name;
   
   G4TessellatedSolid* tessellated = new G4TessellatedSolid(name);
   G4LogicalVolume* volume = new G4LogicalVolume(tessellated,solid_material,name,0,0,0);

   G4cout << "G4STEPRead: Reading solid: " << name << G4endl;

   tessellatedList.push_back(tessellated);
   volumeList.push_back(volume);
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

void G4STEPRead::ReadGeom(const G4String& name) {

   solid_material = new G4Material("Aluminium",13,26.98*g/mole,2.7*g/cm3,kStateSolid);

   G4cout << "G4STEPRead: Reading '" << name << "'..." << G4endl;

   std::ifstream GeomFile(name);
   
   if (!GeomFile) G4Exception("G4STEPRead: ERROR! Can not open file: "+name);

   tessellatedList.clear();
   volumeList.clear();
   std::string line;
   
   while (getline(GeomFile,line)) {
   
      if (line[0] == 'f') tessellatedRead(line); else
      if (line[0] == 'p') facetRead(line);
   }

   G4cout << "G4STEPRead: Reading '" << name << "' done." << G4endl;
}
