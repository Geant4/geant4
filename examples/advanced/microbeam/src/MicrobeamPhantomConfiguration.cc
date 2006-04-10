// -------------------------------------------------------------------
// $Id: MicrobeamPhantomConfiguration.cc,v 1.2 2006-04-10 14:47:32 sincerti Exp $
// -------------------------------------------------------------------

#include "MicrobeamPhantomConfiguration.hh"

G4int MicrobeamPhantomConfiguration::phantomTotalPixels = 0;
G4int MicrobeamPhantomConfiguration::nucleusTotalPixels = 0;
G4int MicrobeamPhantomConfiguration::cytoplasmTotalPixels = 0;
G4float MicrobeamPhantomConfiguration::dx = 0;
G4float MicrobeamPhantomConfiguration::dy = 0;
G4float MicrobeamPhantomConfiguration::dz = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MicrobeamPhantomConfiguration::MicrobeamPhantomConfiguration() {
Initialize();
}

MicrobeamPhantomConfiguration::~MicrobeamPhantomConfiguration()
{
  delete[] voxelThreeVector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int MicrobeamPhantomConfiguration::Initialize() {

  G4int ncols;
  G4float vx, vy, vz, tmp;
  FILE* fMap;
  phantomTotalPixels=0;
  nucleusTotalPixels=0;
  cytoplasmTotalPixels=0;
  dx=0;
  dy=0;
  dz=0;
  
  // READ PHANTOM PARAMETERS
  fMap = fopen("phantom.dat","r");

    ncols = fscanf(fMap,"%i %i %i",&phantomTotalPixels, &nucleusTotalPixels, &cytoplasmTotalPixels);
    ncols = fscanf(fMap,"%f %f %f",&dx, &dy, &dz);
    ncols = fscanf(fMap,"%f %f %f",&tmp, &tmp, &tmp);
    dx = dx * micrometer;
    dy = dy * micrometer;
    dz = dz * micrometer;
    voxelThreeVector = new G4ThreeVector [phantomTotalPixels];
    for (G4int i=0; i<phantomTotalPixels; i++) 
     { 
       ncols = fscanf(fMap,"%f %f %f %f %f %f",&vx, &vy, &vz, &tmp, &tmp, &tmp);
       G4ThreeVector v(vx,vy,vz);
       voxelThreeVector[i] = v;
     }

  fclose(fMap);

  return 0;
}


