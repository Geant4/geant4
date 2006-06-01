// -------------------------------------------------------------------
// $Id: MicrobeamPhantomConfiguration.cc,v 1.3 2006-06-01 22:25:20 sincerti Exp $
// -------------------------------------------------------------------

#include "MicrobeamPhantomConfiguration.hh"

G4int MicrobeamPhantomConfiguration::phantomTotalPixels = 0;
G4int MicrobeamPhantomConfiguration::nucleusTotalPixels = 0;
G4int MicrobeamPhantomConfiguration::cytoplasmTotalPixels = 0;
G4float MicrobeamPhantomConfiguration::dx = 0;
G4float MicrobeamPhantomConfiguration::dy = 0;
G4float MicrobeamPhantomConfiguration::dz = 0;
G4float MicrobeamPhantomConfiguration::nucleusMass = 0;
G4float MicrobeamPhantomConfiguration::cytoplasmMass = 0;

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
  G4float vx, vy, vz, tmp, mat, den, density;
  G4float denCyto1, denCyto2, denCyto3, denNucl1, denNucl2, denNucl3;
  FILE* fMap;
  
  phantomTotalPixels=0;
  nucleusTotalPixels=0;
  cytoplasmTotalPixels=0;
  dx=0;
  dy=0;
  dz=0;
  nucleusMass=0;
  cytoplasmMass=0;

  // READ PHANTOM PARAMETERS
  fMap = fopen("phantom.dat","r");

  ncols = fscanf(fMap,"%i %i %i",&phantomTotalPixels, &nucleusTotalPixels, &cytoplasmTotalPixels);
  ncols = fscanf(fMap,"%f %f %f",&dx, &dy, &dz);
  ncols = fscanf(fMap,"%f %f %f",&tmp, &tmp, &tmp);
  ncols = fscanf(fMap,"%f %f %f",&denCyto1, &denCyto2, &denCyto3);
  ncols = fscanf(fMap,"%f %f %f",&denNucl1, &denNucl2, &denNucl3);
  dx = dx * micrometer;
  dy = dy * micrometer;
  dz = dz * micrometer;
  voxelThreeVector = new G4ThreeVector [phantomTotalPixels];

  for (G4int i=0; i<phantomTotalPixels; i++) 
  { 
    ncols = fscanf(fMap,"%f %f %f %f %f %f",&vx, &vy, &vz, &mat, &den, &tmp);

    if (mat==2) // NUCLEUS
    	{
	  if (den==1) density = denNucl1;
	  if (den==2) density = denNucl2;
	  if (den==3) density = denNucl3;
	  nucleusMass   = nucleusMass   + density * dx * dy * dz ;
    	}

    if (mat==1) // CYTOPLASM
    	{ 
	  if (den==1) density = denCyto1;
	  if (den==2) density = denCyto2;
	  if (den==3) density = denCyto3;
	  cytoplasmMass = cytoplasmMass + density * dx * dy * dz ;
	}
    
    G4ThreeVector v(vx,vy,vz);
    voxelThreeVector[i] = v;
  }

  fclose(fMap);
  
  nucleusMass   = nucleusMass * 1e-6 ;
  cytoplasmMass = cytoplasmMass * 1e-6 ;
  
  return 0;
}


