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
// Rich advanced example for Geant4
// RichTbAnalysisManager.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include "globals.hh"
#include "RichTbGeometryParameters.hh"

G4double PixRowNumSect[NumberOfPadHpdSiPixels];
G4double PixColNumSect[NumberOfPadHpdSiPixels];
G4bool PixelAtSectEdge[NumberOfPadHpdSiPixels];

G4double GetCurAerogelLength(G4int Aerognum ) {
  return 2.0*AgelHalfZ[Aerognum];
}
RichTbPadHpdSiPixPos::RichTbPadHpdSiPixPos(G4int ipixelnum) {
  //coord system such that Origin at the top of a sector;
  // which is near the center of the anode.
  //The X is along the bottom and Y goes radially outward (from top to
  //     bottom).
  


  icurpixel=ipixelnum;
  PadHpdSiPixPosX=PixColNumSect[ipixelnum]*XsizePix;
   if(ipixelnum ==  BigPixelNum ) {
  PadHpdSiPixPosY=PixRowNumSect[ipixelnum]*YsizePix+RowInitPointBigPixel;
  }else {
  PadHpdSiPixPosY=PixRowNumSect[ipixelnum]*YsizePix+RowInitPoint;

  }
}
RichTbPadHpdSiPixPos::~RichTbPadHpdSiPixPos(){ ; }
void InitializeRichTbGeometry(){
  PixRowNumSect[59-1]=0.0;
  PixRowNumSect[70-1]=1.0;
  PixRowNumSect[60-1]=2.0;
  PixRowNumSect[82-1]=3.0;
  PixRowNumSect[47-1]=3.0;
  PixRowNumSect[81-1]=4.0;
  PixRowNumSect[48-1]=4.0;
  
  PixRowNumSect[91-1]=5.0;
  PixRowNumSect[69-1]=5.0;
  PixRowNumSect[38-1]=5.0;
  
  PixRowNumSect[101-1]=6.0;
  PixRowNumSect[80-1]=6.0;
  PixRowNumSect[49-1]=6.0;
  PixRowNumSect[28-1]=6.0;
  
  PixRowNumSect[100-1]=7.0;
  PixRowNumSect[79-1]=7.0;
  PixRowNumSect[50-1]=7.0;
  PixRowNumSect[29-1]=7.0;
  
  PixRowNumSect[99-1]=8.0;
  PixRowNumSect[78-1]=8.0;
  PixRowNumSect[51-1]=8.0;
  PixRowNumSect[30-1]=8.0;
  
  PixRowNumSect[109-1]=9.0;
  PixRowNumSect[90-1]=9.0;
  PixRowNumSect[61-1]=9.0;
  PixRowNumSect[39-1]=9.0;
  PixRowNumSect[20-1]=9.0;
  
  PixRowNumSect[108-1]=10.0;
  PixRowNumSect[89-1]=10.0;
  PixRowNumSect[68-1]=10.0;
  PixRowNumSect[40-1]=10.0;
  PixRowNumSect[21-1]=10.0;
  
  PixRowNumSect[107-1]=11.0;
  PixRowNumSect[88-1]=11.0;
  PixRowNumSect[62-1]=11.0;
  PixRowNumSect[41-1]=11.0;
  PixRowNumSect[22-1]=11.0;
  
  PixRowNumSect[116-1]=12.0;
  PixRowNumSect[98-1]=12.0;
  PixRowNumSect[77-1]=12.0;
  PixRowNumSect[52-1]=12.0;
  PixRowNumSect[31-1]=12.0;
  PixRowNumSect[13-1]=12.0;
  
  PixRowNumSect[115-1]=13.0;
  PixRowNumSect[97-1]=13.0;
  PixRowNumSect[76-1]=13.0;
  PixRowNumSect[53-1]=13.0;
  PixRowNumSect[32-1]=13.0;
  PixRowNumSect[14-1]=13.0;
  
  PixRowNumSect[114-1]=14.0;
  PixRowNumSect[96-1]=14.0;
  PixRowNumSect[75-1]=14.0;
  PixRowNumSect[54-1]=14.0;
  PixRowNumSect[33-1]=14.0;
  PixRowNumSect[15-1]=14.0;
  
  PixRowNumSect[121-1]=15.0;
  PixRowNumSect[106-1]=15.0;
  PixRowNumSect[87-1]=15.0;
  PixRowNumSect[67-1]=15.0;
  PixRowNumSect[42-1]=15.0;
  PixRowNumSect[23-1]=15.0;
  PixRowNumSect[8-1]=15.0;
  
  PixRowNumSect[120-1]=16.0;
  PixRowNumSect[105-1]=16.0;
  PixRowNumSect[86-1]=16.0;
  PixRowNumSect[63-1]=16.0;
  PixRowNumSect[43-1]=16.0;
  PixRowNumSect[24-1]=16.0;
  PixRowNumSect[9-1]=16.0;
  
  PixRowNumSect[119-1]=17.0;
  PixRowNumSect[104-1]=17.0;
  PixRowNumSect[85-1]=17.0;
  PixRowNumSect[66-1]=17.0;
  PixRowNumSect[44-1]=17.0;
  PixRowNumSect[25-1]=17.0;
  PixRowNumSect[10-1]=17.0;
  
  PixRowNumSect[125-1]=18.0;
  PixRowNumSect[113-1]=18.0;
  PixRowNumSect[95-1]=18.0;
  PixRowNumSect[74-1]=18.0;
  PixRowNumSect[55-1]=18.0;
  PixRowNumSect[34-1]=18.0;
  PixRowNumSect[16-1]=18.0;
  PixRowNumSect[4-1]=18.0;
  
  PixRowNumSect[124-1]=19.0;
  PixRowNumSect[112-1]=19.0;
  PixRowNumSect[94-1]=19.0;
  PixRowNumSect[73-1]=19.0;
  PixRowNumSect[56-1]=19.0;
  PixRowNumSect[35-1]=19.0;
  PixRowNumSect[17-1]=19.0;
  PixRowNumSect[5-1]=19.0;
  
  PixRowNumSect[123-1]=20.0;
  PixRowNumSect[111-1]=20.0;
  PixRowNumSect[93-1]=20.0;
  PixRowNumSect[72-1]=20.0;
  PixRowNumSect[57-1]=20.0;
  PixRowNumSect[36-1]=20.0;
  PixRowNumSect[18-1]=20.0;
  PixRowNumSect[6-1]=20.0;
  
  PixRowNumSect[127-1]=21.0;
  PixRowNumSect[118-1]=21.0;
  PixRowNumSect[103-1]=21.0;
  PixRowNumSect[84-1]=21.0;
  PixRowNumSect[65-1]=21.0;
  PixRowNumSect[45-1]=21.0;
  PixRowNumSect[26-1]=21.0;
  PixRowNumSect[11-1]=21.0;
  PixRowNumSect[2-1]=21.0;
  
  PixRowNumSect[126-1]=22.0;
  PixRowNumSect[117-1]=22.0;
  PixRowNumSect[102-1]=22.0;
  PixRowNumSect[83-1]=22.0;
  PixRowNumSect[64-1]=22.0;
  PixRowNumSect[46-1]=22.0;
  PixRowNumSect[27-1]=22.0;
  PixRowNumSect[12-1]=22.0;
  PixRowNumSect[3-1]=22.0;
  
  PixRowNumSect[128-1]=23.0;
  PixRowNumSect[122-1]=23.0;
  PixRowNumSect[110-1]=23.0;
  PixRowNumSect[92-1]=23.0;
  PixRowNumSect[71-1]=23.0;
  PixRowNumSect[58-1]=23.0;
  PixRowNumSect[37-1]=23.0;
  PixRowNumSect[19-1]=23.0;
  PixRowNumSect[7-1]=23.0;
  PixRowNumSect[1-1]=23.0;



  PixColNumSect[59-1]=0.0;
  PixColNumSect[70-1]=0.0;
  PixColNumSect[60-1]=0.0;
  PixColNumSect[82-1]=-0.5;
  PixColNumSect[47-1]=0.5;
  PixColNumSect[81-1]=-0.5;
  PixColNumSect[48-1]=0.5;

  PixColNumSect[91-1]=-1.0;
  PixColNumSect[69-1]=0.0;
  PixColNumSect[38-1]=1.0;

  PixColNumSect[101-1]=-1.5;
  PixColNumSect[80-1]=-0.5;
  PixColNumSect[49-1]=0.5;
  PixColNumSect[28-1]=1.5;

  PixColNumSect[100-1]=-1.5;
  PixColNumSect[79-1]=-0.5;
  PixColNumSect[50-1]=0.5;
  PixColNumSect[29-1]=1.5;

  PixColNumSect[99-1]=-1.5;
  PixColNumSect[78-1]=-0.5;
  PixColNumSect[51-1]=0.5;
  PixColNumSect[30-1]=1.5;

  PixColNumSect[109-1]=-2.0;
  PixColNumSect[90-1]=-1.0;
  PixColNumSect[61-1]=0.0;
  PixColNumSect[39-1]=1.0;
  PixColNumSect[20-1]=2.0;

  PixColNumSect[108-1]=-2.0;
  PixColNumSect[89-1]=-1.0;
  PixColNumSect[68-1]=0.0;
  PixColNumSect[40-1]=1.0;
  PixColNumSect[21-1]=2.0;

  PixColNumSect[107-1]=-2.0;
  PixColNumSect[88-1]=-1.0;
  PixColNumSect[62-1]=0.0;
  PixColNumSect[41-1]=1.0;
  PixColNumSect[22-1]=2.0;

  PixColNumSect[116-1]=-2.5;
  PixColNumSect[98-1]=-1.5;
  PixColNumSect[77-1]=-0.5;
  PixColNumSect[52-1]=0.5;
  PixColNumSect[31-1]=1.5;
  PixColNumSect[13-1]=2.5;

  PixColNumSect[115-1]=-2.5;
  PixColNumSect[97-1]=-1.5;
  PixColNumSect[76-1]=-0.5;
  PixColNumSect[53-1]=0.5;
  PixColNumSect[32-1]=1.5;
  PixColNumSect[14-1]=2.5;

  PixColNumSect[114-1]=-2.5;
  PixColNumSect[96-1]=-1.5;
  PixColNumSect[75-1]=-0.5;
  PixColNumSect[54-1]=0.5;
  PixColNumSect[33-1]=1.5;
  PixColNumSect[15-1]=2.5;

  PixColNumSect[121-1]=-3.0;
  PixColNumSect[106-1]=-2.0;
  PixColNumSect[87-1]=-1.0;
  PixColNumSect[67-1]=0.0;
  PixColNumSect[42-1]=1.0;
  PixColNumSect[23-1]=2.0;
  PixColNumSect[8-1]=3.0;

  PixColNumSect[120-1]=-3.0;
  PixColNumSect[105-1]=-2.0;
  PixColNumSect[86-1]=-1.0;
  PixColNumSect[63-1]=0.0;
  PixColNumSect[43-1]=1.0;
  PixColNumSect[24-1]=2.0;
  PixColNumSect[9-1]=3.0;

  PixColNumSect[119-1]=-3.0;
  PixColNumSect[104-1]=-2.0;
  PixColNumSect[85-1]=-1.0;
  PixColNumSect[66-1]=0.0;
  PixColNumSect[44-1]=1.0;
  PixColNumSect[25-1]=2.0;
  PixColNumSect[10-1]=3.0;

  PixColNumSect[125-1]=-3.5;
  PixColNumSect[113-1]=-2.5;
  PixColNumSect[95-1]=-1.5;
  PixColNumSect[74-1]=-0.5;
  PixColNumSect[55-1]=0.5;
  PixColNumSect[34-1]=1.5;
  PixColNumSect[16-1]=2.5;
  PixColNumSect[4-1]=3.5;

  PixColNumSect[124-1]=-3.5;
  PixColNumSect[112-1]=-2.5;
  PixColNumSect[94-1]=-1.5;
  PixColNumSect[73-1]=-0.5;
  PixColNumSect[56-1]=0.5;
  PixColNumSect[35-1]=1.5;
  PixColNumSect[17-1]=2.5;
  PixColNumSect[5-1]=3.5;

  PixColNumSect[123-1]=-3.5;
  PixColNumSect[111-1]=-2.5;
  PixColNumSect[93-1]=-1.5;
  PixColNumSect[72-1]=-0.5;
  PixColNumSect[57-1]=0.5;
  PixColNumSect[36-1]=1.5;
  PixColNumSect[18-1]=2.5;
  PixColNumSect[6-1]=3.5;

  PixColNumSect[127-1]=-4.0;
  PixColNumSect[118-1]=-3.0;
  PixColNumSect[103-1]=-2.0;
  PixColNumSect[84-1]=-1.0;
  PixColNumSect[65-1]=0.0;
  PixColNumSect[45-1]=1.0;
  PixColNumSect[26-1]=2.0;
  PixColNumSect[11-1]=3.0;
  PixColNumSect[2-1]=4.0;

  PixColNumSect[126-1]=-4.0;
  PixColNumSect[117-1]=-3.0;
  PixColNumSect[102-1]=-2.0;
  PixColNumSect[83-1]=-1.0;
  PixColNumSect[64-1]=0.0;
  PixColNumSect[46-1]=1.0;
  PixColNumSect[27-1]=2.0;
  PixColNumSect[12-1]=3.0;
  PixColNumSect[3-1]=4.0;

  PixColNumSect[128-1]=-4.5;
  PixColNumSect[122-1]=-3.5;
  PixColNumSect[110-1]=-2.5;
  PixColNumSect[92-1]=-1.5;
  PixColNumSect[71-1]=-0.5;
  PixColNumSect[58-1]=0.5;
  PixColNumSect[37-1]=1.5;
  PixColNumSect[19-1]=2.5;
  PixColNumSect[7-1]=3.5;
  PixColNumSect[1-1]=4.5;

for (G4int ipixel=0;ipixel<NumberOfPadHpdSiPixels; ipixel++){
  PixelAtSectEdge[ipixel]=false;
  }
 
PixelAtSectEdge[59-1]=true;
PixelAtSectEdge[70-1]=true;
PixelAtSectEdge[59-1]=true;
PixelAtSectEdge[60-1]=true;
PixelAtSectEdge[82-1]=true;
PixelAtSectEdge[47-1]=true;
PixelAtSectEdge[81-1]=true;
PixelAtSectEdge[48-1]=true;
PixelAtSectEdge[91-1]=true;
PixelAtSectEdge[38-1]=true;
PixelAtSectEdge[101-1]=true;
PixelAtSectEdge[28-1]=true;

PixelAtSectEdge[100-1]=true;
PixelAtSectEdge[29-1]=true;

PixelAtSectEdge[99-1]=true;
PixelAtSectEdge[30-1]=true;

PixelAtSectEdge[109-1]=true;
PixelAtSectEdge[20-1]=true;

PixelAtSectEdge[108-1]=true;
PixelAtSectEdge[21-1]=true;

PixelAtSectEdge[107-1]=true;
PixelAtSectEdge[22-1]=true;

PixelAtSectEdge[116-1]=true;
PixelAtSectEdge[13-1]=true;

PixelAtSectEdge[115-1]=true;
PixelAtSectEdge[14-1]=true;

PixelAtSectEdge[114-1]=true;
PixelAtSectEdge[15-1]=true;

PixelAtSectEdge[121-1]=true;
PixelAtSectEdge[8-1]=true;

PixelAtSectEdge[120-1]=true;
PixelAtSectEdge[9-1]=true;

PixelAtSectEdge[119-1]=true;
PixelAtSectEdge[10-1]=true;

PixelAtSectEdge[125-1]=true;
PixelAtSectEdge[4-1]=true;

PixelAtSectEdge[124-1]=true;
PixelAtSectEdge[5-1]=true;

PixelAtSectEdge[123-1]=true;
PixelAtSectEdge[6-1]=true;

PixelAtSectEdge[127-1]=true;
PixelAtSectEdge[2-1]=true;

PixelAtSectEdge[126-1]=true;
PixelAtSectEdge[3-1]=true;

PixelAtSectEdge[128-1]=true;
PixelAtSectEdge[1-1]=true;
   
   
}







