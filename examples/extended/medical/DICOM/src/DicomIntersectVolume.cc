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
/// \file medical/DICOM/src/DicomIntersectVolume.cc
/// \brief Implementation of the DicomIntersectVolume class
//

#include "DicomIntersectVolume.hh"

#include "G4UIcmdWithAString.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VSolid.hh"
#include "G4tgrSolid.hh"
#include "G4tgbVolume.hh"
#include "G4Material.hh"
#include "G4PhantomParameterisation.hh"
#include "G4PVParameterised.hh"
#include "G4UIcommand.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomIntersectVolume::DicomIntersectVolume()
 : G4UImessenger(),
   fG4VolumeCmd(0),
   fSolid(0),
   fPhantomMinusCorner(),
   fVoxelIsInside(0)
{
  fUserVolumeCmd= new G4UIcmdWithAString("/dicom/intersectWithUserVolume",this);
  fUserVolumeCmd->SetGuidance("Intersects a phantom with a user-defined volume "
     "and outputs the voxels that are totally inside the intersection as"
     " a new phantom file. It must have the parameters: POS_X POS_Y POS_Z "
     "ANG_X ANG_Y ANG_Z SOLID_TYPE SOLID_PARAM_1 (SOLID_PARAM_2 ...)");
  fUserVolumeCmd->SetParameterName("choice",true);
  fUserVolumeCmd->AvailableForStates(G4State_Idle);

  fG4VolumeCmd = new G4UIcmdWithAString("/dicom/intersectWithUserVolume",this);
  fG4VolumeCmd->SetGuidance("Intersects a phantom with a user-defined volume"
     " and outputs the voxels that are totally inside the intersection as "
     " a new phantom file. It must have the parameters: POS_X POS_Y POS_Z "
     "ANG_X ANG_Y ANG_Z SOLID_TYPE SOLID_PARAM_1 (SOLID_PARAM_2 ...)");
  fG4VolumeCmd->SetParameterName("choice",true);
  fG4VolumeCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomIntersectVolume::~DicomIntersectVolume()
{
 if (fUserVolumeCmd) delete fUserVolumeCmd;
 if (fG4VolumeCmd) delete fG4VolumeCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomIntersectVolume::SetNewValue(G4UIcommand * command,
                                             G4String newValues)    
{ 
  G4AffineTransform theVolumeTransform;

  if (command == fUserVolumeCmd)
  {
    std::vector<G4String> params = GetWordsInString( newValues );
    if( params.size() < 8 ) {
      G4Exception("DicomIntersectVolume::SetNewValue",
          " There must be at least 8 parameter: SOLID_TYPE POS_X POS_Y POS_Z "
          "ANG_X ANG_Y ANG_Z SOLID_PARAM_1 (SOLID_PARAM_2 ...)",
                  FatalErrorInArgument,
                  G4String("Number of parameters given = " + 
               G4UIcommand::ConvertToString( G4int(params.size()) )).c_str());
      
    }

    //----- Build G4VSolid 
    BuildUserSolid(params);

    //----- Calculate volume inverse 3D transform
    G4ThreeVector pos = G4ThreeVector( G4UIcommand::ConvertToDouble(params[0]), 
                       G4UIcommand::ConvertToDouble(params[1]),
                       G4UIcommand::ConvertToDouble(params[2]) );
    G4RotationMatrix* rotmat = new G4RotationMatrix;
    std::vector<G4double> angles; 
    rotmat->rotateX( G4UIcommand::ConvertToDouble(params[3]) );
    rotmat->rotateY( G4UIcommand::ConvertToDouble(params[4]) );
    rotmat->rotateY( G4UIcommand::ConvertToDouble(params[5]) );
    theVolumeTransform = G4AffineTransform( rotmat, pos ).Invert();

  }
  else if (command == fG4VolumeCmd)
  {
    std::vector<G4String> params = GetWordsInString( newValues );
    if( params.size() !=1 ) 
      G4Exception("DicomIntersectVolume::SetNewValue",
                  "",
                  FatalErrorInArgument,
                  G4String("Command: "+ command->GetCommandPath() + 
                           "/" + command->GetCommandName() + 
                           " " + newValues + 
                           "  needs 1 argument: VOLUME_NAME").c_str()); 

    //----- Build G4VSolid 
    BuildG4Solid(params);

    //----- Calculate volume inverse 3D transform
    G4VPhysicalVolume* pv = GetPhysicalVolumes( params[0], 1, 1)[0];
    
    theVolumeTransform = 
      G4AffineTransform( pv->GetFrameRotation(), pv->GetFrameTranslation() );
  }

  //----- Calculate relative phantom - volume 3D transform 
  G4PhantomParameterisation* thePhantomParam = GetPhantomParam(true);
  
  G4RotationMatrix* rotph = new G4RotationMatrix(); 
  // assumes the phantom mother is not rotated !!! 
  G4AffineTransform thePhantomTransform( rotph, G4ThreeVector() ); 
  // assumes the phantom mother is not translated !!! 

  G4AffineTransform theTransform = theVolumeTransform*thePhantomTransform;

  G4ThreeVector axisX( 1., 0., 0. );
  G4ThreeVector axisY( 0., 1., 0. );
  G4ThreeVector axisZ( 0., 0., 1. );
  theTransform.ApplyAxisTransform(axisX);
  theTransform.ApplyAxisTransform(axisY);
  theTransform.ApplyAxisTransform(axisZ);
  
  //----- Write phantom header
  G4String thePhantomFileName = "phantom.g4pdcm";
  fout.open(thePhantomFileName);
  std::vector<G4Material*> materials = thePhantomParam->GetMaterials();
  fout << materials.size() << G4endl;
  for( unsigned int ii = 0; ii < materials.size(); ++ii )
  {
    fout << ii << " " << materials[ii]->GetName() << G4endl;
  }

  //----- Loop to pantom voxels
  G4int nx = G4int(thePhantomParam->GetNoVoxelX());
  G4int ny = G4int(thePhantomParam->GetNoVoxelY());
  G4int nz = G4int(thePhantomParam->GetNoVoxelZ());
  G4int nxy = nx*ny;
  fVoxelIsInside = new G4bool[nx*ny*nz];
  G4double voxelHalfWidthX = thePhantomParam->GetVoxelHalfX();
  G4double voxelHalfWidthY = thePhantomParam->GetVoxelHalfY();
  G4double voxelHalfWidthZ = thePhantomParam->GetVoxelHalfZ();

  //----- Write phantom dimensions and limits
  fout << nx << " " << ny << " " << nz << G4endl;
  fout << -voxelHalfWidthX*nx+thePhantomTransform.NetTranslation().x() << " "  
       << voxelHalfWidthX*nx+thePhantomTransform.NetTranslation().x() << G4endl;
  fout << -voxelHalfWidthY*ny+thePhantomTransform.NetTranslation().y() << " " 
       << voxelHalfWidthY*ny+thePhantomTransform.NetTranslation().y() << G4endl;
  fout << -voxelHalfWidthZ*nz+thePhantomTransform.NetTranslation().z() << " " 
       << voxelHalfWidthZ*nz+thePhantomTransform.NetTranslation().z() << G4endl;

  for( G4int iz = 0; iz < nz; ++iz ){
    for( G4int iy = 0; iy < ny; ++iy) {

      G4bool bPrevVoxelInside = true;
      G4bool b1VoxelFoundInside = false;
      G4int firstVoxel = -1;
      G4int lastVoxel = -1;
      for(G4int ix = 0; ix < nx; ++ix ){
        G4ThreeVector voxelCentre( (-nx+ix*2+1)*voxelHalfWidthX, 
                   (-ny+iy*2+1)*voxelHalfWidthY, (-nz+iz*2+1)*voxelHalfWidthZ);
        theTransform.ApplyPointTransform(voxelCentre);
        G4bool bVoxelIsInside = true;
        for( G4int ivx = -1; ivx <= 1; ivx+=2 ) {
          for( G4int ivy = -1; ivy <= 1; ivy+=2 ){
            for( G4int ivz = -1; ivz <= 1; ivz+=2 ) {
              G4ThreeVector voxelPoint = voxelCentre 
                + ivx*voxelHalfWidthX*axisX +
                ivy*voxelHalfWidthY*axisY + ivz*voxelHalfWidthZ*axisZ;
              if( fSolid->Inside( voxelPoint ) == kOutside ) {
                bVoxelIsInside = false;
                break;
              } else {
              }
            }
            if( !bVoxelIsInside ) break;
          }
          if( !bVoxelIsInside ) break;
        }

        G4int copyNo = ix + nx*iy + nxy*iz;
        if( bVoxelIsInside ) {
          if( !bPrevVoxelInside ) {
            G4Exception("DicomIntersectVolume::SetNewValue",
                        "",
                        FatalException,
                 "Volume cannot intersect phantom in discontiguous voxels, "
                 "please use other voxel");
          }
          if( !b1VoxelFoundInside ) {
            firstVoxel = ix;
            b1VoxelFoundInside = true;
          }
          lastVoxel = ix;
          fVoxelIsInside[copyNo] = true;
        } else {
          fVoxelIsInside[copyNo] = false;
        }
           
      }
      fout << firstVoxel << " " << lastVoxel << G4endl;
    }
  }

  //-----  Now write the materials
  for( G4int iz = 0; iz < nz; ++iz ){
    for( G4int iy = 0; iy < ny; ++iy) {
      G4bool b1xFound = false; 
      for(G4int ix = 0; ix < nx; ++ix ){
        size_t copyNo = ix + ny*iy + nxy*iz;
        if( fVoxelIsInside[copyNo] ) {
          fout << thePhantomParam->GetMaterialIndex(copyNo)<< " ";
          b1xFound = true;
        }
      }
      if(b1xFound ) fout << G4endl;
    }
  }

  // Now write densities
  for( G4int iz = 0; iz < nz; ++iz ){
    for( G4int iy = 0; iy < ny; ++iy) {
      G4bool b1xFound = false; 
      for(G4int ix = 0; ix < nx; ++ix ){
        size_t copyNo = ix + ny*iy + nxy*iz;
        if( fVoxelIsInside[copyNo] ) {
          fout <<thePhantomParam->GetMaterial(copyNo)->GetDensity()/g*cm3<< " ";
          b1xFound = true;
        }
      }
      if(b1xFound ) fout << G4endl;
    }
  }

} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomIntersectVolume::BuildUserSolid( std::vector<G4String> params )
{
  for( G4int ii = 0; ii < 6; ++ii ) params.erase( params.begin() ); 
  // take otu position and rotation angles
  params.insert( params.begin(), ":SOLID");
  params.insert( params.begin(), params[1] );
  G4tgrSolid* tgrSolid = new G4tgrSolid(params);
  G4tgbVolume* tgbVolume = new G4tgbVolume();  
  fSolid = tgbVolume->FindOrConstructG4Solid( tgrSolid );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomIntersectVolume::BuildG4Solid( std::vector<G4String> params )
{
  fSolid = GetLogicalVolumes( params[0], 1, 1)[0]->GetSolid();
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PhantomParameterisation*
DicomIntersectVolume::GetPhantomParam(G4bool bMustExist)
{
  G4PhantomParameterisation* paramreg = 0;

  G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance();
  for( auto cite = pvs->cbegin(); cite != pvs->cend(); ++cite )
  {
    if( IsPhantomVolume( *cite ) ) {
      const G4PVParameterised* pvparam = 
        static_cast<const G4PVParameterised*>(*cite);
      G4VPVParameterisation* param = pvparam->GetParameterisation();
      //    if( static_cast<const G4PhantomParameterisation*>(param) ){
      //    if( static_cast<const G4PhantomParameterisation*>(param) ){
      //      G4cout << "G4PhantomParameterisation volume found  "
      // << (*cite)->GetName() << G4endl;
      paramreg = static_cast<G4PhantomParameterisation*>(param);
    }
  }
  
  if( !paramreg && bMustExist ) 
    G4Exception("DicomIntersectVolume::GetPhantomParam",
                "",
                FatalErrorInArgument,
                " No G4PhantomParameterisation found ");
  
  return paramreg;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<G4VPhysicalVolume*> DicomIntersectVolume::GetPhysicalVolumes(
                 const G4String& name, G4bool exists, G4int nVols )
{
  std::vector<G4VPhysicalVolume*> vvolu;
  std::string::size_type ial = name.rfind(":");
  G4String volname = "";
  G4int volcopy = 0;
  if( ial != std::string::npos )
  {
    std::string::size_type ial2 = name.rfind(":",ial-1);
    if( ial2 != std::string::npos )
    {
      G4Exception("DicomIntersectVolume::GetPhysicalVolumes",
                  "",
                  FatalErrorInArgument,
                  G4String("Name corresponds to a touchable " + name).c_str());
    }
    else
    { 
      volname = name.substr( 0, ial );
      volcopy = G4UIcommand::ConvertToInt( name.substr(ial+1, name.length()).c_str() );
    }
  }
  else
  {
    volname = name;
    volcopy = -1;
  }

  G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance();
  for( auto citepv = pvs->cbegin(); citepv != pvs->cend(); ++citepv )
  {
    if( volname == (*citepv)->GetName() 
        && ( (*citepv)->GetCopyNo() == volcopy || -1 == volcopy ) )
    {
      vvolu.push_back( *citepv );
    }
  }
  
  //----- Check that at least one volume was found
  if( vvolu.size() == 0 ) {
    if(exists) {
      G4Exception(" DicomIntersectVolume::GetPhysicalVolumes",
               "",
               FatalErrorInArgument,
               G4String("No physical volume found with name " + name).c_str());
    } else {
      G4cerr << "!!WARNING: DicomIntersectVolume::GetPhysicalVolumes: "<<
        " no physical volume found with name " << name << G4endl;
    }
  }

  if( nVols != -1 && G4int(vvolu.size()) != nVols ) {
    G4Exception("DicomIntersectVolume::GetLogicalVolumes:",
                "Wrong number of physical volumes found",
                FatalErrorInArgument,
                ("Number of physical volumes " + 
         G4UIcommand::ConvertToString(G4int(vvolu.size())) + 
         ", requesting " + G4UIcommand::ConvertToString(nVols)).c_str());
  } 

  return vvolu;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool DicomIntersectVolume::IsPhantomVolume( G4VPhysicalVolume* pv )
{
  EAxis axis;
  G4int nReplicas;
  G4double width,offset;
  G4bool consuming;
  pv->GetReplicationData(axis,nReplicas,width,offset,consuming);
  EVolume type = (consuming) ? kReplica : kParameterised;
  if( type == kParameterised && pv->GetRegularStructureId() == 1 ) {  
    return TRUE;
  } else {
    return FALSE;
  }

} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<G4LogicalVolume*> DicomIntersectVolume::GetLogicalVolumes( 
const G4String& name, G4bool exists, G4int nVols )
{
  //  G4cout << "GetLogicalVolumes " << name << " " << exists << G4endl;
  std::vector<G4LogicalVolume*> vvolu;
  G4int ial = G4int(name.rfind(":"));
  if( ial != -1 ) {
    G4Exception("DicomIntersectVolume::GetLogicalVolumes",
                "",
                FatalErrorInArgument,
                G4String("Name corresponds to a touchable or physcal volume" 
                         + name).c_str());
  }

  G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
  for( auto citelv = lvs->cbegin(); citelv != lvs->cend(); ++citelv )
  {
    if( name == (*citelv)->GetName() ) {
      vvolu.push_back( *citelv );
    }
  }
  
  //----- Check that at least one volume was found
  if( vvolu.size() == 0 ) {
    if(exists) {
      G4Exception("DicomIntersectVolume::GetLogicalVolumes:","",
          FatalErrorInArgument,("no logical volume found with name " 
                                + name).c_str());
    } else {
      G4Exception("DicomIntersectVolume::GetLogicalVolumes:","",
          JustWarning,("no  logical volume found with name " + name).c_str());
    }
  }

  if( nVols != -1 && G4int(vvolu.size()) != nVols ) {
    G4Exception("DicomIntersectVolume::GetLogicalVolumes:",
                "Wrong number of logical volumes found",
                FatalErrorInArgument,
                ("Number of logical volumes " +
         G4UIcommand::ConvertToString(G4int(vvolu.size())) + 
         ", requesting " + G4UIcommand::ConvertToString(nVols)).c_str());
  } 

  return vvolu;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::vector<G4String> DicomIntersectVolume::GetWordsInString( 
                                                const G4String& stemp)
{
  std::vector<G4String> wordlist;

  //---------- Read a line of file:
  //----- Clear wordlist
  G4int ii;
  const char* cstr = stemp.c_str();
  G4int siz = G4int(stemp.length());
  G4int istart = 0;
  G4int nQuotes = 0;
  G4bool lastIsBlank = false;
  G4bool lastIsQuote = false;
  for( ii = 0; ii < siz; ++ii )
  {
    if(cstr[ii] == '\"' ){
      if( lastIsQuote ){
        G4Exception("GmGenUtils:GetWordsFromString","",FatalException,
            ("There cannot be two quotes together " + stemp).c_str() );
      }
      if( nQuotes%2 == 1 ){
        //close word 
        wordlist.push_back( stemp.substr(istart,ii-istart) );
        //        G4cout << "GetWordsInString new word0 " 
    //<< wordlist[wordlist.size()-1] << " istart " << istart 
        //<< " ii " << ii << G4endl;
      } else {
        istart = ii+1;
      }
      ++nQuotes;
      lastIsQuote = true;
      lastIsBlank = false;
    } else if(cstr[ii] == ' ' ){
      //            G4cout << "GetWordsInString blank nQuotes " 
      //<< nQuotes << " lastIsBlank " << lastIsBlank << G4endl;
      if( nQuotes%2 == 0 ){
        if( !lastIsBlank && !lastIsQuote ) {
          wordlist.push_back( stemp.substr(istart,ii-istart) );
          //          G4cout << "GetWordsInString new word1 " 
      //<< wordlist[wordlist.size()-1] << " lastIsBlank " 
          //<< lastIsBlank << G4endl;
        }

        istart = ii+1;
        lastIsQuote = false;
        lastIsBlank = true;
      }
    } else {
      if( ii == siz-1 ) {
        wordlist.push_back( stemp.substr(istart,ii-istart+1) );
        //                G4cout << "GetWordsInString new word2 " 
    //<< wordlist[wordlist.size()-1] << " istart " << istart << G4endl;
      }
      lastIsQuote = false;
      lastIsBlank = false;
    }
  }
  if( nQuotes%2 == 1 ) {
    G4Exception("GmGenUtils:GetWordsFromString","",FatalException, 
        ("unbalanced quotes in line " + stemp).c_str() );
  }

  return wordlist;
}
