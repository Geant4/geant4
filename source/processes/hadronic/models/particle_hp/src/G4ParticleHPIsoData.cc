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
// particle_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
//080901 Avoiding troubles which caused by G4PhysicsVecotor of length 0 by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPIsoData.hh"
#include "G4ParticleHPManager.hh"
#include "G4ParticleHPDataUsed.hh"
#include "G4Neutron.hh"

  //G4bool G4ParticleHPIsoData::Init(G4int A, G4int Z, G4double abun, G4String dirName, G4String aFSType)
G4bool G4ParticleHPIsoData::Init(G4int A, G4int Z, G4int M, G4double abun, G4String dirName, G4String aFSType)
  {
    theChannelData = 0;
    G4double abundance = abun/100.;
    G4String filename;
    G4bool result = true;
    //G4ParticleHPDataUsed aFile = theNames.GetName(A, Z, dirName, aFSType, result);
    G4ParticleHPDataUsed aFile = theNames.GetName(A, Z, M, dirName, aFSType, result);
    filename = aFile.GetName();
//    if(filename=="") return false;
    //std::ifstream theChannel(filename);
    std::istringstream theChannel(filename,std::ios::in);
    G4ParticleHPManager::GetInstance()->GetDataStream(filename,theChannel);

#ifdef G4PHPDEBUG
    if(std::getenv("G4ParticleHPDebug")) G4cout << "G4ParticleHPIsoData::Init = "<< filename <<" "<< A << " " << Z <<G4endl;
#endif
    
    if(Z==1 && (aFile.GetZ()!=Z || std::abs(aFile.GetA()-A)>0.0001) )
    {
      if(std::getenv("G4ParticleHPDebug")) G4cout << "Skipped = "<< filename <<" "<<A<<" "<<Z<<G4endl;
      //080901 TKDB No more necessary below protection, cross sections set to 0 in G4ParticleHPNames
      //And below two lines causes trouble with G4PhysicsVector
      //theChannel.close();
      //return false;
    }
   if(!theChannel) {/*theChannel.close()*/; return false;}
    // accommodating deficiencie of some compilers
    if(theChannel.eof()) {/*theChannel.close()*/; return false;} 
    if(!theChannel) {/*theChannel.close()*/; return false;}
    G4int dummy; 
    theChannel >> dummy >> dummy;
    theChannelData = new G4ParticleHPVector;
    G4int nData;
    theChannel >> nData;
    theChannelData->Init(theChannel, nData, CLHEP::eV, abundance*CLHEP::barn);
//    G4cout << "Channel Data Statistics: "<<theChannelData->GetVectorLength()<<G4endl;
//    G4cout << "Channel data"<<G4endl;
//     G4int hpw;
//     G4cin >> hpw;
//    theChannelData->Dump();
//    theChannel.close();
    return result;
  }
  
  //void G4ParticleHPIsoData::Init(G4int A, G4int Z, G4double abun) //fill PhysicsVector for this Isotope
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
void G4ParticleHPIsoData::Init(G4int A, G4int Z, G4int M,G4double abun, G4ParticleDefinition* projectile, const char* dataDirVariable ) //fill PhysicsVector for this Isotope
  {

  G4String particleName;
  if ( projectile == G4Neutron::Neutron() ) {
     ;
  } else if ( projectile == G4Proton::Proton() ) {
     particleName = "Proton";
  } else if ( projectile == G4Deuteron::Deuteron() ) {
     particleName = "Deuteron";
  } else if ( projectile == G4Triton::Triton() ) {
     particleName = "Triton";
  } else if ( projectile == G4He3::He3() ) {
     particleName = "He3";
  } else if ( projectile == G4Alpha::Alpha() ) {
     particleName = "Alpha";
  } else {
     G4String message("G4ParticleHPInelastic may only be called for neutron, proton, deuteron, triton, He3 or alpha, while it is called for " + projectile->GetParticleName());
     throw G4HadronicException(__FILE__, __LINE__,message.c_str());
  }

  G4String baseName;
  if ( G4FindDataDir( dataDirVariable ) ) {
     baseName = G4FindDataDir( dataDirVariable );
  } else {
     baseName = G4FindDataDir( "G4PARTICLEHPDATA" );
     baseName += "/" + particleName;
  }
    
   // G4String baseName = getenv(dataDirVariable);
    G4String dirName;
    if( projectile == G4Neutron::Neutron() ){
      dirName = baseName+"/Fission";
    //if(Z>89) 
      if(Z>87) //TK Modifed for ENDF VII.0 
	{
	  //Init(A, Z, abun, dirName, "/CrossSection/");
	  Init(A, Z, M, abun, dirName, "/CrossSection");
	}
      else
	{
	  theChannelData = new G4ParticleHPVector;
	}
      theFissionData = theChannelData;
      theChannelData = 0; // fast fix for double delete; revisit later. @@@@@@@

      dirName = baseName+"/Capture";
      //Init(A, Z, abun, dirName, "/CrossSection/");
      Init(A, Z, M, abun, dirName, "/CrossSection");
      theCaptureData = theChannelData;
      theChannelData = 0;

      dirName = baseName+"/Elastic";
      //Init(A, Z, abun, dirName, "/CrossSection/");
      Init(A, Z, M, abun, dirName, "/CrossSection");
      theElasticData = theChannelData;
      theChannelData = 0;
    }

    dirName = baseName+"/Inelastic";
    //Init(A, Z, abun, dirName, "/CrossSection/");
    Init(A, Z, M, abun, dirName, "/CrossSection");
    theInelasticData = theChannelData;
    theChannelData = 0;
    
//    if(theInelasticData!=0) G4cout << "Inelastic Data Statistics: "<<theInelasticData->GetVectorLength()<<G4endl;
//    if(theElasticData!=0) G4cout << "Elastic Data Statistics: "<<theElasticData->GetVectorLength()<<G4endl;
//    if(theCaptureData!=0) G4cout << "Capture Data Statistics: "<<theCaptureData->GetVectorLength()<<G4endl;
//    if(theFissionData!=0) G4cout << "Fission Data Statistics: "<<theFissionData->GetVectorLength()<<G4endl;
//  G4cout << "Inelastic data"<<G4endl;
//  if(theInelasticData!=0) theInelasticData->Dump();
//  G4cout << "Elastic data"<<G4endl;
//  if(theElasticData!=0) theElasticData->Dump();
//  G4cout << "Capture data"<<G4endl;
//  if(theCaptureData!=0) theCaptureData->Dump();
//  G4cout << "Fission data"<<G4endl;
//  if(theFissionData!=0) theFissionData->Dump();

  }
  
  G4String G4ParticleHPIsoData::GetName(G4int A, G4int Z, G4String base, G4String rest)
  {
    G4bool dbool;
    return (theNames.GetName(A, Z, base, rest, dbool)).GetName();
  }

