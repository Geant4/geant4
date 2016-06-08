// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPIsoData.hh"
#include "G4NeutronHPDataUsed.hh"

  G4bool G4NeutronHPIsoData::Init(G4int A, G4int Z, G4double abun, G4String dirName, G4String aFSType)
  {
    theChannelData = NULL;
    G4double abundance = abun/100.;
    G4String filename;
    G4bool result = true;
    G4NeutronHPDataUsed aFile = theNames.GetName(A, Z, dirName, aFSType, result);
    filename = aFile.GetName();
//    if(filename=="") return false;
#ifndef WIN32
    ifstream theChannel(filename);
#else
    ifstream theChannel(filename,ios::in|ios::nocreate);
#endif
    
    if(!theChannel) return false;
    // accommodating deficiencie of some compilers
    if(theChannel.eof()) return false; 
    if(!theChannel) return false;
    G4int count;
    G4int dummy; 
    theChannel >> dummy >> dummy;
    theChannelData = new G4NeutronHPVector;
    G4int nData;
    theChannel >> nData;
    theChannelData->Init(theChannel, nData, eV, abundance*barn);
//    G4cout << "Channel Data Statistics: "<<theChannelData->GetVectorLength()<<endl;
//    G4cout << "Channel data"<<endl;
//     G4int hpw;
//     cin >> hpw;
//    theChannelData->Dump();
    return result;
  }
  
  void G4NeutronHPIsoData::Init(G4int A, G4int Z, G4double abun) //fill PhysicsVector for this Isotope
  {
    G4String dirName;
    if(!getenv("NeutronHPCrossSections")) 
       G4Exception("Please setenv NeutronHPCrossSections to point to the neutron cross-section files.");
    G4String baseName = getenv("NeutronHPCrossSections");
    dirName = baseName+"/Fission";
    if(Z>89) 
    {
      Init(A, Z, abun, dirName, "/CrossSection/");
    }
    else
    {
       theChannelData = new G4NeutronHPVector;
    }
    theFissionData = theChannelData;
    theChannelData = NULL; // fast fix for double delete; revisit later. @@@@@@@
    dirName = baseName+"/Capture";
    Init(A, Z, abun, dirName, "/CrossSection/");
    theCaptureData = theChannelData;
    theChannelData = NULL;
    dirName = baseName+"/Elastic";
    Init(A, Z, abun, dirName, "/CrossSection/");
    theElasticData = theChannelData;
    theChannelData = NULL;
    dirName = baseName+"/Inelastic";
    Init(A, Z, abun, dirName, "/CrossSection/");
    theInelasticData = theChannelData;
    theChannelData = NULL;
    
//    if(theInelasticData!=NULL) G4cout << "Inelastic Data Statistics: "<<theInelasticData->GetVectorLength()<<endl;
//    if(theElasticData!=NULL) G4cout << "Elastic Data Statistics: "<<theElasticData->GetVectorLength()<<endl;
//    if(theCaptureData!=NULL) G4cout << "Capture Data Statistics: "<<theCaptureData->GetVectorLength()<<endl;
//    if(theFissionData!=NULL) G4cout << "Fission Data Statistics: "<<theFissionData->GetVectorLength()<<endl;
//  G4cout << "Inelastic data"<<endl;
//  if(theInelasticData!=NULL) theInelasticData->Dump();
//  G4cout << "Elastic data"<<endl;
//  if(theElasticData!=NULL) theElasticData->Dump();
//  G4cout << "Capture data"<<endl;
//  if(theCaptureData!=NULL) theCaptureData->Dump();
//  G4cout << "Fission data"<<endl;
//  if(theFissionData!=NULL) theFissionData->Dump();

  }
  
  G4String G4NeutronHPIsoData::GetName(G4int A, G4int Z, G4String base, G4String rest)
  {
    G4bool dbool;
    return (theNames.GetName(A, Z, base, rest, dbool)).GetName();
  }
  
