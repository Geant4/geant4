#include "G4tgbRotationMatrixMgr.hh"

#include "G4tgrRotationMatrixFactory.hh"
#include "G4tgrMessenger.hh"


G4tgbRotationMatrixMgr * G4tgbRotationMatrixMgr::theInstance = 0;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4tgbRotationMatrixMgr* G4tgbRotationMatrixMgr::GetInstance()
{
 if( !theInstance ) {
    theInstance = new G4tgbRotationMatrixMgr;
    theInstance->CopyRotMats();
  }
  return theInstance;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4tgbRotationMatrixMgr::~G4tgbRotationMatrixMgr()
{
  mstgbrotm::const_iterator tgbcite;
  for( tgbcite = theTgbRotMats.begin(); tgbcite != theTgbRotMats.end(); tgbcite++) {
    delete (*tgbcite).second;
  }
  theTgbRotMats.clear();

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void G4tgbRotationMatrixMgr::CopyRotMats()
{
  mstgrrotm tgrRotms = G4tgrRotationMatrixFactory::GetInstance()->GetRotMatMap();
  mstgrrotm::iterator cite;
  for( cite = tgrRotms.begin(); cite != tgrRotms.end(); cite++ ){
    G4tgrRotationMatrix* tgr = (*cite).second;
    G4tgbRotationMatrix* tgb = new G4tgbRotationMatrix( tgr );
    theTgbRotMats[tgb->GetName()] = tgb;
  }

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4RotationMatrix* G4tgbRotationMatrixMgr::FindOrBuildG4RotMatrix(const G4String& name)
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbRotationMatrixMgr::FindOrBuildG4RotMatrix " << name << endl;
#endif
  G4RotationMatrix* g4rotm = FindG4RotMatrix( name );
  if( g4rotm == 0 ) {
    G4tgbRotationMatrix* hrotm = FindOrBuildTgbRotMatrix( name );
    // getRotMatrix never returns 0, becuase if it is not found, it will crash
    g4rotm = hrotm->BuildG4RotMatrix();
    //    G4cout << " G4tgbRotationMatrixMgr::FindOrBuildG4RotMatrix rot mat created " << name << " = " << g4rotm << G4endl;
  }
  return g4rotm;

}        


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4RotationMatrix* G4tgbRotationMatrixMgr::FindG4RotMatrix(const G4String& name)
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbRotationMatrixMgr::FindG4RotMatrix " << name << endl;
#endif
  G4RotationMatrix* g4rotm = 0;

  msg4rotm::const_iterator cite = theG4RotMats.find( name );
  if( cite != theG4RotMats.end() ) {
    g4rotm = (*cite).second;
  } 
//  if( name == "NLL" ) {
//    g4rotm = new G4RotationMatrix();
    //---------- Check if corresponding HElement exists
  
  return g4rotm;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4tgbRotationMatrix* G4tgbRotationMatrixMgr::FindOrBuildTgbRotMatrix(const G4String& name)
{
  G4tgbRotationMatrix* rotm = FindTgbRotMatrix( name );

  if( rotm == 0 ) {
    G4Exception(" G4tgbRotationMatrixFactory::FindOrBuildRotMatrix. Rot Matrix " + name + "  not found ");
  }

  return rotm;

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G4tgbRotationMatrix* G4tgbRotationMatrixMgr::FindTgbRotMatrix(const G4String& name)
{
  G4tgbRotationMatrix* rotm = 0;

  mstgbrotm::const_iterator cite = theTgbRotMats.find( name );
  if( cite != theTgbRotMats.end() ) { 
    rotm = (*cite).second;
  }

  return rotm;

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ostream& operator<<(ostream& os , const G4RotationMatrix & rot)
{
  //  os << "( " << rot.xx() << tab << rot.xy() << tab << rot.xz() << " )" << endl;
  //  os << "( " << rot.yx() << tab << rot.yy() << tab << rot.yz() << " )" << endl;
  //  os << "( " << rot.zx() << tab << rot.zy() << tab << rot.zz() << " )" << endl;

  os << "[ " 
     << rot.thetaX()/deg << '\t' << rot.phiX()/deg << '\t'
     << rot.thetaY()/deg << '\t' << rot.phiY()/deg << '\t'
     << rot.thetaZ()/deg << '\t' << rot.phiZ()/deg << " ]"
     << endl;
  return os;
}

