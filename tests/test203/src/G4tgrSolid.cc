#include "G4tgrSolid.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"
#include "G4tgrVolumeMgr.hh"
#include <set>


//-------------------------------------------------------------
G4tgrSolid::G4tgrSolid( const vector<G4String>& wl) 
{
  //---------- set name 
  theName = G4tgrUtils::SubQuotes( wl[1] ); 
  //  G4cout << " G4tgrSolid::G4tgrSolid " << GetName() << G4endl;

  //---------- set solid type
  theType = G4tgrUtils::SubQuotes( wl[2] ); 

  //----- if there are not paremeters 
  uint ii;

  //---------- create only vector<double> of theSolidParams
  FillSolidParams( wl );

  G4tgrVolumeMgr::GetInstance()->RegisterMe( this );

}


//-------------------------------------------------------------
void G4tgrSolid::FillSolidParams( const std::vector<G4String>& wl )
{
  //---- Setting which are angle parameters (for dimensions...)
  std::map< G4String, std::set<G4int> > angleParams;
  std::set<G4int> apar;
  apar.clear(); apar.insert(3);apar.insert(4);
  angleParams["TUBS"] = apar;
  apar.clear(); apar.insert(5);apar.insert(6);
  angleParams["CONS"] = apar;
  apar.clear(); apar.insert(3);apar.insert(4);apar.insert(5);
  angleParams["PARA"] = apar;
  apar.clear(); apar.insert(1);apar.insert(2);apar.insert(6);apar.insert(10);
  angleParams["TRAP"] = apar;
  apar.clear(); apar.insert(2);apar.insert(3);apar.insert(4);apar.insert(5);
  angleParams["SPHERE"] = apar;
  apar.clear(); apar.insert(2);apar.insert(3);apar.insert(4);apar.insert(5);
  angleParams["TORUS"] = apar;
  apar.clear(); apar.insert(0);apar.insert(1);
  angleParams["POLYCONE"] = apar;
  apar.clear(); apar.insert(0);apar.insert(1);
  angleParams["POLYHEDRA"] = apar;
  apar.clear(); apar.insert(2);apar.insert(3);
  angleParams["HYPE"] = apar;
  apar.clear(); apar.insert(0);
  angleParams["TWISTED_BOX"] = apar;
  apar.clear(); apar.insert(0);apar.insert(2);apar.insert(3);apar.insert(10);
  angleParams["TWISTED_TRAP"] = apar;
  apar.clear(); apar.insert(5);
  angleParams["TWISTED_TRD"] = apar;
  apar.clear(); apar.insert(0);apar.insert(4);
  angleParams["TWISTED_TUBS"] = apar;

  vector<double>* vd = new vector<double>;
  theSolidParams.push_back( vd );
  uint noParRead = wl.size()-3;

  G4String solidType = wl[2];
  //--- Default unit (mm) if length, deg if angle
  for(uint ii = 0; ii < noParRead; ii++) {
    G4bool isAngle = 0;
    std::map< G4String, std::set<G4int> >::iterator ite = angleParams.find(solidType);
    if( ite != angleParams.end() ){
      std::set<G4int> apar = (*ite).second;
      if( apar.find(ii) != apar.end() ) {
	isAngle = 1;
	vd->push_back( G4tgrUtils::GetFloat( wl[3+ii], deg ));
#ifdef G4VERBOSE
	if( G4tgrMessenger::GetVerboseLevel() >= 3 ) 
	  G4cout << " G4tgrSolid::FillSolidParams angle param found " << solidType << " " << ii << G4endl;
#endif
      }
    } 
    //-    G4cout << ii<< " G4tgrSolid::FillSolidParam " << wl.size() << G4endl;
    if(!isAngle){
      vd->push_back( G4tgrUtils::GetFloat( wl[3+ii] ) );
    }
  }
}

