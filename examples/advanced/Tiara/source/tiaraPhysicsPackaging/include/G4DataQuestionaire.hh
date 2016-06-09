//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef G4DataQuestionaire_h
#define G4DataQuestionaire_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VPiKBuilder.hh"

#include <vector>

enum G4DataType {no, photon, neutron, radioactive, lowenergy};
class G4DataQuestionaire
{
  public: 
    G4DataQuestionaire(G4DataType t1=no, G4DataType t2=no, G4DataType t3=no, G4DataType t4=no) 
    {
      G4cout <<G4endl<<G4endl;
      G4cout << "You are using the simulation engine reuse library: PACK 2.1"<<G4endl;
      G4cout <<G4endl<<G4endl;
      for(G4int i=0; i<4; i++)
      {
	G4DataType t(no);
	if(i==0) t=t1;
	if(i==1) t=t2;
	if(i==2) t=t3;
	if(i==3) t=t4;
	switch(t)
	{
          case photon:
	    if(!getenv("G4LEVELGAMMADATA") )
	    {
	      G4cout << "Photon-evaporation data are needed."<<G4endl;
	      G4cout << "Please set the environmental variable G4LEVELGAMMADATA"<<G4endl;
	      G4cout << "to point to your PhotonEvaporation directory."<<G4endl;
	      G4cout << "Data are available from the geant4 download page."<<G4endl;
	      G4Exception("Fatal error: Missing mandatory data for this simulation engine");
	    }
	    break;
          case neutron:
	    if(!getenv("NeutronHPCrossSections") )
	    {
	      G4cout << "G4NDL are needed."<<G4endl;
	      G4cout << "Please set the environmental variable NeutronHPCrossSections"<<G4endl;
	      G4cout << "to point to your G4NDL directory."<<G4endl;
	      G4cout << "Data are available from http://cmsdoc.cern.ch/~hpw/G4NDL3.7.tar.gz "<<G4endl;
	      G4cout << "of the geant4 download page."<<G4endl;
	      G4Exception("Fatal error: Missing mandatory data for this simulation engine");
	    }
	    break;
          case radioactive:
	    if(!getenv("G4RADIOACTIVEDATA") )
	    {
	      G4cout << "Radioactive decay data are needed."<<G4endl;
	      G4cout << "Please set the environmental variable G4RADIOACTIVEDATA"<<G4endl;
	      G4cout << "to point to your RadiativeDecay directory."<<G4endl;
	      G4cout << "Data are available from the geant4 download page."<<G4endl;
	      G4Exception("Fatal error: Missing mandatory data for this simulation engine");
	    }
	    break;
          case lowenergy:
	    if(!getenv("G4LEDATA") )
	    {
	      G4cout << "Low energy electromagnetic data are needed."<<G4endl;
	      G4cout << "Please set the environmental variable G4LEDATA"<<G4endl;
	      G4cout << "to point to your G4EMLOW directory."<<G4endl;
	      G4cout << "Data are available from the geant4 download page."<<G4endl;
	      G4Exception("Fatal error: Missing mandatory data for this simulation engine");
	    }
	    break;
	  case no:
	    // all ok
          default:
            if(t1!=no) G4Exception("data type requested is not known to the system");
	}
     }
   }
    ~G4DataQuestionaire() {}
};

// 2002 by J.P. Wellisch

#endif

