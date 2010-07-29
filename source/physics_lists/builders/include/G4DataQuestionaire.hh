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
#ifndef G4DataQuestionaire_h
#define G4DataQuestionaire_h 1

#include "globals.hh"

#include "G4HadronElasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4VPiKBuilder.hh"

enum G4DataType {no, photon, neutron, radioactive, lowenergy, optical, neutronxs};
class G4DataQuestionaire
{
  public: 
    G4DataQuestionaire(G4DataType t1=no, G4DataType t2=no, G4DataType t3=no, 
		       G4DataType t4=no, G4DataType t5=no, G4DataType t6=no) 
    {
      G4cout << G4endl;
      G4cout << "<<< Geant4 Physics List engine packaging library: PACK 5.5"<<G4endl;
      //      G4cout <<G4endl<<G4endl;
      // G4cout << "##### the input "<<t1<<" "<<t2<<" "<<t3<<" "<<t4<<G4endl;
      for(G4int i=0; i<6; ++i)
      {
	G4DataType t(no);
	if(i==0) t=t1;
	if(i==1) t=t2;
	if(i==2) t=t3;
	if(i==3) t=t4;
	if(i==4) t=t5;
	if(i==5) t=t6;
	switch(t)
	{
          case photon:
	    if(!getenv("G4LEVELGAMMADATA") )
	    {
	      G4cout << "Photon-evaporation data are needed."<<G4endl;
	      G4cout << "Please set the environmental variable G4LEVELGAMMADATA"<<G4endl;
	      G4cout << "to point to your PhotonEvaporation directory."<<G4endl;
	      G4cout << "Data are available from the Geant4 download page."<<G4endl;
              G4Exception("G4DataQuestionaire", "007", FatalException,
	                  "Fatal error: Missing mandatory data for this simulation engine");
	    }
	    break;
          case neutron:
	    if(!getenv("G4NEUTRONHPDATA") )
	    {
	      G4cout << "G4NDL are needed."<<G4endl;
	      G4cout << "Please set the environmental variable G4NEUTRONHPDATA"<<G4endl;
	      G4cout << "to point to your G4NDL directory."<<G4endl;
	      G4cout << "Data are available from the Geant4 download page."<<G4endl;
              G4Exception("G4DataQuestionaire", "007", FatalException,
	                  "Fatal error: Missing mandatory data for this simulation engine");
	    }
	    break;
          case radioactive:
	    if(!getenv("G4RADIOACTIVEDATA") )
	    {
	      G4cout << "Radioactive decay data are needed."<<G4endl;
	      G4cout << "Please set the environmental variable G4RADIOACTIVEDATA"<<G4endl;
	      G4cout << "to point to your RadiativeDecay directory."<<G4endl;
	      G4cout << "Data are available from the Geant4 download page."<<G4endl;
              G4Exception("G4DataQuestionaire", "007", FatalException,
	                  "Fatal error: Missing mandatory data for this simulation engine");
	    }
	    break;
          case lowenergy:
	    if(!getenv("G4LEDATA") )
	    {
	      G4cout << "Low energy electromagnetic data are needed."<<G4endl;
	      G4cout << "Please set the environmental variable G4LEDATA"<<G4endl;
	      G4cout << "to point to your G4EMLOW directory."<<G4endl;
	      G4cout << "Data are available from the Geant4 download page."<<G4endl;
              G4Exception("G4DataQuestionaire", "007", FatalException,
	                  "Fatal error: Missing mandatory data for this simulation engine");
	    }
	    break;
          case optical:
	    /*
	    if(!getenv("G4REALSURFACEDATA") )
	    {
	      G4cout << "Data describing surface propeties for optical photons are needed."<<G4endl;
	      G4cout << "Please set the environmental variable G4REALSURFACEDATA"<<G4endl;
	      G4cout << "to point to your RealSurface directory."<<G4endl;
	      G4cout << "Data are available from the Geant4 download page."<<G4endl;
              G4Exception("G4DataQuestionaire", "007", FatalException,
	                  "Fatal error: Missing mandatory data for this simulation engine");
	    }
	    */
	    break;
          case neutronxs:
	    if(!getenv("G4NEUTRONXSDATA") )
	    {
	      G4cout << "G4NEUTRONXS are needed."<<G4endl;
	      G4cout << "Please set the environmental variable G4NEUTRONXSDATA"<<G4endl;
	      G4cout << "to point to your G4NEUTRONXS directory."<<G4endl;
	      G4cout << "Data are available from the Geant4 download page."<<G4endl;
              G4Exception("G4DataQuestionaire", "007", FatalException,
	                  "Fatal error: Missing mandatory data for this simulation engine");
	    }
	    break;
	  case no:
	    // all ok
            break;
	  default:
            if(t!=no) 
	    {
              G4Exception("G4DataQuestionaire", "007", FatalException,
	                  "data type requested is not known to the system");
	    }
	}
     }
   }
    ~G4DataQuestionaire() {}
};

// 2002 by J.P. Wellisch

#endif

