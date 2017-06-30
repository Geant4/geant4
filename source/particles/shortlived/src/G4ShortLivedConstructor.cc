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
// $Id: G4ShortLivedConstructor.cc 102905 2017-03-02 09:50:56Z gcosmo $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//      Add "rho0"                         25  Feb. 2000 H.Kurashige
//      Fix spin/isospin number for quarks 06  Apr. 2001 H.Kurashige
//      update quark mass                  11  Oct. 2006 H.Kurashige
//      update meson/baryon masses         11  Oct. 2006 H.Kurashige

#include "G4ShortLivedConstructor.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"

G4bool G4ShortLivedConstructor::isConstructed = false;

G4ShortLivedConstructor::G4ShortLivedConstructor()
{
}

G4ShortLivedConstructor::~G4ShortLivedConstructor()
{
}


void G4ShortLivedConstructor::ConstructParticle()
{
  if (!isConstructed){
    ConstructQuarks();
    ConstructResonances();
    isConstructed = true;
  }
}

#include "G4DiQuarks.hh"
#include "G4Quarks.hh"
#include "G4Gluons.hh"
void G4ShortLivedConstructor::ConstructQuarks()
{
  G4ParticleDefinition* particle;

  //    Construct Quraks/Gluons as dynamic object
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table 

  // gluon
  particle = new G4Gluons(            
	      "gluon",         0.0*MeV,       0.0*MeV,     0.0*eplus, 
                    2,              -1,             0,          
                    0,               0,             0,             
             "gluons",               0,             0,          21,
		 true,            -1.0,          NULL);
  particle->SetAntiPDGEncoding(21);
  // u-quark
  particle = new G4Quarks(            
	    "u_quark",         2.2*MeV,       0.0*keV,   2./3.*eplus, 
                    1,              +1,             0,          
                    1,              +1,             0,             
             "quarks",               0,             0,           2,
		 true,            -1.0,          NULL);
  // d-quark
  particle = new G4Quarks(            
	    "d_quark",         4.7*MeV,       0.0*keV,  -1./3.*eplus, 
                    1,              +1,             0,          
                    1,              -1,             0,             
             "quarks",               0,             0,           1,
		 true,            -1.0,          NULL);
  // s-quark
  particle = new G4Quarks(            
	    "s_quark",        96.0*MeV,       0.0*keV,  -1./3.*eplus, 
                    1,              +1,             0,          
                    0,               0,             0,             
             "quarks",               0,             0,           3,
		 true,            -1.0,          NULL);
  // c-quark
  particle = new G4Quarks(            
	    "c_quark",        1.27*GeV,       0.0*keV,  +2./3.*eplus, 
                    1,              +1,             0,          
                    0,               0,             0,             
             "quarks",               0,             0,           4,
		 true,            -1.0,          NULL);
  // b-quark
  particle = new G4Quarks(            
	    "b_quark",        4.18*GeV,       0.0*keV,  -1./3.*eplus, 
                    1,              +1,             0,          
                    0,               0,             0,             
             "quarks",               0,             0,           5,
		 true,            -1.0,          NULL);
  // t-quark
  particle = new G4Quarks(            
	    "t_quark",      173.21*GeV,      1.41*GeV,  +2./3.*eplus, 
                    1,              +1,             0,          
                    0,               0,             0,             
             "quarks",               0,             0,           6,
		 true,            -1.0,          NULL);
  // anti u-quark
  particle = new G4Quarks(            
       "anti_u_quark",         2.2*MeV,       0.0*keV,   -2./3.*eplus, 
                    1,              +1,             0,          
                    1,              -1,             0,             
             "quarks",               0,             0,          -2,
		 true,            -1.0,          NULL);
  // anti d-quark
  particle = new G4Quarks(            
       "anti_d_quark",         4.7*MeV,       0.0*keV,   1./3.*eplus, 
                    1,              +1,             0,          
                    1,              +1,             0,             
             "quarks",               0,             0,          -1,
		 true,            -1.0,          NULL);
  // anti s-quark
  particle = new G4Quarks(            
       "anti_s_quark",        96.0*MeV,       0.0*keV,   1./3.*eplus, 
                    1,              +1,             0,          
                    0,               0,             0,             
             "quarks",               0,             0,          -3,
		 true,            -1.0,          NULL);
  // anti c-quark
  particle = new G4Quarks(            
       "anti_c_quark",        1.27*GeV,       0.0*keV,  -2./3.*eplus, 
                    1,              +1,             0,          
                    0,               0,             0,             
             "quarks",               0,             0,          -4,
		 true,            -1.0,          NULL);
  // anti b-quark
  particle = new G4Quarks(            
       "anti_b_quark",        4.18*GeV,       0.0*keV,   1./3.*eplus, 
                    1,              +1,             0,          
                    0,               0,             0,             
             "quarks",               0,             0,          -5,
		 true,            -1.0,          NULL);
  // anti t-quark
  particle = new G4Quarks(            
       "anti_t_quark",      173.21*GeV,      1.41*GeV,  -2./3.*eplus, 
                    1,              +1,             0,          
                    0,               0,             0,             
             "quarks",               0,             0,          -6,
		 true,            -1.0,          NULL);

   // uu1-Diquark
  particle = new G4DiQuarks(            
	"uu1_diquark",         4.6*MeV,       0.0*MeV,   4./3.*eplus, 
                    2,              +1,             0,          
                    2,              +2,             0,             
           "diquarks",               0,             0,        2203,
		 true,            -1.0,          NULL);
    // ud1-Diquark
  particle = new G4DiQuarks(            
	 "ud1_diquark",        7.0*MeV,       0.0*MeV,   1./3.*eplus, 
                    2,              +1,             0,          
                    2,              +0,             0,             
           "diquarks",               0,             0,         2103,
		 true,            -1.0,          NULL);
    // dd1-Diquark
  particle = new G4DiQuarks(            
	"dd1_diquark",         9.6*MeV,       0.0*MeV,   -2./3.*eplus, 
                    2,              +1,             0,          
                    2,              -2,             0,             
           "diquarks",               0,             0,         1103,
		 true,            -1.0,          NULL);
   
   // ud0-Diquark
  particle = new G4DiQuarks(            
	"ud0_diquark",         7.1*MeV,       0.0*MeV,   1./3.*eplus, 
                    0,              +1,             0,          
                    0,              +0,             0,             
           "diquarks",               0,             0,         2101,
		 true,            -1.0,          NULL);
   
   // sd1-Diquark
  particle = new G4DiQuarks(            
	"sd1_diquark",       102.8*MeV,       0.0*MeV,   -2./3.*eplus, 
                    2,              +1,             0,          
                    1,              -1,             0,             
           "diquarks",               0,             0,         3103,
		 true,            -1.0,          NULL);
   
  // su1-Diquark
  particle = new G4DiQuarks(            
	"su1_diquark",       101.4*MeV,       0.0*MeV,   1./3.*eplus, 
                    2,              +1,             0,          
                    1,              +1,             0,             
           "diquarks",               0,             0,         3203,
		 true,            -1.0,          NULL);

    // sd0-Diquark
  particle = new G4DiQuarks(            
	"sd0_diquark",       102.0*MeV,       0.0*MeV,   -2./3.*eplus, 
                    0,              +1,             0,          
                    1,              -1,             0,             
           "diquarks",               0,             0,         3101,
		 true,            -1.0,          NULL);
   
  // su0-Diquark
  particle = new G4DiQuarks(            
	"su0_diquark",       101.4*MeV,       0.0*MeV,   1./3.*eplus, 
                    0,              +1,             0,          
                    1,              +1,             0,             
           "diquarks",               0,             0,         3201,
		 true,            -1.0,          NULL);

   // anti uu1-Diquark
  particle = new G4DiQuarks(            
   "anti_uu1_diquark",         4.6*MeV,       0.0*MeV,  -4./3.*eplus, 
                    2,              +1,             0,          
                    2,              -2,             0,             
           "diquarks",               0,             0,       -2203,
		 true,            -1.0,          NULL);
    // anti ud1-Diquark
  particle = new G4DiQuarks(            
   "anti_ud1_diquark",         7.0*MeV,       0.0*MeV,  -1./3.*eplus, 
                    2,              +1,             0,          
                    2,              +0,             0,             
           "diquarks",               0,             0,        -2103,
		 true,            -1.0,          NULL);
    // anti dd1-Diquark
  particle = new G4DiQuarks(            
   "anti_dd1_diquark",         9.6*MeV,       0.0*MeV,    2./3.*eplus, 
                    2,              +1,             0,          
                    2,              +2,             0,             
           "diquarks",               0,             0,        -1103,
		 true,            -1.0,          NULL);
   
   // anti ud0-Diquark
  particle = new G4DiQuarks(            
   "anti_ud0_diquark",         7.1*MeV,       0.0*MeV,  -1./3.*eplus, 
                    0,              +1,             0,          
                    0,              +0,             0,             
           "diquarks",               0,             0,        -2101,
		 true,            -1.0,          NULL);
   
   // anti  sd1-Diquark
  particle = new G4DiQuarks(            
   "anti_sd1_diquark",       102.8*MeV,       0.0*MeV,    2./3.*eplus, 
                    2,              +1,             0,          
                    1,              +1,             0,             
           "diquarks",               0,             0,        -3103,
		 true,            -1.0,          NULL);
   
  // anti su1-Diquark
  particle = new G4DiQuarks(            
   "anti_su1_diquark",       101.4*MeV,       0.0*MeV,  -1./3.*eplus, 
                    2,              +1,             0,          
                    1,              -1,             0,             
           "diquarks",               0,             0,        -3203,
		 true,            -1.0,          NULL);

    // anti sd0-Diquark
  particle = new G4DiQuarks(            
   "anti_sd0_diquark",       102.0*MeV,       0.0*MeV,    2./3.*eplus, 
                    0,              +1,             0,          
                    1,              +1,             0,             
           "diquarks",               0,             0,        -3101,
		 true,            -1.0,          NULL);
   
  // anti su0-Diquark
  particle = new G4DiQuarks(            
   "anti_su0_diquark",       101.4*MeV,       0.0*MeV,  -1./3.*eplus, 
                    0,              +1,             0,          
                    1,              -1,             0,             
           "diquarks",               0,             0,        -3201,
		 true,            -1.0,          NULL);
    // ss1-Diquark
  particle = new G4DiQuarks(            
	"ss1_diquark",       198.0*MeV,       0.0*MeV,   -2./3.*eplus, 
                    2,              +1,             0,          
                    0,               0,             0,             
           "diquarks",               0,             0,         3303,
		 true,            -1.0,          NULL);
   
    // anti ss1-Diquark
  particle = new G4DiQuarks(            
	"anti_ss1_diquark",  198.0*MeV,       0.0*MeV,    2./3.*eplus, 
                    2,              +1,             0,          
                    0,               0,             0,             
           "diquarks",               0,             0,        -3303,
		 true,            -1.0,          NULL);

  particle = NULL;      
}

#include "G4ExcitedNucleonConstructor.hh"
#include "G4ExcitedDeltaConstructor.hh"
#include "G4ExcitedLambdaConstructor.hh"
#include "G4ExcitedSigmaConstructor.hh"
#include "G4ExcitedXiConstructor.hh"
#include "G4ExcitedMesonConstructor.hh"
void G4ShortLivedConstructor::ConstructResonances()
{
  ConstructBaryons();
  ConstructMesons();

  // N*
  G4ExcitedNucleonConstructor nucleons;
  nucleons.Construct();

  // Delta*
  G4ExcitedDeltaConstructor deltas;
  deltas.Construct();

  // Lambda*
  G4ExcitedLambdaConstructor lamdas;
  lamdas.Construct();

  // Sigma*
  G4ExcitedSigmaConstructor sigmas;
  sigmas.Construct();

  // Xi*
  G4ExcitedXiConstructor xis;
  xis.Construct();
 
  // Mesons
  G4ExcitedMesonConstructor mesons;
  mesons.Construct();

}


#include "G4ExcitedBaryons.hh"
void G4ShortLivedConstructor::ConstructBaryons()
{
  G4DecayTable*   decayTable;
  G4VDecayChannel* mode;
  G4ExcitedBaryons* particle;

  //    Construct Resonace particles as dynamic object
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table 

  // delta baryons
  //  delta(1232)++
  particle = new G4ExcitedBaryons(            
	    "delta++",       1.232*GeV,     120.0*MeV,    +2.0*eplus, 
                    3,              +1,             0,          
                    3,              +3,             0,             
             "baryon",               0,            +1,          2224,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("delta");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of delta++ -> proton + pi+
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("delta++",1.000, 2,
				                        "proton","pi+");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  //  delta(1232)+
  particle = new G4ExcitedBaryons(            
	     "delta+",       1.232*GeV,     120.0*MeV,    +1.0*eplus, 
                    3,              +1,             0,          
                    3,              +1,             0,             
             "baryon",               0,            +1,          2214,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("delta(1232)");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of delta+  -> proton + Gamma
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("delta+", 0.01, 2,
				                        "proton","gamma");
  decayTable->Insert(mode);
  // create decay channel of delta+  -> neutron + pi+
  //                                   parent    BR     #daughters
  // create decay channel of delta+  -> proton + pi0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("delta+", 0.495, 2,
				                        "proton","pi0");
  decayTable->Insert(mode);
  // create decay channel of delta+  -> neutron + pi+
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("delta+", 0.495, 2,
				                        "neutron","pi+");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  //  delta(1232)0
  particle = new G4ExcitedBaryons(            
	     "delta0",       1.232*GeV,     120.0*MeV,    +0.0*eplus, 
                    3,              +1,             0,          
                    3,              -1,             0,             
             "baryon",               0,            +1,          2114,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("delta(1232)");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of delta+  -> neutron + gamma
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("delta0", 0.01, 2,
				                        "neutron","gamma");
  decayTable->Insert(mode);
  // create decay channel of delta+  -> proton + pi-
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("delta0", 0.495, 2,
				                        "proton","pi-");
  decayTable->Insert(mode);
  // create decay channel of delta+  -> neutron + pi0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("delta0", 0.495, 2,
				                        "neutron","pi0");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  //  delta(1232)-
  particle = new G4ExcitedBaryons(            
	     "delta-",       1.232*GeV,     117.0*MeV,    -1.0*eplus, 
                    3,              +1,             0,          
                    3,              -3,             0,             
             "baryon",               0,            +1,          1114,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("delta(1232)");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of delta+  -> neutron + pi-
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("delta-", 1.000, 2,
				                        "neutron","pi-");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);


  ////////////////////////////
  // anti_delta baryons
  //  anti_delta(1232)++
  particle = new G4ExcitedBaryons(            
       "anti_delta++",       1.232*GeV,     120.0*MeV,    -2.0*eplus, 
                    3,              +1,             0,          
                    3,              -3,             0,             
             "baryon",               0,            -1,         -2224,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("delta(1232)");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of delta++ -> anti_proton + pi-
  //                                        parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("anti_delta++",1.000, 2,
				                   "anti_proton","pi-");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  //  anti_delta(1232)+
  particle = new G4ExcitedBaryons(            
	"anti_delta+",       1.232*GeV,     120.0*MeV,    -1.0*eplus, 
                    3,              +1,             0,          
                    3,              -1,             0,             
             "baryon",               0,            -1,         -2214,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("delta(1232)");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of anti_delta+  -> anti_proton + pi0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("anti_delta+", 0.500, 2,
				                   "anti_proton","pi0");
  decayTable->Insert(mode);
  // create decay channel of anti_delta+  -> anti_neutron + pi-
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("anti_delta+", 0.500, 2,
				                  "anti_neutron","pi-");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  //  anti_delta(1232)0
  particle = new G4ExcitedBaryons(            
	"anti_delta0",       1.232*GeV,     120.0*MeV,    +0.0*eplus, 
                    3,              +1,             0,          
                    3,              +1,             0,             
             "baryon",               0,            -1,         -2114,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("delta(1232)");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of anti_delta+  -> anti_proton + pi+
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("anti_delta0", 0.500, 2,
				                        "anti_proton","pi+");
  decayTable->Insert(mode);
  // create decay channel of delta+  -> neutron + pi0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("anti_delta0", 0.500, 2,
				                        "anti_neutron","pi0");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  //  anti_delta(1232)-
  particle = new G4ExcitedBaryons(            
	"anti_delta-",       1.232*GeV,     117.0*MeV,    +1.0*eplus, 
                    3,              +1,             0,          
                    3,              +3,             0,             
             "baryon",               0,            -1,         -1114,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("delta(1232)");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of delta-  -> neutron + pi+
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("anti_delta-", 1.000, 2,
				                        "anti_neutron","pi+");
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);



}
#include "G4ExcitedMesons.hh"
void G4ShortLivedConstructor::ConstructMesons()
{
  G4DecayTable*   decayTable;
  G4VDecayChannel* mode;
  G4ExcitedMesons* particle;

  //    Construct Resonace particles as dynamic object
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table 

  // vector mesons
  // omega
  particle = new G4ExcitedMesons(            
	      "omega",      782.65*MeV,      8.49*MeV,    +0.0*eplus, 
                    2,              -1,            -1,          
                    0,              +0,            -1,             
              "meson",               0,             0,           223,
		false,             0.0,          NULL);
  particle->SetAntiPDGEncoding(223);
  // set sub type
  particle->SetMultipletName("omega");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of omega -> pi+ + pi- + pi0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("omega",0.891, 3,
				                    "pi+","pi-","pi0");
  // add decay table
  decayTable->Insert(mode);

  // create decay channel of omega -> gamma + pi0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("omega",0.0890, 2,
				                    "gamma","pi0");
  // add decay table
  decayTable->Insert(mode);

  // create decay channel of omega -> pi+ + pi- 
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("omega",0.0170, 2,
				                    "pi+","pi-");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  // phi
  particle = new G4ExcitedMesons(            
	        "phi",     1019.46*MeV,     4.266*MeV,    +0.0*eplus, 
                    2,              -1,            -1,          
                    0,              +0,            -1,             
              "meson",               0,             0,           333,
		false,             0.0,          NULL);
  particle->SetAntiPDGEncoding(333);
  // set sub type
  particle->SetMultipletName("phi");
   // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of phi -> kaon+ + kaon-
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("phi",0.492, 2,
				                    "kaon+","kaon-");
  decayTable->Insert(mode);
   // create decay channel of phi -> kaon0S + kaon0L
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("phi",0.340, 2,
				                   "kaon0S","kaon0L");
  // add decay table
  decayTable->Insert(mode);
  // create decay channel of phi -> rho0 + pi0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("phi",0.153, 2,
				                   "rho0","pi0");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);
 
  // rho+
  particle = new G4ExcitedMesons(            
	       "rho+",       775.8*MeV,     150.3*MeV,    +1.0*eplus, 
                    2,              -1,            -1,          
                    2,              +2,            +1,             
              "meson",               0,             0,           213,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("rho");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of rho+ -> pi+ + pi0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("rho+",1.000, 2,
				                    "pi+","pi0");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  // rho-
  particle = new G4ExcitedMesons(            
	       "rho-",       775.8*MeV,     150.3*MeV,    -1.0*eplus, 
                    2,              -1,            -1,          
                    2,              -2,            +1,             
              "meson",               0,             0,          -213,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("rho");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of rho- -> pi- + pi0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("rho-",1.000, 2,
				                    "pi-","pi0");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);
  
  // rho0
  particle = new G4ExcitedMesons(            
               "rho0",       775.26*MeV,    149.1*MeV,         0.0, 
                    2,              -1,            -1,          
                    2,               0,            +1,             
              "meson",               0,             0,         113,
                false,          0.0*ns,          NULL );
  particle->SetAntiPDGEncoding(113);
  // set sub type
  particle->SetMultipletName("rho");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of rho0 -> pi+ + pi-
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("rho0",1.000, 2,
				                    "pi+","pi-");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  // a0(980)+
  particle = new G4ExcitedMesons(            
	   "a0(980)+",       980.0*MeV,      60.0*MeV,    +1.0*eplus, 
                    0,              +1,            +1,          
                    2,              +2,            -1,             
              "meson",               0,             0,       9000211,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("a0(980)");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of a0(980)+ -> eta + pi+
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("a0(980)+",1.000, 2,
				                    "pi+","eta");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  // a0(980)-
  particle = new G4ExcitedMesons(            
	   "a0(980)-",       980.0*MeV,      60.0*MeV,    -1.0*eplus, 
                    0,              +1,            +1,          
                    2,              -2,            -1,             
              "meson",               0,             0,      -9000211,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("a0(980)");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of a0(980)- -> eta + pi-
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("a0(980)-",1.000, 2,
				                    "pi-","eta");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  // a0(980)0
  particle = new G4ExcitedMesons(            
	   "a0(980)0",       980.0*MeV,      75.0*MeV,           0.0,
                    0,              +1,            +1,          
                    2,               0,            -1,             
              "meson",               0,             0,       9000111,
		false,             0.0,          NULL);
  particle->SetAntiPDGEncoding(9000111);
  // set sub type
  particle->SetMultipletName("a0(980)");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of a0(980)0 -> eta + pi0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("a0(980)0",1.000, 2,
				                    "pi0","eta");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  // f0(500) (was f0(500) f0(400-1200))
  particle = new G4ExcitedMesons(            
            "f0(500)",       475.0*MeV,     550.0*MeV,           0.0,
                    0,              +1,            +1,          
                    0,               0,            +1,             
              "meson",               0,             0,       9000221,
		false,             0.0,          NULL);
  particle->SetAntiPDGEncoding(9000221);
  // set sub type
  particle->SetMultipletName("f0(500)");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of f0(500) -> pi + pi
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("f0(500)",1.000, 2,
				                    "pi+","pi-");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);


  // f0(980)
  particle = new G4ExcitedMesons(            
            "f0(980)",       990.0*MeV,      60.0*MeV,           0.0,
                    0,              +1,            +1,          
                    0,               0,            +1,             
              "meson",               0,             0,       9010221,
		false,             0.0,          NULL);
  particle->SetAntiPDGEncoding(9010221);
  // set sub type
  particle->SetMultipletName("f0(980)");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of f0(980) -> pi + pi
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("f0(980)",1.000, 2,
				                    "pi+","pi-");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  // eta(1405)
  particle = new G4ExcitedMesons(            
          "eta(1405)",      1408.8*MeV,      51.0*MeV,           0.0,
                    0,              -1,            +1,          
                    0,               0,            +1,             
              "meson",               0,             0,       9020221,
		false,             0.0,          NULL);
  particle->SetAntiPDGEncoding(9020221);
  // set sub type
  particle->SetMultipletName("eta(1405)");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of eta(1405) -> rho + rho
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("eta(1405)",1.000, 2,
				                   "rho+","rho-");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);
  // f0(1500)
  particle = new G4ExcitedMesons(            
           "f0(1500)",      1504.0*MeV,     109.0*MeV,           0.0,
                    0,              +1,            +1,          
                    0,               0,            +1,             
              "meson",               0,             0,       9030221,
		false,             0.0,          NULL);
  particle->SetAntiPDGEncoding(9030221);
  // set sub type
  particle->SetMultipletName("f0(1500)");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of f0(1500) -> eta + eta
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("f0(1500)",1.000, 2,
				                    "eta","eta");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  // f0(1710)
  particle = new G4ExcitedMesons(            
           "f0(1710)",      1723.0*MeV,     139.0*MeV,           0.0,
                    0,              +1,            +1,          
                    0,               0,            +1,             
              "meson",               0,             0,         10331,
		false,             0.0,          NULL);
  particle->SetAntiPDGEncoding(10331);
  // set sub type
  particle->SetMultipletName("f0(1710)");
  // create decay table
  decayTable =  new G4DecayTable();

  // create decay channel of f0(1710) -> k0 + k0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("f0(1710)",0.40, 2,
				                    "kaon0S","kaon0S");
  // add decay table
  decayTable->Insert(mode);
 
  // create decay channel of f0(1710) -> k+ + k+
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("f0(1710)",0.40, 2,
				                    "kaon+","kaon-");
  // add decay table
  decayTable->Insert(mode);
  
  // create decay channel of f0(1710) -> eta + eta
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("f0(1710)",0.20, 2,
				                    "eta","eta");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);


  // k_star+
  particle = new G4ExcitedMesons(            
	    "k_star+",      891.66*MeV,      50.8*MeV,    +1.0*eplus, 
                    2,              -1,             0,          
                    1,              +1,             0,             
              "meson",               0,             0,           323,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("k_star");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of k_star+ -> kaon+ + pi0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("k_star+",0.500, 2,
				                    "kaon+","pi0");
  // add decay table
  decayTable->Insert(mode);
   // create decay channel of k_star+ -> kaon+ + pi0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("k_star+",0.500, 2,
				                    "kaon0","pi+");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);
  
  // k_star0
  particle = new G4ExcitedMesons(            
	    "k_star0",      895.81*MeV,      47.4*MeV,     0.0*eplus, 
                    2,              -1,             0,          
                    1,              -1,             0,             
              "meson",               0,             0,           313,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("k_star");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of k_star0 -> kaon+ + pi-
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("k_star0",0.500, 2,
				                    "kaon+","pi-");
  // add decay table
  decayTable->Insert(mode);
   // create decay channel of k_star0 -> kaon0 + pi0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("k_star0",0.500, 2,
				                    "kaon0","pi0");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);

  // k_star-
  particle = new G4ExcitedMesons(            
            "k_star-",       891.66*MeV,     50.8*MeV,    -1.0*eplus, 
                    2,              -1,             0,          
                    1,              +1,             0,             
              "meson",               0,             0,          -323,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("k_star");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of k_star- -> kaon- + pi0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("k_star-",0.500, 2,
				                    "kaon-","pi0");
  // add decay table
  decayTable->Insert(mode);
   // create decay channel of k_star- -> anti_kaon0 + pi-
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("k_star-",0.500, 2,
				                    "anti_kaon0","pi-");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);
  
  
  // anti_k_star0
  particle = new G4ExcitedMesons(            
       "anti_k_star0",      895.81*MeV,      47.4*MeV,     0.0*eplus, 
                    2,              -1,             0,          
                    1,              -1,             0,             
              "meson",               0,             0,          -313,
		false,             0.0,          NULL);
  // set sub type
  particle->SetMultipletName("k_star");
  // create decay table
  decayTable =  new G4DecayTable();
  // create decay channel of anti_k_star0 -> kaon- + pi+
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("anti_k_star0",0.500, 2,
				                    "kaon-","pi+");
  // add decay table
  decayTable->Insert(mode);
   // create decay channel of anti_k_star0 -> anti_kaon0 + pi0
  //                                   parent    BR     #daughters
  mode  = new G4PhaseSpaceDecayChannel("anti_k_star0",0.500, 2,
				                    "anti_kaon0","pi0");
  // add decay table
  decayTable->Insert(mode);
  particle->SetDecayTable(decayTable);
  
}



















