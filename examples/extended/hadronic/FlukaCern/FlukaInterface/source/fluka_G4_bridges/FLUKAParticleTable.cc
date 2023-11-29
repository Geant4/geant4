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
#ifdef G4_USE_FLUKA


#include "FLUKAParticleTable.hh"

// G4
#include "G4ParticleTable.hh"


namespace fluka_particle_table {

  static G4bool fIsInitialized = false;

  // FLUKA
  static std::map<G4int, G4String> fFlukaIdtoFlukaName = {
    { -6, "4-HELIUM" },
    { -5, "3-HELIUM" },
    { -4, "TRITON"   },
    { -3, "DEUTERON" },
    { -2, "HEAVYION" },
    { -1, "OPTIPHOT" },
    {  0, "RAY"      },
    {  1, "PROTON"   },
    {  2, "APROTON"  },
    {  3, "ELECTRON" },
    {  4, "POSITRON" },
    {  5, "NEUTRIE"  },
    {  6, "ANEUTRIE" },
    {  7, "PHOTON"   },
    {  8, "NEUTRON"  },
    {  9, "ANEUTRON" },
    { 10, "MUON+"    },
    { 11, "MUON-"    },
    { 12, "KAONLONG" },
    { 13, "PION+"    },
    { 14, "PION-"    },
    { 15, "KAON+"    },
    { 16, "KAON-"    },
    { 17, "LAMBDA"   },
    { 18, "ALAMBDA"  },
    { 19, "KAONSHRT" },
    { 20, "SIGMA-"   },
    { 21, "SIGMA+"   },
    { 22, "SIGMAZER" },
    { 23, "PIZERO"   },
    { 24, "KAONZERO" },
    { 25, "AKAONZER" },
    { 27, "NEUTRIM"  },
    { 28, "ANEUTRIM" },
    { 31, "ASIGMA-"  },
    { 32, "ASIGMAZE" },
    { 33, "ASIGMA+"  },
    { 34, "XSIZERO"  },
    { 35, "AXSIZERO" },
    { 36, "XSI-"     },
    { 37, "AXSI+"    },
    { 38, "OMEGA-"   },
    { 39, "AOMEGA+"  },
    { 40, "WWLOWNEU" },
    { 41, "TAU+"     },
    { 42, "TAU-"     },
    { 43, "NEUTRIT"  },
    { 44, "ANEUTRIT" },
    { 45, "D+"       },
    { 46, "D-"       },
    { 47, "D0"       },
    { 48, "D0BAR"    },
    { 49, "DS+"      },
    { 50, "DS-"      },
    { 51, "LAMBDAC+" },
    { 52, "XSIC+"    },
    { 53, "XSIC0"    },
    { 54, "XSIPC+"   },
    { 55, "XSIPC0"   },
    { 56, "OMEGAC0"  },
    { 57, "ALAMBDC-" },
    { 58, "AXSIC-"   },
    { 59, "AXSIC0"   },
    { 60, "AXSIPC-"  },
    { 61, "AXSIPC0"  },
    { 62, "AOMEGAC0" },
    { 64, "@LASTPAR" },
    { 99, "ISOTOPE"  },
    {201, "ALL-PART" },
    {202, "ALL-CHAR" },
    {203, "ALL-NEUT" },
    {204, "ALL-NEGA" },
    {205, "ALL-POSI" },
    {206, "NUCLEONS" },
    {207, "NUC&PI+-" },
    {208, "ENERGY"   },
    {209, "PIONS+-"  },
    {210, "BEAMPART" },
    {211, "EM-ENRGY" },
    {212, "MUONS"    },
    {213, "E+&E-"    },
    {214, "AP&AN"    },
    {215, "KAONS"    },
    {216, "STRANGE"  },
    {217, "KAONS+-"  },
    {218, "HAD-CHAR" },
    {219, "FISSIONS" },
    {220, "HE-FISS"  },
    {221, "LE-FISS"  },
    {222, "NEU-BALA" },
    {223, "HAD-NEUT" },
    {224, "KAONS0"   },
    {225, "C-MESONS" },
    {226, "C-(A)BAR" },
    {227, "CHARMED"  },
    {228, "DOSE"     },
    {229, "UNB-ENER" },
    {230, "UNB-EMEN" },
    {231, "X-MOMENT" },
    {232, "Y-MOMENT" },
    {233, "Z-MOMENT" },
    {234, "ACTIVITY" },
    {235, "ACTOMASS" },
    {236, "SI1MEVNE" },
    {237, "HADGT20M" },
    {238, "NIEL-DEP" },
    {239, "DPA-SCO"  },
    {240, "DOSE-EQ"  },
    {241, "DOSE-EM"  },
    {242, "NET-CHRG" },
    {243, "DOSEQLET" },
    {244, "RES-NIEL" },
    {245, "DPA-NRES" },
    {246, "LOWENNEU" },
    {247, "NTLOWENE" },
    {248, "ALL-IONS" },
    {249, "HEHAD-EQ" },
    {250, "THNEU-EQ" },
    {251, "RES-NUCL" },
    {252, "DOSE-H2O" },
    {253, "ALPHA-D"  },
    {254, "SQBETA-D" },
    {255, "LGH-IONS" },
    {256, "HVY-IONS" },
    {257, "E+E-GAMM" },
    {258, "ANNIHRST" }
  };

  static std::map<G4int, G4int> fFlukaIdToPdgId = {
    {	1	,	2212	},
    {	2	,	-2212	},
    {	3	,	11	},
    {	4	,	-11	},
    {	5	,	12	},
    {	6	,	-12	},
    {	7	,	22	},
    {	8	,	2112	},
    {	9	,	-2112	},
    {	10	,	-13	},
    {	11	,	13	},
    {	12	,	130	},
    {	13	,	211	},
    {	14	,	-211	},
    {	15	,	321	},
    {	16	,	-321	},
    {	17	,	3122	},
    {	18	,	-3122	},
    {	19	,	310	},
    {	20	,	3112	},
    {	21	,	3222	},
    {	22	,	3212	},
    {	23	,	111	},
    {	24	,	311	},
    {	25	,	-311	},
    {	27	,	14	},
    {	28	,	-14	},
    {	31	,	-3222	},
    {	32	,	-3212	},
    {	33	,	-3112	},
    {	34	,	3322	},
    {	35	,	-3322	},
    {	36	,	3312	},
    {	37	,	-3312	},
    {	38	,	3334	},
    {	39	,	-3334	},
    {	41	,	-15	},
    {	42	,	15	},
    {	43	,	16	},
    {	44	,	-16	},
    {	45	,	411	},
    {	46	,	-411	},
    {	47	,	421	},
    {	48	,	-421	},
    {	49	,	431	},
    {	50	,	-431	},
    {	51	,	4122	},
    {	52	,	4232	},
    {	53	,	4132	},
    {	54	,	4322	},
    {	55	,	4312	},
    {	56	,	4332	},
    {	57	,	-4122	},
    {	58	,	-4232	},
    {	59	,	-4132	},
    {	60	,	-4322	},
    {	61	,	-4312	},
    {	62	,	-4332	}
  };

  // FLUKA -> G4
  static std::unordered_map<G4int, const G4ParticleDefinition*> fFlukaIdToG4Particle;
  static std::map<G4String, const G4ParticleDefinition*> fFlukaNameToG4Particle;
  static std::map<G4String, G4String> fFlukaNameToG4Name = {
    {"4-HELIUM", "alpha"         },
    {"3-HELIUM", "He3"           },
    {"TRITON"  , "triton"        },
    {"DEUTERON", "deuteron"      },
    {"HEAVYION", "ion"           },
    {"OPTIPHOT", "opticalphoton" },
    {"RAY"     , "geantino"      },
    {"PROTON"  , "proton"        },
    {"APROTON" , "anti_proton"   },
    {"ELECTRON", "e-"            },
    {"POSITRON", "e+"            },
    {"NEUTRIE" , "nu_e"          },
    {"ANEUTRIE", "anti_nu_e"     },
    {"PHOTON"  , "gamma"         },
    {"NEUTRON" , "neutron"       },
    {"ANEUTRON", "anti_neutron"  },
    {"MUON+"   , "mu+"           },
    {"MUON-"   , "mu-"           },
    {"KAONLONG", "kaon0L"        },
    {"PION+"   , "pi+"           },
    {"PION-"   , "pi-"           },
    {"KAON+"   , "kaon+"         },
    {"KAON-"   , "kaon-"         },
    {"LAMBDA"  , "lambda"        },
    {"ALAMBDA" , "anti_lambda"   },
    {"KAONSHRT", "kaon0S"        },
    {"SIGMA-"  , "sigma-"        },
    {"SIGMA+"  , "sigma+"        },
    {"SIGMAZER", "sigma0"        },
    {"PIZERO"  , "pi0"           },
    {"KAONZERO", "kaon0"         },
    {"AKAONZER", "anti_kaon0"    },
    {"NEUTRIM" , "nu_mu"         },
    {"ANEUTRIM", "anti_nu_mu"    },
    {"ASIGMA-" , "anti_sigma-"   },
    {"ASIGMAZE", "anti_sigma0"   },
    {"ASIGMA+" , "anti_sigma+"   },
    {"XSIZERO" , "xi0"           },
    {"AXSIZERO", "anti_xi0"      },
    {"XSI-"    , "xi-"           },
    {"AXSI+"   , "anti_xi-"      },
    {"OMEGA-"  , "omega-"        },
    {"AOMEGA+" , "anti_omega-"   },
    {"TAU+"    , "tau+"          },
    {"TAU-"    , "tau-"          },
    {"NEUTRIT" , "nu_tau"        },
    {"ANEUTRIT", "anti_nu_tau"   },
    {"D+"      , "D+"            },
    {"D-"      , "D-"            },
    {"D0"      , "D0"            },
    {"D0BAR"   , "anti_D0"       },
    {"DS+"     , "Ds+"           },
    {"DS-"     , "Ds-"           },
    {"LAMBDAC+", "lambda_c+"     },
    {"XSIC+"   , "xi_c+"         },
    {"XSIC0"   , "xi_c0"         },
    {"OMEGAC0" , "omega_c0"      },
    {"ALAMBDC-", "anti_lambda_c-"},
    {"AXSIC-"  , "anti_xi_c-"    },
    {"AXSIC0"  , "anti_xi_c0"    },
    {"AOMEGAC0", "anti_omega_c0" }
  };
  
  // G4 -> FLUKA
  static std::unordered_map<const G4ParticleDefinition*, G4int> fG4ParticleToFlukaId;
  static std::unordered_map<const G4ParticleDefinition*, G4String> fG4ParticleToFlukaName;

  // ***************************************************************************
  // Initialize all conversion tables.
  // Note that the starting data is the one from FLUKA.
  // ***************************************************************************
  void initialize() {
    if (fIsInitialized) return;
    fIsInitialized = true;

    G4cout << "" << G4endl;
    G4cout << "fluka_particle_table::initialize()" << G4endl;
    
    const auto particleTable = G4ParticleTable::GetParticleTable();

    // Loop on all particles
    for (const auto& particleIt : fFlukaIdtoFlukaName) {

      const G4int flukaId = particleIt.first;
      const G4String& flukaName = particleIt.second;

      const G4String& geantName = fFlukaNameToG4Name[flukaName];
      auto particle = particleTable->FindParticle(geantName);

      // If particle was not found in G4 table (via its name), 
      // try to find it via the PDG id
      if (!particle) {
        const auto& found = fFlukaIdToPdgId.find(flukaId);
        if (found != fFlukaIdToPdgId.end()) {
          const G4int pdgId = found->second;
          particle = particleTable->FindParticle(pdgId);
        }
      }

      if (particle) {
        fFlukaIdToG4Particle.insert(std::make_pair(flukaId, particle));
        fFlukaNameToG4Particle[flukaName] = particle;
        fG4ParticleToFlukaId.insert(std::make_pair(particle, flukaId));
        fG4ParticleToFlukaName.insert(std::make_pair(particle, flukaName));
        //G4cout << "Fluka2Geant: " << flukaId << ' ' << flukaName << ' ' << geantName
        //	<< " particle " << particle->GetParticleName() << G4endl;
      }
      else {
        G4cout << "" << G4endl;
        G4cout << "INFO: fluka_particle_table::initialize():" << G4endl;
        G4cout << "FLUKA particle, or category, named = " << flukaName 
               << " is not found as a particle in G4ParticleTable." 
               << G4endl;
        G4cout << "Should not be needed. If ever needed, an exception will be thrown at run time." << G4endl;
      }
    }
  }


  const G4String& fluka2name(const G4int id) { return fFlukaIdtoFlukaName[id]; }

  // FLUKA -> G4
  const G4ParticleDefinition* fluka2geant(const G4int ij) { return fFlukaIdToG4Particle[ij]; }
  const G4String& fluka2geantName(const G4String& name) { return fFlukaNameToG4Name[name]; }
  const G4ParticleDefinition* fluka2geant(const G4String& name) { return fFlukaNameToG4Particle[name]; }

  // G4 -> FLUKA
  G4int geant2fluka(const G4ParticleDefinition* def) { return fG4ParticleToFlukaId[def]; }
  const G4String& geant2flukaName(const G4ParticleDefinition* def) { return fG4ParticleToFlukaName[def]; }


}


#endif // G4_USE_FLUKA
