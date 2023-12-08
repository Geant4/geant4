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


#include "fluka_interface.hh"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <array>
#include <vector>
#include <string>
#include <unordered_map>
#include <cmath>

// CPP UTILS
#include "string_print.h"

// CPP TO FORTRAN UTILS
#include "types.h"
#include "string_conversion.h"
#include "flush.h"

// FLUKA DEPENDENCIES
#include "fluka_dependencies.hh"
#include "FLUKAParticleTable.hh"
#include "cmcyl.hh"

// G4
#include "Randomize.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4HadFinalState.hh"
#include "G4LorentzVector.hh"

#include "G4IonTable.hh"
#include "G4SystemOfUnits.hh"


#define DEBUG false



namespace fluka_interface {


  static G4bool initializeMaterial = true;
  static std::unordered_map<G4int, G4int> encounteredMaterials;


  // ***************************************************************************
  // Call the necessary FLUKA initialization routines,
  // initialize the necessary FLUKA COMMONs,
  // set physics options as requested.
  // ***************************************************************************
  void initialize(bool activateCoalescence,
                  G4bool activateHeavyFragmentsEvaporation) {
    G4cout << "        !!!!!!!!!!!!! fluka_interface::initialize()" << G4endl;

    // CHECK THAT THE FLUKA-G4 INTERFACE ENV HAS BEEN SOURCED
    if (std::getenv("FLUKA_PATH") == nullptr) {
      G4Exception("fluka_interface",
                  "FLUKA-G4 interface environment",
                  FatalException,
                  "FLUKA-G4 interface environment was not sourced.\n" \
                  "NB: Do not forget to first compile the FLUKA interface itself.\n" \
                  "For example: cd geant4/examples/extended/hadronic/FlukaCern/FlukaInterface/ " \
                  "&& make interface && make env\n" \
                  "FlukaInterface/env_FLUKA_G4_interface.sh then needs to be sourced\n" \
                  "in whichever terminal you want to use the FLUKA interface.\n");
    }

    // FLUKA ZEROING AND INITIALIZATION ROUTINES.
    G4cout << "CALL CMSPPR" << G4endl;
    cmsppr_();
    G4cout << "CALL ZEROIN" << G4endl;
    zeroin_();


    // FLUKA COMMONS INITIALIZATION 
    // Note: obviously, only initlialize the ones needed later on!
    parevt_.lheavy = true;
    parevt_.lgdhpr = true;

    // hadronic interaction GENerator FLaGs 
    genflg_.levpgn = true;
    // hadronic interaction GENerator THResholds
    genthr_.peanct = AINFNT;
    genthr_.peapit = AINFNT;
    genthr_.peakat = AINFNT;
    genthr_.peakbt = AINFNT;
    genthr_.peahyt = AINFNT;
    genthr_.peaaht = AINFNT;
    genthr_.peaant = AINFNT;
    cmcycl_.lout = 7;

    // NB: mmat is dummy here. 
    // All needed commons (including the material common) 
    // are set properly at mmax index on the fly at run time.
    G4int mmat = 3;		
    // FLuKa MATerial properties and atomic data 
    flkmat_.icomp [mmat - 1] = 0;
    flkmat_.rho   [mmat - 1] = ONEONE;

    // PHotoNuClear interaction
    phnccm_.rhphnc[mmat - 1] = ONEONE;
    phnccm_.ifphnc[mmat - 1] = 1111;

    // BEAM properties CoMmon
    G4double txx    = ZERZER;
    G4double tyy    = ZERZER;
    G4double tzz    = ONEONE;
    beamcm_.ubeam  = txx;
    beamcm_.vbeam  = tyy;
    beamcm_.wbeam  = tzz;

	       
    //  PHYSICS CARDS:
    std::array<G4double, 6> what;
    G4String sdum;

    // IMPORTANT
    // ACTIVATE COALESCENCE AND HEAVY FRAGMENTS EVAPORATION HERE!!
    // NB: THEY ARE NOT ACTIVATED BY DEFAULT IN G4!!
    if (activateCoalescence) {
      // ACTIVATE COALESCENCE
      what[0] = 1.;
      what[1] = 0.;
      what[2] = 0.;
      what[3] = 0.;
      what[4] = 0.;
      what[5] = 0.;
      sdum = "COALESCE";	
      G4cout << "CALL PHYCRD, sdum = " << sdum << G4endl;
      phycrd_( what.data(), cpp_to_fortran::convertStringNoSizeNeeded(sdum).data() );
    }
    		
    // IMPORTANT
    // ACTIVATE HEAVY FRAGMENTS EVAPORATION
    // NB: THEY ARE NOT ACTIVATED BY DEFAULT IN G4!! (right?)
    if (activateHeavyFragmentsEvaporation) {
      what[0] = 3.;
      what[1] = 0.;
      what[2] = 0.;
      what[3] = 0.;
      what[4] = 0.;
      what[5] = 0.;
      sdum = "EVAPORAT";	
      G4cout << "CALL PHYCRD, sdum = " << sdum << G4endl;
      phycrd_( what.data(), cpp_to_fortran::convertStringNoSizeNeeded(sdum).data() );
    }

    // IMPORTANT		
    // ACTIVATE CHARMED HADRONS TRANSPORT
    // This prevents charmed hadron decays in FLUKA.
    // Unfortunately, as a result, FLUKA also returns charmed Xi resonances,
    // which do not exist in G4, hence prevents the full use of this option.
    //
    // Could circumvent this issue by activating only the charmed Xi resonances decays in FLUKA,
    // (it is a resonance anyway!!), but this was a hack inside FLUKA code 
    // (not possible from an offical FLUKA release).
    const G4bool activateCharmedHadronTransport = false;
    if (activateCharmedHadronTransport) {
      what[0] = 1.;
      what[1] = 0.;
      what[2] = 0.;
      what[3] = 0.;
      what[4] = 0.;
      what[5] = 0.;
      sdum = "CHARMDEC";	
      G4cout << "CALL PHYCRD, sdum = " << sdum << G4endl;
      phycrd_( what.data(), cpp_to_fortran::convertStringNoSizeNeeded(sdum).data() );
    }


    // INITIALIZES THE FLUKA QUARK MODELS AND EVAPORATION MODULE
    G4cout << "CALL EVVINI" << G4endl;
    sdum = "DUMMY";
    evvini_(what.data(), cpp_to_fortran::convertStringNoSizeNeeded(sdum).data());


    // Needed because the printout buffers from C++ and fortran are distinct!!
    // Hence synchronize them, in order to get a meaningful printout order.
    cpp_to_fortran::flush();


    G4cout << "        !!!!!!!!!!!!! END OF fluka_interface::initialize()" << G4endl;
  }


  // ***************************************************************************
  // Returns FLUKA hadron nucleus inelastic XS.
  //
  // IMPORTANT note: if targetA is not passed as argument, it is set to 0 by default.
  // In that case, FLUKA will loop on the element's natural isotopic composition,
  // (composition defined in FLUKA data).
  // ***************************************************************************
  G4double computeInelasticScatteringXS(const G4DynamicParticle* projectile,
                                      const G4int targetZ,
                                      const G4int targetA) {   

    const G4ParticleDefinition* const projectileDefinition = projectile->GetDefinition();

    const G4bool projectileIsGenericHeavyIon = (projectileDefinition->IsGeneralIon()
                                              && projectileDefinition->GetAtomicNumber() > 2);
    const G4int projectileFLUKAId = (projectileIsGenericHeavyIon ? 
                                   -2
                                   : fluka_particle_table::geant2fluka(projectileDefinition));
    const G4double projectileKineticEnergy = projectile->GetKineticEnergy() / GeV;


    // RENAME INTO FLUKA-WORLD VARIABLES
    G4int kproj = projectileFLUKAId;
    G4double ekproj = projectileKineticEnergy;

    // Transformations needed (extracted from FLUKA code)
    updateFLUKAProjectileId(kproj);
    transformNonSupportedHadrons(kproj, ekproj);

    auto [ekin, pproj] = getKineticEnergyAndMomentum(ekproj, kproj);

    G4int ibtar = targetA;
    G4int ichtar = targetZ;


    // PRINTOUT INTERACTION INFO
#if DEBUG
    G4cout << "        !!!!!!!!!!!!! Call computeInelasticScatteringXS:"
              << " kproj = " << kproj
              << ", ekproj [GeV] = " << ekproj
              << ", pproj = " << pproj
              << ", ibtar = " << ibtar
              << ", ichtar = " << ichtar 
              << G4endl;
#endif

    
    // SET MATERIAL PROPERTIES
    G4int mmat = 3;
    flkmat_.ztar[mmat - 1] = std::abs(ichtar);

    ichtar = std::fabs(ichtar); // 'hole' freezing
    G4double atmss = 0.;
   
    // proton or neutron target:
    if ( ibtar == 1 ) {
      flkmat_.mssnum [mmat - 1] = ibtar;
      atmss  = static_cast<G4double> (ibtar);
    }
    // single isotope:
    else if ( ibtar > 0 ) {
      flkmat_.mssnum [mmat - 1] = ibtar;
      G4double bbtar  = ibtar;
#if DEBUG
      G4cout << "CALL exmsaz" << G4endl;
#endif
      G4int izdumm;
      logical lncmss = false;
      resnuc_.ammtar = bbtar * AMUGEV + EMVGEV * exmsaz_ ( bbtar, flkmat_.ztar [mmat - 1], lncmss, izdumm );
      atmss  = resnuc_.ammtar / AMUGEV;
    }
    // element:
    else {
      flkmat_.mssnum [mmat - 1] = 0;
      atmss  = ZERZER;
      // Loop on the stable isotopes
      for (G4int is = (isotop_.isondx[ichtar-1][0] - 1); is <= (isotop_.isondx[ichtar-1][1] - 1); ++is) {
        G4double bbtar  = isotop_.isomnm [is];
#if DEBUG
        G4cout << "CALL exmsaz" << G4endl;
#endif
        G4int izdumm;
        logical lncmss = false;
        resnuc_.ammtar = bbtar * AMUGEV + EMVGEV * exmsaz_ ( bbtar, flkmat_.ztar [mmat - 1], lncmss, izdumm );
        atmss  = atmss + resnuc_.ammtar / AMUGEV * isotop_.abuiso [is];
      }
    }

    flkmat_.amss[mmat - 1] = atmss;

#if DEBUG
    G4cout << "flkmat_.amss [mmat - 1] = " << flkmat_.amss [mmat - 1] 
              << ", flkmat_.msindx [mmat - 1] = " << flkmat_.msindx [mmat - 1]
              << G4endl;
#endif


    // COMPUTES XS
    G4double sigrea = 0.;

    // electron / positron
    if ( kproj == 3 || kproj == 4 ) {
      sigrea = ZERZER;
    }
    // photon
    else if ( kproj == 7 ) {
#if DEBUG
      G4cout << "CALL PPHCHO" << G4endl;
#endif
      G4double pphnsg = ZERZER;
      pphcho_(ekin, mmat, pphnsg, false);
      sigrea = pphnsg / flkmat_.rho [mmat - 1] * flkmat_.amss [mmat - 1] / AVOGAD * 1.E+27;
    }
    // standard case
    else {
#if DEBUG
      G4cout << "CALL SIGINM" << G4endl;
#endif

      const G4double projectileKineticEnergyPerNucleon = projectileKineticEnergy * (paprop_.ibarch[kproj+6] <= 1 ?
                                                                                  1.
                                                                                  : 1. / paprop_.ibarch[kproj+6]);

      // Assign sgtbcm_.lblsgi: 
      // recognize inelastic XS which should not be set to 0 when kE/n < EKSIG0.
      cksigi_(mmat);      
      
      // VERY IMPORTANT: In FLUKA, all inelastic XS (except the cases below),
      // are set to 0 when kE/n < EKSIG0.
      // See SIGTAB->SGTINL (initialization) and the runtime calls to SGTTOT (from KASKAD).
      if (projectileKineticEnergyPerNucleon < EKSIG0 
          && kproj != 12 && kproj != 19 && kproj != 24 && kproj != 25 // neutron kaons
          && kproj != 9  // antineutron
          && kproj != -3 // deuteron 
          && !sgtbcm_.lblsgi[sgtbcm_.ijsigi[kproj+6] - 1][mmat - 1]
          ) {
        sigrea = 0.;
      }
      else {
        G4double zldum = ZERZER;
        siginm_(kproj, mmat, ekin, pproj, sigrea, zldum);
      }
    }

#if DEBUG
    G4cout << " Adopted reaction sigma :" 
              << cpp_utils::sformat("%#9.4f", sigrea)
              << "[mb]."
              << G4endl;
#endif


    // FLUKA returned XS in mb
    return sigrea * millibarn; 
  }


  // ***************************************************************************
  // Computes the FLUKA hadron nucleus inelastic interaction final state,
  // and set its to the G4HadFinalState object passed as argument. 
  // ***************************************************************************
  void setNuclearInelasticFinalState(G4HadFinalState* const finalState,
                                     const G4HadProjectile& projectile, 
                                     const G4Nucleus& targetNucleus) {

    const G4ParticleDefinition* const projectileDefinition = projectile.GetDefinition();

    const G4bool projectileIsGenericHeavyIon = (projectileDefinition->IsGeneralIon()
                                              && projectileDefinition->GetAtomicNumber() > 2);
    const G4int projectileFLUKAId = (projectileIsGenericHeavyIon ? 
                                   -2
                                   : fluka_particle_table::geant2fluka(projectileDefinition));
    const G4int projectileA = (projectileIsGenericHeavyIon ? 
                             projectileDefinition->GetAtomicMass() 
                             : 0);
    const G4int projectileZ = (projectileIsGenericHeavyIon ? 
                             projectileDefinition->GetAtomicNumber() 
                             : 0);
    const G4double projectileKineticEnergy = projectile.GetKineticEnergy() / GeV;
    const G4int targetA = targetNucleus.GetA_asInt();
    const G4int targetZ = targetNucleus.GetZ_asInt();

    const auto& projectileDirection = projectile.GetMomentumDirection();
    G4double txx    = projectileDirection.x();
    G4double tyy    = projectileDirection.y();
    G4double tzz    = projectileDirection.z();


    // RENAME INTO FLUKA-WORLD VARIABLES
    G4int kproj = projectileFLUKAId;
    G4double ekproj = projectileKineticEnergy;

    const G4int iaproj = projectileA;
    const G4int izproj = projectileZ;		  					   
    G4int ibtar = targetA;
    G4int ichtar = targetZ;


    // PRINTOUT INTERACTION INFO
#if DEBUG
    G4cout << std::setprecision(16) << "        !!!!!!!!!!!!! Call setNuclearInelasticFinalState:"
              << " kproj = " << kproj
              << ", ekproj [GeV] = " << ekproj
      //<< ", pproj = " << pproj
              << ", ibtar = " << ibtar
              << ", ichtar = " << ichtar
              << G4endl;
#endif


    // Some initialization magic directly taken from FLUKA
    if (kproj ==  -1 || kproj ==  -8 ||
        kproj ==  -9 || kproj == -13 || kproj == -14 ||
        kproj == -23 || kproj == -15 || kproj == -16 ||
        kproj == -24 || kproj == -25 || kproj == -11 ||
        kproj ==  -7 ) {
      if ( -kproj != 7 ) {
        parevt_.lpreex = false;
        parevt_.lhlfix = false;
      }
    }
    else if ( kproj >= 94 && kproj < 49999 ) {
      if ( (kproj - 100) != 7 ) {
        parevt_.lpreex = false;
        parevt_.lhlfix = false;
      }
    }
    updateFLUKAProjectileId(kproj);


    // ION CASE
    if ( kproj == -2 ) {
#if DEBUG
      G4cout << "iaproj = " << iaproj << ", izproj = " << izproj << G4endl;
#endif
      beamcm_.ijhion = izproj * 1000 + iaproj;
    }


    // BEAM PROPERTIES COMMON
    beamcm_.ijbeam = kproj;

    parevt_.lparwv = false;

    G4double ekmax = 10000.;
    beamcm_.pbeam  = -ekmax;


    // ION CASE
    if ( kproj == -2 ) {
      beamcm_.ijhion = izproj * 1000 + iaproj;
      beamcm_.ijhion = beamcm_.ijhion * 100  + KXHEAV;
      G4int ionid  = beamcm_.ijhion;
#if DEBUG
      G4cout << "CALL DCDION" << G4endl;
#endif
      dcdion_(ionid);
#if DEBUG
      G4cout << "CALL SETION" << G4endl;
#endif
      setion_(ionid);
      beamcm_.ijbeam = ionid;
      thrscm_.lhvtrn = true;
      parevt_.lheavy = true;
      thrscm_.khvtrn = 6;
    }


    // CHECK WHETHER THE MATERIAL WAS ALREADY ENCOUNTERED
    // - If (targetZ, targetA) has never been encountered: 
    // associate it to a mmat, and initialize the relevant materials COMMONS at mmat index.
    // - If (targetZ, targetA) has already been encountered:
    // just return the associated mmat index.
    // Values from the materials COMMONS are read with mmat index.
    // - IMPORTANT NB: In ANY case, before EACH event, LIKE IS DONE IN FLUKA eventv,
    // mssnum and icriso(2) are always set properly
    G4int mmat = -1;
    const G4int targetLinearizedIndex = targetZ * 100000 + targetA;
    const auto& foundMaterial = encounteredMaterials.find(targetLinearizedIndex);
    if (foundMaterial == encounteredMaterials.end()) {
      initializeMaterial = true;
      mmat = 4 + encounteredMaterials.size();
      encounteredMaterials.insert(std::make_pair(targetLinearizedIndex, mmat));
    }
    else {
      initializeMaterial = false;
      mmat = foundMaterial->second;
    }
#if DEBUG
    G4cout << "targetZ = " << targetZ << G4endl;
    G4cout << "targetA = " << targetA << G4endl;
    G4cout << "initializeMaterial = " << initializeMaterial << G4endl;
    G4cout << "mmat = " << mmat << G4endl;
#endif


    // INITIALIZE MATERIAL IF NEEDED.
    if (initializeMaterial) {
      // This was also placed in initialization (but for dummy mmat = 3)
      // FLuKa MATerial properties and atomic data 
      flkmat_.icomp [mmat - 1] = 0;
      flkmat_.rho   [mmat - 1] = ONEONE;
      // PHotoNuClear interaction
      phnccm_.rhphnc[mmat - 1] = ONEONE;
      phnccm_.ifphnc[mmat - 1] = 1111;


      flkmat_.ztar[mmat - 1] = std::abs(ichtar);

      //  No hole freezing:
      if (ichtar > 0) {
        parevt_.lhlfix = false;
      } 
      // "Hole" freezing:
      else {
        parevt_.lhlfix = true;
        ichtar = -ichtar;
      }


      // SET MATERIAL PROPERTIES
      std::vector<G4int> iaiso;
      std::vector<G4int> iziso;
      std::vector<G4int> iiiso;

      logical lrmsch = false;
      logical lrd1o2 = true;
      logical ltrasp = true;
      
      G4double atmss = 0.;
      //  Proton or neutron target:
      if ( ibtar == 1 ) {
        flkmat_.mssnum [mmat - 1] = ibtar;

        atmss  = static_cast<G4double> (ibtar);

        G4double bbres  = ibtar;
        G4double zzres  = ichtar;
#if DEBUG
        G4cout << "CALL exmsaz" << G4endl;
#endif
        G4int izdum;
        logical lncmss = true;
        resnuc_.amnres = bbres * AMUC12 + EMVGEV * exmsaz_( bbres, zzres, lncmss, izdum );
#if DEBUG
        G4cout << "CALL amnama" << G4endl;
#endif
        resnuc_.ammres = amnama_ ( resnuc_.amnres, ibtar, ichtar );
        G4int niso   = 0;

#if DEBUG
        G4cout << "PROTON / NEUTRON CALL WSTOAP" << G4endl;
#endif			
        wstoap_(niso, iaiso.data(), iziso.data(), iiiso.data(), 
                lrmsch, lrd1o2, ltrasp, 
                0); // 8
#if DEBUG
        cpp_to_fortran::flush();
#endif
      }
      //  Single isotope:
      else if ( ibtar > 0 ) {
        logical setIsotope = false;

        flkmat_.mssnum [mmat - 1] = ibtar;
        G4double bbtar  = ibtar;
#if DEBUG
        G4cout << "CALL exmsaz" << G4endl;
#endif
        G4int izdumm;
        logical lncmss = false;
        resnuc_.ammtar = bbtar * AMUGEV + EMVGEV * exmsaz_ ( bbtar, flkmat_.ztar [mmat - 1], lncmss, izdumm );

        atmss  = resnuc_.ammtar / AMUGEV;

        nucgid_.rhotab [ibtar - 2] = ZERZER;
        G4int niso = 1;
        iaiso.emplace_back(ibtar);
        iziso.emplace_back(ichtar);
        iiiso.resize(1);

        //  Loop on the stable isotopes
        for (G4int is = (isotop_.isondx[ichtar-1][0] - 1); is <= (isotop_.isondx[ichtar-1][1] - 1); ++is) {
          if ( isotop_.isomnm [is] == ibtar ) {
            iiiso [0] = is + 1;
            setIsotope = true;
            break;
          }
        }
        if (!setIsotope) {
          //  Loop on the "quasi" stable isotopes
          for (G4int iq = (isotop_.isqndx[ichtar-1][0] - 1); iq <= (isotop_.isqndx[ichtar-1][1] - 1); ++iq) {
            if ( isotop_.isqmnm [iq] == ibtar ) {
              iiiso [0] = NSTBIS + iq + 1;
              setIsotope = true;
              break;
            }
          }
        }
        if (!setIsotope) {
          iiiso [0] = NTSTIS + 1;
        }

        const G4int isindx = iiiso [0];
        flkmat_.msindx [mmat - 1] = isindx;

        // Initialize nuclear geometry data
        // WARNING: This call is extremely time-consuming
        // (hence the need to call it only the first time the material is encountered).
        if (nucgid_.rhotab [ibtar - 2] < 0.1 ) {
#if DEBUG
          G4cout << "ISOTOPE CALL WSTOAP" << G4endl;
#endif
          wstoap_(niso, iaiso.data(), iziso.data(), iiiso.data(), 
                  lrmsch, lrd1o2, ltrasp, 
                  0);
        }
      }
      // Element:
      else {
        flkmat_.mssnum [mmat - 1] = 0;

        atmss  = ZERZER;

        G4int niso   = 0;
        G4int iwstop = 0;
        // Loop on the stable isotopes
        for (G4int is = (isotop_.isondx[ichtar-1][0] - 1); is <= (isotop_.isondx[ichtar-1][1] - 1); ++is) {
          G4double bbtar  = isotop_.isomnm [is];
#if DEBUG
          G4cout << "CALL exmsaz" << G4endl;
#endif
          G4int izdumm;
          logical lncmss = false;
          resnuc_.ammtar = bbtar * AMUGEV + EMVGEV * exmsaz_ ( bbtar, flkmat_.ztar [mmat - 1], lncmss, izdumm );

          atmss  = atmss + resnuc_.ammtar / AMUGEV * isotop_.abuiso [is];

          if ( isotop_.isomnm [is] > 1 ) {
            nucgid_.rhotab [isotop_.isomnm[is] - 2] = ZERZER;
            niso = niso + 1;
            iaiso.emplace_back(isotop_.isomnm[is]);
            iziso.emplace_back(ichtar);
            iiiso.emplace_back(is + 1);
            if ( nucgid_.rhotab [isotop_.isomnm[is] - 2] < 0.1 ) { iwstop = iwstop + 1; }
          }
        }
        // Initialize nuclear geometry data
        if ( iwstop > 0 ) {
#if DEBUG
          G4cout << "ELEMENT CALL WSTOAP" << G4endl;
#endif				
          wstoap_( niso, iaiso.data(), iziso.data(), iiiso.data(),
                   lrmsch, lrd1o2, ltrasp, 0); // 8
#if DEBUG
          cpp_to_fortran::flush();
#endif
        }

      }


      flkmat_.amss[mmat - 1] = atmss;
    } // end of MATERIAL INITIALIZE


    // NUCLEON DECAY INITIALIZATION
    G4int kkproj = kproj;
    if ( kkproj >= 312500000 ) {
#if DEBUG
      G4cout << "CALL NCDCYI" << G4endl;
#endif
      ncdcyi_(kproj, kkproj, mmat);
#if DEBUG
      G4cout << G4endl;
      G4cout << "   " <<  chpprp_.prname [ndnicm_.ndnitr + 6] << " decay into:" << G4endl;
      G4cout << G4endl;
      for (G4int is = 1; is <= ndnicm_.ncdcsc; ++is) {
        G4cout << "             " <<  chpprp_.prname[ndnicm_.kncdcs[is] + 6] << G4endl;
      }
      G4cout << G4endl;
      G4cout << "CALL NCDCYR" << G4endl;
#endif
      ncdcyr_();
    }
    else if ( kkproj >= 49999 ) {
      const G4double ahelp  = kkproj;
      const G4int iahelp = std::round( ahelp * 1.E-04 ) * 10000;
      kproj = iahelp / 10000;
    }
    else {
      ndnicm_.ndnitr = 0;
      ndnicm_.ncdcsc = 0;
    }

       
    // STATUS SUMMARY
#if DEBUG
    if (parevt_.lpreex) { 
      G4cout << "    Most energetic intranuclear interactions treated explicitly" << G4endl;
    }
    else {
      G4cout << "    No explicit treatment of intranuclear interactions" << G4endl;
    }

    if (parevt_.lhlfix) {
      G4cout << G4endl;
      G4cout << "    Explicitly treated holes frozen " << G4endl;
    }
    if (parevt_.lparwv) {
      G4cout << G4endl;
      G4cout << "    De Broglie w.l. applied " << G4endl;
    }
    if (incpot_.lrhfl1 [currpt_.iptcur - 1]) {
      G4cout << G4endl;
      G4cout << "    Target nucleon center distribution unfolded and used" << G4endl;
    }
    if (incpot_.lrhfl2 [currpt_.iptcur - 1]) {
      G4cout << G4endl;
      G4cout << "    Projectile RMS radius folded with the density distribution" << G4endl;
    }
    if (incpot_.lrhfl3 [currpt_.iptcur - 1]) {
      G4cout << G4endl;
      G4cout << "    Effective range for pion-nuc. G4int. folded in the pion pot." << G4endl;
    }
#endif


    // ELD
    if (paprop_.ichrge [kproj + 6] == 0) { cmelds_.leldss = false; }


    // BEAM AND GEOEMTRY EXTRA SETS
    auto [ekin, pproj] = getKineticEnergyAndMomentum(ekproj, kproj);
#if DEBUG
    G4cout << "pproj = " << pproj << G4endl;
#endif
    beamcm_.pbeam  = pproj;
    // SUMmary COUnters
    sumcou_.weipri = ZERZER;
    caslim_.ncase  = 0;
    // NUClear Geometry Input data
    nucge3_.ievpre = 0;

#if DEBUG
    G4cout << G4endl;
    const G4String particleName = fortran_to_cpp::convertString(chpprp_.prname[kproj + 6], 8);

    G4cout << "  Projectile: " << particleName 
              << " E =" << cpp_utils::sformat("%#9.4f", ekin) << " GeV" 
              << G4endl;
#endif


    G4double wee = ONEONE;
    G4int iev = 1;
		      		      
    cmcycl_.ievt = iev;
    parevt_.levfin = false;


#if DEBUG
    cpp_to_fortran::flush();
#endif
    ///////////////////////////////////////////////////////////////////////////////////
    // IMPORTANT
    // FLUKA HADRONIC EVENT GENERATOR!!
    eventv_(kkproj, pproj, ekin, 
            txx, tyy, tzz, 
            wee,
            mmat);
    ///////////////////////////////////////////////////////////////////////////////////
#if DEBUG
    cpp_to_fortran::flush();
#endif

    if (resnuc_.lfkevt) {
      nucge2_.opacty = AINFNT;
    }


    ///////////////////////////////////////////////////////////////////////////////////
    // IMPORTANT:
    // LOOP ON ALL SECONDARIES, THE HEAVY FRAGMENTS AND THE RESNUC
    // AND ADD THEM TO FINAL STATE
    ///////////////////////////////////////////////////////////////////////////////////


    // SECONDARIES
#if DEBUG
    G4cout << "Number of secondaries = " << genstk_.np << G4endl;
#endif
    for (G4int secondaryIndex = 0; secondaryIndex < genstk_.np; ++secondaryIndex) {

      const G4int flukaId = genstk_.kpart[secondaryIndex];

      // Map FLUKA Id to G4 id.
      const G4ParticleDefinition* const definition = fluka_particle_table::fluka2geant(flukaId);
      if (!definition) {
        G4cout << "ERROR! FLUKA to G4 particle id conversion."
                  << " Could not find FLUKA id = " << flukaId 
                  << G4endl;
      }

      //const G4double momentum = genstk_.plr[secondaryIndex] * GeV;        // GeV/c in FLUKA
      const G4double momentumX = genstk_.cxr[secondaryIndex]; 
      const G4double momentumY = genstk_.cyr[secondaryIndex];
      const G4double momentumZ = genstk_.czr[secondaryIndex];
      const G4double kineticEnergy = genstk_.tki[secondaryIndex] * GeV;   // GeV in FLUKA

#if DEBUG
      G4cout << "FLUKA id = " << flukaId << G4endl;
      G4cout << "kineticEnergy [MeV] = " << kineticEnergy
        //<< ", momentum = " << momentum
                << ", momentumX = " << momentumX
                << ", momentumY = " << momentumY
                << ", momentumZ = " << momentumZ
                << G4endl;
#endif

      const G4ThreeVector momentumDirection = G4ThreeVector(momentumX,
                                                            momentumY,
                                                            momentumZ);

      G4DynamicParticle* const secondaryParticle = new G4DynamicParticle(definition,
                                                                         momentumDirection,
                                                                         kineticEnergy);
      finalState->AddSecondary(secondaryParticle);
    }


    // HEAVY FRAGMENTS
#if DEBUG
    G4cout << G4endl;
    G4cout << "Number of heavy fragments = " << fheavy_.npheav << G4endl;
#endif
    for (G4int fragmentIndex = 0; fragmentIndex < fheavy_.npheav; ++fragmentIndex) {

      const G4int flukaId = -fheavy_.kheavy[fragmentIndex];

      // Map FLUKA Id to G4 id.
      const G4ParticleDefinition* definition = nullptr;
      if (flukaId >= -6 && flukaId <= -3) {
        definition = fluka_particle_table::fluka2geant(flukaId);
#if DEBUG
        G4cout << "fluka_particle_table::fluka2geant(flukaId) = " << fluka_particle_table::fluka2geant(flukaId) << G4endl;
#endif
      }
      else {
        const G4int A = fheavy_.ibheav[-flukaId - 1];
        const G4int Z = fheavy_.icheav[-flukaId - 1];
        definition = G4IonTable::GetIonTable()->GetIon(Z, A);
      }
      if (!definition) {
        G4cout << "ERROR! FLUKA to G4 particle id conversion."
                  << " Could not find FLUKA id = " << flukaId 
                  << G4endl;
      }

      //const G4double momentum = fheavy_.pheavy[fragmentIndex] * GeV;        // GeV/c in FLUKA
      const G4double momentumX = fheavy_.cxheav[fragmentIndex]; 
      const G4double momentumY = fheavy_.cyheav[fragmentIndex];
      const G4double momentumZ = fheavy_.czheav[fragmentIndex];
      const G4double kineticEnergy = fheavy_.tkheav[fragmentIndex] * GeV;   // GeV in FLUKA

#if DEBUG
      G4cout << "FLUKA id = " << flukaId << G4endl;
      G4cout << "kineticEnergy [MeV] = " << kineticEnergy
        //<< ", momentum = " << momentum
                << ", momentumX = " << momentumX
                << ", momentumY = " << momentumY
                << ", momentumZ = " << momentumZ
                << G4endl;
#endif

      const G4ThreeVector momentumDirection = G4ThreeVector(momentumX,
                                                            momentumY,
                                                            momentumZ);

      G4DynamicParticle* const heavyFragment = new G4DynamicParticle(definition,
                                                                     momentumDirection,
                                                                     kineticEnergy);

      finalState->AddSecondary(heavyFragment);			
    }


    // RESNUC
    const G4int residualNucleusA = resnuc_.ibres;
    const G4int residualNucleusZ = resnuc_.icres;

    if (residualNucleusA > 0 && residualNucleusZ != 0) {

      // Map FLUKA Id to G4 id.		
      const G4ParticleDefinition* const definition = G4IonTable::GetIonTable()->GetIon(residualNucleusZ, 
                                                                                       residualNucleusA);
      if (!definition) {
        G4cout << "ERROR! FLUKA to G4 particle id conversion."
                  << " Could not find residual nucleus A = " << residualNucleusA 
                  << ", Z = " << residualNucleusZ
                  << G4endl;
      }

      const G4double momentum = resnuc_.ptres * GeV;        // GeV/c in FLUKA
      const G4double momentumX = resnuc_.pxres * GeV;         // GeV/c in FLUKA
      const G4double momentumY = resnuc_.pyres * GeV;         // GeV/c in FLUKA
      const G4double momentumZ = resnuc_.pzres * GeV;         // GeV/c in FLUKA
      const G4double kineticEnergy = resnuc_.ekres * GeV;     // GeV in FLUKA

#if DEBUG
      G4cout << G4endl;
      G4cout << "Residual nucleus A = " << residualNucleusA 
                << ", Z = " << residualNucleusZ 
                << ", kineticEnergy [MeV] = " << resnuc_.ekres
                << ", momentum [MeV] = " << resnuc_.ptres
                << G4endl;
#endif

      const G4ThreeVector momentumDirection = G4ThreeVector(momentumX / momentum,
                                                            momentumY / momentum,
                                                            momentumZ / momentum);
			
      G4DynamicParticle* const residualNucleus = new G4DynamicParticle(definition,
                                                                       momentumDirection,
                                                                       kineticEnergy);

      finalState->AddSecondary(residualNucleus);
    }

#if DEBUG
    G4cout << G4endl;
#endif


    // IMPORTANT:
    // NOW RESET PARTS OF THE COMMONS (as done in FLUKA)
    genstk_.np = 0;
    genstk_.np0 = 0;
    fheavy_.npheav = 0;
    genstk_.tv = 0.;
    genstk_.tvcms = 0.;
    genstk_.tvrecl = 0.;
    genstk_.tvheav = 0.;
    genstk_.tvbind = 0.;
    resnuc_.ibres  = 0;
    resnuc_.icres  = 0;
  }


  // ***************************************************************************
  // Awful shenanigan from FLUKA to play with particleId,
  // kept exactly as-is.
  // ***************************************************************************
  void updateFLUKAProjectileId(G4int& kproj) {
    if (kproj ==  -1 || kproj ==  -8
        || kproj ==  -9 || kproj == -13 || kproj == -14
        || kproj == -23 || kproj == -15 || kproj == -16
        || kproj == -24 || kproj == -25 || kproj == -11
        || kproj ==  -7
        ) {
      kproj = - kproj;
    }
    else if ( kproj >= 94 && kproj < 49999 ) {
      kproj  = kproj - 100;
    }
  }


  // ***************************************************************************
  // A few particles are converted BEFORE the hadronic interaction,
  // to sth FLUKA hadronic event generator can digest.
  // ***************************************************************************
  void transformNonSupportedHadrons(G4int& kproj, G4double& ekproj) {

    // K0S or K0L -> K0 or K0bar
    if (kproj == 12 || kproj == 19) {
      if (G4UniformRand() > 0.5) { kproj = 24; }
      else { kproj = 25; }
      // NB: No need to modify kinetic energy / momentum here,
      // since the mass is unchanged.
    }
    // Special FLUKA K0 should be a K0
    else if (kproj == 26) {
      kproj = 23;
    }
    else if (kproj == 17
             || kproj == 18
             || kproj == 20
             || kproj == 21
             || kproj == 22
             || kproj == 30
             || kproj == 31
             || kproj == 32
             || kproj == 33
             || kproj == 34
             || kproj == 35
             || kproj == 36
             || kproj == 37
             || kproj == 38
             || kproj == 39) {
      // Boolean used to signal that we are changing particle id.
      // LCHTYP = true;

      // Initial atomic mass (before updating the particle id).
      const G4double amIni = paprop_.am [kproj + 6];

      // change particle: see eventv

      // lambda -> neutron
      if (kproj == 17) {
        kproj = 8;
      }
      // antilambda -> antineutron
      else if (kproj == 18) {
        kproj = 9;
      }
      // sigma- -> neutron
      else if (kproj == 20) {
        kproj = 8;
      }
      // sigma+ -> proton
      else if (kproj == 21) {
        kproj = 1;
      }
      // sigma0 -> neutron
      else if (kproj == 22) {
        kproj = 8;
      }
      // virtual vector meson -> pi0
      else if (kproj == 30) {
        kproj = 23;
      }
      // antisigma- -> antiproton
      else if (kproj == 31) {
        kproj = 2;
      }
      // antiSigma0 -> antineutron
      else if (kproj == 32) {
        kproj = 9;
      }
      // antiSigma+ -> antineutron
      else if (kproj == 33) {
        kproj = 9;
      }
      // xsi0 -> neutron
      else if (kproj == 34) {
        kproj = 8;
      }
      // antixsi0 -> antineutron
      else if (kproj == 35) {
        kproj = 9;
      }
      // xsi- -> proton
      else if (kproj == 36) {
        kproj = 1;
      }
      // antiXsi+ -> antiproton
      else if (kproj == 37) {
        kproj = 2;
      }
      // omega- -> proton
      else if (kproj == 38) {
        kproj = 1;
      }
      // antiomega+ -> antiproton
      else if (kproj == 39) {
        kproj = 2;
      }

      // IMPORTANT
      // Update kinetic energy.
      ekproj += amIni - paprop_.am [kproj + 6];

      // NB: Momentum is computed later on, with getKineticEnergyAndMomentum
    }
  }


  // ***************************************************************************
  // Kinetic energy and momentum utility from FLUKA.
  // ***************************************************************************
  std::pair<G4double, G4double> getKineticEnergyAndMomentum(const G4double ekproj, const G4int kproj) {
    const G4double ekin = ekproj * (kproj != -2 
                                  ? 1.
                                  : paprop_.am [-2 + 6] / AMUC12);
    const G4double pproj = std::sqrt(ekin * ( ekin + TWOTWO * paprop_.am [kproj + 6]));
    return {ekin, pproj};
  }	

} // namespace fluka_interface


#endif // G4_USE_FLUKA
