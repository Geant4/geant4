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
// *                                                                  *
// * Parts of this code which have been  developed by Abdel-Waged     *
// * et al under contract (31-465) to the King Abdul-Aziz City for    *
// * Science and Technology (KACST), the National Centre of           *
// * Mathematics and Physics (NCMP), Saudi Arabia.                    *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file G4UrQMD1_3Model.cc
/// \brief Implementation of the G4UrQMD1_3Model class

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:          G4UrQMD1_3Model.hh
//
// Version:          0.B
// Date:           25/01/12
// Authors:        Kh. Abdel-Waged and Nuha Felemban
// Revised by:      V.V. Uzhinskii
//                  SPONSERED BY
// Customer:        KAUST/NCMP
// Contract:        31-465
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
#ifdef G4_USE_URQMD

#  include "G4UrQMD1_3Model.hh"

#  include "G4UrQMD1_3Interface.hh"
//-------------------------------
#  include "G4CollisionOutput.hh"
#  include "G4DynamicParticle.hh"
#  include "G4HadronicException.hh"
#  include "G4IonTable.hh"
#  include "G4LorentzRotation.hh"
#  include "G4Nucleus.hh"
#  include "G4ParticleDefinition.hh"
#  include "G4ParticleTable.hh"
#  include "G4PhysicalConstants.hh"
#  include "G4SystemOfUnits.hh"
#  include "G4Track.hh"
#  include "G4V3DNucleus.hh"
#  include "Randomize.hh"
#  include "globals.hh"

// AND->
#  include "G4Version.hh"
// AND<-
//----------------new_anti
#  include "G4AntiAlpha.hh"
#  include "G4AntiDeuteron.hh"
#  include "G4AntiHe3.hh"
#  include "G4AntiTriton.hh"
//---------------------------
#  include <fstream>
#  include <string>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
G4UrQMD1_3Model::G4UrQMD1_3Model(const G4String& nam)
  : G4VIntraNuclearTransportModel(nam), verbose(0)
{
  if (verbose > 3) {
    G4cout << " >>> G4UrQMD1_3Model default constructor" << G4endl;
  }

  //
  // Set the minimum and maximum range for the UrQMD model

  //  SetMinEnergy(100.0*MeV);
  //  SetMaxEnergy(200.0*GeV);

  //

  //
  WelcomeMessage();
  //
  CurrentEvent = 0;
  //

  InitialiseDataTables();

  //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Destructor
//
G4UrQMD1_3Model::~G4UrQMD1_3Model() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ReactionProductVector* G4UrQMD1_3Model::Propagate(G4KineticTrackVector*, G4V3DNucleus*)
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// ApplyYourself
//
// Member function to process an event, and get information about the products.

G4HadFinalState* G4UrQMD1_3Model::ApplyYourself(const G4HadProjectile& theTrack,
                                                G4Nucleus& theTarget)
{
  // anti_new
  //   ------------------define anti_light_nucleus
  const G4ParticleDefinition* anti_deu = G4AntiDeuteron::AntiDeuteron();

  const G4ParticleDefinition* anti_he3 = G4AntiHe3::AntiHe3();

  const G4ParticleDefinition* anti_tri = G4AntiTriton::AntiTriton();

  const G4ParticleDefinition* anti_alp = G4AntiAlpha::AntiAlpha();
  //---------------------------------------------------
  //
  // The secondaries will be returned in G4HadFinalState &theResult -
  // initialise this.  The original track will always be discontinued and
  // secondaries followed.
  //
  theResult.Clear();
  theResult.SetStatusChange(stopAndKill);

  G4DynamicParticle* cascadeParticle = 0;
  //
  //
  // Get relevant information about the projectile and target (A, Z, energy/nuc,
  // momentum, etc).
  //

  const G4ParticleDefinition* definitionP = theTrack.GetDefinition();
  const G4double AP = definitionP->GetBaryonNumber();
  const G4double ZP = definitionP->GetPDGCharge();
  G4double AT = theTarget.GetA_asInt();  // mass number A (was GetN(); NOT GetN_asInt!)
  G4double ZT = theTarget.GetZ_asInt();  // charge Z
  //  -----------------------------------------------
  G4int id = definitionP->GetPDGEncoding();  // get particle encoding
  // ------------------------------------------------
  G4int AP1 = G4lrint(AP);
  G4int ZP1 = G4lrint(ZP);
  G4int AT1 = G4lrint(AT);
  G4int ZT1 = G4lrint(ZT);
  //  G4cout<<"------ap1--=="<<AP1<<"---zp1---=="<<ZP1<<"---id-=="<<id<< G4endl;
  //
  // ****************************************************************************
  // The following is the parameters necessary to initiate Uinit() and UrQMD()
  // ----------------------------------------------------------------------------
  urqmdparams_.u_sptar = 0;  //! 0= normal proj/target, 1=special proj/tar
  urqmdparams_.u_spproj = 1;  // projectile is a special particle

  // new_anti

  if (AP1 > 1 || definitionP == anti_deu || definitionP == anti_he3 || definitionP == anti_tri
      || definitionP == anti_alp)
  {
    urqmdparams_.u_ap = AP1;
    urqmdparams_.u_zp = ZP1;

    urqmdparams_.u_spproj = 0;
  }
  else if (id == 2212) {  //! proton
    urqmdparams_.u_ap = 1;
    urqmdparams_.u_zp = 1;
  }
  else if (id == -2212) {  //! anti-proton
    urqmdparams_.u_ap = -1;
    urqmdparams_.u_zp = -1;
  }
  else if (id == 2112) {  //! neutron
    urqmdparams_.u_ap = 1;
    urqmdparams_.u_zp = -1;
  }
  else if (id == -2112) {  //! anti-neutron
    urqmdparams_.u_ap = -1;
    urqmdparams_.u_zp = 1;
  }
  else if (id == 211) {  //! pi+
    urqmdparams_.u_ap = 101;
    urqmdparams_.u_zp = 2;
  }
  else if (id == -211) {  //! pi-
    urqmdparams_.u_ap = 101;
    urqmdparams_.u_zp = -2;
  }
  else if (id == 321) {  //! K+
    urqmdparams_.u_ap = 106;
    urqmdparams_.u_zp = 1;
  }
  else if (id == -321) {  //! K-
    urqmdparams_.u_ap = -106;
    urqmdparams_.u_zp = -1;
  }
  else if (id == 130 || id == 310) {  //  ! K0
    urqmdparams_.u_ap = 106;
    urqmdparams_.u_zp = -1;
  }
  else if (id == -130 || id == -310) {  // ! K0bar
    urqmdparams_.u_ap = -106;
    urqmdparams_.u_zp = 1;
  }
  else {
    G4cout << " Sorry, No definition for particle for UrQMD::" << id << "found" << G4endl;

    // AND->
#  if G4VERSION_NUMBER >= 950
    // New signature (9.5) for G4Exception
    // Using G4HadronicException
    throw G4HadronicException(__FILE__, __LINE__, "Sorry, no definition for particle for UrQMD");
#  else
    G4Exception(" ");
#  endif
    // AND<-
  }  // end if id
  //-------------------------------------------------------

  urqmdparams_.u_at = AT1;  // Target identified
  urqmdparams_.u_zt = ZT1;
  //----------------------------------------------------
  //  identify Energy
  //
  G4ThreeVector Pbefore = theTrack.Get4Momentum().vect();
  G4double T = theTrack.GetKineticEnergy();
  G4double E = theTrack.GetTotalEnergy();
  G4double TotalEbefore = E * AP1 + theTarget.AtomicMass(AT1, ZT1) + theTarget.GetEnergyDeposit();
  //    -----------------------------------------------------------------

  if (AP1 > 1) {
    urqmdparams_.u_elab = T / (AP1 * GeV);  // Units are GeV/nuc for UrQMD

    E = E / AP1;  // Units are GeV/nuc
  }
  else {
    urqmdparams_.u_elab = T / GeV;  // units are GeV

    TotalEbefore = E + theTarget.AtomicMass(AT1, ZT1) + theTarget.GetEnergyDeposit();
  }

  //------------------------------------------------------------
  // identify impact parameter (negative => UrQMD samples b, area-weighted, up to
  // |u_imp|). Use the interaction radius R_proj + R_target = 1.1*(Ap^1/3 + At^1/3),
  // matching UrQMD's own definition nucrad(Ap)+nucrad(At) in init.f. The previous
  // value used only the target radius, which omits the projectile size: negligible
  // for a hadron (Ap=1) but about half the radius for a heavy ion, forcing
  // over-central collisions.
  urqmdparams_.u_imp =
    -(1.1 * (std::pow(G4double(AP1), (1. / 3.)) + std::pow(G4double(AT1), (1. / 3.))));
  // units are in fm for UrQMD;
  //------------------------------------------------------------
  ///////////////////////// initialise/////////////////////
  // uinit() must run before every UrQMD() call, as the reference UrqmdEvent
  // driver does: it re-samples the projectile/target nucleons (positions and
  // Fermi momenta) via cascinit() for THIS event. Previously it ran only once
  // (guarded by CurrentEvent), which froze the nucleus for the whole run and
  // made re-sampled interactions identical, so had012 warnings could never be
  // recovered and escalated to a fatal had006. uinit()'s own internal "ifirst"
  // guard keeps the expensive one-time table loading (loadwtab, potentials,
  // norm_init) to the first call, so repeating it per event only redoes the
  // cheap per-event nucleus sampling.
  G4int io = 0;
  if (CurrentEvent == 0) {
    G4cout << "\n creation of table, wait-------\n" << G4endl;
  }

  uinit_(&io);

  if (CurrentEvent == 0) {
    G4cout << "\n end to create  table " << G4endl;
    CurrentEvent = 1;
  }
  ////////////////////////////////////////////////////////

  // #ifdef debug_G4UrQMD1_3Model

  G4cout << "UrQMDModel running-------------" << G4endl;

  urqmd_();

  // #endif

  // G4cout <<"Number of produced particles:  " <<sys_.npart<<G4endl;

  G4int nProduced = sys_.npart;  // no of produced particles
  if (nProduced < 2) {
    G4cout << "===============Warning================" << G4endl;
    G4cout << "======================================" << G4endl;

    G4cout << "Number of produced particles is very low:  " << sys_.npart << G4endl;
    G4cout << "============================================" << G4endl;

// AND->
#  if G4VERSION_NUMBER >= 950
    // New signature (9.5) for G4Exception
    // Using G4HadronicException instead of base class
    throw G4HadronicException(__FILE__, __LINE__, "Number of produced particle is very low");
#  else
    G4Exception(" ");  // stop
#  endif
    // AND<-
  }
  else {
    for (G4int i = 0; i < nProduced; i++) {
      G4int pid = pdgid_(&isys_.ityp[i], &isys_.iso3[i]);

      // Particle is a final state secondary and not a nucleus.
      // Determine what this secondary particle is, and if valid, load dynamic
      // parameters.
      //

      G4ParticleDefinition* pd = G4ParticleTable::GetParticleTable()->FindParticle(pid);

      if (pd) {
        // The secondary's physical momentum includes the Fermi motion that UrQMD
        // keeps in the separate ffermp array: init.f moves the Fermi momentum out of
        // px into ffermp (setting the grid px), and coload.f uses px+ffermp as the
        // real momentum in collisions. UrQMD's stored energy p0, however, is on-shell
        // for px ALONE (p0 = sqrt(px^2+fmass^2)), so the original code paired momentum
        // px+ffermp with energy et = p0 and produced an OFF-SHELL particle
        // (E^2 != p^2 + m^2). Keep the physical momentum px+ffermp and give it the
        // matching on-shell energy sqrt((px+ffermp)^2 + fmass^2), using UrQMD's
        // dynamical mass fmass.
        G4double pxf = coor_.px[i] + ffermi_.ffermpx[i];
        G4double pyf = coor_.py[i] + ffermi_.ffermpy[i];
        G4double pzf = coor_.pz[i] + ffermi_.ffermpz[i];
        G4double px = pxf * GeV;
        // units are in MeV/c for G4
        G4double py = pyf * GeV;
        G4double pz = pzf * GeV;

        G4double et = std::sqrt(pxf * pxf + pyf * pyf + pzf * pzf
                                + coor_.fmass[i] * coor_.fmass[i])
                      * GeV;

        //    ------------------------------Use only "Lorentz vector"----------

        G4LorentzVector lorenzvec = G4LorentzVector(px, py, pz, et);

        cascadeParticle = new G4DynamicParticle(pd, lorenzvec);  //

        theResult.AddSecondary(cascadeParticle);

        //======================================================================

      }  // if
    }  // for

  }  // if warning

  //=======================================================================
  if (verbose >= 3) {
    //
    G4double TotalEafter = 0.0;
    G4ThreeVector TotalPafter;
    G4double charge = 0.0;
    G4int baryon = 0;
    G4int nSecondaries = theResult.GetNumberOfSecondaries();

    for (G4int j = 0; j < nSecondaries; j++) {
      TotalEafter += theResult.GetSecondary(j)->GetParticle()->GetTotalEnergy();

      TotalPafter += theResult.GetSecondary(j)->GetParticle()->GetMomentum();

      G4ParticleDefinition* pd = theResult.GetSecondary(j)->GetParticle()->GetDefinition();

      charge += pd->GetPDGCharge();
      baryon += pd->GetBaryonNumber();

    }  // for secondaries

    G4cout << "----------------------------------------"
           << "----------------------------------------" << G4endl;
    G4cout << "Total energy before collision  = " << TotalEbefore  /// MeV
           << " MeV" << G4endl;
    G4cout << "Total energy after collision    = " << TotalEafter  // MeV
           << " MeV" << G4endl;

    G4cout << "----------------------------------------" << G4endl;

    G4cout << "Total momentum before collision = " << Pbefore  // MeV
           << " MeV/c" << G4endl;
    G4cout << "Total momentum after collision  = " << TotalPafter  // MeV
           << " MeV/c" << G4endl;
    G4cout << "----------------------------------------" << G4endl;

    if (verbose >= 4) {
      G4cout << "Total charge before collision  = " << (ZP + ZT) * eplus << G4endl;
      G4cout << "Total charge after collision    = " << charge << G4endl;

      G4cout << "----------------------------------------" << G4endl;

      G4cout << "Total baryon number before collision = " << AP + AT << G4endl;
      G4cout << "Total baryon number after collision  = " << baryon << G4endl;
      G4cout << "----------------------------------------" << G4endl;

    }  // if verbose4

    G4cout << "----------------------------------------"
           << "----------------------------------------" << G4endl;

  }  // if verbose3

  return &theResult;
}  // G4hadfinal

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// WelcomeMessage
//
void G4UrQMD1_3Model::WelcomeMessage() const
{
  G4cout << G4endl;
  G4cout << " *****************************************************************" << G4endl;
  G4cout << " Interface to        G4UrQMD_1.3                      activated" << G4endl;
  G4cout << " Version number : 00.00.0B          File date : 25/01/12" << G4endl;
  G4cout << " (Interface written by Kh. Abdel-Waged et al. for the KACST/NCMP)" << G4endl;
  G4cout << G4endl;
  G4cout << " *****************************************************************" << G4endl;
  G4cout << G4endl;

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4UrQMD1_3Model::InitialiseDataTables()
{
  //
  //
  // The next line is to make sure the block data statements are
  // executed.
  //

  g4urqmdblockdata_();

  // Seed the UrQMD random generator from Geant4's random engine, so that UrQMD
  // results depend on Geant4's seed (/random/setSeeds) and are reproducible.
  // The old code used a hardcoded constant (1097569630), which ignored Geant4's
  // engine and made every run identical; the commented-out -time_() variant
  // seeded from the wall clock, which is not reproducible either. sseed() needs
  // a positive seed (it negates it internally to prime the ran2 generator).
  G4int ranseed = 1 + static_cast<G4int>(G4UniformRand() * 2147483646.0);

  G4cout << "\n UrQMD seed (from Geant4 engine):  " << ranseed << G4endl;

  sseed_(&ranseed);

  loginit_();
}

#endif  // G4_USE_URQMD
