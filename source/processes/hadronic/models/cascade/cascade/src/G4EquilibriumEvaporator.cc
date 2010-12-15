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
// $Id: G4EquilibriumEvaporator.cc,v 1.49 2010-12-15 07:41:03 gunter Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100308  M. Kelsey -- Bug fix for setting masses of evaporating nuclei
// 20100319  M. Kelsey -- Use new generateWithRandomAngles for theta,phi stuff;
//		eliminate some unnecessary std::pow()
// 20100319  M. Kelsey -- Bug fix in new GetBindingEnergy() use right after
//		goodRemnant() -- Q1 should be outside call.
// 20100412  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100413  M. Kelsey -- Pass buffers to paraMaker[Truncated]
// 20100419  M. Kelsey -- Handle fission output list via const-ref
// 20100517  M. Kelsey -- Use G4CascadeInterpolator for QFF
// 20100520  M. Kelsey -- Inherit from common base class, make other colliders
//		simple data members.  Rename timeToBigBang() to override
//		base explosion().
// 20100617  M. Kelsey -- Remove "RUN" preprocessor flag and all "#else" code,
//		pass verbosity to colliders.
// 20100620  M. Kelsey -- Use local "bindingEnergy()" function to call through.
// 20100701  M. Kelsey -- Don't need to add excitation to nuclear mass; compute
//		new excitation energies properly (mass differences)
// 20100702  M. Kelsey -- Simplify if-cascades, indentation
// 20100712  M. Kelsey -- Add conservation checking
// 20100714  M. Kelsey -- Move conservation checking to base class.  Use
//		_generated_ evaporate energy (S) to adjust EEXS directly,
//		and test for S < EEXS before any particle generation; shift
//		nucleus momentum (PEX) by evaporate momentum directly
// 20100719  M. Kelsey -- Remove duplicative EESX_new calculation.
// 20100923  M. Kelsey -- Migrate to integer A and Z

#include "G4EquilibriumEvaporator.hh"
#include "G4BigBanger.hh"
#include "G4CascadeInterpolator.hh"
#include "G4CollisionOutput.hh"
#include "G4Fissioner.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4LorentzConvertor.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"

using namespace G4InuclSpecialFunctions;


G4EquilibriumEvaporator::G4EquilibriumEvaporator()
  : G4CascadeColliderBase("G4EquilibriumEvaporator") {}

G4EquilibriumEvaporator::~G4EquilibriumEvaporator() {}


void G4EquilibriumEvaporator::collide(G4InuclParticle* /*bullet*/,
				      G4InuclParticle* target,
				      G4CollisionOutput& output) {
  if (verboseLevel) {
    G4cout << " >>> G4EquilibriumEvaporator::collide" << G4endl;
  }

  // Sanity check
  G4InuclNuclei* nuclei_target = dynamic_cast<G4InuclNuclei*>(target);
  if (!nuclei_target) {
    G4cerr << " EquilibriumEvaporator -> target is not nuclei " << G4endl;    
    return;
  }

  if (verboseLevel > 1) {
    G4cout << " evaporating target: " << G4endl;
    target->printParticle();
  }

  theFissioner.setVerboseLevel(verboseLevel);
  theBigBanger.setVerboseLevel(verboseLevel);

  // simple implementation of the equilibium evaporation a la Dostrowski
  const G4double huge_num = 50.0;
  const G4double small = -50.0;
  const G4double prob_cut_off = 1.0e-15;
  const G4double Q1[6] = { 0.0, 0.0, 2.23, 8.49, 7.72, 28.3 };
  const G4int AN[6] = { 1, 1, 2, 3, 3, 4 };
  const G4int Q[6] =  { 0, 1, 1, 1, 2, 2 };
  const G4double G[6] = { 2.0, 2.0, 6.0, 6.0, 6.0, 4.0 };
  const G4double BE = 0.0063;
  const G4double fisssion_cut = 1000.0;
  const G4double cut_off_energy = 0.1;

  const G4double BF = 0.0242;
  const G4int itry_max = 1000;
  const G4int itry_global_max = 1000;
  const G4double small_ekin = 1.0e-6;
  const G4int itry_gam_max = 100;

  G4double W[8], u[6], V[6], TM[6];
  G4int A1[6], Z1[6];

  G4int A = nuclei_target->getA();
  G4int Z = nuclei_target->getZ();
  G4LorentzVector PEX = nuclei_target->getMomentum();
  G4double EEXS = nuclei_target->getExitationEnergy();
  
  if (verboseLevel > 3) G4cout << " after noeq: eexs " << EEXS << G4endl;

  G4InuclElementaryParticle dummy(small_ekin, 1);
  G4LorentzConvertor toTheNucleiSystemRestFrame;
  //*** toTheNucleiSystemRestFrame.setVerbose(verboseLevel);
  toTheNucleiSystemRestFrame.setBullet(dummy);

  G4LorentzVector ppout;
  
  // See if fragment should just be dispersed
  if (explosion(A, Z, EEXS)) {
    if (verboseLevel > 1) G4cout << " big bang in eql start " << G4endl;
    theBigBanger.collide(0, target, output);

    validateOutput(0, target, output);		// Check energy conservation
    return;
  }

  // If nucleus is in ground state, no evaporation
  if (EEXS < cut_off_energy) {
    if (verboseLevel > 1) G4cout << " no energy for evaporation" << G4endl;
    output.addOutgoingNucleus(*nuclei_target);

    validateOutput(0, target, output);		// Check energy conservation
    return;
  }

  // Initialize evaporation attempts
  G4double coul_coeff = (A >= 100.0) ? 1.4 : 1.2;
   
  G4LorentzVector pin = PEX;	// Save original target for testing
    
  G4bool try_again = true;  
  G4bool fission_open = true;
  G4int itry_global = 0;
    
  // Buffer for parameter sets
  std::pair<std::vector<G4double>, std::vector<G4double> > parms;
    
  while (try_again && itry_global < itry_global_max) {
    itry_global++;

    // Set rest frame of current (recoiling) nucleus
    toTheNucleiSystemRestFrame.setTarget(PEX);
    toTheNucleiSystemRestFrame.toTheTargetRestFrame();
      
    if (verboseLevel > 2) {
      G4double nuc_mass = G4InuclNuclei::getNucleiMass(A, Z, EEXS);
      G4cout << " A " << A << " Z " << Z << " mass " << nuc_mass
	     << " EEXS " << EEXS << G4endl;
    }
      
    if (explosion(A, Z, EEXS)) { 			// big bang
      if (verboseLevel > 2) 
	G4cout << " big bang in eql step " << itry_global << G4endl;
	
      G4InuclNuclei nuclei(PEX, A, Z, EEXS, 6);        
      theBigBanger.collide(0, &nuclei, output);

      validateOutput(0, target, output);	// Check energy conservation
      return;	
    } 

    if (EEXS < cut_off_energy) {	// Evaporation not possible
      if (verboseLevel > 2)
	G4cout << " no energy for evaporation in eql step " << itry_global 
	       << G4endl;

      try_again = false;
      break;
    }

    // Normal evaporation chain
    G4double E0 = getE0(A);
    G4double parlev = getPARLEVDEN(A, Z);
    G4double u1 = parlev * A;

    paraMaker(Z, parms);
    const std::vector<G4double>& AK = parms.first;
    const std::vector<G4double>& CPA = parms.second;

    G4double DM0 = bindingEnergy(A,Z);
    G4int i(0);

    for (i = 0; i < 6; i++) {
      A1[i] = A - AN[i];
      Z1[i] = Z - Q[i];
      u[i] = parlev * A1[i];
      TM[i] = -0.1;

      if (goodRemnant(A1[i], Z1[i])) {
	G4double QB = DM0 - bindingEnergy(A1[i],Z1[i]) - Q1[i];
	V[i] = coul_coeff * Z * Q[i] * AK[i] / (1.0 + EEXS / E0) /
	  (G4cbrt(A1[i]) + G4cbrt(AN[i]));
	TM[i] = EEXS - QB - V[i] * A / A1[i];  
      };
    }; 
      
    G4double ue = 2.0 * std::sqrt(u1 * EEXS);
    G4double prob_sum = 0.0;

    W[0] = 0.0;
    if (TM[0] > cut_off_energy) {
      G4double AL = getAL(A);
      W[0] = BE * G4cbrt(A1[0]*A1[0]) * G[0] * AL;
      G4double TM1 = 2.0 * std::sqrt(u[0] * TM[0]) - ue;

      if (TM1 > huge_num) TM1 = huge_num;
      else if (TM1 < small) TM1 = small;

      W[0] *= std::exp(TM1);
      prob_sum += W[0];
    }
      
    for (i = 1; i < 6; i++) {
      W[i] = 0.0;
      if (TM[i] > cut_off_energy) {
	W[i] = BE * G4cbrt(A1[i]*A1[i]) * G[i] * (1.0 + CPA[i]);
	G4double TM1 = 2.0 * std::sqrt(u[i] * TM[i]) - ue;

	if (TM1 > huge_num) TM1 = huge_num;
	else if (TM1 < small) TM1 = small;

	W[i] *= std::exp(TM1);
	prob_sum += W[i];
      }
    }

    // fisson part
    W[6] = 0.0;
    if (A >= 100.0 && fission_open) {
      G4double X2 = Z * Z / A;
      G4double X1 = 1.0 - 2.0 * Z / A; 
      G4double X = 0.019316 * X2 / (1.0 - 1.79 * X1 * X1);
      G4double EF = EEXS - getQF(X, X2, A, Z, EEXS);
	  
      if (EF > 0.0) {
	G4double AF = u1 * getAF(X, A, Z, EEXS);
	G4double TM1 = 2.0 * std::sqrt(AF * EF) - ue;

	if (TM1 > huge_num) TM1 = huge_num;
	else if (TM1 < small) TM1 = small;

	W[6] = BF * std::exp(TM1);
	if (W[6] > fisssion_cut*W[0]) W[6] = fisssion_cut*W[0]; 	     

	prob_sum += W[6];
      }
    } 

    // again time to decide what next
    if (verboseLevel > 2){
      G4cout << " wn " << W[0] << " wp " << W[1] << " wd " << W[2] << G4endl
	     << " wh3 " << W[3] << " wt " << W[4] << " whe4 " << W[5] << G4endl
	     << " wfi " << W[6] << G4endl;
    }

    G4int icase = -1;

    if (prob_sum < prob_cut_off) { 		// photon emission chain
      G4double UCR0 = 2.5 + 150.0 / A;
      G4double T00 = 1.0 / (std::sqrt(u1 / UCR0) - 1.25 / UCR0);
      G4int itry_gam = 0;

      while (EEXS > cut_off_energy && try_again) {
	itry_gam++;
	G4int itry = 0;
	G4double T04 = 4.0 * T00;
	G4double FMAX;

	if (T04 < EEXS) {
	  FMAX = (T04*T04*T04*T04) * std::exp((EEXS - T04) / T00);
	} else {
	  FMAX = EEXS*EEXS*EEXS*EEXS;
	}; 

	G4double S(0);
	while (itry < itry_max) {
	  itry++;
	  S = EEXS * inuclRndm();
	  G4double X1 = (S*S*S*S) * std::exp((EEXS - S) / T00);

	  if (X1 > FMAX * inuclRndm()) break;
	};

	if (itry == itry_max) {		// Maximum attempts exceeded
	  try_again = false;
	  break;
	}

	if (verboseLevel > 2) G4cout << " photon escape ?" << G4endl;

	if (S < EEXS) {		// Valid evaporate
	  S /= GeV;				// Convert to Bertini units

	  G4double pmod = S;
	  G4LorentzVector mom = generateWithRandomAngles(pmod, 0.);

	  // Push evaporated particle into current rest frame
	  mom = toTheNucleiSystemRestFrame.backToTheLab(mom);

	  if (verboseLevel > 3) {
	    G4cout << " nucleus   px " << PEX.px() << " py " << PEX.py()
		   << " pz " << PEX.pz() << " E " << PEX.e() << G4endl
		   << " evaporate px " << mom.px() << " py " << mom.py()
		   << " pz " << mom.pz() << " E " << mom.e() << G4endl;
	  }

	  PEX -= mom;			// Remaining four-momentum
	  EEXS -= S*GeV;		// New excitation energy (in MeV)

	  G4InuclElementaryParticle particle(mom, 10, 6);
	  output.addOutgoingParticle(particle);
	  
	  if (verboseLevel > 3) particle.printParticle();
	  
	  ppout += mom;
	} else {
	  if (itry_gam == itry_gam_max) try_again = false;
	}
      }		// while (EEXS > cut_off
      try_again = false;
    } else {			// if (prob_sum < prob_cut_off)
      G4double SL = prob_sum * inuclRndm();
      G4double S1 = 0.0;

      for (G4int i = 0; i < 7; i++) {	// Select evaporation scenario
	S1 += W[i];
	if (SL <= S1) {
	  icase = i;
	  break;
	};
      };

      if (icase < 6) { // particle or light nuclei escape
	if (verboseLevel > 2) 
	  G4cout << " particle/light-ion escape ?" << G4endl;

	G4double uc = 2.0 * std::sqrt(u[icase] * TM[icase]);
	G4double ur = (uc > huge_num ? std::exp(huge_num) : std::exp(uc));
	G4double d1 = 1.0 / ur;
	G4double d2 = 1.0 / (ur - 1.0);	    
	G4int itry1 = 0;
	G4bool bad = true;

	while (itry1 < itry_max && bad) {
	  itry1++; 
	  G4int itry = 0;
	  G4double EPR = -1.0;
	  G4double S = 0.0;

	  while (itry < itry_max && EPR < 0.0) {
	    itry++;
	    G4double uu = uc + std::log((1.0 - d1) * inuclRndm() + d2);
	    S = 0.5 * (uc * uc - uu * uu) / u[icase];
	    EPR = TM[icase] - S * A / (A - 1.0) + V[icase];
	  }; 
	    
	  if (EPR > 0.0 && S > V[icase]) { // real escape
	    if (verboseLevel > 2)
	      G4cout << " escape itry1 " << itry1 << " icase "
		     << icase << " S (MeV) " << S << G4endl;

	    if (S < EEXS) {		// Valid evaporate
	      S /= GeV;				// Convert to Bertini units

	      if (icase < 2) { 	// particle escape
		G4int ptype = 2 - icase;
		if (verboseLevel > 2)
		  G4cout << " particle " << ptype << " escape" << G4endl;
		
		G4InuclElementaryParticle particle(ptype);
		particle.setModel(6);
		
		// generate particle momentum
		G4double mass = particle.getMass();
		G4double pmod = std::sqrt((2.0 * mass + S) * S);
		G4LorentzVector mom = generateWithRandomAngles(pmod, mass);

		// Push evaporated particle into current rest frame
		mom = toTheNucleiSystemRestFrame.backToTheLab(mom);

		if (verboseLevel > 2) {
		  G4cout << " nucleus   px " << PEX.px() << " py " << PEX.py()
			 << " pz " << PEX.pz() << " E " << PEX.e() << G4endl
			 << " evaporate px " << mom.px() << " py " << mom.py()
			 << " pz " << mom.pz() << " E " << mom.e() << G4endl;
		}
		
		// New excitation energy depends on residual nuclear state
		G4double mass_new =
		  G4InuclNuclei::getNucleiMass(A1[icase],Z1[icase]);

		G4double EEXS_new = ((PEX-mom).m() - mass_new)*GeV;
		if (EEXS_new < 0.0) continue;	// Sanity check for new nucleus

		PEX -= mom;		// Remaining four-momentum
		EEXS = EEXS_new;

		A = A1[icase];
		Z = Z1[icase]; 	      

		particle.setMomentum(mom);
		output.addOutgoingParticle(particle);
		if (verboseLevel > 3) particle.printParticle();

		ppout += mom;
		bad = false;
	      } else {	// if (icase < 2)
		if (verboseLevel > 2) {
		  G4cout << " nucleus A " << AN[icase] << " Z " << Q[icase]
			 << " escape icase " << icase << G4endl;
		}
		
		G4InuclNuclei nuclei(AN[icase], Q[icase]);
		nuclei.setModel(6);
		G4double mass = nuclei.getMass();
		// generate particle momentum
		G4double pmod = std::sqrt((2.0 * mass + S) * S);
		G4LorentzVector mom = generateWithRandomAngles(pmod,mass);
		
		// Push evaporated particle into current rest frame
		mom = toTheNucleiSystemRestFrame.backToTheLab(mom);

		if (verboseLevel > 2) {
		  G4cout << " nucleus   px " << PEX.px() << " py " << PEX.py()
			 << " pz " << PEX.pz() << " E " << PEX.e() << G4endl
			 << " evaporate px " << mom.px() << " py " << mom.py()
			 << " pz " << mom.pz() << " E " << mom.e() << G4endl;
		}
		
		// New excitation energy depends on residual nuclear state
		G4double mass_new =
		  G4InuclNuclei::getNucleiMass(A1[icase],Z1[icase]);

		G4double EEXS_new = ((PEX-mom).m() - mass_new)*GeV;
		if (EEXS_new < 0.0) continue;	// Sanity check for new nucleus

		PEX -= mom;		// Remaining four-momentum
		EEXS = EEXS_new;

		A = A1[icase];
		Z = Z1[icase];

		nuclei.setMomentum(mom);
		output.addOutgoingNucleus(nuclei);
		if (verboseLevel > 3) nuclei.printParticle();

		ppout += mom;
		bad = false;
	      }		// if (icase < 2)
	    }		// if (S < EEXS)
	  }		// if (EPR > 0.0 ...
	}		// while (itry1 ...

	if (itry1 == itry_max || bad) try_again = false;
      } else { 	// if (icase < 6)
	G4InuclNuclei nuclei(A, Z, EEXS);        
	nuclei.setModel(6);

	if (verboseLevel > 2) {
	  G4cout << " fission: A " << A << " Z " << Z << " eexs " << EEXS
		 << " Wn " << W[0] << " Wf " << W[6] << G4endl;
	}

	// Catch fission output separately for verification
	G4CollisionOutput foutput;
	theFissioner.collide(0, &nuclei, foutput);

	if (foutput.getOutgoingNuclei().size() == 2) { // fission o'k
	  if (verboseLevel > 2) G4cout << " fission done in eql" << G4endl;

	  // Move fission fragments to lab frame for processing
	  foutput.boostToLabFrame(toTheNucleiSystemRestFrame);

	  // Now evaporate the fission fragments individually
	  G4bool prevDoChecks = doConservationChecks;	// Turn off checking
	  setConservationChecks(false);

	  std::vector<G4InuclNuclei> nuclea = foutput.getOutgoingNuclei();
	  this->collide(0, &nuclea[0], output);
	  this->collide(0, &nuclea[1], output);

	  setConservationChecks(prevDoChecks);	// Restore previous flag value
	  validateOutput(0, target, output);	// Check energy conservation
	  return;
	} else { // fission forbidden now
	  fission_open = false;
	}
      }		// End of fission case
    }		// if (prob_sum < prob_cut_off)
  }		// while (try_again

  // this time it's final nuclei

  if (itry_global == itry_global_max) {
    if (verboseLevel > 3) {
      G4cout << " ! itry_global " << itry_global_max << G4endl;
    }
  }

  G4LorentzVector pnuc = pin - ppout;

  G4InuclNuclei nuclei(pnuc, A, Z, EEXS, 6);

  /***** THIS SHOULD NOT BE NECESSARY IF EEXS WAS COMPUTED RIGHT
  pnuc = nuclei.getMomentum(); 
  G4double eout = pnuc.e() + ppout.e();  
  G4double eex_real = 1000.0 * (pin.e() - eout);        
  nuclei.setExitationEnergy(eex_real);
  *****/

  if (verboseLevel > 3) {
    G4cout << " remaining nucleus " << G4endl;
    nuclei.printParticle();
  }

  output.addOutgoingNucleus(nuclei);

  validateOutput(0, target, output);		// Check energy conservation
  return;
}		     

G4bool G4EquilibriumEvaporator::explosion(G4int a, 
					  G4int z, 
					  G4double e) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4EquilibriumEvaporator::explosion" << G4endl;
  }

  const G4double be_cut = 3.0;

  // Different criteria from base class, since nucleus more "agitated"
  G4bool bigb = (!(a >= 12 && z >= 0 && z < 3*(a-z)) &&
		 (e >= be_cut * bindingEnergy(a,z))
		 );

  return bigb;
}

G4bool G4EquilibriumEvaporator::goodRemnant(G4int a, 
					    G4int z) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4EquilibriumEvaporator::goodRemnant" << G4endl;
  }

  return a > 1 && z > 0 && a > z;
}

G4double G4EquilibriumEvaporator::getQF(G4double x, 
					G4double x2, 
					G4int a,
					G4int /*z*/, 
					G4double ) const {
  if (verboseLevel > 3) {
    G4cout << " >>> G4EquilibriumEvaporator::getQF" << G4endl;
  }
  
  static const G4double QFREP[72] = {  
    //     TL201 *     *   *    *
    //      1    2     3   4    5
    22.5, 22.0, 21.0, 21.0, 20.0,
    //     BI209 BI207 PO210 AT213 *    TH234
    //      6     7    8     9     10   11
    20.6, 20.6, 18.6, 15.8, 13.5, 6.5,
    //     TH233 TH232 TH231 TH230 TX229 PA233 PA232 PA231 PA230 U240
    //     12    13    14    15    16    17    18    19    20    21
    6.65, 6.22, 6.27, 6.5,  6.7,  6.2,  6.25, 5.9,  6.1,  5.75,
    //     U239 U238 U237  U236 U235 U234 U233 U232 U231
    //     22   23   24    25   26   27   28   29   30
    6.46, 5.7, 6.28, 5.8, 6.15, 5.6, 5.8, 5.2, 5.8,
    //     NP238 NP237 NP236 NP235 PU245 NP234  PU244 NP233
    //     31    32    33    34    35    36     37    38
    6.2 , 5.9 , 5.9,  6.0,  5.8,  5.7,   5.4,  5.4,
    //     PU242 PU241 PU240 PU239 PU238 AM247 PU236 AM245 AM244 AM243
    //     39    40    41    42    43    44    45    46    47    48
    5.6,  6.1,  5.57, 6.3,  5.5,  5.8,  4.7,  6.2,  6.4,  6.2,
    //     AM242 AM241 AM240 CM250 AM239 CM249 CM248 CM247 CM246
    //     49    50    51    52    53    54    55    56    57
    6.5,  6.2,  6.5,  5.3,  6.4,  5.7,  5.7,  6.2,  5.7,
    //     CM245 CM244 CM243 CM242 CM241 BK250 CM240
    //     58    59    60    61    62    63    64
    6.3,  5.8,  6.7,  5.8,  6.6,  6.1,  4.3,
    //     BK249 CF252 CF250 CF248 CF246 ES254 ES253 FM254
    //     65    66    67    68    69    70    71    72
    6.2,  3.8,  5.6,  4.0,  4.0,  4.2,  4.2,  3.5 };
     
  static const G4double XREP[72] = {
    //      1      2     3      4      5
    0.6761, 0.677, 0.6788, 0.6803, 0.685,
    //      6     7     8     9     10     11
    0.6889, 0.6914, 0.6991, 0.7068, 0.725, 0.7391,
    //     12  13    14   15   16    17  18    19    20    21
    0.74, 0.741, 0.742, 0.743, 0.744, 0.7509, 0.752, 0.7531, 0.7543, 0.7548,
    //     22    23    24
    0.7557, 0.7566, 0.7576,
    //      25     26   27    28    29   30   31    32     33    34
    0.7587, 0.7597, 0.7608, 0.762, 0.7632, 0.7644, 0.7675, 0.7686, 0.7697, 0.7709,
    //      35    36    37    38    39   40    41
    0.7714, 0.7721, 0.7723, 0.7733, 0.7743, 0.7753, 0.7764,
    //      42    43    44    45    46    47    48   49
    0.7775, 0.7786, 0.7801, 0.781, 0.7821, 0.7831, 0.7842, 0.7852,
    //     50     51    52    53    54    55    56    57    58
    0.7864, 0.7875, 0.7880, 0.7887, 0.7889, 0.7899, 0.7909, 0.7919, 0.7930,
    //      59    60    61    62    63    64
    0.7941, 0.7953, 0.7965, 0.7977, 0.7987, 0.7989,
    //      65    66    67    68    69    70    71    72
    0.7997, 0.8075, 0.8097, 0.8119, 0.8143, 0.8164, 0.8174, 0.8274 };

  const G4double G0 = 20.4;
  const G4double XMIN = 0.6761;
  const G4double XMAX = 0.8274;

  G4double QFF = 0.0;

  if (x < XMIN || x > XMAX) {
    G4double X1 = 1.0 - 0.02 * x2;
    G4double FX = (0.73 + (3.33 * X1 - 0.66) * X1) * (X1*X1*X1);

    QFF = G0 * FX * G4cbrt(a*a);
  } else {
    static G4CascadeInterpolator<72> interp(XREP);	// Only need one!
    QFF = interp.interpolate(x, QFREP);
  }

  if (QFF < 0.0) QFF = 0.0;

  return QFF; 
}

G4double G4EquilibriumEvaporator::getAF(G4double , 
					G4int /*a*/, 
					G4int /*z*/, 
					G4double e) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4EquilibriumEvaporator::getAF" << G4endl;
  }

  // ugly parameterisation to fit the experimental fission cs for Hg - Bi nuclei
  G4double AF = 1.285 * (1.0 - e / 1100.0);

  if(AF < 1.06) AF = 1.06;
  // if(AF < 1.08) AF = 1.08;

  return AF;
}	

G4double G4EquilibriumEvaporator::getPARLEVDEN(G4int /*a*/, 
					       G4int /*z*/) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4EquilibriumEvaporator::getPARLEVDEN" << G4endl;
  }

  const G4double par = 0.125;

  return par;
}

G4double G4EquilibriumEvaporator::getE0(G4int /*a*/) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4EquilibriumEvaporator::getE0" << G4endl;
  }

  const G4double e0 = 200.0;   

  return e0;   
}
