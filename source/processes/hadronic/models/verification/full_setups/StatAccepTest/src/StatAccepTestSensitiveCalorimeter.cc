#include "StatAccepTestSensitiveCalorimeter.hh"
#include "StatAccepTestCalorimeterHit.hh"
#include "G4SDManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "StatAccepTestAnalysis.hh"


StatAccepTestSensitiveCalorimeter::
StatAccepTestSensitiveCalorimeter( const G4String name, const G4int numberOfReplicasIn )
  : G4VSensitiveDetector( name ), calCollection( 0 ), 
    numberOfReplicas( numberOfReplicasIn ) {
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    positionID.push_back( -1 );
  }
  G4String HCname;
  collectionName.insert( HCname="calCollection" );
}


StatAccepTestSensitiveCalorimeter::~StatAccepTestSensitiveCalorimeter() {}


void StatAccepTestSensitiveCalorimeter::Initialize( G4HCofThisEvent* ) {
  calCollection = 
    new StatAccepTestCalorimeterHitsCollection( SensitiveDetectorName,collectionName[0] );
  positionID.clear();
  for ( int i = 0; i < numberOfReplicas; i++ ) {
    positionID.push_back( -1 );
  }
}


G4bool StatAccepTestSensitiveCalorimeter::
ProcessHits( G4Step* aStep, G4TouchableHistory* ) {

  G4double edep = aStep->GetTotalEnergyDeposit() * aStep->GetTrack()->GetWeight();
  // Multiply the energy deposit with the weight of the track,
  // to allow the use of biasing.

  if ( edep < 0.001*eV ) return false;

  // -----------------------------------------------------------------
  //                        BIRKS' LAW
  //                        ==========
  //
  // 6-Apr-2006 : In the case of Scintillator as active medium, we can
  //              describe the quenching effects with the Birks' law, 
  //              using the expression and the coefficients taken from 
  //              the paper NIM 80 (1970) 239-244 for the organic 
  //              scintillator NE-102:
  //                                     S*dE/dr 
  //                    dL/dr = -----------------------------------
  //                               1 + C1*(dE/dr) + C2*(dE/dr)^2 
  //              with:
  //                    C1 = 1.29 x 10^-2  g*cm^-2*MeV^-1 
  //                    C2 = 9.59 x 10^-6  g^2*cm^-4*MeV^-2 
  //              These are the same values used by ATLAS TileCal 
  //              and CMS HCAL (and also the default in Geant3).
  //              You can try different values for these parameters,
  //              to have an idea on the uncertainties due to them,
  //              by uncommenting one of the lines below.
  //              To get the "dE/dr" that appears in the formula,
  //              which has the dimensions
  //                    [ dE/dr ] = MeV * cm^2 / g
  //              we have to divide the energy deposit in MeV by the
  //              product of the step length (in cm) and the density
  //              of the scintillator:
  //                    rho_scintillator = 1.032  g/cm3 .
  //              Of course, in our case we use only the denominator
  //              of the Birks' formula above, because we do not have
  //              digitization, i.e. we only have energy deposit but
  //              not conversion to photons.
  //              If you do not want to apply Birks' law, simply set
  //              isBirksOn = false.
  //
  // 26-Aug-2006: In the case of LAr as active medium, we can also
  //              take into consideration the electron recombination 
  //              in liquid argon with a Birks' law, using the 
  //              expression taken from the paper: NIM A 523 (2004) 275-286  
  //                                      Q0
  //                    Q  = A * -------------------------
  //                              1 + (k/E_field)*(dE/dx)  
  //              with:
  //                    A  adimensional constant (fitteed on the data);
  //                    k  [kV/cm] is constant (fitted on the data);
  //                    E_field  [kV/cm] is the electric field.
  //              In the case of ATLAS LAr :
  //                    E_field = 10 kV/cm
  //              so it is a much higher field than the one discussed
  //              in the paper ( 0.1 < E_field < 1.0 kV/cm ) so we 
  //              cannot use directly the values of A and k reported 
  //              there. However, they discuss also the results reported
  //              by a previous analysis (Scalettar et al., 
  //              Phys Rev. A 25 (1882) 2419), which was done for higher
  //              fields, up to 10 kV/cm. In this analysis, they discuss
  //              the recombination effect as a function of the electric
  //              field only, and not as a function of the particle
  //              stopping power dE/dx, using the phenomenological
  //              formula:             Q0
  //                         Q = ----------------
  //                              1 + k_E/E_field
  //              where  k_E = 0.53 kV/cm .
  //              Now, we could use our formula above, where dE/dx
  //              appears, by identifying:
  //                    A = 1 
  //                    k = k_E/(dE/dx)
  //              where as  dE/dx  we have to use the average stopping
  //              power for the case considered in that analysis, i.e.
  //              364 keV electrons (from a 113Sn conversion electron
  //              source). This correspond, in liquid argon:
  //                    dE/dx = 2.4 MeV/(g/cm^2)
  //              which gives:
  //                    k = 0.22 (kV/cm)((g/cm^2)/MeV)
  //              In conclusion, to study the electron recombination
  //              for ATLAS LAr we can use the same formula as for
  //              the quenching in scintillator, with the following
  //              coefficients:
  //                    C1 = k/E_field = 0.22/10,0 = 0.022 (g/cm^2)/MeV
  //                    C2 = 0.0
  //              and with the following density:
  //                    rho_lAr = 1.396  g/cm^3
  //
  //              Update (24-Oct-2006) : the coefficient C1 used 
  //              originally by P.Sala in ATLAS:  
  //                    C1 = 0.005  gr/(MeV*cm^2) 
  //
  // 14-Nov-2006: Birks should not be applied in the case of gamma
  //              energy depositions (which happens only for the 
  //              photoelectric process), because in this case the
  //              step length is related to the photoelectric cross
  //              section, and not to the distance in which the
  //              energy is actually deposited, that is what is
  //              needed in order to determine dE/dx which enters
  //              in the Birks' law.
  //              Similarly, for neutron energy depositions (which
  //              have been introduced in Geant4 8.1 as an effective
  //              way to describe the elastic nuclei recoils below
  //              a certain kinetic threshold, set by default to 
  //              100 keV), the step length is related to the neutron
  //              elastic cross-section, and not to the real ionization
  //              range of the recoiled nuclei, which should instead
  //              be considered for the dE/dx in the Birks' law.
  //              In the case of neutron energy depositions, the 
  //              most correct way to apply the Birks quench would
  //              be to eliminate them by setting the kinetic 
  //              threshold to 0.0 (variable  "edepLimit"  in the 
  //              file  "G4HadronElasticPhysics.cc"  in 
  //               "geant4/physics_lists/hadronic/Packaging/src/" ),
  //              so that the recoiled nuclei tracks are always 
  //              generated explicitly. This, of course, costs in
  //              term of CPU performance.
  //
  //--------------------------------------------------------------------
  bool isBirksOn = false;         //***LOOKHERE***
  bool isScintillator = true;     //              true = Quenching in Scintillator 
                                  //              false = Recombination in LAr.
  if ( isBirksOn ) {

    // --- Quenching in Scintillator ---
    double C1_scintillator = 1.29e-2;                    // [g*cm^-2*MeV^-1]  
    double C2_scintillator = 9.59e-6;                    // [g^2*cm^-4*MeV^-2] 
    // Uncomment one line below if you want to try different
    // values for the scintillator Birks' coefficient C1, with C2 = 0.
    //C1_scintillator = 1.31e-2; C2_scintillator = 0.0;  // Standard Birks law.
    //C1_scintillator = 8.5e-3;  C2_scintillator = 0.0;  // Zeus SCSN38, lower limit (-35%)
    //C1_scintillator = 1.59e-2; C2_scintillator = 0.0;  // Pilot B (same paper), upper limit (+21%)
    double rho_scintillator = 1.032;                     // [g/cm^3]

    // --- Recombination in LiquidArgon --
    double C1_lAr = 0.022;                // my estimation  [g*cm^-2*MeV^-1]   
    //double C1_lAr = 0.005;                // Fluka value    [g*cm^-2*MeV^-1]  
    double C2_lAr = 0.0;
    double rho_lAr = 1.396;                              // [g/cm^3]

    double C1, C2, rho;
    if ( isScintillator ) {
      C1 = C1_scintillator;
      C2 = C2_scintillator;
      rho = rho_scintillator;
    } else {
      C1 = C1_lAr;
      C2 = C2_lAr;
      rho = rho_lAr;
    }

    static bool isFirstCall = true;
    if ( isFirstCall ) {
      G4cout << "\t *** StatAccepTestSensitiveCalorimeter::ProcessHits *** " << G4endl
	     << "\t isBirksOn = " << isBirksOn << G4endl
	     << "\t isScintillator = " << isScintillator << G4endl
	     << "\t C1 = " << C1 << G4endl
             << "\t C2 = " << C2 << G4endl
	     << "\t rho = " << rho << G4endl
             << "\t ******************* " << G4endl; 
      isFirstCall = false;
    }

    double stepLength_cm = aStep->GetStepLength() / cm ;  // [cm] 

    // Do not apply Birks for gamma deposits (corresponding to 
    // the photoelectric process). 
    if ( stepLength_cm > 1.0e-8 ) { 
      double dedx = edep / ( rho * stepLength_cm );       // [MeV*cm^2/g]
      double birksFactor = 1.0 / ( 1.0 + C1*dedx + C2*dedx*dedx );
      //G4cout << " birksFactor=" << birksFactor << G4endl; //***DEBUG***

      if ( aStep->GetTrack()->GetDefinition() == G4Gamma::GammaDefinition() ) { 
	birksFactor = 1.0;
	static bool isFirstWarningGamma = true;
	if ( isFirstWarningGamma ) {
	  G4cout << " ***BIRKS WARNING*** : BIRKS NOT APPLIED TO GAMMA ENERGY DEPOSITIONS!" << G4endl;
	  isFirstWarningGamma = false;	  
        }
      } 

      static bool isFirstWarningNeutron = true;
      if ( isFirstWarningNeutron &&
           aStep->GetTrack()->GetDefinition() == G4Neutron::NeutronDefinition() ) {
	G4cout << " ***BIRKS WARNING*** : BIRKS INCORRECTLY APPLIED TO NEUTRON ENERGY DEPOSITIONS!" << G4endl
	       << "                       BETTER TO SET edepLimit = 0.0 " << G4endl;
	isFirstWarningNeutron = false;	  
      }

      //if ( birksFactor > 0.0  &&  1.0/birksFactor > 100.0 ) {  //***DEBUG***
      //	G4cout << "Birks quenching = " << 1.0/birksFactor << "\t"
      //	       << aStep->GetTrack()->GetDefinition()->GetParticleName() 
      //	       << "\t stepLength=" << aStep->GetStepLength() << " mm" 
      //	       << G4endl
      //	       << "\t Ekin=" << aStep->GetTrack()->GetKineticEnergy() 
      //	       << " MeV ; edep=" << aStep->GetTotalEnergyDeposit() 
      //	       << " MeV ; process=";
      //	if ( aStep->GetPostStepPoint()->GetProcessDefinedStep() ) {
      //	  G4cout << aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
      //	}
      //  G4cout << G4endl; 
      //}

      edep = edep * birksFactor;

    }

  }
  // -----------------------------------------------------------------

  G4int replica = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber( 1 );

  G4double radius = std::sqrt( aStep->GetPreStepPoint()->GetPosition().x() *
			       aStep->GetPreStepPoint()->GetPosition().x() +
			       aStep->GetPreStepPoint()->GetPosition().y() *
			       aStep->GetPreStepPoint()->GetPosition().y() );  

  G4int particlePDG = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();

  //G4cout << " StatAccepTestSensitiveCalorimeter::ProcessHits : DEBUG Info" << G4endl
  //       << "\t z = " << aStep->GetPreStepPoint()->GetPosition().z()/mm 
  //       << " mm ;  replica = " << replica  << " ;  radius = " << radius/mm 
  //       << " mm ;  particlePDG - " << particlePDG << G4endl; //***DEBUG***

  if ( replica >= numberOfReplicas ) {
    G4cout << " StatAccepTestSensitiveCalorimeter::ProcessHits: ***ERROR*** "
           << " \t replica=" << replica 
	   << "  >=  numberOfReplicas=" << numberOfReplicas << G4endl;
    replica = numberOfReplicas - 1;  // To avoid a crash
  }

  if ( positionID[ replica ] == -1 ) {
    StatAccepTestCalorimeterHit* calHit = new StatAccepTestCalorimeterHit();
    calHit->SetEdep( edep );
    calHit->SetLayer( replica );

    G4int iPos = calCollection->insert( calHit );
    positionID[ replica ] = iPos - 1;
    //G4cout << " StatAccepTestSensitiveCalorimeter::ProcessHits: DEBUG Info " << G4endl
    //       << "\t New Calorimeter Hit : Energy = " 
    //       << edep / MeV << " MeV " << G4endl; //***DEBUG*** 
  } else { 
    (*calCollection)[ positionID[ replica ] ]->AddEdep( edep );
    //G4cout << " StatAccepTestSensitiveCalorimeter::ProcessHits : DEBUG Info " 
    //       << G4endl << "\t Old Calorimter Hit : Added energy = " 
    //       << " \t Added energy = " << edep / MeV << " MeV " << G4endl; //***DEBUG*** 
  }

  // We need to fill the shower profile information at this point,
  // because in the EventAction only the hit collection is available
  // and the hits do not have information about the radius of each
  // energy deposition that contributed (and also about the PDG code
  // of the particle that made such energy deposition).
  StatAccepTestAnalysis* analysis = StatAccepTestAnalysis::getInstance();
  if ( analysis ) {
    analysis->fillShowerProfile( replica, radius, edep, particlePDG );
  }

  return true;
}


void StatAccepTestSensitiveCalorimeter::EndOfEvent( G4HCofThisEvent* HCE ) {
  static G4int HCID = -1;
  if ( HCID < 0 ) { 
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID( collectionName[0] );
  }
  HCE->AddHitsCollection( HCID, calCollection );
  if ( ! calCollection->entries() ) {
    // G4cout << " StatAccepTestSensitiveCalorimeter::EndOfEvent : DEBUG Info " << G4endl
    //        << "\t NO calorimeter hit! " << G4endl; //***DEBUG***
  }   
}


void StatAccepTestSensitiveCalorimeter::clear() {} 


void StatAccepTestSensitiveCalorimeter::DrawAll() {} 


void StatAccepTestSensitiveCalorimeter::PrintAll() {} 

