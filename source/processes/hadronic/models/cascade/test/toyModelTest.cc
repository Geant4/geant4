#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4Proton.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4Fancy3DNucleus.hh"
//#include "G4Cascade.hh"

int main() {
  G4cout << G4endl << "       ... intra-nuclear cascade toy model ..." << G4endl;
	 
  // Model based on: 
  // Methods in Subnuclear Physics, Vol IV, Part 3 pp. 137-144,
  // Proceedings of the International School of Elelmantary Particle Physics, 
  // Herceg-Novi, Yugoslavia, 1968   
  // Ed. M. Mikolic

  // Set constanst and create terminology
  G4double cutOff = 20 * MeV;
  G4double fm = 1.0 * pow(10.0, -15);
  G4double oneThird = 0.3333333333333;
  enum particleType {proton, neutron};
  // Eenum channelType {pionProduction, };
  enum collisionType {elastic, inelastic};

  vector<G4KineticTrack*>  particleVector;   

  // Set test parameters
  G4int        numberOfCascades       = 3;
  G4int        targetNucleusN         = 10;
  G4int        targetNucleusZ         = 10;
  particleType incidentParticleType   = proton;
  G4double     incidentParticleEnergy = 1000.0 * MeV;                         
  G4Fancy3DNucleus      theCarbon;
  G4double              theA = 12.0;
  G4double              theZ =  6.0;
  theCarbon.Init(theA, theZ);

  G4cout << G4endl << " Nucleus mass number = " << theCarbon.GetMassNumber() << G4endl
	 << " Nucleus mass = " << theCarbon.GetMass() << G4endl
	 << " Nucleus charge = " << theCarbon.GetCharge() << G4endl;
 
  G4cout << " projectile                 : " << incidentParticleType             << G4endl;
  G4cout << " projectile energy          : " << incidentParticleEnergy << " MeV" << G4endl;
  G4cout << " neutrons in target nucleus : " << targetNucleusN                   << G4endl;
  G4cout << " protons  in target nucleus : " << targetNucleusZ                   << G4endl;

  for (G4int cascadeIndex = 1; cascadeIndex <= numberOfCascades; cascadeIndex++){

    G4cout << G4endl << "........................ cascade " << cascadeIndex << 
                        " .............................." << G4endl;

    // Select the impact parameter b
    G4double targetNucleusA = targetNucleusN + targetNucleusZ;
    G4double radius0        = 1.0 * fm;                                       // radius of hydrogen atom nucleus
    G4double bMax           = radius0 * fm / pow(targetNucleusA, oneThird) ;  // radius of atom with a = n + z
    G4double b              = bMax * sqrt(G4UniformRand());
    G4double z              = b;                                              // x, y, z coordinates of impact point
    G4double x              = - sqrt(sqr(bMax) - sqr(z));
    G4double y              = 0;
    //    G4cout << radius0 << " " << bMax << endl;
    G4cout << "impact point (x, y, z) : "  << "\t" << x <<", " << "\t" << y << ", " << "\t" << z  << G4endl;

    // *) Find the interaction distance b
    G4double xSecNeutron = 1 * millibarn;
    G4double xSecProton  = 1 * millibarn;
    G4double lambda      = 4 * pi * targetNucleusA * pow(radius0, 3) /
      (3 * (targetNucleusZ * xSecProton + targetNucleusN * xSecNeutron));
    b = - lambda * log(G4UniformRand());

    G4cout << "interaction distance   : " << "\t" << b << G4endl;

    // Choose the sruck particle (n/p)
    particleType targetParticle;
    if (G4UniformRand() < targetNucleusZ * xSecProton / 
                        ( targetNucleusZ * xSecProton + targetNucleusN * xSecNeutron))
      targetParticle = proton;
    else 
      targetParticle = neutron;

    G4cout << "target particle        : " << "\t" << targetParticle << G4endl;

    // Choose the struck nucleon momentum
    G4double qMax = (300.0 / radius0) * pow(10, -15) * MeV;
    G4double q    = qMax * sqrt(G4UniformRand());

    G4cout << "target  momentum       : " << "\t" << q << G4endl;

    // Choose the direction of the struck nucleon
    G4double mu  = 1 - 2 * G4UniformRand();
    G4double phi = 2 * pi * G4UniformRand();
    G4cout << "target direction       : " << "\t" << mu << ", " << "\t" << phi << G4endl;

    // Choose the nature of the collision (elastic/inelastic)
    collisionType collision;
    if (G4UniformRand() < 0.8) // replace with more realistic 
      collision = elastic;
    else
      collision = inelastic;
    G4cout << "collision type         : " << "\t" << collision << G4endl;

    if (collision == elastic) {
      // Choose the scattering angle theta
      G4double thetaCos = 2 *G4UniformRand() - 1;                           // dummy
      G4cout << " scattering angle      : " << "\t" << thetaCos << G4endl;
    }

    if (collision == inelastic) {

      G4ParticleDefinition* theProtonDefinition = G4Proton::ProtonDefinition();
      G4double              theProtonFormationTime = 0.0*ns;
      G4ThreeVector         theProtonInitialPosition(0.0 * cm, 0.0 * cm, 0.0 * cm);
      G4LorentzVector       theProton4Momentum(0.0 * MeV, 0.0 * MeV, 1.0 * MeV, 938.27 * MeV);
      G4KineticTrack        aProton(theProtonDefinition, 
				    theProtonFormationTime,
				    theProtonInitialPosition,
				    theProton4Momentum);

      G4KineticTrack     *xx = new G4KineticTrack(theProtonDefinition, 
						  theProtonFormationTime,
						  theProtonInitialPosition,
						  theProton4Momentum);
                          
      G4Fancy3DNucleus      theCarbon;
      G4double              theA = 12.0;
      G4double              theZ =  6.0;
      theCarbon.Init(theA, theZ);
   
      particleVector.push_back(xx);
      G4cout << " pion production       : " << G4endl;

    }

    // Choose the azimuthal angle of rotation
    G4double alpha = 2 * pi * G4UniformRand();
    G4cout << " azimuthal angle       : " << "\t" << alpha << G4endl;

    // Transform the energy back to the laboratory system (lab)
    G4double energy =  50 * MeV* G4UniformRand();                            //dummy

    // For all items in particleVectos
    // Pauli exclusion principle
    G4double M = 0.5 ;                                                       //dummy
    G4double fermiEnergy = sqr(qMax)/(2 * M);
    if (energy <  fermiEnergy)  {                                            // state occupied
      G4cout << "pauli blocking " << G4endl;
      // Add article to traced list  
    }          

    // Transform the allowed collision products to the laboratory system

    // Follow the products until the escape the nucleus or energy is below cut-off
    if (energy < cutOff) {
      G4cout << "energy below cut-off   "<< G4endl;
      // Add energy to atom exitation
      // exitation += particle.energy 
      // update atom quantum numbers
      // if particle = proton, a += 1, z += 1 
      // if particle = neutron, a += 1, n += 1 
    }
  }
   
}    
/*
       ... intra-nuclear cascade toy model ...

 Nucleus mass number = 12
 Nucleus mass = 11174.9
 Nucleus charge = 6
 projectile                 : 0
 projectile energy          : 1000 MeV
 neutrons in target nucleus : 10
 protons  in target nucleus : 10

........................ cascade 1 ..............................
impact point (x, y, z) :        -2.82281e-31,   0,      2.36724e-31
interaction distance   :        3.54704e-20
target particle        :        1
target  momentum       :        152.393
target direction       :        -0.713471,      2.86607
collision type         :        0
 scattering angle      :        0.159015
 azimuthal angle       :        4.14197
pauli plocking

........................ cascade 2 ..............................
impact point (x, y, z) :        -2.82281e-31,   0,      2.36724e-31
interaction distance   :        3.54704e-20
target particle        :        1
target  momentum       :        152.393
target direction       :        -0.713471,      2.86607
collision type         :        0
 scattering angle      :        0.159015
 azimuthal angle       :        4.14197
pauli plocking
energy below cut-off  

*/




