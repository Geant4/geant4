// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// --------------------------------------------------------------
//      GEANT 4 class implementation file --- Copyright CERN 1998
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      2nd December 1995, G.Cosmo
//      ------------ G4MuonMinusCaptureAtRest physics process --------
//                   by Vladimir Ivanchenko
//                     E-mail: Vladimir.Ivantchenko@cern.ch
//                            April 2000
// **************************************************************
//-----------------------------------------------------------------------------

#include "G4MuMinusCaptureCascade.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// constructor
G4MuMinusCaptureCascade::G4MuMinusCaptureCascade()
{ 
  theElectron = G4Electron::Electron();
  theGamma = G4Gamma::Gamma();
  Emass = theElectron->GetPDGMass();
  MuMass = G4MuonMinus::MuonMinus()->GetPDGMass();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// destructor
G4MuMinusCaptureCascade::~G4MuMinusCaptureCascade()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuMinusCaptureCascade::GetKShellEnergy(G4double Z)
{
 
  // Calculate the Energy of K Mesoatom Level for this Element using
  // the Energy of Hydrogen Atom taken into account finite size of the
  // nucleus (V.Ivanchenko)
   const size_t ListK = 27;
   const G4double ListZK[ListK] = {
      2.,  4.,  6.,  8., 11., 14., 17., 18., 21., 24.,
     26., 29., 32., 38., 40., 41., 44., 49., 53., 55.,
     60., 65., 70., 75., 81., 85., 92.};
   const G4double ListKEnergy[ListK] = {
     0.011, 0.043, 0.098, 0.173, 0.326,
     0.524, 0.765, 0.853, 1.146, 1.472,
     1.708, 2.081, 2.475, 3.323, 3.627, 
     3.779, 4.237, 5.016, 5.647, 5.966,
     6.793, 7.602, 8.421, 9.249, 10.222,
    10.923,11.984};

  // Energy with finit size corrections
   G4double KEnergy = GetLinApprox(ListK,ListZK,ListKEnergy,Z);

 
  return KEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  G4double G4MuMinusCaptureCascade::GetLinApprox(const size_t N, 
                                               const G4double X[], 
                                               const G4double Y[], 
                                               G4double Xuser)
{
  G4double Yuser = 0.0;
  G4int i;

  if(N < 1) return Yuser;

  else if(Xuser < X[0])   Yuser = Y[0];

  else if(Xuser > X[N-1]) Yuser = Y[N-1];

  else {
    for (i = 1; i < N - 1; i++){
      if(Xuser < X[i]) {break;} 
    }    

    Yuser = X[i] - X[i-1];
    if(Yuser != 0.0){
      Yuser = Y[i-1] + (Y[i] - Y[i-1]) * (Xuser - X[i-1]) / Yuser;
    }
  }
  return Yuser;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector G4MuMinusCaptureCascade::GetRandomVec()
{
   //
   // generate uniform vector
   //

   G4double Theta = (2.0 * G4UniformRand() - 1.0) * pi ;
   G4double Phi  = twopi * G4UniformRand() ;
   G4double sinTheta = sin(Theta);
   G4double dirx = sinTheta * cos(Phi); 
   G4double diry = sinTheta * sin(Phi); 
   G4double dirz = cos(Theta); 

   return G4ThreeVector(dirx, diry, dirz);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuMinusCaptureCascade::AddNewParticle(G4ParticleDefinition* aParticle,
                                             G4ThreeVector Momentum,
                                             G4double mass,
                                             G4int* nParticle,
                                             G4GHEKinematicsVector* Cascade)
{

  // Store particle in the HEK vector and increment counter
   
   Cascade[*nParticle].SetZero();
   Cascade[*nParticle].SetMass( mass );
   Cascade[*nParticle].SetMomentumAndUpdate(Momentum.x(), Momentum.y(), Momentum.z());
   Cascade[*nParticle].SetParticleDef( aParticle );
   (*nParticle)++;

   return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4MuMinusCaptureCascade::DoCascade(const G4double Z, const G4double massA, 
                                               G4GHEKinematicsVector* Cascade)
{
  // Inicialization - cascade start from 14th level
  // N.C.Mukhopadhyay Phy. Rep. 30 (1977) 1.
  G4int nPart = 0;
  G4double EnergyLevel[14];

  G4double mass = MuMass * massA / (MuMass + massA) ;

  const G4double KEnergy = 13.6 * eV * Z * Z * mass/ electron_mass_c2;

  EnergyLevel[0] = GetKShellEnergy(Z);
  for( G4int i = 2; i < 15; i++ ) {
    EnergyLevel[i-1] = KEnergy / (i*i) ;
  }

  G4int nElec  = G4int(Z);
  G4int nAuger = 1;
  G4int nLevel = 13;
  G4double DeltaE;
  G4double pGamma = Z*Z*Z*Z;

  // Capture on 14-th level
  G4double ptot = sqrt(EnergyLevel[13]*(EnergyLevel[13] + 2.0*Emass));
  G4ThreeVector moment = ptot * GetRandomVec();

  AddNewParticle(theElectron,moment,Emass,&nPart,Cascade);

  // Emit new photon or electron
  // Simplified model for probabilities
  // N.C.Mukhopadhyay Phy. Rep. 30 (1977) 1.
  do {

    // case of Auger electrons
    if((nAuger < nElec) && ((pGamma + 10000.0) * G4UniformRand() < 10000.0) ) {
        nAuger++;
        DeltaE = EnergyLevel[nLevel-1] - EnergyLevel[nLevel];
        nLevel--;

        ptot = sqrt(DeltaE * (DeltaE + 2.0*Emass));
        moment = ptot * GetRandomVec();

        AddNewParticle(theElectron, moment, Emass, &nPart, Cascade);

    } else {

      // Case of photon cascade, probabilities from
      // C.S.Wu and L.Wilets, Ann. Rev. Nuclear Sci. 19 (1969) 527.

        G4double var = (10.0 + G4double(nLevel - 1) ) * G4UniformRand();
        G4int iLevel = nLevel - 1 ;
        if(var > 10.0) iLevel -= G4int(var-10.0) + 1;
        if( iLevel < 0 ) iLevel = 0;
        DeltaE = EnergyLevel[iLevel] - EnergyLevel[nLevel];
        nLevel = iLevel;
        moment = DeltaE * GetRandomVec();
        AddNewParticle(theGamma, moment, 0.0, &nPart, Cascade);
    }

  } while( nLevel > 0 );

  return nPart;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuMinusCaptureCascade::DoBoundMuonMinusDecay(G4double Z, G4double massA, 
                                                    G4int* nCascade, 
                                                    G4GHEKinematicsVector* Cascade)
{
  // Simulation on Decay of mu- on a K-shell of the muonic atom

  G4double Energy, r, x;
  G4double xmax = ( 1.0 + Emass*Emass/ (MuMass*MuMass) );
  G4double KEnergy = GetKShellEnergy(Z);

  // Calculate electron energy
  do {
    do {
      x = xmax*G4UniformRand();
    } while (G4UniformRand() < (3.0 - 2.0*x)*x*x );
    Energy = x*MuMass*0.5 - Emass - KEnergy;
  } while (Energy < 0.0);

   //
   // generate uniform vector
   //
  G4double ptot = sqrt(Energy * (Energy + 2.0*Emass));
  G4ThreeVector moment = ptot * GetRandomVec();

  AddNewParticle(theElectron, moment, Emass, nCascade, Cascade);
  
  // Calculate rest frame parameters of 2 neutrinos
  G4double E = MuMass*( 1.0 - x*0.5 );
  G4double P = sqrt( MuMass*MuMass*x*x*0.25 - Emass*Emass );

  if(P >= E) {P = E;}
  G4double ecm = 0.5 * sqrt( E*E - P*P );

   //
   // generate uniform vector
   //

  moment *= -P / (ptot * E);
  G4ThreeVector p1 = ecm * GetRandomVec();

  // Create Neutrinos
  G4LorentzVector N1 = G4LorentzVector(p1,ecm);
  N1.boost(moment);


  AddNewParticle(G4AntiNeutrinoE::AntiNeutrinoE(),G4ThreeVector(N1),0.0,nCascade,Cascade);

  G4LorentzVector N2 = G4LorentzVector(-p1,ecm);
  N2.boost(moment);

  AddNewParticle(G4NeutrinoMu::NeutrinoMu(),G4ThreeVector(N2),0.0,nCascade,Cascade);

  return;
}







