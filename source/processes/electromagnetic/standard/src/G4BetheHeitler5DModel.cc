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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4BetheHeitler5DModel.cc
//
// Authors:
// Igor Semeniouk and Denis Bernard,
// LLR, Ecole polytechnique & CNRS/IN2P3, 91128 Palaiseau, France
//
// Acknowledgement of the support of the French National Research Agency
// (ANR-13-BS05-0002).
//
// Reference: Nucl. Instrum. Meth. A 899 (2018) 85 (arXiv:1802.08253 [hep-ph])
//
// Class Description:
//
// Generates the conversion of a high-energy photon to an e+e- pair, either in the field of an 
// atomic electron (triplet) or nucleus (nuclear).
// Samples the five-dimensional (5D) differential cross-section analytical expression:
// . Non polarized conversion:
//   H.A. Bethe, W. Heitler, Proc. R. Soc. Lond. Ser. A 146 (1934) 83.
// . Polarized conversion:
//   T. H. Berlin and L. Madansky, Phys. Rev. 78 (1950) 623,
//   M. M. May, Phys. Rev. 84 (1951) 265,
//   J. M. Jauch and F. Rohrlich, The theory of photons and electrons, 1976.
//
// All the above expressions are named "Bethe-Heitler" here.
//
// Bethe & Heitler, put in Feynman diagram parlance, compute only the two dominant diagrams of
// the first order Born development, which is an excellent approximation for nuclear conversion
// and for high-energy triplet conversion.
//
// Only the linear polarisation of the incoming photon takes part in these expressions.
// The circular polarisation of the incoming photon does not (take part) and no polarisation 
// is transfered to the final leptons.
//
// In case conversion takes place in the field of an isolated nucleus or electron, the bare 
// Bethe-Heitler expression is used.
//
// In case the nucleus or the electron are part of an atom, the screening of the target field 
// by the other electrons of the atom is described by a simple form factor, function of q2:
// . nuclear: N.F. Mott, H.S.W. Massey, The Theory of Atomic Collisions, 1934.
// . triplet: J.A. Wheeler and W.E. Lamb, Phys. Rev. 55 (1939) 858.
//
// The nuclear form factor that affects the probability of very large-q2 events, is not considered.
//
// In principle the code is valid from threshold, that is from 2 * m_e c^2 for nuclear and from 
// 4 * m_e c^2 for triplet, up to infinity, while in pratice the divergence of the differential 
// cross section at small q2 and, at high-energy, at small polar angle, make it break down at 
// some point that depends on machine precision.
//
// Very-high-energy (above a few tens of TeV) LPM suppression effects in the normalized differential
// cross-section are not considered.
//
// The 5D differential cross section is sampled without any high-energy nor small 
// angle approximation(s).
// The generation is strictly energy-momentum conserving when all particles in the final state 
// are taken into account, that is, including the recoiling target.
// (In contrast with the BH expressions taken at face values, for which the electron energy is
// taken to be EMinus = GammaEnergy - EPlus)
//
// Tests include the examination of 1D distributions: see TestEm15
//
// Total cross sections are not computed (we inherit from other classes).
// We just convert a photon on a target when asked to do so.
//
// Pure nuclear, pure triplet and 1/Z triplet/nuclear mixture can be generated. 
//
// -------------------------------------------------------------------

#include "G4BetheHeitler5DModel.hh"
#include "G4EmParameters.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4IonTable.hh"
#include "G4NucleiProperties.hh"

#include "Randomize.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Pow.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4BetheHeitler5DModel::G4BetheHeitler5DModel(const G4ParticleDefinition* pd,
                                             const G4String& nam)
  : G4BetheHeitlerModel(pd, nam), fVerbose(1), fConversionType(0), iraw(false)
{
  SetLowEnergyLimit(2*CLHEP::electron_mass_c2);  
  theIonTable = G4IonTable::GetIonTable();
  // Verbosity levels: ( Can redefine as needed, but some consideration )
  // 0 = nothing
  // > 2 print results
  // > 3 print rejection warning from transformation (fix bug from gammaray .. )
  // > 4 print photon direction & polarisation
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4BetheHeitler5DModel::~G4BetheHeitler5DModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4BetheHeitler5DModel::Initialise(const G4ParticleDefinition* part,
				       const G4DataVector& vec)
{
  G4BetheHeitlerModel::Initialise(part, vec);

  G4EmParameters* theManager = G4EmParameters::Instance();
  // place to initialise model parameters
  fVerbose = theManager->Verbose();
  fConversionType  = theManager->GetConversionType();
  //////////////////////////////////////////////////////////////
  // iraw :
  //      true  : isolated electron or nucleus.
  //      false : inside atom -> screening form factor
  iraw = theManager->OnIsolated();
  // G4cout << "BH5DModel::Initialise verbose " << fVerbose
  // 	 << " isolated " << iraw << " ctype "<< fConversionType << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


//
// Converting from pair coordinate
//
void
G4BetheHeitler5DModel::BoostG4LorentzVector(const G4LorentzVector& p,
                                            const G4LorentzVector& q,
                                            G4LorentzVector& res) const
{
// p : 4-vector  which will be boosted
// q : 4-vector of new origin in the old coordinates
  const G4double pq = p.x()*q.x() + p.y()*q.y() + p.z()*q.z();
  const G4double qq = q.x()*q.x() + q.y()*q.y() + q.z()*q.z();
  const G4double mass2  = q.t()*q.t()-qq;
  if ( mass2 > 0.0 ) {
    const G4double mass = std::sqrt(q.t()*q.t()-qq);
    const G4double lf = ((q.t()-mass)*pq/qq+p.t())/mass;
    res.set( (p.x()+q.x()*lf), (p.y()+q.y()*lf), (p.z()+q.z()*lf),
	     ((p.t()*q.t()+pq)/mass) );
  } else {
    res = p;
    if ( fVerbose > 3 ) {
      G4cout << "G4BetheHeitler5DModel::BoostG4LorentzVector Warning point not converted"
	     << G4endl << "secondary in " << p
	     << G4endl << "Pair1 " << q << G4endl;
    }
  }
}

// assuming that q.x=q.y=0.0
void 
G4BetheHeitler5DModel::BoostG4LorentzVector(const G4LorentzVector& p,
                                            const G4double qz,
                                            const G4double qt,
                                            const G4double lffac,
                                            const G4double imass, 
                                            G4LorentzVector& res) const
{
// p : 4-vector  which will be boosted
// q : 4-vector of new origin in the old coordinates
  const G4double pq = p.z()*qz;
  const G4double lf = (lffac*pq+p.t())*imass;
  res.setZ(p.z()+qz*lf);
  res.setT((p.t()*qt+pq)*imass);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4BetheHeitler5DModel::MaxDiffCrossSection(const G4double* par, 
                                                    G4double Z, 
                                                    G4double e, 
                                                    G4double loge) const
{
  const G4double Q = e/par[9];
  return par[0] * G4Exp((par[2]+loge*par[4])*loge)
         / (par[1]+ G4Exp(par[3]*loge)+G4Exp(par[5]*loge))
         * (1+par[7]*G4Exp(par[8]*G4Log(Z))*Q/(1+Q));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4BetheHeitler5DModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                         const G4MaterialCutsCouple* couple,
                                         const G4DynamicParticle* aDynamicGamma,
                                         G4double, G4double)
{
  // MeV
  static const G4double ElectronMass   = CLHEP::electron_mass_c2;
  static const G4double ElectronMass2  = ElectronMass*ElectronMass;
  static const G4double alpha0         = CLHEP::fine_structure_const;
  // mm 
  static const G4double r0             = CLHEP::classic_electr_radius;
  // mbarn
  static const G4double r02            = r0*r0*1.e+25;
  static const G4double twoPi          = CLHEP::twopi;
  static const G4double factor         = alpha0 * r02 / (twoPi*twoPi);
  //  static const G4double factor1        = pow((6.0 * pi),(1.0/3.0))/(8.*alpha0*ElectronMass);
  static const G4double factor1        = 2.66134007899/(8.*alpha0*ElectronMass);
  //
  static const G4double PairInvMassMin = 2.*ElectronMass;
  //
  static const G4double nu[10] = { 0.0227436, 0.0582046, 3.0322675,  2.8275065, 
                           -0.0034004, 1.1212766, 1.8989468, 68.3492750,  
                            0.0211186, 14.4 };
  static const G4double tr[10] = { 0.0332350, 4.3942537, 2.8515925,  2.6351695, 
                           -0.0031510, 1.5737305, 1.8104647, 20.6434021,
                           -0.0272586, 28.9};
  // 
  static const G4double para[3][2] = { {11., -16.},{-1.17, -2.95},{-2., -0.5} };
  //
  static const G4double correctionIndex = 1.4;
  //
  const G4double GammaEnergy  = aDynamicGamma->GetKineticEnergy();
  const G4double GammaEnergy2 = GammaEnergy*GammaEnergy;
  // Will not be true tot cross section = 0
  if ( GammaEnergy <= 2.0*ElectronMass) { return; }
  //
  const G4ParticleMomentum GammaDirection = aDynamicGamma->GetMomentumDirection();
  G4ThreeVector GammaPolarization = aDynamicGamma->GetPolarization();

  // The protection polarization perpendicular to the direction vector,
  // as it done in G4LivermorePolarizedGammaConversionModel,
  // assuming Direction is unitary vector
  //  (projection to plane) p_proj = p - (p o d)/(d o d) x d
  if ( GammaPolarization.howOrthogonal(GammaDirection) != 0) {
    GammaPolarization -= GammaPolarization.dot(GammaDirection) * GammaDirection;
  }
  // End of Protection
  //
  const G4double GammaPolarizationMag = GammaPolarization.mag();
  //////////////////////////////////////////////////////////////
  // target element
  // select randomly one element constituting the material
  const G4Element* anElement  = SelectRandomAtom(couple, fTheGamma, GammaEnergy);
  // Atomic number
  const G4int Z       = anElement->GetZasInt();
  const G4int A       = SelectIsotopeNumber(anElement);
  const G4double iZ13 = 1./anElement->GetIonisation()->GetZ3();
  const G4double targetMass = G4NucleiProperties::GetNuclearMass(A, Z);

  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();

  // itriplet : true -- triplet, false -- nuclear.
  G4bool itriplet = false;
  if (fConversionType == 1) {
    itriplet = false;
  } else if (fConversionType == 2) {
    itriplet = true;
    if ( GammaEnergy <= 4.0*ElectronMass ) return;
  } else if ( GammaEnergy > 4.0*ElectronMass ) {
    // choose triplet or nuclear from a triplet/nuclear=1/Z
    // total cross section ratio.
    // approximate at low energies !
    if(rndmEngine->flat()*(Z+1) < 1.)  {
      itriplet = true;
    }
  }
  //
  const G4double RecoilMass  = itriplet ? ElectronMass : targetMass;
  const G4double RecoilMass2 = RecoilMass*RecoilMass;
  const G4double sCMS        = 2.*RecoilMass*GammaEnergy + RecoilMass2;
  const G4double sCMSPlusRM2 = sCMS + RecoilMass2;
  const G4double sqrts       = std::sqrt(sCMS);
  const G4double isqrts2     = 1./(2.*sqrts);
  //
  const G4double PairInvMassMax   = sqrts-RecoilMass;
  const G4double PairInvMassRange = PairInvMassMax/PairInvMassMin;
  const G4double lnPairInvMassRange = G4Log(PairInvMassRange);

  // initial state. Defines z axis of "0" frame as along photon propagation.
  // create 4-vectors: gamma0 + target0 and CMS=gamma0+target0
  // Since CMS(0., 0., GammaEnergy, GammaEnergy+RecoilMass) set some constants 
  // for the special boost that makes use of the form of CMS 4-vector 
  const G4double CMSqz    = GammaEnergy;
  const G4double CMSt     = GammaEnergy+RecoilMass;
  const G4double iCMSmass = 1./std::sqrt(RecoilMass*(RecoilMass+2.*GammaEnergy));
  const G4double CMSfact  = (CMSt-1./iCMSmass)/(CMSqz*CMSqz);
  // maximum value of pdf
  const G4double EffectiveZ = iraw ? 0.5 : Z;
  const G4double Threshold  = itriplet ? 4.*ElectronMass : 2.*ElectronMass;
  const G4double AvailableEnergy    = GammaEnergy - Threshold;
  const G4double LogAvailableEnergy = G4Log(AvailableEnergy);
  //
  const G4double MaxDiffCross = itriplet
    ? MaxDiffCrossSection(tr, EffectiveZ, AvailableEnergy, LogAvailableEnergy)
    : MaxDiffCrossSection(nu, EffectiveZ, AvailableEnergy, LogAvailableEnergy);
  //
  // 50% safety marging factor
  const G4double ymax = 1.5 * MaxDiffCross;
  // x1 bounds
  const G4double xu1 =   (LogAvailableEnergy > para[2][0])
                       ? para[0][0] + para[1][0]*LogAvailableEnergy
                       : para[0][0] + para[2][0]*para[1][0];
  const G4double xl1 =   (LogAvailableEnergy > para[2][1])
                       ? para[0][1] + para[1][1]*LogAvailableEnergy
                       : para[0][1] + para[2][1]*para[1][1];
  //
  G4LorentzVector Recoil0;
  G4LorentzVector Positron0;
  G4LorentzVector Electron0;
  G4double pdf    = 0.;
  // START Sampling
  do {
    ////////////////////////////////////////////////// 
    // pdf  pow(x,c) with c = 1.4
    // integral y = pow(x,(c+1))/(c+1) @ x = 1 =>  y = 1 /(1+c)
    // invCdf exp( log(y /* *( c + 1.0 )/ (c + 1.0 ) */ ) /( c + 1.0) )
    ////////////////////////////////////////////////// 
    const G4double X1 =
      G4Exp(G4Log(rndmEngine->flat())/(correctionIndex + 1.0));
 
    const G4double x0       = G4Exp(xl1 + (xu1 - xl1)*rndmEngine->flat());
    const G4double dum0     = 1./(1.+x0);
    const G4double cosTheta = (x0-1.)*dum0;
    const G4double sinTheta = std::sqrt(4.*x0)*dum0;

    const G4double PairInvMass  = PairInvMassMin*G4Exp(X1*X1*lnPairInvMassRange);

    G4double rndmv3[3];
    rndmEngine->flatArray(3, rndmv3);
    //--------------------------------------------------------------------------
    // const G4double ThetaLept    = pi*rndmv3[0];
    // const G4double cosThetaLept = std::cos(ThetaLept);
    // const G4double sinThetaLept = std::sin(ThetaLept);
    // 
    // const G4double PhiLept      = twoPi*rndmv3[1]-pi;
    // const G4double cosPhiLept   = std::cos(PhiLept);
    // const G4double sinPhiLept   = std::sin(PhiLept);
    //
    // const G4double Phi          = twoPi*rndmv3[2]-pi;
    // const G4double cosPhi       = std::cos(Phi);
    // const G4double sinPhi       = std::sin(Phi);
    //---------------------------------------------------------------------------
    // cos and sin theta-lepton
    const G4double cosThetaLept = std::cos(pi*rndmv3[0]);
    // sin(ThetaLept) is always in [0,+1] if ThetaLept is in [0,pi]
    const G4double sinThetaLept = std::sqrt((1.-cosThetaLept)*(1.+cosThetaLept)); 
    // cos and sin phi-lepton
    const G4double cosPhiLept   = std::cos(twoPi*rndmv3[1]-pi);
    // sin(PhiLept) is in [-1,0] if PhiLept in [-pi,0) and
    //              is in [0,+1] if PhiLept in [0,+pi]
    const G4double sinPhiLept   = std::copysign(std::sqrt((1.-cosPhiLept)*(1.+cosPhiLept)),rndmv3[1]-0.5);
    // cos and sin phi
    const G4double cosPhi       = std::cos(twoPi*rndmv3[2]-pi);
    const G4double sinPhi        = std::copysign(std::sqrt((1.-cosPhi)*(1.+cosPhi)),rndmv3[2]-0.5);

    // frames:
    // 0 : the laboratory Lorentz frame, axes along photon direction and polarisation
    // 1 : the center-of-mass Lorentz frame
    // 2 : the pair Lorentz frame
    // 3 : the laboratory Lorentz frame, Geant4 axes definition

    // in the center-of-mass frame
    const G4double RecEnergyCMS  = (sCMSPlusRM2-PairInvMass*PairInvMass)*isqrts2;
    const G4double LeptonEnergy2 = PairInvMass*0.5;
    // Denis ** correction
    //    const G4double thePRecoil    = std::sqrt( (RecEnergyCMS-RecoilMass)
    //                                         *(RecEnergyCMS+RecoilMass));
    const G4double ap1 = 2.0*GammaEnergy*RecoilMass -
      PairInvMass*PairInvMass + 2.0*PairInvMass*RecoilMass;
    const G4double bp1 = 2.0*GammaEnergy*RecoilMass -
      PairInvMass*PairInvMass - 2.0*PairInvMass*RecoilMass;
    
    if (bp1 <= 0.0 ) {
      if ( fVerbose > 3 ) {
	G4cout
	  << "G4BetheHeitler5DModel::SampleSecondaries Warning bp1 "
	  << bp1 << "point rejected" <<  G4endl
	  << "GammaEnergy " << GammaEnergy << G4endl
	  << "PairInvMass " << PairInvMass << G4endl;
    }
      pdf = -1.0; // force next iteration
      continue;
    }
    const G4double thePRecoil = std::sqrt(ap1 * bp1) * isqrts2;
    // Denis ** correction

    // back to the center-of-mass frame
    const G4LorentzVector Recoil1( thePRecoil*sinTheta*cosPhi,
			     thePRecoil*sinTheta*sinPhi,
			     thePRecoil*cosTheta,
			     RecEnergyCMS);

    const G4LorentzVector Pair1(-Recoil1.x(),
			  -Recoil1.y(),
			  -Recoil1.z(),
			  sqrts-RecEnergyCMS);

    // in the pair frame
    const G4double thePLepton    = std::sqrt( (LeptonEnergy2-ElectronMass)
                                             *(LeptonEnergy2+ElectronMass));
    
    const G4LorentzVector Positron2( thePLepton*sinThetaLept*cosPhiLept,
			       thePLepton*sinThetaLept*sinPhiLept,
			       thePLepton*cosThetaLept,
			       LeptonEnergy2);

    const G4LorentzVector Electron2(-Positron2.x(),
				    -Positron2.y(),
				    -Positron2.z(),
				    LeptonEnergy2);


    // Normalisation of final state phase space:
    // Section 47 of Particle Data Group, Chin. Phys. C, 40, 100001 (2016)
    const G4double Norme = Recoil1.vect().mag() * Positron2.vect().mag();
    //
    G4LorentzVector Positron1;
    G4LorentzVector Electron1;
    BoostG4LorentzVector(Positron2, Pair1, Positron1);
    BoostG4LorentzVector(Electron2, Pair1, Electron1);
    
    // back to the laboratory frame (make use of the CMS(0,0,Eg,Eg+RM)) form
    Recoil0.setX(Recoil1.x());
    Recoil0.setY(Recoil1.y());
    BoostG4LorentzVector(Recoil1  , CMSqz, CMSt, CMSfact, iCMSmass, Recoil0);

    Positron0.setX(Positron1.x());
    Positron0.setY(Positron1.y());
    BoostG4LorentzVector(Positron1, CMSqz, CMSt, CMSfact, iCMSmass, Positron0);

    Electron0.setX(Electron1.x());
    Electron0.setY(Electron1.y());    
    BoostG4LorentzVector(Electron1,  CMSqz, CMSt, CMSfact, iCMSmass, Electron0);

    // Jacobian factors
    const G4double Jacob0 = x0*dum0*dum0;
    const G4double Jacob1 = 2.*X1*lnPairInvMassRange*PairInvMass;
    const G4double Jacob2 = std::abs(sinThetaLept);

    const G4double EPlus = Positron0.t();
    const G4double PPlus = Positron0.vect().mag();
    const G4double sinThetaPlus = Positron0.vect().perp()/PPlus;
    const G4double cosThetaPlus = Positron0.vect().cosTheta();

    const G4double pPX  = Positron0.x();
    const G4double pPY  = Positron0.y();
    const G4double dum1 = 1./std::sqrt( pPX*pPX + pPY*pPY );
    const G4double cosPhiPlus = pPX*dum1;
    const G4double sinPhiPlus = pPY*dum1;

    // denominators:
    // the two cancelling leading terms for forward emission at high energy, removed
    const G4double elMassCTP = ElectronMass*cosThetaPlus;
    const G4double ePlusSTP  = EPlus*sinThetaPlus; 
    const G4double DPlus     = (elMassCTP*elMassCTP + ePlusSTP*ePlusSTP)
                              /(EPlus + PPlus*cosThetaPlus);

    const G4double EMinus = Electron0.t();
    const G4double PMinus = Electron0.vect().mag();
    const G4double sinThetaMinus = Electron0.vect().perp()/PMinus;
    const G4double cosThetaMinus = Electron0.vect().cosTheta();

    const G4double ePX  = Electron0.x();
    const G4double ePY  = Electron0.y();
    const G4double dum2 = 1./std::sqrt( ePX*ePX + ePY*ePY ); 
    const G4double cosPhiMinus =  ePX*dum2;
    const G4double sinPhiMinus =  ePY*dum2;

    const G4double elMassCTM = ElectronMass*cosThetaMinus;
    const G4double eMinSTM   = EMinus*sinThetaMinus; 
    const G4double DMinus    = (elMassCTM*elMassCTM + eMinSTM*eMinSTM) 
                              /(EMinus + PMinus*cosThetaMinus);

    // cos(phiMinus-PhiPlus)
    const G4double cosdPhi = cosPhiPlus*cosPhiMinus + sinPhiPlus*sinPhiMinus;
    const G4double PRec    = Recoil0.vect().mag();
    const G4double q2      = PRec*PRec;
    const G4double BigPhi  = -ElectronMass2 / (GammaEnergy*GammaEnergy2 * q2*q2);

    G4double FormFactor = 1.;
    if (!iraw) {
      if (itriplet) {
	const G4double qun = factor1*iZ13*iZ13;
	const G4double nun = qun * PRec;
	if (nun < 1.) {
          FormFactor =  (nun < 0.01) ? (13.8-55.4*std::sqrt(nun))*nun 
                                     : std::sqrt(1-(nun-1)*(nun-1)); 
	} // else FormFactor = 1 by default
      } else {
        const G4double dum3 = 217.*PRec*iZ13;
	const G4double AFF  = 1./(1. + dum3*dum3);
	FormFactor = (1.-AFF)*(1-AFF);
      }
    } // else FormFactor = 1 by default

    G4double betheheitler;
    if (GammaPolarizationMag==0.) {
      const G4double pPlusSTP   = PPlus*sinThetaPlus;
      const G4double pMinusSTM  = PMinus*sinThetaMinus;
      const G4double pPlusSTPperDP  = pPlusSTP/DPlus;
      const G4double pMinusSTMperDM = pMinusSTM/DMinus;
      const G4double dunpol = BigPhi*( 
                  pPlusSTPperDP *pPlusSTPperDP *(4.*EMinus*EMinus-q2)
                + pMinusSTMperDM*pMinusSTMperDM*(4.*EPlus*EPlus - q2)
                + 2.*pPlusSTPperDP*pMinusSTMperDM*cosdPhi
                    *(4.*EPlus*EMinus + q2 - 2.*GammaEnergy2)
                - 2.*GammaEnergy2*(pPlusSTP*pPlusSTP+pMinusSTM*pMinusSTM)/(DMinus*DPlus));
      betheheitler = dunpol * factor;
    } else {
      const G4double pPlusSTP  = PPlus*sinThetaPlus;
      const G4double pMinusSTM = PMinus*sinThetaMinus;
      const G4double pPlusSTPCPPperDP  = pPlusSTP*cosPhiPlus/DPlus;
      const G4double pMinusSTMCPMperDM = pMinusSTM*cosPhiMinus/DMinus;
      const G4double caa = 2.*(EPlus*pMinusSTMCPMperDM+EMinus*pPlusSTPCPPperDP);
      const G4double cbb = pMinusSTMCPMperDM-pPlusSTPCPPperDP;
      const G4double ccc = (pPlusSTP*pPlusSTP + pMinusSTM*pMinusSTM 
                          +2.*pPlusSTP*pMinusSTM*cosdPhi)/ (DMinus*DPlus);
      const G4double dtot= 2.*BigPhi*( caa*caa - q2*cbb*cbb - GammaEnergy2*ccc);
      betheheitler = dtot * factor;
    }
    //
    const G4double cross =  Norme * Jacob0 * Jacob1 * Jacob2 * betheheitler 
                          * FormFactor * RecoilMass / sqrts;
    pdf = cross * (xu1 - xl1) / G4Exp(correctionIndex*G4Log(X1)); // cond1;
  } while ( pdf < ymax * rndmEngine->flat() ); 
  // END of Sampling

  if ( fVerbose > 2 ) {
    G4double recul = std::sqrt(Recoil0.x()*Recoil0.x()+Recoil0.y()*Recoil0.y()
                              +Recoil0.z()*Recoil0.z());
    G4cout << "BetheHeitler5DModel GammaEnergy= " << GammaEnergy
	   << " PDF= " <<  pdf << " ymax= " << ymax 
           << " recul= " << recul << G4endl;
  }
  
  // back to Geant4 system

  if ( fVerbose > 4 ) {
    G4cout << "BetheHeitler5DModel GammaDirection " << GammaDirection << G4endl;
    G4cout << "BetheHeitler5DModel GammaPolarization " << GammaPolarization << G4endl;
    G4cout << "BetheHeitler5DModel GammaEnergy " << GammaEnergy << G4endl;
    G4cout << "BetheHeitler5DModel Conv "
	   << (itriplet ? "triplet" : "nucl") << G4endl;
  }

  if (GammaPolarizationMag == 0.0) {
    // set polarization axis orthohonal to direction
    GammaPolarization = GammaDirection.orthogonal().unit();
  } else {
    // GammaPolarization not a unit vector
    GammaPolarization /= GammaPolarizationMag;
  }

  // The unit norm vector that is orthogonal to the two others
  G4ThreeVector yGrec = GammaDirection.cross(GammaPolarization);
  // rotation
  G4ThreeVector Rot =  Recoil0.x()*GammaPolarization + Recoil0.y()*yGrec 
    + Recoil0.z()*GammaDirection;
  Recoil0.setVect(Rot); 
  Rot =  Positron0.x()*GammaPolarization + Positron0.y()*yGrec 
    + Positron0.z()*GammaDirection;
  Positron0.setVect(Rot);
  Rot =  Electron0.x()*GammaPolarization + Electron0.y()*yGrec 
    + Electron0.z()*GammaDirection;
  Electron0.setVect(Rot);
  //
  if ( fVerbose > 2 ) {
    G4cout << "BetheHeitler5DModel Recoil0 " << Recoil0.x() << " " << Recoil0.y() << " " << Recoil0.z() 
	   << " " << Recoil0.t() << " " << G4endl;
    G4cout << "BetheHeitler5DModel Positron0 " << Positron0.x() << " " << Positron0.y() << " " 
	   << Positron0.z() << " " << Positron0.t() << " " << G4endl;
    G4cout << "BetheHeitler5DModel Electron0 " << Electron0.x() << " " << Electron0.y() << " " 
	   << Electron0.z() << " " << Electron0.t() << " " << G4endl;
  }
  
  // Create secondaries
  
  // electron
  G4DynamicParticle* aParticle1 = new G4DynamicParticle(fTheElectron,Electron0);
  // positron
  G4DynamicParticle* aParticle2 = new G4DynamicParticle(fThePositron,Positron0);
  // create G4DynamicParticle object for the particle3 ( recoil )
  G4ParticleDefinition* RecoilPart;
  if (itriplet) {
    // triplet
    RecoilPart = fTheElectron;
  } else{
    RecoilPart = theIonTable->GetIon(Z, A, 0);
  }
  G4DynamicParticle* aParticle3 = new G4DynamicParticle(RecoilPart,Recoil0);
  
  // Fill output vector
  fvect->push_back(aParticle1);
  fvect->push_back(aParticle2);
  fvect->push_back(aParticle3);

  // kill incident photon
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
