///////////////////////////////////////////////////////////////////////////////
//
// MODULE:        G4GeneralParticleSource.hh
//
// Version:       1.1
// Date:          18/10/00
// Author:        C Ferguson, F Lei & P Truscott
// Organisation:  University of Southampton / DERA
// Customer:      ESA/ESTEC
//
///////////////////////////////////////////////////////////////////////////////
//
// DESCRIPTION
// -----------
// The General Particle Source is designed to extend the functionality of the
// G4ParticleGun class. It is designed to allow specification of input
// particles in terms of position, direction (or angular) and energy
// distributions.  This class is derived from G4VPrimaryGenerator.
//
///////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//
// Version 1.0, 28 February 2000, C Ferguson, Created.
//
// Version 1.1, 18 October 2000, Modified to inherit from G4VPrimaryGenerator.
// New name at the request of M. Asai.
//
///////////////////////////////////////////////////////////////////////////////
//
// MEMBER FUNCTIONS
// ----------------
//
// G4GeneralParticleSource ()
//    Constructor: Initializes variables and instantiates the 
//                 Messenger and Navigator classes
//
// ~G4GeneralParticleSource ()
//    Destructor:  deletes Messenger and prints out run information.
//
// void GeneratePrimaryVertex(G4Event *evt)
//    Generate the particles initial parameters.
//
// void SetPosDisType(G4String)
//    Allows user to choose Point, Plane, Surface or Volume source
//    position distributions.
//
// void SetPosDisShape(G4String)
//    Allows the user to choose the particular shape they wish for the
//    position distribution. Choices are Square, Circle, Ellipse, Rectangle,
//    Sphere, Ellipsoid, Cylinder, Parallelepiped.
//
// void SetCentreCoords(G4ThreeVector)
//    Sets the co-ordinates of the centre of the position distribution.
//
// void SetPosRot1(G4ThreeVector)
//    Used to specify the co-ordinate system for the position distribution
//    along with SetPosRot2. SetPosRot1 sets the vector x' and need not be
//    a unit vector.
//
// void SetPosRot2(G4ThreeVector)
//    Used in connection with SetPosRot1. This sets a vector in the plane
//    x'y'. By a series of cross products x', y', z' are generated. Again
//    need not be a unit vector.
// 
// void SetHalfX(G4double)
//    Sets the half length in x.
//
// void SetHalfY(G4double)
//    Sets the half length in y.
//
// void SetHalfZ(G4double)
//    Sets the half length in z.
//
// void SetRadius(G4double)
//    Sets the radius where appropriate for source distribution shapes.
//
// void SetRadius0(G4double)
//    Sets the inner radius where appropriate for source distribution shapes.
//
// void SetParAlpha(G4double)
//    Sets the angle Alpha in the Parallelepiped shapes.
//
// void SetParTheta(G4double)
//    Sets the angle Theta in the Parallelepiped shapes.
//
// void SetParPhi(G4double)
//    Sets the angle Phi in the Parallelepiped shapes.
//
// void ConfineSourceToVolume(G4String)
//    Used to confine the start positions to a particular volume.
//
// void GenerateRotationMatrices()
//    This is used to calculate the cross products and determine the
//    vectors x', y', z' for the position distribution.
//
// void GeneratePointSource()
//    Generates a point source.
//
// void GeneratePointsInPlane()
//    Generates starting positions confined to a plane source.
//
// void GeneratePointsOnSurface()
//    Generates starting positions confined to the surface of a shape.
//
// void GeneratePointsInVolume()
//    Generates starting positions within the volume of a shape.
//
// G4bool IsSourceConfined()
//    Checks the source is confined to the requested volume. 
//
// Angular Distribution Methods:
//
// void SetAngDistType(G4String)
//    Used to set the type of angular distribution wanted. Arguments
//    are iso, cos and user for isotropic, cosine-law and user-defined
//    respectively.
//
// void DefineAngRefAxes(G4String, G4ThreeVector)
//    DefineAngRefAxes is used in a similar way as SetPosRot to
//    define vectors, one x' and one in the plane x'y', to create
//    a rotated set of axes for the angular distribution.
//
// void SetMinTheta(G4double)
//    Sets the minimum value for the angle theta.
//
// void SetMinPhi(G4double)
//    Sets the minimum value for phi.  
//
// void SetMaxTheta(G4double)
//    Sets the maximum value for theta.
//
// void SetMaxPhi(G4double)
//    Sets the maximum value for phi.
//
// void UserDefAngTheta(G4ThreeVector)
//    This method allows the user to define a histogram in Theta.
//
// void UserDefAngPhi(G4ThreeVector)
//    This method allows the user to define a histogram in phi.
//
// void GenerateIsotropicFlux()
//    This method generates momentum vectors for particles according
//    to an isotropic distribution.
//
// void GenerateCosineLawFlux()
//    This method generates momentum vectors for particles according
//    to a cosine-law distribution.
//
// void GenerateUserDefFlux()
//    Controls generation of momentum vectors according to user-defined
//    distributions.
//
// G4double GenerateUserDefTheta()
//    Generates the theta angle according to a user-defined distribution.
//
// G4double GenerateUserDefPhi()
//    Generates phi according to a user-defined distribution.
//
// void SetUserWRTSurface(G4bool)
//    Allows user to have user-defined spectra either with respect to the
//    co-ordinate system (default) or with respect to the surface normal.
//
// Energy Distribution methods:
//
// void SetEnergyDisType(G4String)
//    Allows the user to choose the energy distribution type. The arguments
//    are Mono (mono-energetic), Lin (linear), Pow (power-law), Exp 
//    (exponential), Brem (bremsstrahlung), BBody (black-body), Cdg
//    (cosmic diffuse gamma-ray), User (user-defined), Arb (arbitrary
//    point-wise), Epn (energy per nucleon).
//
// void SetEmin(G4double)
//    Sets the minimum energy.
//
// void SetEmax(G4double)
//    Sets the maximum energy.
//
// void SetMonoEnergy(G4double)
//    Sets energy for mono-energetic distribution.
//
// void SetAlpha(G4double)
//    Sets alpha for a power-law distribution.
//
// void SetTemp(G4double)
//    Sets Temperature for a Brem or BBody distributions.
//
// void SetEzero(G4double)
//    Sets Ezero for an exponential distribution.
//
// void SetGradient(G4double)
//    Sets gradient for a linear distribution.
//
// void SetInterCept(G4double)
//    Sets intercept for a linear distribution.
//
// void UserEnergyHisto(G4ThreeVector)
//    Allows user to defined a histogram for the energy distribution.
//
// void ArbEnergyHisto(G4ThreeVector)
//    Allows the user to define an Arbitrary set of points for the
//    energy distribution.
//
// void EpnEnergyHisto(G4ThreeVector)
//    Allows the user to define an Energy per nucleon histogram.
//
// void Calculate()
//    Controls the calculation of Integral PDF for the Cdg and BBody
//    distributions.
//
// void CalculateCdgSpectrum()
//    Calculates the integral PDF for the Cdg distribution.
//
// void CalculateBbodySpectrum()
//    Calculates the Integral PDF for the Bbody distribution.
//
// void InputEnergySpectra(G4bool)
//    Allows the user to choose between momentum and energy histograms
//    for user-defined histograms and arbitrary point-wise spectr.
//    The default is true (energy).
//
// void InputDifferentialSpectra(G4bool)
//    Allows the user to choose between integral and differential 
//    distributions when using the arbitrary point-wise option.
//
// void ArbInterpolate(G4String)
//    ArbInterpolate allows the user to specify the type of function to
//    interpolate the Arbitrary points spectrum with.
//
// void LinearInterpolation()
//    Interpolates arbitrary points with a series of line segments.
//
// void LogInterpolation()
//    Interpolates arbitrary points with a series of power-laws.
//
// void ExpInterpolation()
//    Interpolates arbitrary points with a series of exponentials.
//
// void SplineInterpolation()
//    Interpolates arbitrary points using cubic splines.
//
// void GenerateMonoEnergetic()
//    Generates a mono-energetic source.
//
// void GenerateLinearEnergies()
//    Generates particle energies according to a linear distribution.
//
// void GeneratePowEnergies()
//    Generates particle energies according to a power-law distribution.
//
// void GenerateExpEnergies()
//    Generates particle energies according to an exponential distribution.
//
// void GenerateBremEnergies()
//    Generates particle energies according to a bremsstrahlung distribution.
//
// void GenerateBbodyEnergies()
//    Generates particle energies according to a black-body distribution.
//
// void GenerateCdgEnergies()
//    Generates particle energies according to a Cdg distribution.
//
// void GenUserHistEnergies()
//    Generates particle energies according to a user-defined distribution.
//
// void GenEpnHistEnergies()
//    Generates particle energies according to a energy per nucleon
//    distribution.
//
// void GenArbPointEnergies()
//    Generates particle energies according to an arbitrary point-wise
//    spectrum.
//
// void ConvertEPNToEnergy()
//    Converts energy per nucleon histograms to energy histograms.
//
// Biasing Methods:
//
// void SetXBias(G4ThreeVector)
//    Allows the user to re-distribute the random
//    numbers used to generate x co-ordinates.
//
// void SetYBias(G4ThreeVector)
//    Allows the user to re-distribute the random
//    numbers used to generate y co-ordinates.
//
// void SetZBias(G4ThreeVector)
//    Allows the user to re-distribute the random
//    numbers used to generate z co-ordinates.
//
// void SetThetaBias(G4ThreeVector)
//    Allows the user to re-distribute the random
//    numbers used to generate values of theta.
//
// void SetPhiBias(G4ThreeVector)
//    Allows the user to re-distribute the random
//    numbers used to generate values of phi.
//
// void SetEnergyBias(G4ThreeVector)
//    Allows the user to re-distribute the random
//    numbers used to generate the energies.
//
// G4double GenRandX()
//    Generates the random number for x, with or without biasing.
//
// G4double GenRandY()
//    Generates the random number for y, with or without biasing.
//
// G4double GenRandZ()
//    Generates the random number for z, with or without biasing.
//
// G4double GenRandTheta()
//    Generates the random number for theta, with or without biasing.
//
// G4double GenRandPhi()
//    Generates the random number for phi, with or without biasing.
//
// G4double GenRandEnergy()
//    Generates the random number for energy, with or without biasing.
//
// void SetVerbosity(G4int)
//    Sets the verbosity level.
//
///////////////////////////////////////////////////////////////////////////////
//
#include "G4VPrimaryGenerator.hh"
#include "G4Navigator.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleDefinition.hh"
#include "G4DataInterpolation.hh"

#include "G4GeneralParticleSourceMessenger.hh"

class G4GeneralParticleSource : public G4VPrimaryGenerator
{
public:
  G4GeneralParticleSource (); 
  ~G4GeneralParticleSource ();
  void GeneratePrimaryVertex(G4Event *evt);

  // methods to create source position dist.
  void SetPosDisType(G4String); // Point, Plane, Surface, Volume
  inline G4String GetPosDisType()
  { return SourcePosType; }
  void SetPosDisShape(G4String);
  inline G4String GetPosDisShape()
  { return Shape; }
  // SetPosDisShape - Square, Circle, Annulus, Ellipse, Rectangle, Sphere,
  // Ellipsoid, Cylinder, Right (parallelepiped).
  void SetCentreCoords(G4ThreeVector);
  inline G4ThreeVector GetCentreCoords()
  { return CentreCoords; }
  void SetPosRot1(G4ThreeVector); 
  void SetPosRot2(G4ThreeVector); 
  void SetHalfX(G4double);
  inline G4double GetHalfX()
  { return halfx; }
  void SetHalfY(G4double);
  inline G4double GetHalfY()
  { return halfy; }
  void SetHalfZ(G4double);
  inline G4double GetHalfZ()
  { return halfz; }
  void SetRadius(G4double);
  inline G4double GetRadius()
  { return Radius; }
  void SetRadius0(G4double);
  void SetParAlpha(G4double);
  void SetParTheta(G4double);
  void SetParPhi(G4double);
  void ConfineSourceToVolume(G4String);
  void GenerateRotationMatrices();
  // the following routines generate the source position
  void GeneratePointSource();
  void GeneratePointsInPlane();
  void GeneratePointsOnSurface();
  void GeneratePointsInVolume();

  G4bool IsSourceConfined();

  // Angular Distribution Methods
  void SetAngDistType(G4String);
  void DefineAngRefAxes(G4String, G4ThreeVector);
  void SetMinTheta(G4double);
  void SetMinPhi(G4double);
  void SetMaxTheta(G4double);
  void SetMaxPhi(G4double);
  void UserDefAngTheta(G4ThreeVector);
  void UserDefAngPhi(G4ThreeVector);
  // These methods generate the momentum vectors for the particles.
  void GenerateIsotropicFlux();
  void GenerateCosineLawFlux();
  void GenerateUserDefFlux();
  G4double GenerateUserDefTheta();
  G4double GenerateUserDefPhi();
  void SetUserWRTSurface(G4bool);

  // Energy Distribution methods
  void SetEnergyDisType(G4String);
  inline G4String GetEnergyDisType()
  {return EnergyDisType;}
  void SetEmin(G4double);
  inline G4double GetEmin()
  {return Emin;}
  inline G4double GetArbEmin()
  {return ArbEmin;}
  void SetEmax(G4double);
  inline G4double GetEmax()
  {return Emax;}
  inline G4double GetArbEmax()
  {return ArbEmax;}
  void SetMonoEnergy(G4double);
  void SetAlpha(G4double);
  void SetTemp(G4double);
  void SetEzero(G4double);
  void SetGradient(G4double);
  void SetInterCept(G4double);
  void UserEnergyHisto(G4ThreeVector);
  void ArbEnergyHisto(G4ThreeVector);
  void EpnEnergyHisto(G4ThreeVector);
  void Calculate();
  void CalculateCdgSpectrum(); 
  void CalculateBbodySpectrum();
  void InputEnergySpectra(G4bool);
  void InputDifferentialSpectra(G4bool);
  void ArbInterpolate(G4String);
  inline G4String GetIntType()
  {return IntType;}
  void LinearInterpolation();
  void LogInterpolation();
  void ExpInterpolation();
  void SplineInterpolation();
  // The following methods generate energies according to the spectral
  // parameters defined above.
  void GenerateMonoEnergetic();
  void GenerateLinearEnergies();
  void GeneratePowEnergies();
  void GenerateExpEnergies();
  void GenerateBremEnergies();
  void GenerateBbodyEnergies();
  void GenerateCdgEnergies();
  void GenUserHistEnergies();
  void GenEpnHistEnergies();
  void GenArbPointEnergies();
  // converts energy per nucleon to energy.
  void ConvertEPNToEnergy();

  // Biasing Methods 
  void SetXBias(G4ThreeVector);
  void SetYBias(G4ThreeVector);
  void SetZBias(G4ThreeVector);
  void SetThetaBias(G4ThreeVector);
  void SetPhiBias(G4ThreeVector);
  void SetEnergyBias(G4ThreeVector);
  G4double GenRandX();
  G4double GenRandY();
  G4double GenRandZ();
  G4double GenRandTheta();
  G4double GenRandPhi();
  G4double GenRandEnergy();

  // Set the verbosity level.
  void SetVerbosity(G4int);

  // Set the particle species
  void SetParticleDefinition (G4ParticleDefinition * aParticleDefinition);
  inline G4ParticleDefinition * GetParticleDefinition ()
  {return particle_definition;}

  // SR1.3 - allowing user to define an isotope by A,Z,energy.
  //  void  SetNucleus(Nucleus theIon1); // Sets the isotope.
  //inline Nucleus GetNucleus() {return theIon;} // Returns the isotope.

  inline void SetParticleCharge(G4double aCharge)
  { particle_charge = aCharge; }

  // Set polarization
  inline void SetParticlePolarization (G4ThreeVector aVal)
  {particle_polarization = aVal;}
  inline G4ThreeVector GetParticlePolarization ()
  {return particle_polarization;}

  // Set Time.
  inline void SetParticleTime(G4double aTime)
  { particle_time = aTime; }
  inline G4double GetParticleTime()
  { return particle_time; }

  inline void SetNumberOfParticles(G4int i)
  { NumberOfParticlesToBeGenerated = i; }
  inline G4int GetNumberOfParticles()
  { return NumberOfParticlesToBeGenerated; }

  inline G4ThreeVector GetParticlePosition()
  { return particle_position;}
  
  inline G4ThreeVector GetParticleMomentumDirection()
  { return particle_momentum_direction;}

  inline G4double GetTheta()
  { return Theta;}
  
  inline G4double GetPhi()
  { return Phi;}

  inline G4double GetParticleEnergy()
  {return particle_energy;}

private:

  // Position distribution Variables
  G4String SourcePosType; //Point,Plane,Surface,Volume
  G4String Shape; //Circle,Square,Rectangle etc..
  G4double halfx, halfy, halfz; //half lengths
  G4double Radius; //Radius for circles or spheres
  G4double Radius0; // The inner radius of an annulus
  G4ThreeVector CentreCoords; // Coords of centre of input shape
  G4ThreeVector Rotx, Roty, Rotz; // Unit vectors defining rotation matrix
  G4double ParAlpha, ParTheta, ParPhi; //Angle for Right Parallellepipeds
  G4bool Confine; //If true confines source distribution to VolName
  G4String VolName;
  G4ThreeVector SideRefVec1,SideRefVec2,SideRefVec3; //Side rotation matrices

  // Angular distribution variables.
  G4String AngDistType; // String to hold Ang dist type iso, cos, user
  G4ThreeVector AngRef1, AngRef2, AngRef3; // Reference axes for ang dist
  G4double MinTheta, MaxTheta, MinPhi, MaxPhi; // min/max theta/phi
  G4double Theta, Phi; // Store these for use with DEBUG
  G4bool IPDFThetaExist, IPDFPhiExist; // tell whether IPDF histos exist
  G4PhysicsOrderedFreeVector UDefThetaH; // Theta histo data
  G4PhysicsOrderedFreeVector IPDFThetaH; //Cumulative Theta histogram.
  G4PhysicsOrderedFreeVector UDefPhiH; // Phi histo bins
  G4PhysicsOrderedFreeVector IPDFPhiH; // Cumulative phi histogram.
  G4String UserDistType; //String to hold user distributions
  G4bool UserWRTSurface; // G4bool to tell whether user wants distribution wrt
                       // surface normals or co-ordinate system

  // Energy Distribution variables
  G4String EnergyDisType; // energy dis type Variable  - Mono,Lin,Exp,etc
  G4double MonoEnergy; //Mono-energteic energy
  G4double Emin, Emax; // emin and emax
  G4double alpha, Ezero, Temp; // alpha (pow), E0 (exp) and Temp (bbody,brem)
  G4double grad, cept; // gradient and intercept for linear spectra
  G4bool EnergySpec; // true - energy spectra, false - momentum spectra
  G4bool DiffSpec;  // true - differential spec, false integral spec
  G4bool ApplyRig; // false no rigidity cutoff, true then apply one
  G4double ERig; // energy of rigidity cutoff
  G4PhysicsOrderedFreeVector UDefEnergyH; // energy hist data
  G4PhysicsOrderedFreeVector IPDFEnergyH;  
  G4bool IPDFEnergyExist, IPDFArbExist, Epnflag;
  G4PhysicsOrderedFreeVector ArbEnergyH; // Arb x,y histogram
  G4PhysicsOrderedFreeVector IPDFArbEnergyH; // IPDF for Arb
  G4PhysicsOrderedFreeVector EpnEnergyH; 
  G4double CDGhist[3]; // cumulative histo for cdg
  G4double BBHist[10001], Bbody_x[10001];
  G4String IntType; // Interpolation type
  G4double Arb_grad[256], Arb_cept[256]; // grad and cept for 256 segments
  G4double Arb_alpha[256], Arb_Const[256]; // alpha and constants
  G4double Arb_ezero[256]; // ezero
  G4double ArbEmin, ArbEmax; // Emin and Emax for the whole arb distribution
  //use primarily for debug.

  // Bias variables
  G4bool XBias, IPDFXBias;
  G4PhysicsOrderedFreeVector XBiasH;
  G4PhysicsOrderedFreeVector IPDFXBiasH;
  G4bool YBias, IPDFYBias;
  G4PhysicsOrderedFreeVector YBiasH;
  G4PhysicsOrderedFreeVector IPDFYBiasH;
  G4bool ZBias, IPDFZBias;
  G4PhysicsOrderedFreeVector ZBiasH;
  G4PhysicsOrderedFreeVector IPDFZBiasH;
  G4bool ThetaBias, IPDFThetaBias;
  G4PhysicsOrderedFreeVector ThetaBiasH;
  G4PhysicsOrderedFreeVector IPDFThetaBiasH;
  G4bool PhiBias, IPDFPhiBias;
  G4PhysicsOrderedFreeVector PhiBiasH;
  G4PhysicsOrderedFreeVector IPDFPhiBiasH;
  G4bool EnergyBias, IPDFEnergyBias;
  G4PhysicsOrderedFreeVector EnergyBiasH;
  G4PhysicsOrderedFreeVector IPDFEnergyBiasH;
  G4double bweights[6], bweight; //record x,y,z,theta,phi,energy weights

  // Other particle properties
  G4int                  NumberOfParticlesToBeGenerated;
  G4ParticleDefinition * particle_definition;
  G4ParticleMomentum     particle_momentum_direction;
  G4double               particle_energy;
  G4double               particle_charge;
  G4ThreeVector          particle_position;
  G4double               particle_time;
  G4ThreeVector          particle_polarization;

  // Verbosity
  G4int verbosityLevel;

private:

  G4DataInterpolation *SplineInt; // holds Spline stuff

  G4GeneralParticleSourceMessenger *theMessenger;
  G4Navigator *gNavigator;

  //  Nucleus theIon;
  
};








