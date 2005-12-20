//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4DNATest.cc,v 1.13 2005-12-20 13:53:56 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "globals.hh"
#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include <memory>
#include <cstdlib>

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ProcessManager.hh"

#include "G4EnergyLossTables.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4ForceCondition.hh"

#include "G4LowEnergyBremsstrahlung.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4eIonisation.hh"
#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4ComptonScattering.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4RunManager.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4GRSVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4Step.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"

#include "G4UnitsTable.hh"

#include "AIDA/IManagedObject.h"
#include "AIDA/IAnalysisFactory.h"
#include "AIDA/ITreeFactory.h"
#include "AIDA/ITree.h"
#include "AIDA/IHistogramFactory.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"
#include "AIDA/IHistogram3D.h"
#include "AIDA/ITupleFactory.h"
#include "AIDA/ITuple.h"

// DNA
#include "G4DNAGenericIonsManager.hh"
#include "G4DNAElectronElasticBrenner.hh"
#include "G4DNAElectronElasticEmfietzoglou.hh"
#include "G4DNAProtonExcitation.hh"
#include "G4DNAHeliumExcitation.hh"
#include "G4DNAAlphaPlusExcitation.hh"
#include "G4DNAAlphaPlusPlusExcitation.hh"
#include "G4DNAElectronExcitation.hh"
#include "G4DNAProtonRuddIonization.hh"
#include "G4DNAProtonChargeDecrease.hh"
#include "G4DNAHydrogenChargeIncrease.hh"
#include "G4DNAHydrogenRuddIonization.hh"
#include "G4DNAProtonBornExcitation.hh"
#include "G4DNAElectronBornExcitation.hh"
#include "G4DNAAlphaPlusChargeDecrease.hh"
#include "G4DNAAlphaPlusChargeIncrease.hh"
#include "G4DNAAlphaPlusPlusChargeDecrease.hh"
#include "G4DNAHeliumChargeIncrease.hh"

//! \brief Options structure
struct Options
{
 //! \brief Mean free path test
 bool meanFreePathTest;
 //! \brief Post step do it test
 bool postStepDoItTest;
 //! \brief Post step do it test
 bool randomEnergy;
 //! \brief Output file name
 const char *outputFileName;
 //! \brief Material name
 const char *material;
 //! \brief Process name
 const char *process;
 //! \brief Particle name
 const char *particle;
 //! \brief Minimum energy
 G4double minEnergy;
 //! \brief Maximum energy
 G4double maxEnergy;
 //! \brief Number of energy step
 G4int nEnergySteps;
 //! \brief Number of interactions
 G4int nIterations;
};

//! \brief Default output file name
struct Options defaultOptions = { false, false, false, "G4DNATest.hbook", "Water", "G4DNAElectronElasticBrenner", "electron", 7.5*eV, 200*eV, 300, 1 };

//! \brief Creates some materials
void CreateMaterials(void)
{
 G4Element * H       = new G4Element ("Hydrogen", "H",   1.,   1.01*g/mole);
 G4Element * O       = new G4Element ("Oxygen",   "O",   8.,  16.00*g/mole);
 G4Element * C       = new G4Element ("Carbon",   "C",   6.,  12.00*g/mole);
 G4Element * Cs      = new G4Element ("Cesium",   "Cs", 55., 132.905*g/mole);
 G4Element * I       = new G4Element ("Iodine",   "I",  53., 126.9044*g/mole);

 G4Material * Si     = new G4Material("Silicon",   14., 28.055*g/mole,  2.33*g/cm3);
 G4Material * Fe     = new G4Material("Iron",      26.,  55.85*g/mole,  7.87*g/cm3);
 G4Material * Cu     = new G4Material("Copper",    29.,  63.55*g/mole,  8.96*g/cm3);
 G4Material * W      = new G4Material("Tungsten",  74., 183.85*g/mole, 19.30*g/cm3);
 G4Material * Pb     = new G4Material("Lead",      82., 207.19*g/mole, 11.35*g/cm3);
 G4Material * U      = new G4Material("Uranium",   92., 238.03*g/mole, 18.95*g/cm3);
 G4Material * maO    = new G4Material("Oxygen",     8.,  16.00*g/mole,   1.1*g/cm3);
 G4Material * water  = new G4Material ("Water",  1.*g/cm3,     2);
 water->AddElement(H, 2);
 water->AddElement(O, 1);

 G4Material* ethane = new G4Material ("Ethane", 0.4241*g/cm3, 2);
 ethane->AddElement(H, 6);
 ethane->AddElement(C, 2);

 G4Material* csI    = new G4Material ("CsI",    4.53*g/cm3,   2);
 csI->AddElement(Cs, 1);
 csI->AddElement(I, 1);
 
 // This is needed to suppress some warnings. These lines can be deleted;
 Si->GetTemperature();
 Fe->GetTemperature();
 Cu->GetTemperature();
 W->GetTemperature();
 Pb->GetTemperature();
 U->GetTemperature();
 maO->GetTemperature();
 water->GetTemperature();
 ethane->GetTemperature();
 csI->GetTemperature();
}

//! \brief Process the options arguments
//! \param argc Number of arguments
//! \param argv Pointer to the arguments
//! \param options Structure to fill-in
void processOptions(int argc, char ** argv, struct Options * options)
{
 options->meanFreePathTest = defaultOptions.meanFreePathTest;
 options->postStepDoItTest = defaultOptions.meanFreePathTest;
 options->randomEnergy     = defaultOptions.randomEnergy;
 options->outputFileName   = defaultOptions.outputFileName;
 options->material         = defaultOptions.material;
 options->process          = defaultOptions.process;
 options->particle         = defaultOptions.particle;
 options->minEnergy        = defaultOptions.minEnergy;
 options->maxEnergy        = defaultOptions.maxEnergy;
 options->nEnergySteps     = defaultOptions.nEnergySteps;
 options->nIterations      = defaultOptions.nIterations;
 
 int i(1);
 
 while (i<argc)
 {
  if (argv[i][0]=='-' && argv[i][2]==0)
  {
   switch(argv[i][1])
   {
    case 'h':
    case '?':
     G4cout << argv[0] << " [-h|-?] [-a] [-b] [-r] [-o <file name>] [-m <material name>] [-p <process name>] [-P <particle name>] [-e <min energy in eV>] [-E <max energy in eV>] [-s <energy steps>] [-n <iterations>] " << G4endl
            << G4endl
            << "-h|-?     Shows this help" << G4endl
            << "-a        Enables mean free path test" << G4endl
            << "-b        Enables post step do it test" << G4endl
            << "-r        Energy is choosen at random within the range" << G4endl
            << "-o <arg>  Set the output file name (default: \"" << defaultOptions.outputFileName << "\")" << G4endl
            << "-m <arg>  Set the material (default: \"" << defaultOptions.material << "\")" << G4endl
            << "-p <arg>  Set the process (default: \"" << defaultOptions.process << "\")" << G4endl
            << "-P <arg>  Set the incoming particle (default: \"" << defaultOptions.particle << "\")" << G4endl
            << "-e <arg>  Set the low energy range in eV (default: " << defaultOptions.minEnergy/eV << " eV)" << G4endl
            << "-E <arg>  Set the high energy range in eV (default: " << defaultOptions.maxEnergy/eV << " eV)" << G4endl
            << "-s <arg>  Set the energy range step (default: " << defaultOptions.nEnergySteps << ")"<< G4endl
            << "-n <arg>  Set the number of iterations for the post step do it (default: " << defaultOptions.nIterations << ")" << G4endl;
     exit(0);
     break;
     
    case 'a':
     options->meanFreePathTest=true;
     break;
     
    case 'b':
     options->postStepDoItTest=true;
     break;
     
    case 'r':
     options->randomEnergy=true;
     break;
     
    case 'o':
     i++;
     if (i<argc)
     {
      options->outputFileName = argv[i];
      break;
     }

    case 'm':
     i++;
     if (i<argc)
     {
      options->material       = argv[i];
      break;
     }
     
    case 'p':
     i++;
     if (i<argc)
     {
      options->process        = argv[i];
      break;
     }

    case 'P':
     i++;
     if (i<argc)
     {
      options->particle       = argv[i];
      break;
     }

    case 'e':
     i++;
     if (i<argc)
     {
      options->minEnergy      = std::atof(argv[i])*eV;
      if (options->minEnergy <= 0.)
      {
       G4cout << argv[0] << ": Energy must be > 0." << G4endl;
       exit(-1);
      }

      break;
     }

    case 'E':
     i++;
     if (i<argc)
     {
      options->maxEnergy      = std::atof(argv[i])*eV;
      if (options->maxEnergy <= 0.)
      {
       G4cout << argv[0] << ": Energy must be > 0." << G4endl;
       exit(-1);
      }

      break;
     }

    case 's':
     i++;
     if (i<argc)
     {
      options->nEnergySteps   = atoi(argv[i]);
      if (options->nEnergySteps <= 1)
      {
       G4cout << argv[0] << ": Expected at least two steps." << G4endl;
       exit(-1);
      }
      
      break;
     }
     
     G4cout << argv[0] << ": Expected one more parameter in " << argv[i] << " option. Use -h option for help." << G4endl;
     exit(-1);
     
    case 'n':
     i++;
     if (i<argc)
     {
      options->nIterations    = atoi(argv[i]);
      if (options->nIterations <= 0)
      {
       G4cout << argv[0] << ": Expected at least one iteration." << G4endl;
       exit(-1);
      }
      
      break;
     }
     
     G4cout << argv[0] << ": Expected one more parameter in " << argv[i] << " option. Use -h option for help." << G4endl;
     exit(-1);
     
    default:
     G4cout << argv[0] << ": Unknown " << argv[i] << " option. Use -h option for help." << G4endl;
     exit(-1);
   }
  }
  else
  {
   G4cout << argv[0] << ": Bad arguments. Use -h option for help." << G4endl;
   exit(-1);
  }
  
  i++;
 }
 
 if (options->minEnergy >= options->maxEnergy)
 {
  G4cout << argv[0] << ": Mininum energy is higher than maximum energy" << G4endl;
  exit(-1);
 }

 G4cout << "Mean free path test:      ";
 if (options->meanFreePathTest)
  G4cout << "On";
 else
  G4cout << "Off";
 G4cout << G4endl << "Post step do it test:     ";
 if (options->postStepDoItTest)
  G4cout << "On";
 else
  G4cout << "Off";
 G4cout << G4endl << "Random energy generation: ";
 if (options->randomEnergy)
  G4cout << "On";
 else
  G4cout << "Off";
 G4cout << G4endl << "Output file:              " << options->outputFileName << G4endl;
 G4cout << "Material:                 " << options->material << G4endl;
 G4cout << "Process:                  " << options->process << G4endl;
 G4cout << "Min energy:               " << options->minEnergy/eV << " eV" << G4endl;
 G4cout << "Max energy:               " << options->maxEnergy/eV << " eV" << G4endl;
 G4cout << "N energy steps:           " << options->nEnergySteps << G4endl;
 G4cout << "N iterations:             " << options->nIterations << G4endl;
}

//! \brief Return the selected material
//! \param options Options for the material choice
//! \return The material
G4Material * GetSelectedMaterial(const struct Options & options)
{
 const G4MaterialTable* theMaterialTable=G4Material::GetMaterialTable();

 G4int i(G4Material::GetNumberOfMaterials());
 
 while (i>0)
 {
  i--;
  
  if ((*theMaterialTable)[i]->GetName()==options.material)
   return (*theMaterialTable)[i];
 }
 
 i=G4Material::GetNumberOfMaterials();
 
 G4cout << "Available materials are: " << G4endl;
 while (i>0)
 {
  i--;
  G4cout << (*theMaterialTable)[i]->GetName();

  if (i>0)
   G4cout << ", ";
 }
 
 G4cout << G4endl;
 
 exit(-2);
 return 0;
}

//! \brief Creates the geometry
//! \param options Options for the material choice
//! \return The world volume
G4PVPlacement * CreateGeometry(const struct Options & options)
{
 G4Box* theFrame = new G4Box ("Frame", 1*mm, 1*mm, 1*mm);
  
 G4LogicalVolume* logicalFrame = new G4LogicalVolume(theFrame, GetSelectedMaterial(options), "LFrame", 0, 0, 0);
  
 G4PVPlacement * placement = new G4PVPlacement(0, G4ThreeVector(), "PFrame", logicalFrame, 0, false, 0);
   
 G4cout << "[OK] Geometry built" << G4endl;
 return placement;
}

//! \brief Get process from options
//! \param options Options for the process choice
//! \return The choosen process
G4VLowEnergyTestableDiscreteProcess * GetSelectedProcess(const struct Options & options)
{
 static G4VLowEnergyTestableDiscreteProcess ** processes=0;
 if (!processes)
 {
  G4DNAGenericIonsManager * genericIonsManager;
  genericIonsManager=G4DNAGenericIonsManager::Instance();
  
  processes=new G4VLowEnergyTestableDiscreteProcess * [18];
  processes[0]=new G4DNAElectronElasticBrenner;
  processes[1]=new G4DNAElectronElasticEmfietzoglou;
  processes[2]=new G4DNAProtonExcitation;
  processes[3]=new G4DNAHeliumExcitation;
  processes[4]=new G4DNAAlphaPlusExcitation;
  processes[5]=new G4DNAAlphaPlusPlusExcitation;
  processes[6]=new G4DNAElectronExcitation;
  processes[7]=new G4DNAProtonRuddIonization;
  processes[8]=new G4DNAProtonChargeDecrease;
  processes[9]=new G4DNAHydrogenChargeIncrease;
  processes[10]=new G4DNAHydrogenRuddIonization;
  processes[11]=new G4DNAProtonBornExcitation;
  processes[12]=new G4DNAElectronBornExcitation;
  processes[13]=new G4DNAAlphaPlusChargeDecrease;
  processes[14]=new G4DNAAlphaPlusChargeIncrease;
  processes[15]=new G4DNAAlphaPlusPlusChargeDecrease;
  processes[16]=new G4DNAHeliumChargeIncrease;
  processes[17]=0;
 }
 
 unsigned long i(0);
 while (processes[i])
 {
  if (processes[i]->GetProcessName()==options.process)
   return processes[i];
   
  i++;
 }
 
 G4cout << "Available processes are: " << G4endl;
 i=0;
 while (processes[i])
 {
  G4cout << processes[i]->GetProcessName();
  i++;
  
  if (processes[i])
   G4cout << ", ";
 }
 
 G4cout << G4endl;
 
 exit(-2);
 return 0;
}

G4ParticleDefinition * const * ParticleList()
{
 static G4ParticleDefinition * * particles=0;

 if (!particles)
 {
  G4DNAGenericIonsManager * genericIonsManager;
  genericIonsManager=G4DNAGenericIonsManager::Instance();
  
  particles=new G4ParticleDefinition * [10];
  particles[0]=G4Electron::Electron();
  particles[1]=G4Positron::Positron();
  particles[2]=G4Gamma::Gamma();
  particles[3]=G4AntiProton::AntiProton();
  particles[4]=G4Proton::Proton();
  particles[5]=genericIonsManager->GetIon("alpha++");
  particles[6]=genericIonsManager->GetIon("alpha+");
  particles[7]=genericIonsManager->GetIon("helium");
  particles[8]=genericIonsManager->GetIon("hydrogen");
  particles[9]=0;
 }
 
 return particles;
}

G4int GetParticleIndex(G4ParticleDefinition * const particle)
{
 G4int i(0);
 G4ParticleDefinition * const * particles(ParticleList());

 while (*particles)
 {
  if ((*particles)==particle)
   return i;
   
  i++;
  particles++;
 }
 
 return -1;
}

G4ParticleDefinition * GetSelectedParticle(const struct Options & options)
{
 G4ParticleDefinition * const * particles(ParticleList());
 while (*particles)
 {
  if ((*particles)->GetParticleName()==options.particle)
   return (*particles);
  
  particles++;
 }
 
 G4cout << "Available particles are: " << G4endl;
 particles=ParticleList();

 while (*particles)
 {
  G4cout << (*particles)->GetParticleName();
  
  particles++;

  if (*particles)
   G4cout << ", ";
 }
 
 G4cout << G4endl;
 
 exit(-2);
 return 0;
}

//! \brief Setup processes
//! \param options Options for the process choice
void SetPhysics(const struct Options & options, G4VPhysicalVolume * world)
{
 G4ProductionCutsTable * cutsTable(G4ProductionCutsTable::GetProductionCutsTable());
 G4ProductionCuts * cuts(cutsTable->GetDefaultProductionCuts());

 G4ParticleDefinition * const * particles(ParticleList());
 while (*particles)
 {
  cuts->SetProductionCut(1*micrometer, *particles);
  particles++;
 }

 cutsTable->UpdateCoupleTable(world);
 G4cout << "[OK] Cuts are defined " << G4endl;

 G4ParticleDefinition * particle(GetSelectedParticle(options));
 G4VProcess * eProcess=GetSelectedProcess(options);
 if (! (eProcess->IsApplicable(*particle)))
 {
  G4cout<< "Process " << eProcess->GetProcessName() << " is not applicable to " << particle->GetParticleName() << G4endl;
  exit(0);
  return;
 }
 
 G4cout<< "[OK] Process " << eProcess->GetProcessName() << " is applicable to " << particle->GetParticleName() << G4endl;
  
 G4ProcessManager * gProcessManager(new G4ProcessManager(particle));
 particle->SetProcessManager(gProcessManager);
 gProcessManager->AddDiscreteProcess(eProcess);

 G4cout << "[OK] Processes are defined " << G4endl;
 

 G4cout << "[OK] Building physics tables" << G4endl;
 eProcess->BuildPhysicsTable(* particle);

 G4cout << "[OK] Physics tables built" << G4endl;
}

//! \brief Generates the step
//! \param options Options related to the track generation
//! \return The generated track
G4Step * GenerateStep(const struct Options & options)
{
 G4ThreeVector momentumDirection;
 
 momentumDirection.setRThetaPhi(1., std::acos(2*G4UniformRand()-1.), twopi * G4UniformRand());

 G4double lnEnergyMin=std::log(options.minEnergy);
 G4double lnEnergyMax=std::log(options.maxEnergy); 
 G4DynamicParticle * dynamicElectron(new G4DynamicParticle(GetSelectedParticle(options), momentumDirection, std::exp(lnEnergyMin+(lnEnergyMax-lnEnergyMin)*G4UniformRand())));
  
 G4Track * aTrack(new G4Track(dynamicElectron, 0., G4ThreeVector(0., 0., 0.)));

 G4Step* aStep(new G4Step());  
 aStep->SetTrack(aTrack);
 aTrack->SetStep(aStep);

 G4Material * material(GetSelectedMaterial(options));
 G4ProductionCutsTable * cutsTable(G4ProductionCutsTable::GetProductionCutsTable());
 const G4MaterialCutsCouple * theCouple(cutsTable->GetMaterialCutsCouple(material, cutsTable->GetDefaultProductionCuts()));

 G4StepPoint * aPoint(new G4StepPoint());
 aPoint->SetPosition(G4ThreeVector(0., 0., 0.));
 aPoint->SetMaterial(material);
 aPoint->SetMaterialCutsCouple(theCouple);
 aPoint->SetSafety(10000.*cm);

 aStep->SetPreStepPoint(aPoint);  
 
 return aStep;
}

void ProgressBar(G4int remainingIterations)
{
 static time_t startingTime;
 static time_t nextDumpTime;
 static G4int startingIteration(0);
 time_t now;
 
 if (remainingIterations==0)
 {
  startingIteration=0;
 }
 else if (startingIteration==0)
 {
  startingTime=time(0);
  nextDumpTime=startingTime+3;
  startingIteration=remainingIterations;
 }
 else
 {
  now=time(0);
  if (now>nextDumpTime)
  {
   nextDumpTime=now+10;
   G4double time;
   G4double perc;
   
   time=std::floor(static_cast<G4double>(now-startingTime)/static_cast<G4double>(startingIteration-remainingIterations)*static_cast<G4double>(remainingIterations)+0.5);
   perc=std::floor(static_cast<G4double>(remainingIterations)/static_cast<G4double>(startingIteration)*200.+.5)/2.;
   
   G4cout << "  " << perc << " % Remaining time: " << time << " s        \r";
   G4cout.flush();
  }
 }
}

//! \brief Test the mean free path table
//! \param tupleFactory The tuple factory
//! \param options Options related to the mean free path test
void MeanFreePathTest(AIDA::ITupleFactory * tupleFactory, const struct Options & options)
{
 AIDA::ITuple* iTuple = tupleFactory->create("1", "Mean Free Path Ntuple", "double k, log_k, mfp, log_mfp, cpu_time");
 
 G4double energy(options.minEnergy);
 G4double stpEnergy(std::pow(options.maxEnergy/energy, 1./static_cast<G4double>(options.nEnergySteps-1)));
 G4int step(options.nEnergySteps);
 
 G4ForceCondition condition;
 G4VLowEnergyTestableDiscreteProcess * process(GetSelectedProcess(options));
 
 G4double mfp;
 clock_t time;
 
 ProgressBar(0);
 while (step>0)
 {
  G4Step * aStep(GenerateStep(options));
  G4Track * aTrack(aStep->GetTrack());

  if (!options.randomEnergy)
  {
   aTrack->SetKineticEnergy(energy);
   energy*=stpEnergy;
  }
  ProgressBar(step);
  step--;

  time=clock();
  mfp=process->DumpMeanFreePath(*aTrack, 1.*mm, &condition)/m;
  time=clock()-time;

  iTuple->fill(iTuple->findColumn("k"), aTrack->GetKineticEnergy()/eV);
  iTuple->fill(iTuple->findColumn("log_k"), std::log10(aTrack->GetKineticEnergy()/eV));
  iTuple->fill(iTuple->findColumn("mfp"), mfp);
  iTuple->fill(iTuple->findColumn("log_mfp"), std::log10(mfp));
  iTuple->fill(iTuple->findColumn("cpu_time"), static_cast<G4double>(time)/static_cast<G4double>(CLOCKS_PER_SEC));
  iTuple->addRow();
 
  delete aTrack;
  delete aStep;
 }
}

//! \brief Test the post step do it
//! \param tupleFactory The tuple factory
//! \param options Options related to the post step do it test
void PostStepDoItTest(AIDA::ITupleFactory * tupleFactory, const struct Options & options)
{
 AIDA::ITuple* iTuple = tupleFactory->create("2", "Post Step Do It Test", "double iteration, step, in_k, log_in_k, in_theta, in_phi, in_pol_theta, in_pol_phi, e_deposit, log_e_deposit, trk_status, out_k, log_out_k, out_theta, out_phi, out_pol_theta, out_pol_phi, n_secondaries, sec0_type, sec0_k, log_sec0_k, sec0_theta, sec0_phi, sec0_pol_theta, sec0_pol_phi, sec1_type, sec1_k, log_sec1_k, sec1_theta, sec1_phi, sec1_pol_theta, sec1_pol_phi, sec2_type, sec2_k, log_sec2_k, sec2_theta, sec2_phi, sec2_pol_theta, sec2_pol_phi, sec3_type, sec3_k, log_sec3_k, sec3_theta, sec3_phi, sec3_pol_theta, sec3_pol_phi, cpu_time");
 
 G4double energy(options.minEnergy);
 G4double stpEnergy(std::pow(options.maxEnergy/energy, 1./static_cast<G4double>(options.nEnergySteps-1)));
 G4int step(options.nEnergySteps);
 
 G4VLowEnergyTestableDiscreteProcess * process(GetSelectedProcess(options));
 clock_t time;
 
 char strBuffer[20];
 
 ProgressBar(0);
 while (step>0)
 {
  G4Step * aStep(GenerateStep(options));
  G4Track * aTrack(aStep->GetTrack());
  const G4DynamicParticle * aParticle(aTrack->GetDynamicParticle());
  G4ThreeVector vector;
   
  if (!options.randomEnergy)
  {
   aTrack->SetKineticEnergy(energy);
   energy*=stpEnergy;
  }
  ProgressBar(step);
  step--;

  G4int iteration(options.nIterations);
  
  while (iteration>0)
  {
   iteration--;
   
   aStep->SetStepLength(1*micrometer);
      
   iTuple->fill(iTuple->findColumn("iteration"), iteration);
   iTuple->fill(iTuple->findColumn("step"), aStep->GetStepLength()/m);

   iTuple->fill(iTuple->findColumn("in_k"), aParticle->GetKineticEnergy()/eV);
   iTuple->fill(iTuple->findColumn("log_in_k"), std::log10(aParticle->GetKineticEnergy()/eV));
   vector=aParticle->GetMomentumDirection();
   iTuple->fill(iTuple->findColumn("in_theta"), vector.theta());
   iTuple->fill(iTuple->findColumn("in_phi"), vector.phi());
   vector=aParticle->GetPolarization();
   iTuple->fill(iTuple->findColumn("in_pol_theta"), vector.theta());
   iTuple->fill(iTuple->findColumn("in_pol_phi"), vector.phi());
   
   time=clock();
   G4ParticleChange * particleChange(dynamic_cast<G4ParticleChange *>(process->PostStepDoIt(*aTrack, *aStep)));
   time=clock()-time;
   
   aTrack->SetKineticEnergy(particleChange->GetEnergy());
   aTrack->SetMomentumDirection(*particleChange->GetMomentumDirection());
   aTrack->SetPolarization(*particleChange->GetPolarization());
   
   iTuple->fill(iTuple->findColumn("e_deposit"), particleChange->GetLocalEnergyDeposit()/eV);
   iTuple->fill(iTuple->findColumn("log_e_deposit"), std::log10(particleChange->GetLocalEnergyDeposit()/eV));
   iTuple->fill(iTuple->findColumn("trk_status"), particleChange->GetTrackStatus());
   
   iTuple->fill(iTuple->findColumn("out_k"), aParticle->GetKineticEnergy()/eV);
   iTuple->fill(iTuple->findColumn("log_out_k"), std::log10(aParticle->GetKineticEnergy()/eV));
   vector=aParticle->GetMomentumDirection();
   iTuple->fill(iTuple->findColumn("out_theta"), vector.theta());
   iTuple->fill(iTuple->findColumn("out_phi"), vector.phi());
   vector=aParticle->GetPolarization();
   iTuple->fill(iTuple->findColumn("out_pol_theta"), vector.theta());
   iTuple->fill(iTuple->findColumn("out_pol_phi"), vector.phi());

   G4int n(particleChange->GetNumberOfSecondaries());
   iTuple->fill(iTuple->findColumn("n_secondaries"), n);
   
   while (n>0)
   {
    n--;
    G4Track * aSecTrack(particleChange->GetSecondary(n));
    const G4DynamicParticle * aSecParticle(aSecTrack->GetDynamicParticle());
    
    sprintf(strBuffer, "sec%d_type", n);
    iTuple->fill(iTuple->findColumn(strBuffer), GetParticleIndex(aSecParticle->GetDefinition()));
    sprintf(strBuffer, "sec%d_k", n);
    iTuple->fill(iTuple->findColumn(strBuffer), aSecParticle->GetKineticEnergy()/eV);
    sprintf(strBuffer, "log_sec%d_k", n);
    iTuple->fill(iTuple->findColumn(strBuffer), std::log10(aSecParticle->GetKineticEnergy()/eV));
    vector=aSecParticle->GetMomentumDirection();
    sprintf(strBuffer, "sec%d_theta", n);
    iTuple->fill(iTuple->findColumn(strBuffer), vector.theta());
    sprintf(strBuffer, "sec%d_phi", n);
    iTuple->fill(iTuple->findColumn(strBuffer), vector.phi());
    vector=aSecParticle->GetPolarization();
    sprintf(strBuffer, "sec%d_pol_theta", n);
    iTuple->fill(iTuple->findColumn(strBuffer), vector.theta());
    sprintf(strBuffer, "sec%d_pol_phi", n);
    iTuple->fill(iTuple->findColumn(strBuffer), vector.phi());
    
    delete aSecTrack;
   }
   
   iTuple->fill(iTuple->findColumn("cpu_time"), static_cast<G4double>(time)/static_cast<G4double>(CLOCKS_PER_SEC));
   iTuple->addRow();

   particleChange->Clear();
  }

  delete aTrack;
  delete aStep;
 }
}

//! \brief Main function
//! \param argc Number of arguments
//! \param argv Pointer to the arguments
//! \return The exit value 
int main(int argc, char ** argv)
{
 struct Options options;
 processOptions(argc, argv, &options);
 
 CreateMaterials();
 
 GetSelectedProcess(options);
 GetSelectedMaterial(options);
 
 G4RunManager* rm = new G4RunManager();
 rm->GeometryHasBeenModified();
 G4VPhysicalVolume * world(CreateGeometry(options));
 rm->DefineWorldVolume(world);
 G4cout << "[OK] World is defined " << G4endl;
 
 SetPhysics(options, world);
 
 if (!(options.meanFreePathTest || options.postStepDoItTest))
 {
  G4cout << "[OK] Program completed" << G4endl;
  return 0;
 }
 
 // HBOOK initialization 
 AIDA::IAnalysisFactory * analysisFactory(AIDA_createAnalysisFactory());
 AIDA::ITreeFactory * treeFactory(analysisFactory->createTreeFactory());
 AIDA::ITree * tree(treeFactory->create(options.outputFileName, "hbook", false, true));
 G4cout << "[OK] Tree store: " << tree->storeName() << G4endl;
 
 AIDA::ITupleFactory * tupleFactory(analysisFactory->createTupleFactory(*tree));
 
 // Mean free path test
 if (options.meanFreePathTest)
 {
  G4cout << "[OK] Mean free path test started" << G4endl;
  MeanFreePathTest(tupleFactory, options);
  G4cout << "[OK] Mean free path test completed" << G4endl;
 }
 
 // Post step do it test
 if (options.postStepDoItTest)
 {
  G4cout << "[OK] Post step do it test started" << G4endl;
  PostStepDoItTest(tupleFactory, options);
  G4cout << "[OK] Post step do it test completed" << G4endl;
 }
 
 G4cout << "[OK] Storing analysis data" << G4endl;
 tree->commit();
 tree->close();
 
 G4cout << "[OK] Deleting analysis data" << G4endl;
 delete tupleFactory;
 delete tree;
 delete treeFactory;
 delete analysisFactory;

 G4cout << "[OK] Program completed" << G4endl;
 return 0;
}
