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
// INCL++ intra-nuclear cascade model
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
// INCL++ revision: v5.0_rc3
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLParticleType.hh"
#include "G4INCLConfig.hh"
#include "G4INCLParticleType.hh"
#include "G4INCLParticleTable.hh"

#ifdef HAS_BOOST_PROGRAM_OPTIONS
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "G4INCLCascade.hh"
#include "G4INCLLogger.hh"
#endif

namespace G4INCL {

  Config::Config() {}

  Config::Config(G4int A, G4int Z, G4INCL::ParticleType proj, G4double projectileE)
    :verbosity(1), inputFileName(""),
     title("INCL default run title"), outputFileRoot("incl-output"),
     logFileName("-"), // Log to the standard output (uses G4cout if in G4)
     nShots(1000),              // Irrelevant in Geant4
     targetA(A), targetZ(Z), // Dummy target. G4INCLXXFactory will
     // change these.
     naturalTarget(false),
     projectileString("projectile"), projectileType(proj),
     projectileKineticEnergy(projectileE),
     verboseEvent(-1),
     randomSeed1(666), randomSeed2(777), // In Geant4 these seeds
					 // aren't actually used
     pauliString("strict-statistical"), pauliType(StrictStatisticalPauli),
     CDPP(true),
     coulombString("non-relativistic"), coulombType(NonRelativisticCoulomb),
     potentialString("isospin-energy"), potentialType(IsospinEnergyPotential),
     pionPotential(true),
     localEnergyBBString("first-collision"), localEnergyBBType(FirstCollisionLocalEnergy),
     localEnergyPiString("first-collision"), localEnergyPiType(FirstCollisionLocalEnergy),
     clusterAlgorithmString("G4intercomparison"), clusterAlgorithmType(IntercomparisonClusterAlgorithm),
     clusterMaxMass(8),
     backToSpectator(true),
     backToSpectatorThreshold(18.)
  {
    //    std::cout << echo() << std::endl;
  }

  // NOT used in Geant4 mode
#ifdef HAS_BOOST_PROGRAM_OPTIONS
  Config::Config(G4int argc, char *argv[], G4bool isFullRun) : naturalTarget(false) {
    const std::string suggestHelpMsg("You might want to run `INCLCascade -h' to get a help message.\n");

    try {

      // Hidden options
      boost::program_options::options_description hiddenOptDesc("Hidden options");
      hiddenOptDesc.add_options()
        ("input-file", boost::program_options::value<std::string>(&inputFileName), "input file")
        ;

      // Generic options
      std::stringstream verbosityDescription;
      verbosityDescription << "set verbosity level:\n"
        << "  0: \tquiet, suppress all output messages\n"
        << "  " << InfoMsg << ": \tminimal logging\n"
        << "  " << FatalMsg << ": \tlog fatal error messages as well\n"
        << "  " << ErrorMsg << ": \tlog error messages as well\n"
        << "  " << WarningMsg << ": \tlog warning messages as well\n"
        << "  " << DebugMsg << ": \tlog debug messages as well\n"
        << "  " << DataBlockMsg << ": \tlog data-block messages as well";

      boost::program_options::options_description genericOptDesc("Generic options");
      genericOptDesc.add_options()
        ("help,h", "produce this help message")
        ("version", "prG4int version string and exit")
        ("verbosity,v", boost::program_options::value<G4int>(&verbosity)->default_value(4), verbosityDescription.str().c_str())
        ;

      // Run-specific options
      boost::program_options::options_description runOptDesc("Run options");
      runOptDesc.add_options()
        ("title,t", boost::program_options::value<std::string>(&title)->default_value("INCL default run title"), "run title")
        ("output,o", boost::program_options::value<std::string>(&outputFileRoot)->default_value("incl-output"), "root for generating output file names. Suffixes (.root, .out, etc.) will be appended to this root. Defaults to the input file name")
        ("logfile,l", boost::program_options::value<std::string>(&logFileName), "log file name. Defaults to `<output_root>.log'. Use `-' if you want to redirect logging to stdout")
        ("number-shots,N", boost::program_options::value<G4int>(&nShots), "* number of shots")
        ("target-mass,A", boost::program_options::value<G4int>(&targetA), "mass number of the target. If omitted, natural target composition is assumed.")
        ("target-charge,Z", boost::program_options::value<G4int>(&targetZ), "* charge number of the target")
        ("projectile,p", boost::program_options::value<std::string>(&projectileString), "* projectile name:\n"
         "  \tproton, p\n"
         "  \tneutron, n\n"
         "  \tpi+, piplus, pion+, pionplus\n"
         "  \tpi0, pizero, pion0, pionzero\n"
         "  \tpi-, piminus, pion-, pionminus")
        ("energy,E", boost::program_options::value<G4float>(&projectileKineticEnergy), "* kinetic energy of the projectile, in MeV")
        ("verbose-event", boost::program_options::value<G4int>(&verboseEvent)->default_value(-1), "request verbose logging for the specified event only")
        ("random-seed-1", boost::program_options::value<G4int>(&randomSeed1)->default_value(666), "first seed for the random-number generator")
        ("random-seed-2", boost::program_options::value<G4int>(&randomSeed2)->default_value(777), "second seed for the random-number generator")
        ("ablav3p-cxx-datafile-path", boost::program_options::value<std::string>(&ablav3pCxxDataFilePath)->default_value("./ablaxx/data/G4ABLA3.0/"))
        ("abla07-datafile-path", boost::program_options::value<std::string>(&abla07DataFilePath)->default_value("./smop/tables/"))
        ;

      // Physics options
      boost::program_options::options_description physicsOptDesc("Physics options");
      physicsOptDesc.add_options()
        ("pauli", boost::program_options::value<std::string>(&pauliString)->default_value("strict-statistical"), "Pauli-blocking algorithm:\n"
         "  \tstrict-statistical (default)\n"
         "  \tstrict\n"
         "  \tstatistical\n"
         "  \tglobal\n"
         "  \tnone")
        ("cdpp", boost::program_options::value<G4bool>(&CDPP)->default_value(true), "whether to apply CDPP after collisions:\n  \ttrue, 1 (default)\n  \tfalse, 0")
        ("coulomb", boost::program_options::value<std::string>(&coulombString)->default_value("non-relativistic"), "Coulomb-distortion algorithm:\n  \tnon-relativistic (default)\n  \tnone")
        ("potential", boost::program_options::value<std::string>(&potentialString)->default_value("isospin-energy"), "nucleon potential:\n  \tisospin-energy (default)\n  \tisospin\n  \tconstant")
        ("pion-potential", boost::program_options::value<G4bool>(&pionPotential)->default_value("true"), "whether to use a pion potential:\n  \ttrue, 1 (default)\n  \tfalse, 0")
        ("local-energy-BB", boost::program_options::value<std::string>(&localEnergyBBString)->default_value("first-collision"), "local energy in baryon-baryon collisions:\n  \talways\n  \tfirst-collision (default)\n  \tnever")
        ("local-energy-pi", boost::program_options::value<std::string>(&localEnergyPiString)->default_value("first-collision"), "local energy in pi-N collisions and in delta decays:\n  \talways\n  \tfirst-collision (default)\n  \tnever")
        ("de-excitation", boost::program_options::value<std::string>(&deExcitationString)->default_value("none"), "which de-excitation model to use:\n  \tnone (default)\n  \tABLAv3p\n  \tABLA07")
        ("cluster-algorithm", boost::program_options::value<std::string>(&clusterAlgorithmString)->default_value("G4intercomparison"), "clustering algorithm for production of composites:\n  \tG4intercomparison (default)\n  \tnone")
        ("cluster-max-mass", boost::program_options::value<G4int>(&clusterMaxMass)->default_value(8), "maximum mass of produced composites:\n  \tminimum 2\n  \tmaximum 12")
        ("back-to-spectator", boost::program_options::value<G4bool>(&backToSpectator)->default_value("true"), "whether to use back-to-spectator:\n  \ttrue, 1 (default)\n  \tfalse, 0")
        ("back-to-spectator-threshold", boost::program_options::value<G4float>(&backToSpectatorThreshold)->default_value(18.), "threshold for back-to-spectator, in MeV, measured from the Fermi energy")
        ;

      // Select options allowed on the command line
      boost::program_options::options_description cmdLineOptions;
      cmdLineOptions.add(hiddenOptDesc).add(genericOptDesc).add(runOptDesc).add(physicsOptDesc);

      // Select options allowed in config files
      boost::program_options::options_description configFileOptions;
      configFileOptions.add(runOptDesc).add(physicsOptDesc);

      // Select visible options
      boost::program_options::options_description visibleOptions;
      visibleOptions.add(genericOptDesc).add(runOptDesc).add(physicsOptDesc);

      // Declare input-file as a positional option (if we just provide a file
      // name on the command line, it should be G4interpreted as an input-file
      // option).
      boost::program_options::positional_options_description p;
      p.add("input-file", 1);

      // Disable guessing of option names
      G4int cmdstyle =
        boost::program_options::command_line_style::default_style &
        ~boost::program_options::command_line_style::allow_guessing;

      // Result of the option processing
      boost::program_options::variables_map variablesMap;
      boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
          style(cmdstyle).
          options(cmdLineOptions).positional(p).run(), variablesMap);
      boost::program_options::notify(variablesMap);

      // If an input file was specified, merge the options with the command-line
      // options.
      if(variablesMap.count("input-file")) {
        std::ifstream inputFileStream(inputFileName.c_str());
        if(!inputFileStream) {
          std::cerr << "Cannot open input file: " << inputFileName << std::endl;
          std::exit(EXIT_FAILURE);
        } else {
          // Merge options from the input file
          boost::program_options::store(boost::program_options::parse_config_file(inputFileStream, configFileOptions), variablesMap);
          boost::program_options::notify(variablesMap);
          if(!variablesMap.count("output")) outputFileRoot = inputFileName;
        }
        inputFileStream.close();
      }

      // Process the options from the user-specific config file ~/.inclxxrc
      if(getenv("HOME")) { // Check if we can find the home directory
        std::string homeDirectory = getenv("HOME");
        std::string configFileName = homeDirectory + "/.inclxxrc";
        std::ifstream configFileStream(configFileName.c_str());
        std::cout << "Reading config file " << configFileName << std::endl;
        if(!configFileStream) {
          std::cerr << "INCL++ config file " << configFileName
            << " not found. Continuing the run regardless."
            << std::endl;
        } else {
          boost::program_options::store(boost::program_options::parse_config_file(configFileStream, configFileOptions), variablesMap);
          boost::program_options::notify(variablesMap);
        }
        configFileStream.close();
      } else {
        std::cerr << "Could not determine the user's home directory. "
          << "Are you running Linux, Unix or BSD?"<< std::endl;
        std::exit(EXIT_FAILURE);
      }

      /* *******************
       * Process the options
       * *******************/

      // -h/--help: prG4int the help message and exit successfully
      if(variablesMap.count("help")) {
        std::cout << "Usage: INCLCascade [options] <input_file>" << std::endl;
        std::cout << std::endl << "Options marked with a * are compulsory, i.e. they must be provided either on\nthe command line or in the input file." << std::endl;
        std::cout << visibleOptions << std::endl;
        std::exit(EXIT_SUCCESS);
      }

      // --version: prG4int the version string and exit successfully
      if(variablesMap.count("version")) {
        std::cout <<"INCL++ version " << getVersionID() << std::endl;
        std::exit(EXIT_SUCCESS);
      }

      // Check if the required options are present
      if(isFullRun) {
        std::string missingOption("");
        if(!variablesMap.count("number-shots"))
          missingOption = "number-shots";
        else if(!variablesMap.count("target-charge"))
          missingOption = "target-charge";
        else if(!variablesMap.count("projectile"))
          missingOption = "projectile";
        else if(!variablesMap.count("energy"))
          missingOption = "energy";
        if(!missingOption.empty()) {
          std::cerr << "Required option " << missingOption << " is missing." << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      } else {
        std::cout <<"Not performing a full run. This had better be a test..." << std::endl;
      }

      // -A/--target-mass: if not specified, assume natural target
      if(!variablesMap.count("target-mass"))
        naturalTarget=true;

      // --pauli
      if(variablesMap.count("pauli")) {
        std::string pauliNorm = pauliString;
        std::transform(pauliNorm.begin(), pauliNorm.end(), pauliNorm.begin(), ::tolower);
        if(pauliNorm=="statistical")
          pauliType = StatisticalPauli;
        else if(pauliNorm=="strict")
          pauliType = StrictPauli;
        else if(pauliNorm=="strict-statistical")
          pauliType = StrictStatisticalPauli;
        else if(pauliNorm=="global")
          pauliType = GlobalPauli;
        else if(pauliNorm=="none")
          pauliType = NoPauli;
        else {
          std::cerr << "Unrecognized Pauli-blocking algorithm. Must be one of:" << std::endl
            << "  strict-statistical (default)" << std::endl
            << "  strict" << std::endl
            << "  statistical" << std::endl
            << "  global" << std::endl
            << "  none" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      }

      // --coulomb
      if(variablesMap.count("coulomb")) {
        std::string coulombNorm = coulombString;
        std::transform(coulombNorm.begin(), coulombNorm.end(), coulombNorm.begin(), ::tolower);
        if(coulombNorm=="non-relativistic")
          coulombType = NonRelativisticCoulomb;
        else if(coulombNorm=="none")
          coulombType = NoCoulomb;
        else {
          std::cerr << "Unrecognized Coulomb-distortion algorithm. Must be one of:" << std::endl
            << "  non-relativistic (default)" << std::endl
            << "  none" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      }

      // --potential
      if(variablesMap.count("potential")) {
        std::string potentialNorm = potentialString;
        std::transform(potentialNorm.begin(), potentialNorm.end(), potentialNorm.begin(), ::tolower);
        if(potentialNorm=="isospin-energy") {
          potentialType = IsospinEnergyPotential;
        } else if(potentialNorm=="isospin")
          potentialType = IsospinPotential;
        else if(potentialNorm=="constant")
          potentialType = ConstantPotential;
        else {
          std::cerr << "Unrecognized potential type. Must be one of:" << std::endl
            << "  isospin-energy (default)" << std::endl
            << "  isospin" << std::endl
            << "  constant" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      }

      // --local-energy-BB
      if(variablesMap.count("local-energy-BB")) {
        std::string localEnergyBBNorm = localEnergyBBString;
        std::transform(localEnergyBBNorm.begin(), localEnergyBBNorm.end(), localEnergyBBNorm.begin(), ::tolower);
        if(localEnergyBBNorm=="always") {
          localEnergyBBType = AlwaysLocalEnergy;
        } else if(localEnergyBBNorm=="first-collision")
          localEnergyBBType = FirstCollisionLocalEnergy;
        else if(localEnergyBBNorm=="never")
          localEnergyBBType = NeverLocalEnergy;
        else {
          std::cerr << "Unrecognized local-energy-BB type. Must be one of:" << std::endl
            << "  always" << std::endl
            << "  first-collision (default)" << std::endl
            << "  never" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      }

      // --local-energy-pi
      if(variablesMap.count("local-energy-pi")) {
        std::string localEnergyPiNorm = localEnergyPiString;
        std::transform(localEnergyPiNorm.begin(), localEnergyPiNorm.end(), localEnergyPiNorm.begin(), ::tolower);
        if(localEnergyPiNorm=="always") {
          localEnergyPiType = AlwaysLocalEnergy;
        } else if(localEnergyPiNorm=="first-collision")
          localEnergyPiType = FirstCollisionLocalEnergy;
        else if(localEnergyPiNorm=="never")
          localEnergyPiType = NeverLocalEnergy;
        else {
          std::cerr << "Unrecognized local-energy-pi type. Must be one of:" << std::endl
            << "  always" << std::endl
            << "  first-collision (default)" << std::endl
            << "  never" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      }

      // --de-excitation
      if(variablesMap.count("de-excitation")) {
        std::string deExcitationNorm = deExcitationString;
        std::transform(deExcitationNorm.begin(),
            deExcitationNorm.end(),
            deExcitationNorm.begin(), ::tolower);
        if(deExcitationNorm=="none")
          deExcitationType = DeExcitationNone;
        else if(deExcitationNorm=="ablav3p")
          deExcitationType = DeExcitationABLAv3p;
        else if(deExcitationNorm=="abla07")
          deExcitationType = DeExcitationABLA07;
        else {
          std::cerr << "Unrecognized de-excitation model. "
            << "Must be one of:" << std::endl
            << "  none (default)" << std::endl
            << "  ABLAv3p" << std::endl
            << "  ABLA07" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      } else {
        deExcitationType = DeExcitationNone;
      }

      // --cluster-algorithm
      if(variablesMap.count("cluster-algorithm")) {
        std::string clusterAlgorithmNorm = clusterAlgorithmString;
        std::transform(clusterAlgorithmNorm.begin(),
            clusterAlgorithmNorm.end(),
            clusterAlgorithmNorm.begin(), ::tolower);
        if(clusterAlgorithmNorm=="none")
          clusterAlgorithmType = NoClusterAlgorithm;
        else if(clusterAlgorithmNorm=="G4intercomparison")
          clusterAlgorithmType = IntercomparisonClusterAlgorithm;
        else {
          std::cerr << "Unrecognized cluster algorithm. "
            << "Must be one of:" << std::endl
            << "  G4intercomparison (default)" << std::endl
            << "  none" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      } else {
        clusterAlgorithmType = IntercomparisonClusterAlgorithm;
      }

      // --cluster-max-mass
      if(variablesMap.count("cluster-max-mass") && clusterMaxMass < 2 && clusterMaxMass > 12) {
        std::cerr << "Maximum cluster mass outside the allowed range. Must be between 2 and 12 (included)"
          << std::endl
          << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
      }

      // -l/--logfile
      if(!variablesMap.count("logfile"))
        logFileName = outputFileRoot + ".log";

    }
    catch(std::exception& e)
    {
      std::cerr << e.what() << "\n";
      std::cerr << suggestHelpMsg;
      std::exit(EXIT_FAILURE);
    }

    // Process the rest of the options here:

    // Set projectileType:
    projectileType = ParticleTable::getParticleType(projectileString);
    if(projectileType == G4INCL::UnknownParticle && isFullRun) {
      std::cerr << "Error: unrecognized particle type " << projectileString << std::endl;
      std::cerr << suggestHelpMsg;
      std::exit(EXIT_FAILURE);
    }
  }
#else
    Config::Config(G4int /*argc*/, char ** /*argv*/, G4bool /*isFullRun*/) :
      verbosity(1), inputFileName(""),
      title("INCL default run title"), outputFileRoot("incl-output"),
      logFileName("incl-output.log"),
      nShots(1000),
      targetA(208), targetZ(82),
      naturalTarget(false),
      projectileString("proton"), projectileType(G4INCL::Proton), projectileKineticEnergy(1000.0),
      verboseEvent(-1),
      randomSeed1(666), randomSeed2(777),
      pauliString("strict-statistical"), pauliType(StrictStatisticalPauli),
      CDPP(true),
      coulombString("non-relativistic"), coulombType(NonRelativisticCoulomb),
      potentialString("isospin-energy"), potentialType(IsospinEnergyPotential),
      pionPotential(true),
      localEnergyBBString("first-collision"), localEnergyBBType(FirstCollisionLocalEnergy),
      localEnergyPiString("first-collision"), localEnergyPiType(FirstCollisionLocalEnergy),
      clusterAlgorithmString("G4intercomparison"), clusterAlgorithmType(IntercomparisonClusterAlgorithm),
      clusterMaxMass(8),
      backToSpectator(true), backToSpectatorThreshold(18.)
    {}
#endif

  Config::~Config()
  {}

  std::string Config::summary() {
    std::stringstream message;
    message << "INCL++ version " << getVersionID() << std::endl;
    message << "Projectile: " << ParticleTable::getName(projectileType) << std::endl;
    message << "  energy = " << projectileKineticEnergy << std::endl;
    message << "Target: A = " << targetA << " Z = " << targetZ << std::endl;
    message << "Number of shots = " << nShots << std::endl;
    return message.str();
  }

  std::string const Config::echo() const {
    std::stringstream ss;
    ss << std::boolalpha;
    ss << "###########################" << std::endl
      << "### Start of input echo ###" << std::endl
      << "###########################" << std::endl << std::endl
      << " # You may re-use this snippet of the log file as an input file!" << std::endl
      << " # Options marked with a * are compulsory." << std::endl
      << std::endl
      << "# Run options" << std::endl
      << "title = " << title << "\t# run title" << std::endl
      << "output = " << outputFileRoot << "\t# root for generating output file names. Suffixes (.root, .out, etc.) will be appended to this root. Defaults to the input file name" << std::endl
      << "logfile = " << logFileName << "\t# log file name. Defaults to <output_root>.log. Use `-' if you want to redirect logging to stdout" << std::endl
      << "number-shots = " << nShots << "\t# * number of shots" << std::endl
      << "ablav3p-cxx-datafile-path = " << ablav3pCxxDataFilePath << std::endl
      << "abla07-datafile-path = " << abla07DataFilePath << std::endl
      << std::endl << "# Projectile and target definitions" << std::endl
      << "target-mass = " << targetA << "\t# mass number of the target. If omitted, natural target composition is assumed" << std::endl
      << "target-charge = " << targetZ << "\t# * charge number of the target" << std::endl
      << "projectile = " << projectileString << "\t# * projectile name (proton, neutron, pi+, pi0, pi-...)" << std::endl
      << "energy = " << projectileKineticEnergy << "\t# * kinetic energy of the projectile, in MeV" << std::endl
      << std::endl << "# Physics options " << std::endl
      << "pauli = " << pauliString << "\t# Pauli-blocking algorithm. Must be one of: strict-statistical (default), strict, statistical, global, none" << std::endl
      << "cdpp = " << CDPP << "\t# whether to apply CDPP after collisions" << std::endl
      << "coulomb = " << coulombString << "\t# Coulomb-distortion algorithm. Must be one of: non-relativistic (default), none" << std::endl
      << "potential = " << potentialString << "\t# nucleon potential. Must be one of: isospin-energy (default), isospin, constant" << std::endl
      << "pion-potential = " << pionPotential << "\t# whether to use a pion potential" << std::endl
      << "local-energy-BB = " << localEnergyBBString << "\t# local energy in baryon-baryon collisions. Must be one of: always, first-collision (default), never" << std::endl
      << "local-energy-pi = " << localEnergyPiString << "\t# local energy in pi-N collisions and in delta decays. Must be one of: always, first-collision (default), never" << std::endl
      << "de-excitation = " << deExcitationString << "\t # which de-excitation model to use. Must be one of: none (default), ABLAv3p, ABLA07" << std::endl
      << "cluster-algorithm = " << clusterAlgorithmString << "\t# clustering algorithm for production of composites. Must be one of: G4intercomparison (default), none" << std::endl
      << "cluster-max-mass = " << clusterMaxMass << "\t# maximum mass of produced composites. Must be between 2 and 12 (included)" << std::endl
      << "back-to-spectator = " << backToSpectator << "\t# whether to use back-to-spectator" << std::endl
      << "back-to-spectator-threshold = " << backToSpectatorThreshold << "\t# threshold for back-to-spectator, in MeV, measured from the Fermi energy" << std::endl
      << std::endl << "# Technical options " << std::endl
      << "verbosity = " << verbosity << "\t# from 0 (quiet) to 10 (most verbose)" << std::endl
      << "verbose-event = " << verboseEvent << "\t# request verbose logging for the specified event only" << std::endl
      << "random-seed-1 = " << randomSeed1 << "\t# first seed for the random-number generator" << std::endl
      << "random-seed-2 = " << randomSeed2 << "\t# second seed for the random-number generator" << std::endl
      << std::endl << "#########################" << std::endl
      << "### End of input echo ###" << std::endl
      << "#########################" << std::endl;

    return ss.str();
  }

}
