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
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLParticleType.hh"
#include "G4INCLConfig.hh"
#include "G4INCLParticleSpecies.hh"
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

  Config::Config()
  {
    init();
  }

  Config::Config(G4int /*A*/, G4int /*Z*/, G4INCL::ParticleSpecies proj, G4double projectileE)
  {
    init();
    projectileSpecies = proj;
    projectileKineticEnergy = projectileE;
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
        ("impact-parameter", boost::program_options::value<G4double>(&impactParameter)->default_value(-1.), "impact parameter")
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
        ("version", "print version string and exit")
        ;

      // Run-specific options
      boost::program_options::options_description runOptDesc("Run options");
      runOptDesc.add_options()
        ("title", boost::program_options::value<std::string>(&title)->default_value("INCL default run title"), "run title")
        ("output,o", boost::program_options::value<std::string>(&outputFileRoot), "root for generating output file names. Suffixes (.root, .out, etc.) will be appended to this root. Defaults to the input file name, if given; otherwise, defaults to a string composed of the explicitly specified options")
        ("logfile,l", boost::program_options::value<std::string>(&logFileName), "log file name. Defaults to `<output_root>.log'. Use `-' if you want to redirect logging to stdout")
        ("number-shots,N", boost::program_options::value<G4int>(&nShots), "* number of shots")
        ("target,t", boost::program_options::value<std::string>(&targetString), "* target nuclide. Can be specified as Fe56, 56Fe, Fe-56, 56-Fe, Fe_56, 56_Fe or Fe. If the mass number is omitted, natural target composition is assumed.")
        ("projectile,p", boost::program_options::value<std::string>(&projectileString), "* projectile name:\n"
         "  \tproton, p\n"
         "  \tneutron, n\n"
         "  \tpi+, piplus, pion+, pionplus\n"
         "  \tpi0, pizero, pion0, pionzero\n"
         "  \tpi-, piminus, pion-, pionminus\n"
         "  \td, t, a, deuteron, triton, alpha\n"
         "  \tHe-4, He4, 4He (and so on)\n")
        ("energy,E", boost::program_options::value<G4float>(&projectileKineticEnergy), "* total kinetic energy of the projectile, in MeV")
        ("verbose-event", boost::program_options::value<G4int>(&verboseEvent)->default_value(-1), "request verbose logging for the specified event only")
        ("random-seed-1", boost::program_options::value<G4int>(&randomSeed1)->default_value(666), "first seed for the random-number generator")
        ("random-seed-2", boost::program_options::value<G4int>(&randomSeed2)->default_value(777), "second seed for the random-number generator")
        ("inclxx-datafile-path", boost::program_options::value<std::string>(&INCLXXDataFilePath)->default_value("./data/"))
#ifdef INCL_DEEXCITATION_ABLAXX
        ("ablav3p-cxx-datafile-path", boost::program_options::value<std::string>(&ablav3pCxxDataFilePath)->default_value("./de-excitation/ablaxx/data/G4ABLA3.0/"))
#endif
#ifdef INCL_DEEXCITATION_ABLA07
        ("abla07-datafile-path", boost::program_options::value<std::string>(&abla07DataFilePath)->default_value("./de-excitation/abla07/upstream/tables/"))
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
        ("geminixx-datafile-path", boost::program_options::value<std::string>(&geminixxDataFilePath)->default_value("./de-excitation/geminixx/upstream/"))
#endif
        ("verbosity,v", boost::program_options::value<G4int>(&verbosity)->default_value(4), verbosityDescription.str().c_str())
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
        ("potential", boost::program_options::value<std::string>(&potentialString)->default_value("isospin-energy"), "nucleon potential:\n  \tisospin-energy-smooth\n  \tisospin-energy (default)\n  \tisospin\n  \tconstant")
        ("pion-potential", boost::program_options::value<G4bool>(&pionPotential)->default_value("true"), "whether to use a pion potential:\n  \ttrue, 1 (default)\n  \tfalse, 0")
        ("local-energy-BB", boost::program_options::value<std::string>(&localEnergyBBString)->default_value("first-collision"), "local energy in baryon-baryon collisions:\n  \talways\n  \tfirst-collision (default)\n  \tnever")
        ("local-energy-pi", boost::program_options::value<std::string>(&localEnergyPiString)->default_value("first-collision"), "local energy in pi-N collisions and in delta decays:\n  \talways\n  \tfirst-collision (default)\n  \tnever")
        ("de-excitation", boost::program_options::value<std::string>(&deExcitationString)->default_value("none"), "which de-excitation model to use:"
         "\n  \tnone (default)"
#ifdef INCL_DEEXCITATION_ABLAXX
         "\n  \tABLAv3p"
#endif
#ifdef INCL_DEEXCITATION_ABLA07
         "\n  \tABLA07"
#endif
#ifdef INCL_DEEXCITATION_SMM
         "\n  \tSMM"
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
         "\n  \tGEMINIXX"
#endif
         )
        ("cluster-algorithm", boost::program_options::value<std::string>(&clusterAlgorithmString)->default_value("intercomparison"), "clustering algorithm for production of composites:\n  \tintercomparison (default)\n  \tnone")
        ("cluster-max-mass", boost::program_options::value<G4int>(&clusterMaxMass)->default_value(8), "maximum mass of produced composites:\n  \tminimum 2\n  \tmaximum 12")
        ("back-to-spectator", boost::program_options::value<G4bool>(&backToSpectator)->default_value("true"), "whether to use back-to-spectator:\n  \ttrue, 1 (default)\n  \tfalse, 0")
        ("use-real-masses", boost::program_options::value<G4bool>(&useRealMasses)->default_value("true"), "whether to use real masses for the outgoing particle energies:\n  \ttrue, 1 (default)\n  \tfalse, 0")
        ("separation-energies", boost::program_options::value<std::string>(&separationEnergyString)->default_value("INCL"), "how to assign the separation energies of the INCL nucleus:\n  \tINCL (default)\n  \treal\n  \treal-light")
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
      // name on the command line, it should be interpreted as an input-file
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
          boost::program_options::parsed_options parsedOptions = boost::program_options::parse_config_file(inputFileStream, configFileOptions, true);

          // Make sure that the unhandled options are all "*-datafile-path"
          std::vector<std::string> unhandledOptions =
            boost::program_options::collect_unrecognized(parsedOptions.options, boost::program_options::exclude_positional);
          G4bool ignoreNext = false;
          const std::string match = "-datafile-path";
          for(std::vector<std::string>::const_iterator i=unhandledOptions.begin(); i!=unhandledOptions.end(); ++i) {
            if(ignoreNext) {
              ignoreNext=false;
              continue;
            }
            if(i->rfind(match) == i->length()-match.length()) {
              std::cerr << "Ignoring unrecognized option " << *i << std::endl;
              ignoreNext = true;
            } else {
              std::cerr << "Error: unrecognized option " << *i << std::endl;
              std::cerr << suggestHelpMsg;
              std::exit(EXIT_FAILURE);
            }
          }

          // Store the option values in the variablesMap
          boost::program_options::store(parsedOptions, variablesMap);
          boost::program_options::notify(variablesMap);
        }
        inputFileStream.close();
      }

      // Process the options from the user-specific config file ~/.inclxxrc
      std::string configFileName;
      const char * const configFileVar = getenv("INCLXXRC");
      if(configFileVar)
        configFileName = configFileVar;
      else {
        const char * const homeDirectoryPointer = getenv("HOME");
        if(homeDirectoryPointer) { // Check if we can find the home directory
          std::string homeDirectory(homeDirectoryPointer);
          configFileName = homeDirectory + "/.inclxxrc";
        } else {
          std::cerr << "Could not determine the user's home directory. "
            << "Are you running Linux, Unix or BSD?"<< std::endl;
          std::exit(EXIT_FAILURE);
        }
      }

      std::ifstream configFileStream(configFileName.c_str());
      std::cout << "Reading config file " << configFileName << std::endl;
      if(!configFileStream) {
        std::cerr << "INCL++ config file " << configFileName
          << " not found. Continuing the run regardless."
          << std::endl;
      } else {
        // Merge options from the input file
        boost::program_options::parsed_options parsedOptions = boost::program_options::parse_config_file(configFileStream, configFileOptions, true);
        boost::program_options::store(parsedOptions, variablesMap);

        // Make sure that the unhandled options are all "*-datafile-path"
        std::vector<std::string> unhandledOptions =
          boost::program_options::collect_unrecognized(parsedOptions.options, boost::program_options::exclude_positional);
        G4bool ignoreNext = false;
        const std::string match = "-datafile-path";
        for(std::vector<std::string>::const_iterator i=unhandledOptions.begin(); i!=unhandledOptions.end(); ++i) {
          if(ignoreNext) {
            ignoreNext=false;
            continue;
          }
          if(i->rfind(match) == i->length()-match.length()) {
            std::cerr << "Ignoring unrecognized option " << *i << std::endl;
            ignoreNext = true;
          } else {
            std::cerr << "Error: unrecognized option " << *i << std::endl;
            std::cerr << suggestHelpMsg;
            std::exit(EXIT_FAILURE);
          }
        }

        // Store the option values in the variablesMap
        boost::program_options::store(parsedOptions, variablesMap);
        boost::program_options::notify(variablesMap);
      }
      configFileStream.close();

      /* *******************
       * Process the options
       * *******************/

      // -h/--help: print the help message and exit successfully
      if(variablesMap.count("help")) {
        std::cout << "Usage: INCLCascade [options] <input_file>" << std::endl;
        std::cout << std::endl << "Options marked with a * are compulsory, i.e. they must be provided either on\nthe command line or in the input file." << std::endl;
        std::cout << visibleOptions << std::endl;
        std::exit(EXIT_SUCCESS);
      }

      // --version: print the version string and exit successfully
      if(variablesMap.count("version")) {
        std::cout <<"INCL++ version " << getVersionID() << std::endl;
        std::exit(EXIT_SUCCESS);
      }

      // Check if the required options are present
      if(isFullRun) {
        std::string missingOption("");
        if(!variablesMap.count("number-shots"))
          missingOption = "number-shots";
        else if(!variablesMap.count("target"))
          missingOption = "target";
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

      // -p/--projectile: projectile species
      projectileSpecies = ParticleSpecies(projectileString);
      if(projectileSpecies.theType == G4INCL::UnknownParticle && isFullRun) {
        std::cerr << "Error: unrecognized particle type " << projectileString << std::endl;
        std::cerr << suggestHelpMsg;
        std::exit(EXIT_FAILURE);
      }

      // -t/--target: target species
      if(variablesMap.count("target")) {
        targetSpecies = ParticleSpecies(targetString);
        if(targetSpecies.theType!=Composite) {
          std::cerr << "Unrecognized target. You specified: " << targetString << std::endl
            << "  The target nuclide must be specified in one of the following forms:" << std::endl
            << "    Fe56, 56Fe, Fe-56, 56-Fe, Fe_56, 56_Fe, Fe" << std::endl
            << "  You can also use IUPAC element names (such as Uuh)." << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
        if(targetSpecies.theA==0)
          naturalTarget = true;
      }

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
            << "  non-relativistic-heavy-ion (default)" << std::endl
            << "  non-relativistic" << std::endl
            << "  none" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      }

      // --potential
      if(variablesMap.count("potential")) {
        std::string potentialNorm = potentialString;
        std::transform(potentialNorm.begin(), potentialNorm.end(), potentialNorm.begin(), ::tolower);
        if(potentialNorm=="isospin-energy-smooth") {
          potentialType = IsospinEnergySmoothPotential;
        } else if(potentialNorm=="isospin-energy") {
          potentialType = IsospinEnergyPotential;
        } else if(potentialNorm=="isospin")
          potentialType = IsospinPotential;
        else if(potentialNorm=="constant")
          potentialType = ConstantPotential;
        else {
          std::cerr << "Unrecognized potential type. Must be one of:" << std::endl
            << "  isospin-energy-smooth" << std::endl
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
            << "  first-collision" << std::endl
            << "  never (default)" << std::endl;
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
#ifdef INCL_DEEXCITATION_ABLAXX
        else if(deExcitationNorm=="ablav3p")
          deExcitationType = DeExcitationABLAv3p;
#endif
#ifdef INCL_DEEXCITATION_ABLA07
        else if(deExcitationNorm=="abla07")
          deExcitationType = DeExcitationABLA07;
#endif
#ifdef INCL_DEEXCITATION_SMM
        else if(deExcitationNorm=="smm")
          deExcitationType = DeExcitationSMM;
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
        else if(deExcitationNorm=="geminixx")
          deExcitationType = DeExcitationGEMINIXX;
#endif
        else {
          std::cerr << "Unrecognized de-excitation model. "
            << "Must be one of:" << std::endl
            << "  none (default)" << std::endl
#ifdef INCL_DEEXCITATION_ABLAXX
            << "  ABLAv3p" << std::endl
#endif
#ifdef INCL_DEEXCITATION_ABLA07
            << "  ABLA07" << std::endl
#endif
#ifdef INCL_DEEXCITATION_SMM
            << "  SMM" << std::endl
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
            << "  GEMINIXX" << std::endl
#endif
            ;
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
        else if(clusterAlgorithmNorm=="intercomparison")
          clusterAlgorithmType = IntercomparisonClusterAlgorithm;
        else {
          std::cerr << "Unrecognized cluster algorithm. "
            << "Must be one of:" << std::endl
            << "  intercomparison (default)" << std::endl
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

      // --separation-energies
      if(variablesMap.count("separation-energies")) {
        std::string separationEnergyNorm = separationEnergyString;
        std::transform(separationEnergyNorm.begin(),
            separationEnergyNorm.end(),
            separationEnergyNorm.begin(), ::tolower);
        if(separationEnergyNorm=="incl")
          separationEnergyType = INCLSeparationEnergy;
        else if(separationEnergyNorm=="real")
          separationEnergyType = RealSeparationEnergy;
        else if(separationEnergyNorm=="real-light")
          separationEnergyType = RealForLightSeparationEnergy;
        else {
          std::cerr << "Unrecognized separation-energies option. "
            << "Must be one of:" << std::endl
            << "  INCL (default)" << std::endl
            << "  real" << std::endl
            << "  real-light" << std::endl;
          std::cerr << suggestHelpMsg;
          std::exit(EXIT_FAILURE);
        }
      } else {
        separationEnergyType = INCLSeparationEnergy;
      }

      // --output: construct a reasonable output file root if not specified
      if(!variablesMap.count("output") && isFullRun) {
        // If an input file was specified, use its name as the output file root
        if(variablesMap.count("input-file"))
          outputFileRoot = inputFileName;
        else {
          std::stringstream outputFileRootStream;
          outputFileRootStream.precision(0);
          outputFileRootStream.setf(std::ios::fixed, std::ios::floatfield);
          outputFileRootStream <<
            ParticleTable::getShortName(projectileSpecies) << "_" <<
            ParticleTable::getShortName(targetSpecies) << "_" <<
            projectileKineticEnergy;

          // Append suffixes to the output file root for each explicitly specified CLI option
          typedef boost::program_options::variables_map::const_iterator BPOVMIter;
          for(BPOVMIter i=variablesMap.begin(); i!=variablesMap.end(); ++i) {
            std::string const &name = i->first;
            // Only process CLI options
            if(name!="projectile"
                && name!="target"
                && name!="energy"
                && name!="number-shots"
                && name!="random-seed-1"
                && name!="random-seed-2"
                && name!="inclxx-datafile-path"
#ifdef INCL_DEEXCITATION_ABLA07
                && name!="abla07-datafile-path"
#endif
#ifdef INCL_DEEXCITATION_ABLAXX
                && name!="ablav3p-cxx-datafile-path"
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
                && name!="geminixx-datafile-path"
#endif
                ) {
              boost::program_options::variable_value v = i->second;
              if(!v.defaulted()) {
                const std::type_info &type = v.value().type();
                if(type==typeid(std::string))
                  outputFileRootStream << "_" << name << "=" << v.as<std::string>();
                else if(type==typeid(G4float))
                  outputFileRootStream << "_" << name << "=" << v.as<G4float>();
                else if(type==typeid(G4int))
                  outputFileRootStream << "_" << name << "=" << v.as<G4int>();
                else if(type==typeid(G4bool))
                  outputFileRootStream << "_" << name << "=" << v.as<G4bool>();
              }
            }
          }

          outputFileRoot = outputFileRootStream.str();
        }
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

  }
#else
    Config::Config(G4int /*argc*/, char * /*argv*/ [], G4bool /*isFullRun*/)
    {
      init();
    }
#endif

  Config::~Config()
  {}

  void Config::init() {
      verbosity = 1;
      inputFileName = "";
      title = "INCL default run title";
      nShots = 1000;
      naturalTarget = false;
      projectileString = "proton";
      projectileSpecies = G4INCL::Proton;
      projectileKineticEnergy = 1000.0;
      verboseEvent = -1;
      randomSeed1 = 666;
      randomSeed2 = 777;
      pauliString = "strict-statistical";
      pauliType = StrictStatisticalPauli;
      CDPP = true;
      coulombString = "non-relativistic";
      coulombType = NonRelativisticCoulomb;
      potentialString = "isospin-energy";
      potentialType = IsospinEnergyPotential;
      pionPotential = true;
      localEnergyBBString = "first-collision";
      localEnergyBBType = FirstCollisionLocalEnergy;
      localEnergyPiString = "first-collision";
      localEnergyPiType = FirstCollisionLocalEnergy;
      deExcitationString = "none";
      deExcitationType = DeExcitationNone;
      clusterAlgorithmString = "intercomparison";
      clusterAlgorithmType = IntercomparisonClusterAlgorithm;
      clusterMaxMass = 8;
      backToSpectator = true;
      useRealMasses = true;
      impactParameter = -1.;
      separationEnergyString = "INCL";
      separationEnergyType = INCLSeparationEnergy;
  }

  std::string Config::summary() {
    std::stringstream message;
    message << "INCL++ version " << getVersionID() << std::endl;
    if(projectileSpecies.theType != Composite)
      message << "Projectile: " << ParticleTable::getName(projectileSpecies) << std::endl;
    else
      message << "Projectile: composite, A=" << projectileSpecies.theA << ", Z=" << projectileSpecies.theZ << std::endl;
    message << "  energy = " << projectileKineticEnergy << std::endl;
    if(targetSpecies.theA>0)
      message << "Target: A = " << targetSpecies.theA << " Z = " << targetSpecies.theZ << std::endl;
    else
      message << "Target: natural isotopic composition, Z = " << targetSpecies.theZ << std::endl;
    message << "Number of requested shots = " << nShots << std::endl;
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
      << "output = " << outputFileRoot << "\t# root for generating output file names. Suffixes (.root, .out, etc.) will be appended to this root. Defaults to the input file name, if given; otherwise, defaults to a string composed of the explicitly specified options" << std::endl
      << "logfile = " << logFileName << "\t# log file name. Defaults to `<output_root>.log'. Use `-' if you want to redirect logging to stdout" << std::endl
      << "number-shots = " << nShots << "\t# * number of shots" << std::endl
      << "inclxx-datafile-path = " << INCLXXDataFilePath << std::endl
#ifdef INCL_DEEXCITATION_ABLAXX
      << "ablav3p-cxx-datafile-path = " << ablav3pCxxDataFilePath << std::endl
#endif
#ifdef INCL_DEEXCITATION_ABLA07
      << "abla07-datafile-path = " << abla07DataFilePath << std::endl
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
      << "geminixx-datafile-path = " << geminixxDataFilePath << std::endl
#endif
      << std::endl << "# Projectile and target definitions" << std::endl
      << "target = " << targetString << "\t# * target nuclide. Can be specified as Fe56, 56Fe, Fe-56, 56-Fe, Fe_56, 56_Fe or Fe. If the mass number is omitted, natural target composition is assumed." << std::endl
      << "         " << "# the target nuclide was parsed as Z=" << targetSpecies.theZ;
    if(targetSpecies.theA>0)
      ss << ", A=" << targetSpecies.theA;
    else
      ss << ", natural target";
    ss << std::endl
      << "projectile = " << projectileString << "\t# * projectile name (proton, neutron, pi+, pi0, pi-, d, t, a, He-4...)" << std::endl
      << "         " << "# the projectile nuclide was parsed as Z=" << projectileSpecies.theZ << ", A=" << projectileSpecies.theA << std::endl
      << "energy = " << projectileKineticEnergy << "\t# * total kinetic energy of the projectile, in MeV" << std::endl
      << std::endl << "# Physics options " << std::endl
      << "pauli = " << pauliString << "\t# Pauli-blocking algorithm. Must be one of: strict-statistical (default), strict, statistical, global, none" << std::endl
      << "cdpp = " << CDPP << "\t# whether to apply CDPP after collisions" << std::endl
      << "coulomb = " << coulombString << "\t# Coulomb-distortion algorithm. Must be one of: non-relativistic (default), none" << std::endl
      << "potential = " << potentialString << "\t# nucleon potential. Must be one of: isospin-energy-smooth, isospin-energy (default), isospin, constant" << std::endl
      << "pion-potential = " << pionPotential << "\t# whether to use a pion potential" << std::endl
      << "local-energy-BB = " << localEnergyBBString << "\t# local energy in baryon-baryon collisions. Must be one of: always, first-collision (default), never" << std::endl
      << "local-energy-pi = " << localEnergyPiString << "\t# local energy in pi-N collisions and in delta decays. Must be one of: always, first-collision (default), never" << std::endl
      << "de-excitation = " << deExcitationString << "\t # which de-excitation model to use. Must be one of:"
      " none (default)"
#ifdef INCL_DEEXCITATION_ABLAXX
      ", ABLAv3p"
#endif
#ifdef INCL_DEEXCITATION_ABLA07
      ", ABLA07"
#endif
#ifdef INCL_DEEXCITATION_SMM
      ", SMM"
#endif
#ifdef INCL_DEEXCITATION_GEMINIXX
      ", GEMINIXX"
#endif
      << std::endl
      << "cluster-algorithm = " << clusterAlgorithmString << "\t# clustering algorithm for production of composites. Must be one of: intercomparison (default), none" << std::endl
      << "cluster-max-mass = " << clusterMaxMass << "\t# maximum mass of produced composites. Must be between 2 and 12 (included)" << std::endl
      << "back-to-spectator = " << backToSpectator << "\t# whether to use back-to-spectator" << std::endl
      << "use-real-masses = " << useRealMasses << "\t# whether to use real masses for the outgoing particle energies" << std::endl
      << "separation-energies = " << separationEnergyString << "\t# how to assign the separation energies of the INCL nucleus. Must be one of: INCL (default), real, real-light" << std::endl
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
