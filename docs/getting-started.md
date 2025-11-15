# Getting Started with Geant4

## Overview

Geant4 is a C++ toolkit for simulating the passage of particles through matter. This guide will help you get started with building and using Geant4.

## System Requirements

- **CMake**: Version 3.16 or higher (up to 3.27)
- **C++ Compiler**: Supporting C++17 or later
- **Operating Systems**: Linux, macOS, Windows

## Installation

Geant4 uses CMake as its build system and requires an out-of-source build.

### Basic Build Steps

1. **Download the source code**
   ```bash
   git clone https://github.com/Geant4/geant4.git
   cd geant4
   ```

2. **Create a build directory**
   ```bash
   mkdir build
   cd build
   ```

3. **Configure with CMake**
   ```bash
   cmake ..
   ```

4. **Build**
   ```bash
   make -j$(nproc)
   ```

5. **Install** (optional)
   ```bash
   make install
   ```

### Important Notes

- **Out-of-source builds are required**: Geant4 will not allow in-source builds
- The build system enforces this to keep the source tree clean
- See CMakeLists.txt:6-14 for the build enforcement logic

## Version Information

Current version: **11.4.0**

Version is defined in CMakeLists.txt:32-35:
- Major: 11
- Minor: 4
- Patch: 0

## Next Steps

- Review the [Architecture](./architecture) to understand Geant4's design
- Explore the [Source Modules](./reference/source-modules) for detailed component information
- Check the [official installation guide](http://cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html) for advanced configuration options

## Additional Resources

- **User Documentation**: [http://cern.ch/geant4/support/user_documentation](http://cern.ch/geant4/support/user_documentation)
- **User Forum**: [http://cern.ch/geant4-forum](http://cern.ch/geant4-forum)
- **Examples**: Check the `examples/` directory in the source tree

## Getting Help

If you have questions or issues:

1. **Search the Forum**: The [User Forum](http://cern.ch/geant4-forum) is fully searchable
2. **Ask Questions**: Post in the relevant category on the forum
3. **Report Bugs**: Use the [Bugzilla Tracker](http://bugzilla-geant4.kek.jp)
