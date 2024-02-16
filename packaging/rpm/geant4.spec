%global geant4_version 11.2.1

%global NEUTRONHPDATA G4NDL.4.7
%global LEDATA G4EMLOW.8.5
%global LEVELGAMMADATA G4PhotonEvaporation.5.7
%global RADIOACTIVEDATA G4RadioactiveDecay.5.6
%global PARTICLEXSDATA G4PARTICLEXS.4.0
%global PIIDATA G4PII.1.3
%global REALSURFACEDATA G4RealSurface.2.2
%global SAIDXSDATA G4SAIDDATA.2.0
%global ABLADATA G4ABLA.3.3
%global INCLDATA G4INCL.1.2
%global ENSDFSTATEDATA G4ENSDFSTATE.2.3
%global TENDLDATA G4TENDL.1.4

Name: geant4
Version: %{geant4_version}
Release: 1%{?dist}
Summary: Toolkit for the simulation of the passage of particles through matter.
License: Geant4
URL: https://geant4.web.cern.ch
Source0:  https://geant4-data.web.cern.ch/releases/%{name}-v%{version}.tar.gz
Source1:  https://geant4-data.web.cern.ch/datasets/%{NEUTRONHPDATA}.tar.gz
Source2:  https://geant4-data.web.cern.ch/datasets/%{LEDATA}.tar.gz
Source3:  https://geant4-data.web.cern.ch/datasets/%{LEVELGAMMADATA}.tar.gz
Source4:  https://geant4-data.web.cern.ch/datasets/%{RADIOACTIVEDATA}.tar.gz
Source5:  https://geant4-data.web.cern.ch/datasets/%{PARTICLEXSDATA}.tar.gz
Source6:  https://geant4-data.web.cern.ch/datasets/%{PIIDATA}.tar.gz
Source7:  https://geant4-data.web.cern.ch/datasets/%{REALSURFACEDATA}.tar.gz
Source8:  https://geant4-data.web.cern.ch/datasets/%{SAIDXSDATA}.tar.gz
Source9:  https://geant4-data.web.cern.ch/datasets/%{ABLADATA}.tar.gz
Source10: https://geant4-data.web.cern.ch/datasets/%{INCLDATA}.tar.gz
Source11: https://geant4-data.web.cern.ch/datasets/%{ENSDFSTATEDATA}.tar.gz
Source12: https://geant4-data.web.cern.ch/datasets/%{TENDLDATA}.tar.gz

%undefine __cmake_in_source_build
%undefine __cmake3_in_source_build

%bcond_with vtk
%bcond_without examples
%bcond_without qt5
%bcond_without threads
%bcond_without trajectories

%if %{?rhel}%{!?rhel:9} >= 9
%bcond_without tbb
%bcond hdf5 %{without threads}
%bcond inventor %[%{?fedora:1}%{!?fedora:0}]
%endif

%if %{?rhel}%{!?rhel:0} == 7
BuildRequires:  cmake3 >= 3.16
BuildRequires:  devtoolset-7-toolchain
%else
BuildRequires:  cmake >= 3.16
BuildRequires:  gcc-c++
%endif

BuildRequires: make
BuildRequires: expat-devel
BuildRequires: freetype-devel
BuildRequires: motif-devel
BuildRequires: xerces-c-devel
BuildRequires: zlib-devel

%if %{with hdf5}
BuildRequires: hdf5-devel
%endif

%if %{with inventor}
BuildRequires: SoQt-devel
%endif

%if %{with qt5}
BuildRequires: qt5-qt3d-devel
BuildRequires: qt5-qtbase-devel
%endif

%if %{with tbb}
BuildRequires: tbb-devel
%endif

%if %{with vtk}
BuildRequires: vtk-devel
BuildRequires: java-latest-openjdk-devel
%endif

Requires: %{name}-data = %{version}-%{release}
Requires: %{name}-libs%{?_isa} = %{version}-%{release}

%description
Geant simulates the passage of subatomic particles through matter,
for instance, particle detectors. Geant 3 simulations are performed
by linking Fortran code supplied by the user with the Geant libraries,
then running the resulting executable. This package includes gxint,
the script used to perform this linking step. Geant 4 is a complete
rewrite in C++ with addition of other modern features and detectors.

%package libs
Summary: Geant4 libraries

%description libs
This package contains Geant4 libraries used by simulation applications.

%package devel
Summary: Development files for Geant4 (CMake modules and header files)
Provides:  %{name}-libs-devel = %{version}-%{release}
Provides:  %{name}-libs-devel%{?_isa} = %{version}-%{release}
Obsoletes: %{name}-libs-devel < %{version}-%{release}
Requires:  %{name}-libs%{?_isa} = %{version}-%{release}

%if %{?rhel}%{!?rhel:0} == 7
Requires:  cmake3 >= 3.16
Requires:  devtoolset-7-toolchain
%else
Requires:  cmake >= 3.16
Requires:  gcc-c++
%endif

Requires: make
Requires: expat-devel
Requires: freetype-devel
Requires: motif-devel
Requires: xerces-c-devel
Requires: zlib-devel

%if %{with hdf5}
Requires: hdf5-devel
%endif

%if %{with inventor}
Requires: SoQt-devel
%endif

%if %{with qt5}
Requires: qt5-qt3d-devel
Requires: qt5-qtbase-devel
%endif

%if %{with tbb}
Requires: tbb-devel
%endif

%if %{with vtk}
Requires: vtk-devel
Requires: java-latest-openjdk-devel
%endif

%description devel
Geant4 development components such as CMake modules and header files.

%package data
BuildArch: noarch
Summary: Geant4 data files required for physics models

%description data
Geant4 data files required for physics models

%if %{with examples}
%package examples
BuildArch: noarch
Summary: Geant4 user examples

%description examples
Geant4 user examples
%endif

%prep
%setup -q -n %{name}-v%{version}

%build
%if %{?rhel}%{!?rhel:0} == 7
. /opt/rh/devtoolset-7/enable
%endif

%cmake3 \
  -DCMAKE_INSTALL_PREFIX:PATH=%{_prefix} \
  -DCMAKE_INSTALL_DATADIR:PATH=%{_datadir}/%{name} \
  -DGEANT4_BUILD_BUILTIN_BACKTRACE:BOOL=OFF \
  -DGEANT4_BUILD_MULTITHREADED:BOOL=%{with threads} \
  -DGEANT4_BUILD_STORE_TRAJECTORY:BOOL=%{with trajectories} \
  -DGEANT4_BUILD_TLS_MODEL:STRING=global-dynamic \
  -DGEANT4_BUILD_VERBOSE_CODE:BOOL=OFF \
  -DGEANT4_ENABLE_TESTING:BOOL=OFF \
  -DGEANT4_INSTALL_DATA:BOOL=OFF \
  -DGEANT4_INSTALL_EXAMPLES:BOOL=%{with examples} \
  -DGEANT4_INSTALL_PACKAGE_CACHE:BOOL=OFF \
  -DGEANT4_USE_FREETYPE:BOOL=ON \
  -DGEANT4_USE_G3TOG4:BOOL=ON \
  -DGEANT4_USE_GDML:BOOL=ON \
  -DGEANT4_USE_HDF5:BOOL=%{with hdf5} \
  -DGEANT4_USE_INVENTOR:BOOL=OFF \
  -DGEANT4_USE_INVENTOR_QT:BOOL=%{with inventor} \
  -DGEANT4_USE_OPENGL_X11:BOOL=ON \
  -DGEANT4_USE_QT:BOOL=%{with qt5} \
  -DGEANT4_USE_RAYTRACER_X11:BOOL=ON \
  -DGEANT4_USE_SMARTSTACK:BOOL=OFF \
  -DGEANT4_USE_SYSTEM_CLHEP:BOOL=OFF \
  -DGEANT4_USE_SYSTEM_EXPAT:BOOL=ON \
  -DGEANT4_USE_SYSTEM_PTL:BOOL=OFF \
  -DGEANT4_USE_SYSTEM_ZLIB:BOOL=ON \
  -DGEANT4_USE_TBB:BOOL=%{with tbb} \
  -DGEANT4_USE_TIMEMORY:BOOL=OFF \
  -DGEANT4_USE_VTK:BOOL=%{with vtk} \
  -DGEANT4_USE_XM:BOOL=ON

%cmake3_build

%install
%if %{?rhel}%{!?rhel:0} == 7
. /opt/rh/devtoolset-7/enable
%endif
%cmake3_install
rm -f %{buildroot}/%{_bindir}/geant4.sh
rm -f %{buildroot}/%{_bindir}/geant4.csh

%if %{?rhel}%{!?rhel:0} == 7
find %{buildroot}/%{_datadir}/%{name}/examples -name '*.py' -delete
%endif

mkdir -p %{buildroot}/%{_datadir}/%{name}/data
tar xzf %{SOURCE1} -C %{buildroot}/%{_datadir}/%{name}/data
tar xzf %{SOURCE2} -C %{buildroot}/%{_datadir}/%{name}/data
tar xzf %{SOURCE3} -C %{buildroot}/%{_datadir}/%{name}/data
tar xzf %{SOURCE4} -C %{buildroot}/%{_datadir}/%{name}/data
tar xzf %{SOURCE5} -C %{buildroot}/%{_datadir}/%{name}/data
tar xzf %{SOURCE6} -C %{buildroot}/%{_datadir}/%{name}/data
tar xzf %{SOURCE7} -C %{buildroot}/%{_datadir}/%{name}/data
tar xzf %{SOURCE8} -C %{buildroot}/%{_datadir}/%{name}/data
tar xzf %{SOURCE9} -C %{buildroot}/%{_datadir}/%{name}/data
tar xzf %{SOURCE10} -C %{buildroot}/%{_datadir}/%{name}/data
tar xzf %{SOURCE11} -C %{buildroot}/%{_datadir}/%{name}/data
tar xzf %{SOURCE12} -C %{buildroot}/%{_datadir}/%{name}/data

%files
# Empty

%files data
%{_datadir}/%{name}/data
%{_datadir}/%{name}/fonts

%files devel
%{_bindir}/geant4-config
%{_datadir}/%{name}/geant4make
%{_includedir}/Geant4/*
%{_libdir}/cmake/Geant4/*
%{_libdir}/pkgconfig/*

%files libs
%license LICENSE
%{_datadir}/%{name}/tools.license
%{_libdir}/Geant4*
%{_libdir}/lib*.so*

%if %{with examples}
%files examples
%{_datadir}/%{name}/examples
%endif

%changelog
* Wed Dec 12 2023 Guilherme Amadio <amadio@cern.ch> - 11.2.0
- Update to version 11.2.0

* Wed Jun 28 2023 Guilherme Amadio <amadio@cern.ch> - 11.1.2
- Update to version 11.1.2

* Wed May 18 2022 Guilherme Amadio <amadio@cern.ch> - 11.0.1
- Initial version
