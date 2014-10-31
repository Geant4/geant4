#!/usr/bin/env bash
#
# usage:
# source generate_upload_script_pbar_regre.sh geant4-10.00-ref08 internal
#

g4version=${1}
status=${2}

# g4version=geant4-10.00-ref07
# status=internal

testname=test48

script=${testname}-upload-pbar-regre.xml

printf "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" >> ${script}
printf "\n" >> ${script}
printf "<testslist>\n" >> ${script}
printf "\n" >> ${script}

printf "<tests>\n" >> ${script}

printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}

printf "<observable>Charged Pion Momentun; Multiplicity of Secondary Pions</observable>\n" >> ${script}
printf "<secondary>pi+,pi-,pi0</secondary>\n" >> ${script}

printf "<reaction>annihilation pbar on H</reaction>\n" >> ${script}

printf "<beam>antiproton</beam>\n" >> ${script}
printf "<beammomentum>0.0 MeV/c</beammomentum>\n" >> ${script}
printf "<beamenergy></beamenergy>\n" >> ${script}
printf "<target>Hydrogen</target>\n" >> ${script}
printf "<status>${status}</status>\n" >> ${script}
printf "<score>\n" >> ${script}
printf "<type>expert</type>\n" >> ${script}
printf "<value>passed</value>\n" >> ${script}
printf "</score>\n" >> ${script}

printf "<imageblob>pbar-H-FTF-regre.gif</imageblob>\n" >> ${script}

printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>FTF</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Regression testing</name>\n" >> ${script}
printf "<value>g4-9.6-p03,g4-10.00-p02,g4-10.00-ref07,${g4version}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}
printf "</tests>\n" >> ${script}
printf "\n" >> ${script}

printf "</testslist>\n" >> ${script}
