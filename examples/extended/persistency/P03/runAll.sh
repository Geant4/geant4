#!/bin/sh

ls sg4ge*.txt | grep -v cutsPerRegion | grep -v _SD | grep -v parallel | 
while read filename ; do
  type1="`echo "$filename" | cut -d "_" -f2`"
  type="`echo "$type1" | cut -d "." -f1`"
  echo "File: $filename $type"
  cp $filename g4geom.txt
 $G4INSTALL/bin/Linux-g++/exTextGeom run.mac | tee zz_$type.out 
 mv g4_00.wrl g4_$type.wrl
sleep 5
done

sed -i s/"new ExTGDetectorConstruction"/"new ExTGDetectorConstructionWithSD"/g exTextGeom.cc
make;  $G4INSTALL/bin/Linux-g++/exTextGeom run.mac | tee zz_WithSD.out 
sed -i s/"new ExTGDetectorConstructionWithSD"/"new ExTGDetectorConstructionWithCpp"/g exTextGeom.cc
make;  $G4INSTALL/bin/Linux-g++/exTextGeom run.mac | tee zz_WithCpp.out 
sed -i s/"new ExTGDetectorConstructionWithCpp"/"new ExTGDetectorConstructionWithCuts"/g exTextGeom.cc
make;  $G4INSTALL/bin/Linux-g++/exTextGeom run.mac | tee zz_WithCuts.out 
sed -i s/"new ExTGDetectorConstructionWithCuts"/"new ExTGDetectorConstructionWithParallel"/g exTextGeom.cc
sed -i s/"new ExTGPhysicsList"/"new ExTGPhysicsListWithParallel"/g exTextGeom.cc
make;  $G4INSTALL/bin/Linux-g++/exTextGeom run.mac | tee zz_WithParallel.out 

sed -i s/"new ExTGDetectorConstructionWithParallel"/"new ExTGDetectorConstruction"/g exTextGeom.cc
sed -i s/"new ExTGPhysicsListWithParallel"/"new ExTGPhysicsList"/g exTextGeom.cc
make

exit
