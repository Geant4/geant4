#!/bin/sh

ls g4ge*.txt | grep -v cutsPerRegion | grep -v _SD |
while read filename ; do
  type1="`echo "$filename" | cut -d "_" -f2`"
  type="`echo "$type1" | cut -d "." -f1`"
  echo "File: $filename $type"
  cp $filename g4geom.txt
 $G4INSTALL/bin/Linux-g++/textGeom run.mac | tee zz_$type.out 
 mv g4_00.wrl g4_$type.wrl
sleep 5
done

sed -i s/"new ExTGDetectorConstruction"/"new ExTGDetectorConstructionWithSD"/g textGeom.cc
make;  $G4INSTALL/bin/Linux-g++/textGeom run.mac | tee zz_WithSD.out 
sed -i s/"new ExTGDetectorConstructionWithSD"/"new ExTGDetectorConstructionWithCpp"/g textGeom.cc
make;  $G4INSTALL/bin/Linux-g++/textGeom run.mac | tee zz_WithCpp.out 
sed -i s/"new ExTGDetectorConstructionWithCpp"/"new ExTGDetectorConstructionWithCuts"/g textGeom.cc
make;  $G4INSTALL/bin/Linux-g++/textGeom run.mac | tee zz_WithCuts.out 
sed -i s/"new ExTGDetectorConstructionWithCuts"/"new ExTGDetectorConstruction"/g textGeom.cc
make

exit
