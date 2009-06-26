#/bin/csh
#================================================
#     Macro file for hadr00 run over all Physics Lists
#     26.06.2009 V.Ivanchneko
#================================================

rm -f $1.out
$G4BIN/$G4SYSTEM/hadr00 $2.in $1 >& $1.out 

#
