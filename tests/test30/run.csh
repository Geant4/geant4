source setup_aida.csh

cd $TARGET
setenv G4TARGET data
rm ${G4MY}/${G4TARGET}
gmake
${G4MY}/${G4TARGET} exfor.dat
echo "Data are built"

if ( -r exfor_1.dat ) then
  setenv G4TARGET data_1
  rm ${G4MY}/${G4TARGET}
  gmake
  ${G4MY}/${G4TARGET} exfor_1.dat
  echo "Data_1 are built"
endif


cd ../

gmake
rm $TARGET/res.out
echo "Start of run"

$G4MY/test30 $TARGET/run.mac  

mv par.paw   $TARGET/
mv kin.paw   $TARGET/
mv chips.paw $TARGET/
mv prec.paw  $TARGET/

rm $TARGET/*.old
rm $TARGET/*~

echo "Is done!"

