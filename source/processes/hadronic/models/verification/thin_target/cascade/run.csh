source setup_aida.csh

gmake
rm $TARGET/res.out
echo "Start of run"

$G4MY/test30 $TARGET/run.mac  > /dev/null

mv par.paw   $TARGET/
mv kin.paw   $TARGET/
mv chips.paw $TARGET/
mv prec.paw  $TARGET/

rm $TARGET/*.old
rm $TARGET/*~

echo "Is done!"

