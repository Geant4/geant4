
$G4MY/test30 $TARGET/run.mac  >& /dev/null

mv par.paw   $TARGET/
mv kin.paw   $TARGET/
mv chips.paw $TARGET/
mv prec.paw  $TARGET/

echo $TARGET " is done!"

