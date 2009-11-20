#/bin/csh

foreach tPart (p he4 c12)
    ${G4BIN}/${G4SYSTEM}/reader_test44 $tPart $1
end
