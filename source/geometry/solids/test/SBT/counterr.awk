BEGIN { nn = 0 }
{ if (substr($1,1,1) != "%") nn += 1 }
END { print "Total errors printed: ", nn; exit(nn) }
