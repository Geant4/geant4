BEGIN { nn = 0 }
{ if (substr($1,1,5) == "VOXEL") nn += 1 }
END { print "Total voxel errors: ", nn; exit(nn) }
