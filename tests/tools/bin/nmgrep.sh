#
# Use in conjunction with find to find an entry point in libs :
# Example :
#  UNIX> find $G4LIB/$G4SYSTEM -name "*.a" -exec ./nmgrep.sh {} name \;
#
nm $1 | grep $2 && echo $1
