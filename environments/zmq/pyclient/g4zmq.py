#/usr/bin/env python
# ==================================================================
# python client for zmq-geant4
#
# *** Python3 ***
# ==================================================================
import zmq

# global vars
charset = 'utf-8'
context = zmq.Context(1)
socket = context.socket(zmq.REQ)

# ===============================================================
def connect(endpoint="tcp://127.0.0.1:5555") :
  socket.setsockopt(zmq.LINGER, 0)
  socket.connect(endpoint)
  ping()

def ping() :
  socket.send(b"@@ping")
  poller = zmq.Poller()
  poller.register(socket, zmq.POLLIN)
  if poller.poll(1*1000) :
    output = socket.recv()
    print ("@@ G4ZMQ server connected.")
  else :
    raise ConnectionError("*** connection timeout")

def debug(dflag) :
  if dflag :
    socket.send(b"@@debug")
  else :
    socket.send(b"@@nodebug")
  output = socket.recv()
  print(output.decode(charset))

def apply(command) :
  cmd_str= command.encode(charset)
  socket.send(cmd_str)
  output = socket.recv()
  print(output.decode(charset))

def execute(command) :
  apply(command)

def pwd() :
  socket.send(b"pwd")
  output = socket.recv()
  print(output.decode(charset))

def cwd() :
  socket.send(b"cwd")
  output = socket.recv()
  print(output.decode(charset))

def cd(dir="") :
  cmd = "cd " + dir
  socket.send(cmd.encode(charset))
  output = socket.recv()
  print(output.decode(charset))

def lc(target="") :
  cmd = "lc " + target
  socket.send(cmd.encode(charset))
  output = socket.recv()
  print(output.decode(charset))

def ls(target="") :
  cmd = "ls " + target
  socket.send(cmd.encode(charset))
  output = socket.recv()
  print(output.decode(charset))

def help(target="") :
  if ( target == "" ) :
    raise SyntaxWarning("*** no command specified.")
  else :
    cmd = "help " + target
    socket.send(cmd.encode(charset))
    output = socket.recv()
    print(output.decode(charset))

def history() :
  socket.send(b"history")
  output = socket.recv()
  print(output.decode(charset))

def beamOn(nevent) :
  cmd = "/run/beamOn " + str(nevent)
  apply(cmd)

def getvalue(target="") :
  cmd = "?" + target
  socket.send(cmd.encode(charset))
  output = socket.recv()
  print(output.decode(charset))

def exit() :
  socket.send(b"exit")
  output = socket.recv()
  print(output.decode(charset))

# ===============================================================
if __name__ == '__main__' :
  pass
