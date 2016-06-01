// GAG (Geant4 Adaptive GUI)
// Requires : JDK 1.1.3 (or later) + swing-1.0.2
// Version : alpha  1998/Mar/31
// GEANT4 GUI Group, Toshiaki Kodama
// 1998 July 6 GEANT4 Beta-01

package GAG;

import java.io.*;     

class GAGpipe {
  private Process ps;
  private BufferedReader inpS, errS;
  private PrintWriter outS;
  private String errMsg;
  private boolean isErr;

  GAGpipe(String[] execFile){
    errMsg = null; isErr = false;
    try {
      ps = Runtime.getRuntime().exec(execFile);
    } catch(IOException e){
      isErr = true;
      errMsg = e.getMessage();
      return;
    }
    inpS = new BufferedReader(new InputStreamReader(ps.getInputStream()));
    errS = new BufferedReader(new InputStreamReader(ps.getErrorStream()));
    outS = new PrintWriter(new DataOutputStream(ps.getOutputStream()));
  }
  String readStdLine(){
    String line = null;
    try{
      line=inpS.readLine();
    }catch(IOException e){
      isErr = true;
      errMsg = e.getMessage();
    }
    return line;
  }
  String readErrLine(){
    String line = null;
    try{
      line=errS.readLine();
    }catch(IOException e){
      isErr = true;
      errMsg = e.getMessage();
      System.err.println(errMsg);
    }
    return line;
  }
  synchronized boolean writeLine(String st){
    outS.println(st);
    return outS.checkError();
  }
  synchronized void processKill(){
    ps.destroy();
  }
  boolean isError(){
    return isErr;
  }
  String getErrorMsg(){
    if (!isErr) return "";
    return errMsg;
  }
}
