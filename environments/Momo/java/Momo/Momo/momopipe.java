//package MOMO;

import java.io.*;     

class momopipe {
  private Process ps;
  private BufferedReader inpS, errS;
  private PrintWriter outS;
  private String errMsg;
  private boolean isErr;

  momopipe(String execFile){
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

