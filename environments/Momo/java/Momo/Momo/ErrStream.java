import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.lang.Runtime.*;
import java.io.*;
import com.sun.java.swing.text.*;

class ErrStream extends Thread implements ActionListener{
  private MOMOpipe pipe;
  private Frame errFrame;
  private TextArea errArea;
  private Button close;

  ErrStream(MOMOpipe pipe){    
    this.pipe = pipe;
    errFrame = new Frame("Momo TextEditor(standerd output)");
    errFrame.setLayout(new BorderLayout());
    errFrame.add("Center", errArea = new TextArea());
    errFrame.add("South", close = new Button("close"));
    close.addActionListener(this);
    errFrame.setSize(500,200);
  }
  public void run(){
    String sline, eline;
    sline = pipe.readStdLine();
    errArea.append("Welcome to GGE\n\n");
    errArea.append("Making myDetector...\n\n");
    while((sline=pipe.readStdLine())!=null){
      errArea.append(sline+"\n");
      if (!errFrame.isVisible()){ 
        errFrame.setVisible(true);
      }
    }
    eline = pipe.readStdLine();
    errArea.append("\n Compilation ended\n\n");
    errArea.append("\n std error outputs if any\n\n");
    while((eline=pipe.readErrLine())!=null){
      errArea.append(eline+"\n");
      if (!errFrame.isVisible()){
        errFrame.setVisible(true);
      }
    }
  }
  public void actionPerformed(ActionEvent ae){
    errFrame.setVisible(false);
    errArea.setText("");
  }
  void close(){
    errFrame.dispose();
    stop();
  }
}


















