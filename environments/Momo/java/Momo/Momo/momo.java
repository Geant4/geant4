// 1998 Sep 20 the first start of Java Momo
// makeMydetector has  a pipe to display compilation messages 
// But it don't have good finding path to any files.

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.lang.Runtime.*;
import java.io.*;
import com.sun.java.swing.text.*;
import java.lang.System.*;


public class momo extends JFrame implements ActionListener{
//   editorFrame ef;
   private JButton gagButton, ggeButton, makeButton, newshellButton,pluginButton;
   JMenuBar mbar;
   JMenu jm;
   JMenuItem xtermitem, netitem, exititem;
   JPanel momopanel;
   private MOMOpipe pipe;
   private ErrStream errStream;

   public momo(){
//      ef = new editorFrame();

      momopanel = new JPanel();
      getContentPane().add(momopanel);
      gagButton = new JButton("GAG");
      gagButton.setFont(new Font("Serif", Font.BOLD, 14));
      ggeButton = new JButton("GGE");
      ggeButton.setFont(new Font("Serif", Font.BOLD, 14));
      makeButton = new JButton("make myDetector");
      makeButton.setFont(new Font("Serif", Font.BOLD, 14));
      newshellButton = new JButton("new shell");
      newshellButton.setFont(new Font("Serif", Font.BOLD, 14));
      pluginButton = new JButton("Plugins");
      pluginButton.setFont(new Font("Serif", Font.BOLD, 14));

      mbar = new JMenuBar();
      setJMenuBar(mbar);         
      jm = new JMenu("File");
      mbar.add(jm);
      //      jm.add(xtermitem = new JMenuItem("xterm"));
      jm.add(netitem = new JMenuItem("Netscape"));
      jm.addSeparator();
      jm.add(exititem = new JMenuItem("exit"));
     
      momopanel.setLayout(new GridLayout(1,3));
      momopanel.add(gagButton);
      momopanel.add(ggeButton);
      momopanel.add(makeButton);
      momopanel.add(newshellButton);
      momopanel.add(pluginButton);

      gagButton.setForeground(Color.black);
      gagButton.addActionListener(this);
      ggeButton.setForeground(Color.black);
      ggeButton.addActionListener(this);
      makeButton.setForeground(Color.black);
      makeButton.addActionListener(this);
      newshellButton.setForeground(Color.black);
      newshellButton.addActionListener(this);
      pluginButton.setEnabled(false);
      pluginButton.addActionListener(this);
      //      xtermitem.addActionListener(this);
      netitem.addActionListener(this);
      exititem.addActionListener(this);
  }

  public void actionPerformed(ActionEvent ae){
     if(ae.getSource() == exititem){
     System.exit(0); 
     }else if(ae.getSource() == newshellButton){
       try{
         String cmd="xterm";
         Runtime r = Runtime.getRuntime();
         Process p = r.exec(cmd);
       }catch(IOException e){
         System.out.println(e.toString());
       }
     }else if(ae.getSource() == netitem){
       try{
         String cmdNets="netscape";
         Runtime rNets = Runtime.getRuntime();
         Process pNets = rNets.exec(cmdNets);
       }catch(IOException ie){
         System.out.println(ie.toString());
       }
     }else if(ae.getSource() == ggeButton){
       try{
         System.out.println("GGE Started");
         String cmdGGE="java gge";
         Runtime rGGE = Runtime.getRuntime();
         Process pGGE = rGGE.exec(cmdGGE);
       }catch(IOException ie){
         System.out.println(ie.toString());
       }
     }else if(ae.getSource() == gagButton){
       try{
         System.out.println("GAG Started");
         String cmdGAG="java gag";
         Runtime rGAG = Runtime.getRuntime();
         Process pGAG = rGAG.exec(cmdGAG);
        }catch(IOException ie){
         System.out.println(ie.toString());
       }
     }else if(ae.getSource() == makeButton){
//       try{
         String userHome = System.getProperty("user.home");
         String cmdMAKE=" tcsh ./makeMyDetector"; // path is the problem
         System.out.println("Make Started");
         if((pipe = new MOMOpipe(cmdMAKE)) == null){
	   System.out.println("MOMOpipe error");}
	   else {
         errStream = new ErrStream(pipe);
         errStream.start();
	   }
//         ef.setVisible(true);
//         String cmdMAKE="makemyDetector";
//         String cmdMAKE="ls";
//         Runtime rMAKE = Runtime.getRuntime();
//         Process pMAKE = rMAKE.exec(cmdMAKE);
//       }catch(IOException ie){
//         System.out.println(ie.toString());
//       }
     }
  }
  public static void main(String args[]){
     momo m = new momo();
     m.setTitle("MOMO");
     m.pack();
     m.setVisible(true);
  }
}
































