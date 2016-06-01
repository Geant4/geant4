import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;
import com.sun.java.swing.preview.filechooser.*;
import ExampleFileFilter;

import com.sun.java.swing.text.*;
import com.sun.java.swing.preview.*;
import com.sun.java.swing.event.*;
import com.sun.java.swing.*;

class DAWNSourceFrame extends JFrame implements ActionListener {
  private JMenuItem close, save, cut, copy, paste;
  private String fileName;
  public JTextArea dawnEditor;
//  protected FileDialog fileDialog;
  private JFileChooser fileChooser;
  DAWNSourceFrame(){
    super(" DAWN Text Editor");
    getContentPane().setLayout(new BorderLayout());
    dawnEditor = new JTextArea();
    //editor.setFont(new Font("Courier", Font.PLAIN, 12));
    JScrollPane scroller = new JScrollPane();
    //scroller.setVerticalScrollBarPolicy(JScrollPane.SCROLLBARS_ALWAYS);
    JViewport port = scroller.getViewport();
    port.add(dawnEditor);
    port.setBackingStoreEnabled(false);
    setJMenuBar( createMenubar() );
    getContentPane().add("Center", scroller);
    pack();
    setSize(350, 350);
//    show();
  }

    JMenuBar createMenubar(){
    JMenu file, edit;
    JMenuBar mb = new JMenuBar();
    file = new JMenu("File");
    file.add(close = new JMenuItem("Close"));
    file.add(save = new JMenuItem("Save"));

    edit = new JMenu("Edit"); 
    edit.add(cut = new JMenuItem("Cut"));
    edit.add(copy = new JMenuItem("Copy"));
    edit.add(paste = new JMenuItem("Paste"));
    mb.add(file);
    mb.add(edit);

    close.addActionListener(this);
    save.addActionListener(this);
    cut.addActionListener(this);
    copy.addActionListener(this);
    paste.addActionListener(this);
    return mb;
  }

  public static void main(String[] args) {
    new DAWNSourceFrame();
  }
  public void actionPerformed(ActionEvent e) {
    Object o = e.getSource();
    if (o == save){saveAction(); return;}
    if (o == close){setVisible(false);}
    if (o == cut){dawnEditor.cut(); return;}
    if (o == copy){dawnEditor.copy(); return;}
    if (o == paste){dawnEditor.paste(); return;}
  }
  void saveAction(){
    String saveText;
    if (fileChooser == null) {
      fileChooser = new JFileChooser(".");
    }
    fileChooser.setDialogTitle("Save C++ Source Code");
    fileChooser.setDialogType(1);
    ExampleFileFilter primFile = new ExampleFileFilter("prim","PRIM Source Code");
    fileChooser.setFileFilter(primFile);
    fileChooser.showDialog(this, null);
    fileChooser.getSelectedFile();
    
    File f = fileChooser.getSelectedFile();
    fileName = f.getPath();
    try{
      FileOutputStream fo = new FileOutputStream(fileName);
      PrintWriter outf = new PrintWriter(new DataOutputStream(fo));
      outf.println(dawnEditor.getText());
      outf.close();
    }catch(IOException e){ System.err.println(e.getMessage()); }
  } 
}













