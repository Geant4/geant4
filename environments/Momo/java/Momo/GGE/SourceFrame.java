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

class SourceFrame extends JFrame implements ActionListener {
  private JMenuItem close, save, cut, copy, paste;
  private String fileName;
  public JTextArea editor;
//  protected FileDialog fileDialog;
  private JFileChooser fileChooser;
  SourceFrame(){
    super("C++ Text Editor");
    getContentPane().setLayout(new BorderLayout());
    editor = new JTextArea();
    //editor.setFont(new Font("Courier", Font.PLAIN, 12));
    JScrollPane scroller = new JScrollPane();
    //scroller.setVerticalScrollBarPolicy(JScrollPane.SCROLLBARS_ALWAYS);
    JViewport port = scroller.getViewport();
    port.add(editor);
    port.setBackingStoreEnabled(false);
    setJMenuBar( createMenubar() );
    getContentPane().add("Center", scroller);
    pack();
    setSize(700, 650);
//    show();
  }

    JMenuBar createMenubar(){
    JMenu file, edit;
    JMenuBar mb = new JMenuBar();
    file = new JMenu("File");
//    file.add(newf = new JMenuItem("New"));
//    file.add(open = new JMenuItem("Open"));
    file.add(save = new JMenuItem("Save"));
    file.add(close = new JMenuItem("Close"));
    edit = new JMenu("Edit"); 
    edit.add(cut = new JMenuItem("Cut"));
    edit.add(copy = new JMenuItem("Copy"));
    edit.add(paste = new JMenuItem("Paste"));
    mb.add(file);
    mb.add(edit);
//    newf.addActionListener(this);
    close.addActionListener(this);
    save.addActionListener(this);
    cut.addActionListener(this);
    copy.addActionListener(this);
    paste.addActionListener(this);
    return mb;
  }

  public static void main(String[] args) {
    new SourceFrame();
  }
  public void actionPerformed(ActionEvent e) {
    Object o = e.getSource();
    if (o == close){setVisible(false); return;}
    if (o == save){saveAction(); return;}
    if (o == cut){editor.cut(); return;}
    if (o == copy){editor.copy(); return;}
    if (o == paste){editor.paste(); return;}
  }
  void saveAction(){
    String saveText;
    if (fileChooser == null) {
      fileChooser = new JFileChooser(".");
    }
    fileChooser.setDialogTitle("Save C++ Source Code");
    fileChooser.setDialogType(1);
    ExampleFileFilter ccFile = new ExampleFileFilter("cc","C++ Source Code");
    fileChooser.setFileFilter(ccFile);
    fileChooser.showDialog(this, null);
    fileChooser.getSelectedFile();
    
    File f = fileChooser.getSelectedFile();
    fileName = f.getPath();
    try{
      FileOutputStream fo = new FileOutputStream(fileName);
      PrintWriter outf = new PrintWriter(new DataOutputStream(fo));
      outf.println(editor.getText());
      outf.close();
    }catch(IOException e){ System.err.println(e.getMessage()); }
  } 
}













