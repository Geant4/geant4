// jdk-1.1.3 (or later) + swing-1.0.2
// This is not required in GAG

import com.sun.java.swing.*;

import java.awt.*;
import java.net.URL;
import java.net.MalformedURLException;
import java.io.IOException;

public class HtmlFrame extends JFrame {
  public HtmlFrame() {
    // setBackground(Color.white);
    setLayout(new BorderLayout());
    try {
      URL url = new URL("http://erlt1.naruto-u.ac.jp/~geant4/README_Momo.html");
      JEditorPane html = new JEditorPane(url);
      JScrollPane scroller = new JScrollPane(html);
      add("Center", scroller);
    } catch (MalformedURLException e) {
      System.out.println("Malformed URL: " + e);
    } catch (IOException e) {
      System.out.println("IOException: " + e);
    }
    setSize(400, 400);
    show();
  }
  public static void main(String arg[]){
    new HtmlFrame();
  }
}
