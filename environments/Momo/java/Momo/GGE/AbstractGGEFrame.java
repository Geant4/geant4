import java.awt.event.*;
import com.sun.java.swing.*;
import java.io.Serializable;

abstract class AbstractGGEFrame extends JFrame implements WindowListener, Serializable {
  AbstractGGEFrame(String title){
    super(title);
    addWindowListener(this);
  }
  public void windowDeiconified(WindowEvent event) {}
  public void windowIconified(WindowEvent event) {}
  public void windowActivated(WindowEvent event) {}
  public void windowDeactivated(WindowEvent event) {}
  public void windowOpened(WindowEvent event) {}
  public void windowClosed(WindowEvent event) {}
  public void windowClosing(WindowEvent event) {
    System.exit(0);
  }
}
