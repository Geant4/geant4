/*
 *
 */

import java.awt.*;
import java.awt.event.*;

import com.sun.java.swing.*;
import com.sun.java.swing.event.*;

abstract class SolidDialog extends JDialog {
  SolidDialog(Frame parent, String title){
    super(parent, title);
  }
  abstract void editStop();
  abstract SolidItem getValues();
}
