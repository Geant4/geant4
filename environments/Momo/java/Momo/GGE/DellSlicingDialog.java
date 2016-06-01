
//example
//8.5

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

public class DellSlicingDialog extends JDialog implements ActionListener{
   private JButton b1, b2;
   private JLabel label1;
   private SlicingAxialFrame slicingFrame;
   public DellSlicingDialog(SlicingAxialFrame slicingFrame){
      super(slicingFrame,"Delete Slicing Axial",false);
      this.slicingFrame = slicingFrame;
      JPanel createSlicingPanel = new JPanel();
       createSlicingPanel.setLayout(new BorderLayout());
       b1 = new JButton("OK");
       b2 = new JButton("CANCEL");
       b1.addActionListener(this);
       b2.addActionListener(this);

        label1 = new JLabel("   Click the pName to delete");
        label1.setFont(new Font("Serif",Font.BOLD, 13));

       createSlicingPanel.add(label1, BorderLayout.NORTH);
       createSlicingPanel.add(b1, BorderLayout.CENTER);
       createSlicingPanel.add(b2, BorderLayout.EAST);
      getContentPane().add(createSlicingPanel, BorderLayout.CENTER);
//      pack();
      setSize(250, 80);

   }
   public void actionPerformed(ActionEvent e){
      if(e.getSource() == b1){
        slicingFrame.dellSlicing();
//        this.setVisible(false); 
      }else if(e.getSource() == b2){
        this.setVisible(false);       
        slicingFrame.dellSlicingCloseAct();
      }
   }
 
}






