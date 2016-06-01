
//example
//8.5

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

public class CreateSlicingDialog extends JDialog implements ActionListener{
   private JButton b1, b2, b3;
   private JLabel label1, label2;
   private SlicingAxialFrame slicingFrame;
   private InsertSlicingDialog insertSlicingDialog;
   public CreateSlicingDialog(SlicingAxialFrame slicingFrame){
      super(slicingFrame,"Create Slicing Axial Smmetric",false);
      this.slicingFrame = slicingFrame;
      JPanel createSlicingPanel = new JPanel();
       createSlicingPanel.setLayout(new BorderLayout());
       b1 = new JButton("APPEND");
       b2 = new JButton("INSERT");
       b3 = new JButton("END");
       b1.addActionListener(this);
       b2.addActionListener(this);
       b3.addActionListener(this);
       JPanel labelPanel = new JPanel();
        labelPanel.setLayout(new BorderLayout());
        label1 = new JLabel("  Choose Append or Insert");
        label1.setFont(new Font("Serif",Font.BOLD, 12));
//        label2 = new JLabel(" ");
//        label2.setFont(new Font("Serif",Font.BOLD, 12));
        labelPanel.add(label1, BorderLayout.NORTH);
//        labelPanel.add(label2, BorderLayout.CENTER);
       createSlicingPanel.add(labelPanel, BorderLayout.NORTH);
       createSlicingPanel.add(b3, BorderLayout.EAST);
       createSlicingPanel.add(b2, BorderLayout.CENTER);
       createSlicingPanel.add(b1, BorderLayout.WEST);
      getContentPane().add(createSlicingPanel, BorderLayout.CENTER);
//      pack();
      setSize(250, 100);

   }
   public void actionPerformed(ActionEvent e){
      if(e.getSource() == b1){
        slicingFrame.appendSlicing(); 
      }else if(e.getSource() == b2){
        insertSlicingDialog = new InsertSlicingDialog(slicingFrame);
        insertSlicingDialog.setVisible(true);
      }else if(e.getSource() == b3){
        this.setVisible(false);      
        slicingFrame.createSlicingCloseAct();
      } 
   }    
 
}





