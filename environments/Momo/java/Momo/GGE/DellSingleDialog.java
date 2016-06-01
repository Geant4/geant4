//Tetsuya Yamada
//example
//8.21

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

public class DellSingleDialog extends JDialog implements ActionListener{
   private JButton b1, b2;
   private JLabel label1;
   private SingleVolume singleFrame;
   public DellSingleDialog(SingleVolume singleFrame){
      super(singleFrame,"Delete Single Positioned Volume",false);
      this.singleFrame = singleFrame;
      JPanel createSinglePanel = new JPanel();
       createSinglePanel.setLayout(new BorderLayout());
       b1 = new JButton("OK");
       b2 = new JButton("CANCEL");
       b1.addActionListener(this);
       b2.addActionListener(this);

        label1 = new JLabel("   Click the pName to delete");
        label1.setFont(new Font("Serif",Font.BOLD, 13));

       createSinglePanel.add(label1, BorderLayout.NORTH);
       createSinglePanel.add(b1, BorderLayout.CENTER);
       createSinglePanel.add(b2, BorderLayout.EAST);
      getContentPane().add(createSinglePanel, BorderLayout.CENTER);
//      pack();
      setSize(250, 80);

   }
   public void actionPerformed(ActionEvent e){
      if(e.getSource() == b1){
        singleFrame.dellSingle();
//        this.setVisible(false); 
      }else if(e.getSource() == b2){
        this.setVisible(false);       
        singleFrame.dellSingleCloseAct();
      }
   }
 
}










