//Tetsuya Yamada
//GGE 
//example
//8.21

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

public class InsertSingleDialog extends JDialog implements ActionListener{
   private JButton b1, b2;
   private JLabel label1, label2;
   private SingleVolume singleFrame;
   public InsertSingleDialog(SingleVolume singleFrame){
      super(singleFrame,"Insert Single Positioned Volume",false);
      this.singleFrame = singleFrame;
      JPanel createSinglePanel = new JPanel();
       createSinglePanel.setLayout(new BorderLayout());

       b1 = new JButton("OK");
       b2 = new JButton("CANCEL");
       b1.addActionListener(this);
       b2.addActionListener(this);
       JPanel labelPanel = new JPanel();
        labelPanel.setLayout(new BorderLayout());
        label1 = new JLabel("   Click the pName after which ");
        label1.setFont(new Font("Serif",Font.BOLD, 12));
        label2 = new JLabel("   you are going to insert");
        label2.setFont(new Font("Serif",Font.BOLD, 13));
        labelPanel.add(label1, BorderLayout.NORTH);
        labelPanel.add(label2, BorderLayout.CENTER);
       createSinglePanel.add(labelPanel, BorderLayout.NORTH);
       createSinglePanel.add(b1, BorderLayout.CENTER);
       createSinglePanel.add(b2, BorderLayout.EAST);
      getContentPane().add(createSinglePanel, BorderLayout.CENTER);
//      pack();
      setSize(250, 100);

   }
   public void actionPerformed(ActionEvent e){
      if(e.getSource() == b1){
        singleFrame.insertSingle();
//        this.setVisible(false); 
      }else if(e.getSource() == b2){
        this.setVisible(false);
        singleFrame.insertSingleCloseAct();;       
      }
   }
 
}









