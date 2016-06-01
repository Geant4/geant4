/*
 * Copyright (c) 1997 John Jensen. All rights reserved.
 *
 * This software is FREE FOR COMMERCIAL AND NON-COMMERCIAL USE,
 * provided the following condition is met.
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for any purpose and without fee is hereby granted,
 * provided that any copy or derivative of this software or documentation
 * retaining the name "John Jensen" also retains this condition and the
 * following disclaimer.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * CopyrightVersion 1.0
 */

package TextEdit;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.awt.event.*;

class TextCanvas extends Canvas implements MouseMotionListener, MouseListener, ComponentListener, AdjustmentListener, KeyListener, FocusListener {
  final public int NOACTION = 0;			// a null note action
  final private int EDGE = 10;

  private Scrollbar	horiz;				// control
  private Scrollbar	vert;				// control

  private boolean		eactive,hactive;	// is there currently a range?
  private boolean		dirty;				// document dirty flag
  private boolean		mouseDown;
	
  private int			sx,sy,ny;			// current scroll positions
  private int			line,column,pix;	// the current column
  private int			eline,ecolumn,epix;	// the current end column
  private int			widest,wpix;		// the longest line
  private int			oldlines;			// repaint to old line count (if longer)
  private int			fontHeight;			// save time
  private int			fontDescent;		// save time
  private int			tabSize = 4;		// preference

  private String		lineSeparator;
	
  private Font		font;				// preference
  private FontMetrics fontMetrics;		// save time
  private Frame		parent;				// freindly object
  private TextCursor	textCursor;			// freindly object
  private EditMan		editMan;			// freindly object
  private Dimension	dimension;			// save time
  private Vector		lines;				// line buffer
  private TextJournal	journal;			// journal for "undo"
  private TextScroller textScroller = null;

  public TextCanvas(Frame p, Scrollbar h, Scrollbar v) {
    super();

    setBackground(Color.white);
    setForeground(Color.black);

    parent = p;

    lineSeparator = System.getProperty("line.separator");

    horiz = h;
    horiz.addAdjustmentListener(this);
    horiz.addKeyListener(this);
		
    vert = v;
    vert.addAdjustmentListener(this);
    vert.addKeyListener(this);

    lines = new Vector();

    journal = new TextJournal(this,lines,lineSeparator);

    addComponentListener(this);
    addMouseListener(this);
    addKeyListener(this);
    addFocusListener(this);

    clear();
    validate();
    requestFocus();

    textCursor = new TextCursor(this);
    textCursor.start();
  }

  //Toshiaki added
  public void writeString(String str){
    paste(str);
    journal.clear();
  }
  public void writeLine(String str){
    pause_cursor();
    textCursor.suspend();
    if (hactive) copy(true);
    String s = (String)lines.elementAt(line);
    lines.setElementAt(s + str, line);
    line++; column = 0;
    lines.insertElementAt("", line);
    if (editMan != null) {
      editMan.updateCopyItems(false);
      editMan.updateUndoItems(false, false);
    }
    journal.clear();
    if (dimension == null) dimension = getSize();
    shiftVert(line);
    shiftHoriz(EDGE);
    repaint();
    textCursor.resume();
    dirty = true;
  }

  private boolean gotFocus = false;

  public void focusGained(FocusEvent e) {
    gotFocus = true;
    release_cursor();
  }

	private void release_cursor()
	{
		if ((textCursor != null) && !eactive && gotFocus)
			textCursor.release_cursor(true);
	}

	public void focusLost(FocusEvent e)
	{
		pause_cursor();
		gotFocus = false;
	}

	private void pause_cursor()
	{
		if ((textCursor != null) && !eactive && gotFocus)
			textCursor.pause_cursor(true);
	}

	public void setEditMan(EditMan em)
	{
		editMan = em;
	}

	public void clear()
	{
		clear_area(null);
		hactive = eactive = dirty = false;
		if (editMan != null)
		{
			editMan.updateCopyItems(false);
			editMan.updateUndoItems(false, false);
		}
		oldlines = wpix = widest = line = column = sx = sy = 0;
		pix = EDGE;
		horiz.setValue(0);
		vert.setValue(0);
		lines.removeAllElements();
		journal.clear();
	}

	public void validate()
	{
		if (lines.size() == 0)
			appendLine(new String(""));
		redoCanvas();
	}

	public void setFont(String fs, int size)
	{
		font = new Font(fs, Font.PLAIN, size);
		Graphics g = getGraphics();
		
		if (g != null)
		{
			clear_area(g);
			updateFonts(g);
			resizeLines();
			redoCanvas();
			repaint();
			g.dispose();
		}
		else
			fontMetrics = null;
	}

	private void clear_area(Graphics g)
	{
		boolean foobar = false;
		
		pause_cursor();

		if (dimension == null)
			return;
			
		if (g == null)
		{
			g = getGraphics();
			foobar = true;
		}
		
		if (g != null)
			g.clearRect(0,0,dimension.width,dimension.height);

		if (foobar && (g != null))
			g.dispose();
	}
			
	public void setTab(int tab)
	{
		clear_area(null);
		tabSize = tab;
		resizeLines();
		repaint();
	}

	public void undo(boolean undo_it)
	{
		TextPosition tp;

		pause_cursor();

		if (undo_it)
			tp = journal.undo();
		else
			tp = journal.redo();

		if (tp == null)
			release_cursor();
		else
		{
			hactive = eactive = false;
			line = eline = tp.line;
			column = ecolumn = tp.column;
			pix = epix = pix_at(line,column);
			shiftVert(line);
			shiftHoriz(pix);
			repaint();
		}
	}

	public void updateUndoItems(boolean have_undo, boolean have_redo)
	{
		if (editMan != null)
			editMan.updateUndoItems(have_undo, have_redo);
	}
	
	private void updateFonts(Graphics g)
	{
		boolean foobar = false;
		
		if (g == null)
		{
			g = getGraphics();
			foobar = true;
		}
		
		if (g != null)
		{
			g.setFont(font);
			fontMetrics = g.getFontMetrics(font);
			fontHeight = fontMetrics.getHeight();
			fontDescent = fontMetrics.getDescent();
			
			if (foobar)
				g.dispose();
		}
	}

	public void print(String name)
	{
		PrintMan pm = new PrintMan(parent, name, lines, font, tabSize);

		if (pm != null)
			pm.start();
	}

	public boolean find( String pattern )
	{
		int i,max,ct,start,from;
		String whole;

		pause_cursor();

		max = lines.size();
		i = line;
		ct = 0;

		while (ct++ <= max)
		{
			whole = (String)lines.elementAt(i);

			if (( i == line ) && ( ct == 1 ))
			{
				if (eactive)
					from = column + 1;
				else
					from = column;
			}
			else
				from = 0;

			start = whole.indexOf(pattern,from);

			if (start >= 0)
			{
				line = eline = i;
				column = start;
				pix = stringLength(whole.substring(0,column)) + EDGE;
				ecolumn = start + pattern.length();
				epix = stringLength(whole.substring(0,ecolumn)) + EDGE;
				eactive = true;
				setup_h();
				save_h();
				if (editMan != null)
					editMan.updateCopyItems(eactive);
				shiftVert(line);
				shiftHoriz(pix);
				shiftHoriz(epix);
				repaint();
				return true;
			}

			i++;
			if (i >= max)
				i = 0;
		}

		release_cursor();
		return false;
	}

	public synchronized void mousePressed(MouseEvent e)
	{
		Graphics g;

		requestFocus();

		if (mouseDown)
			return;

		mouseDown = true;
		
		pause_cursor();
		g = getGraphics();
		updateFonts(g);

		clickPosition( e.getX(), e.getY() );

		if ((cline < 0) || (cline >= lines.size()))
			return;

		eline = cline;
		ecolumn = ccolumn;
		epix = cpix;

		int clickCount = e.getClickCount();

		if (clickCount == 3)
		{
			eline = line = cline;
			column = 0;
			pix = EDGE;
			String whole = (String)lines.elementAt(cline);
			ecolumn = whole.length();
			epix = stringLength(whole) + EDGE;
			eactive = true;
			setup_h();
			save_h();
			if (editMan != null)
				editMan.updateCopyItems(eactive);
			repaint();
		}
		else
		if (clickCount == 2)
		{
			char c;
			eline = line = cline;
			String whole = (String)lines.elementAt(cline);
			column = ecolumn = ccolumn;
			while (true)
			{
				if (column == 0)
					break;

				if (!Character.isLetterOrDigit(whole.charAt(column-1)))
					break;

				column--;
			}
			pix = stringLength(whole.substring(0,column)) + EDGE;
			int max = whole.length();
			while (true)
			{
				if (ecolumn >= max)
					break;

				if (!Character.isLetterOrDigit(whole.charAt(ecolumn)))
					break;

				ecolumn++;
			}
			epix = stringLength(whole.substring(0,ecolumn)) + EDGE;
			eactive = true;
			setup_h();
			save_h();
			if (editMan != null)
				editMan.updateCopyItems(eactive);
			repaint();
		}
		else
		if (e.isShiftDown())
		{
			eline = cline;
			ecolumn = ccolumn;
			epix = cpix;
			eactive = true;
			setup_h();
			save_h();
			if (editMan != null)
				editMan.updateCopyItems(eactive);
			repaint();
		}
		else
		{
			if (hactive)
				flip_h( g, oline, opix, oeline, oepix );
			eline = line = cline;
			ecolumn = column = ccolumn;
			epix = pix = cpix;
			hactive = false;
			eactive = true;
			if (editMan != null)
				editMan.updateCopyItems(eactive);
			addMouseMotionListener(this);
		}
		
		g.dispose();
	}
	
	private int opix, oepix, oline, oeline;
	private int hpix, hepix, hline, heline, hcolumn, hecolumn;
	private int lastx,lasty;

	public synchronized void mouseDragged(MouseEvent e)
	{
		Graphics g;
		
		g = getGraphics();
		updateFonts(g);
		
		if (e != null)
		{
			lastx = e.getX();
			lasty = e.getY();
		}

		clickPosition( lastx, lasty );

		eline = cline;
		ecolumn = ccolumn;
		epix = cpix;

		eactive = ((eline != line) || (ecolumn != column));
		if (eactive || hactive)
		{
			setup_h();
			
			if ((eline < sy) || 
				(eline >= sy+ny) ||
				(epix < sx) ||
				(epix >= dimension.width))
			{
				if (textScroller == null)
				{
					textScroller = new TextScroller(this);
					textScroller.start();
				}
				shiftVert(eline);
				shiftHoriz(epix);
				setup_h();
				save_h();
				for (int i = 0; i < ny; i++)
					drawLine(g,sy+i);
			}
			else
			{
				if (textScroller != null)
				{
					textScroller.stop();
					textScroller = null;
				}
				if (hactive)
				{
					if ((hline < oline) || ((hline == oline) && (hpix < opix)))
						flip_h( g, hline, hpix, oline, opix );
					if ((hline > oline) || ((hline == oline) && (hpix > opix)))
						flip_h( g, oline, opix, hline, hpix );
					if ((heline < oeline) || ((heline == oeline) && (hepix < oepix)))
						flip_h( g, heline, hepix, oeline, oepix );
					if ((heline > oeline) || ((heline == oeline) && (hepix > oepix)))
						flip_h( g, oeline, oepix, heline, hepix );
				}
				else
					flip_h( g, hline, hpix, heline, hepix );
				
				hactive = eactive;
				save_h();
			}
		}	
		g.dispose();
	}

	private int cline, ccolumn, cpix;

	void clickPosition( int x, int y )
	{
		cline = (y / fontHeight) + sy;

		if (cline < 0)
		{
			cline = ccolumn = cpix = 0;
			return;
		}
		
		if (cline >= lines.size())
		{
			cline = lines.size() - 1;
			String whole = (String)lines.elementAt(cline);
			ccolumn = whole.length();
			cpix = stringLength(whole) + EDGE;
			return;
		}

		x += sx;

		ccolumn = 0;
		cpix = 0;

		if (x > EDGE)
		{
			String whole = (String)lines.elementAt(cline);
			String part;

			ccolumn = whole.length();
			
			if (( ccolumn > 0 ) && ( x < (cpix = stringLength(whole) + EDGE) ))
			{
				if (ccolumn > 10)
				{
					int delta = ccolumn / 5;

					do {
						ccolumn -= delta;
						if (ccolumn < 0)
						{
							ccolumn = 0;
							part = null;
						}
						else
							part = whole.substring(0,ccolumn);
					} while (( ccolumn > 0 ) && ( x < (cpix = stringLength(part) + EDGE) ));

					ccolumn += delta;
				}

				while (( ccolumn > 0 ) && ( x < (cpix = stringLength(whole) + EDGE) ))
					whole = whole.substring(0,--ccolumn);
			}
		}
	}

	private void flip_h(Graphics g, int sline, int spix, int eline, int epix)
	{
		int i;
		int bx,ex,by;
		
		g.setXORMode(Color.white);
		
		if (eline >= sy + ny)
		{
			eline = sy + ny - 1;
			epix = 5000;
		}
	
		for (i = sline; i <= eline; i++)
		{
			by = (i - sy) * fontHeight;
			bx = 0;
			ex = 5000;
			if (i == sline)
			{
				bx = spix - sx;
			}
			if (i == eline)
			{
				ex = epix - sx;
				if (sline == eline)
					ex -= bx;
			}
			g.fillRect(bx,by,ex,fontHeight);
		}
	}
		
	private void setup_h()
	{
		if ((line < eline) || ((line == eline) && (column <= ecolumn)))
		{
			hline = line;
			hcolumn = column;
			hpix = pix;
			heline = eline;
			hecolumn = ecolumn;
			hepix = epix;
		}
		else
		{
			hline = eline;
			hcolumn = ecolumn;
			hpix = epix;
			heline = line;
			hecolumn = column;
			hepix = pix;
		}
	}

	private void save_h()
	{
		opix = hpix;
		oepix = hepix;
		oline = hline;
		oeline = heline;
		hactive = eactive;
	}
		
	public synchronized void mouseReleased(MouseEvent e)
	{
		if (!mouseDown)
			return;

		mouseDown = false;
		
		if (textScroller != null)
		{
			textScroller.stop();
			textScroller = null;
		}
		removeMouseMotionListener(this);
		eactive = ((eline != line) || (ecolumn != column));
		release_cursor();
		if (editMan != null)
			editMan.updateCopyItems(eactive);
	}
	
	public void mouseEntered(MouseEvent e) {}
	public void mouseExited(MouseEvent e) {}
	public void mouseClicked(MouseEvent e) {}
	public void mouseMoved(MouseEvent e) {}
	public void keyTyped(KeyEvent e) {}
	public void keyReleased(KeyEvent e) {}
	
	private final int myMask = KeyEvent.CTRL_MASK | KeyEvent.META_MASK | KeyEvent.ALT_MASK;

    public void keyPressed(KeyEvent e)
	{
		int		keyCode = e.getKeyCode();
		char	keyChar;
		int		max;
		String	s;

		if (!gotFocus)
			return;
		
		if ((e.getModifiers() & myMask) != 0)
			return;

		if ((keyCode == e.VK_UP) && (line == 0))
			return;

		pause_cursor();

		switch (keyCode) {
		case e.VK_UP:
			line--;
			cursorAdjust();
			return;
		case e.VK_DOWN:
			line++;
			cursorAdjust();
			return;
		case e.VK_RIGHT:
			column++;
			cursorAdjust();
			return;
		case e.VK_LEFT:
			column--;
			cursorAdjust();
			return;
		default:
			keyChar = e.getKeyChar();
			if (keyChar != KeyEvent.CHAR_UNDEFINED)
			{
				if (keyCode == KeyEvent.VK_ENTER)
				{
					if (hactive)
						copy(true);
					journal.split_line(line,column);
					shiftVert(++line);
					column = 0;
					pix = EDGE;
					shiftHoriz(0);
					dirty = true;
					repaint();
					return;
				}
				else
				if (keyCode == KeyEvent.VK_BACK_SPACE)
				{
					if (hactive)
					{
						copy(true);
						return;
					}
					if (column > 0)
					{
						journal.delete_prev_char(line,column);
						column--;
						redoLine(line,column+1);
						return;
					}
					else
					if (line > 0)
					{
						line--;
						s = (String)lines.elementAt(line);
						column = s.length();
						shiftHoriz(stringLength(s) + EDGE);
						journal.join_line(line,column);
						shiftVert(line);
						dirty = true;
						repaint();
						return;
					}
				}
				else
				if (keyCode == KeyEvent.VK_DELETE)
				{
					if (hactive)
					{
						copy(true);
						return;
					}
					s = (String)lines.elementAt(line);
					max = s.length();
					if (column < max)
					{
						journal.delete_next_char(line,column);
						redoLine(line,column+1);
						return;
					}
					else
					if (line+1 < lines.size())
					{
						journal.join_line(line,column);
						dirty = true;
						repaint();
						return;
					}
				}
				else
				if (keyCode != KeyEvent.VK_ESCAPE &&
					keyChar != '\n' &&
					keyChar != '\b')
				{
					if (hactive)
						copy(true);
					journal.insert_char(line,column,keyChar);
					column++;
					redoLine(line,column);
					release_cursor();
					return;
				}
			}
		}
		release_cursor();
	}
	
	public void paste(String s)
	{
		int oldline;
		int oldcolumn;
		TextPosition tp;
		
		pause_cursor();

		if (hactive)
			copy(true);

		oldline = line;
		oldcolumn = column;

		tp = journal.insert_section(line,column,s);
		line = tp.line;
		column = tp.column;

		pix = pix_at( line, column );
		shiftVert(line);
		shiftHoriz(pix);
		
		if (oldline != line)
		{
			repaint();
		}
		else
		{
			redoLine(oldline,oldcolumn+1);
//			release_cursor();
		}

		dirty = true;
	}

	int pix_at(int line, int column)
	{
		if (column == 0)
			return EDGE;
		else
			return stringLength(((String)lines.elementAt(line)).substring(0, column)) + EDGE;
	}

	public String copy(boolean cut)
	{
		String s;
		int oldline;
		int oldcolumn;
		boolean paint = false;

		if (!hactive)
			return null;

		if (hline == heline)
		{
			oldline = hline;
			oldcolumn = hcolumn;
		}
		else
			paint = true;
		
		s = journal.delete_section(hline, hcolumn, heline, hecolumn, cut);

		if (cut)
		{
			dirty = true;
			line = hline;
			column = hcolumn;
			hactive = eactive = false;

			if (shiftVert(line))
				paint = true;

			pix = pix_at(line,column);
			if (shiftHoriz(pix))
				paint = true;

			if (paint)
			{
				repaint();
			}
			else
			{
				redoLine(line,column+1);
//				release_cursor();
			}

			if (editMan != null)
				editMan.updateCopyItems(eactive);
		}

		return s;
	}
	
	private void cursorAdjust()
	{
		int min, max;
		boolean paint = false;
		
		if (eactive)
		{
			paint = true;
			hactive = eactive = false;
			if (editMan != null)
				editMan.updateCopyItems(eactive);
		}

		// straighten up lines first

		if (line < 0)
			line = 0;
		else
		{
			max = lines.size();

			if (line >= max )
				line = max - 1;
		}

		// straighten up columns

		if (column < 0)
			column = 0;
		else
		{
			max = ((String)lines.elementAt(line)).length();
			if (column > max)
				column = max;
		}

		// off page?

		if (shiftVert(line))
			paint = true;

		if (shiftHoriz(stringLength(((String)lines.elementAt(line)).substring(0, column)) + EDGE))
			paint = true;

		if (paint)
			repaint();
		else
			release_cursor();
	}

	private boolean shiftVert( int line )
	{
		if (line < sy)
		{
			sy = line;
			if (sy < 0)
				sy = 0;
			vert.setValue(sy);
			return true;
		}

		if (line >= sy+ny)
		{
			sy = line - ny + 1;
			vert.setValue(sy);
			return true;
		}

		return false;
	}

	private boolean shiftHoriz( int x )
	{
		if (x <= sx)
		{
			sx = x - (dimension.width / 5);

			if (sx < 0)
				sx = 0;

			horiz.setValue(sx);
			return true;
		}

		if (x >= sx + dimension.width)
		{
			sx = x - dimension.width + (dimension.width / 5);
			horiz.setValue(sx);
			return true;
		}

		return false;
	}


	private void redoLine(int line, int column)
	{
		String s;
		int x,y,height,width;

		s = (String)lines.elementAt(line);

		if (longest(s,line))
			redoControls(horiz.getValue(), vert.getValue(),false);
		dirty = true;
		updateFonts(null);
		height = fontHeight;
		y = (line - sy) * height;
		String	sub = s.substring(0, column-1);
		x = stringLength(sub) + EDGE;
		if (shiftHoriz(x))
		{
			repaint();
		}
		else
		{
			x = x - sx - EDGE;
			width = dimension.width - x;
			repaint(10, x, y, width, height);
		}

	}
	
	public void componentHidden(ComponentEvent e) {}
	public void componentShown(ComponentEvent e) {}
	public void componentMoved(ComponentEvent e) {}

	public void componentResized(ComponentEvent e)
	{
		redoControls(horiz.getValue(), vert.getValue(),true);
	}

	private void redoCanvas()
	{
		redoControls(horiz.getValue(), vert.getValue(),false);
	}

	private void redoControls(int h, int v, boolean clear)
	{
		int maxline, maxpix, bubble;

		Graphics g = getGraphics();

		if (g == null)
			return;

		updateFonts(g);

		if (g != null)
		{
			dimension = getSize();

			if (clear)
				clear_area(g);

			maxpix = wpix;
			bubble = dimension.width;

			if ( bubble > maxpix)
				maxpix = bubble;

			if ( h > maxpix )
				h = maxpix;

			if ( maxpix > 0 )
				maxpix += 12;

			horiz.setValues(h, bubble, 0, maxpix);

			maxline = lines.size();

			if ( maxline <= 0 )
				maxline = 10;

			bubble = ny = dimension.height / fontHeight;

			if (bubble > maxline)
				bubble = maxline;

			if (v > maxline)
				v = maxline;

			vert.setValues(v, bubble, 0, maxline);
			repaint();

			g.dispose();
		}
	}

	private void appendLine(String s)
	{
	        int count = lines.size();
		lines.insertElementAt( s, count );
		longest(s,count);	
	}

	private boolean longest(String s, int line)
	{
		int len = stringLength(s);
		if (len > wpix)
		{
			wpix = (len * 5) / 4;
			widest = line;
			return true;
		}
		return false;
	}

	private void resizeLines()
	{
		if (widest < lines.size())
		{
			int len = stringLength((String)lines.elementAt(widest));
			if (len > wpix)
			{
				wpix = len;
				redoControls(horiz.getValue(), vert.getValue(),false);
			}
		}

		if ((hline >= 0) && (hline < lines.size()))
			opix = hpix = pix_at(hline,hcolumn);

		if ((heline >= 0) && (heline < lines.size()))
			oepix = hepix = pix_at(heline,hecolumn);		
	}

	private int stringLength(String s)
	{
		if (fontMetrics == null)
			return -1;

		return fontMetrics.stringWidth(detabbed(s));
	}

	private String detabbed(String s)
	{
		if (s.indexOf('\t') < 0)
			return s;

		char c;
		String t = new String("");
		int j=0;
		int tabs;
		int max = s.length();

		for (int i=0; i<max; i++)
		{
			c = s.charAt(i);
			if (c == '\t')
			{
				tabs = tabSize - (j % tabSize);
				j += tabs;
				while (tabs-- > 0) t = t + ' ';
			}
			else
			{
				t = t + c;
				j++;
			}
		}

		return t;
	}

	public boolean isDirty()
	{
		return dirty;
	}

	public void getCursorPos(Graphics g, Rectangle r)
	{
		updateFonts(g);
		r.width = 2;
		r.height = fontHeight;
		r.y = (line - sy) * r.height;

		if (column > 0)
		{
			String s = ((String)lines.elementAt(line)).substring(0, column);
			r.x = stringLength(s) + EDGE - sx;
		}
		else
			r.x = EDGE - sx;
	}

	public void adjustmentValueChanged(AdjustmentEvent e)
	{
		pause_cursor();
		if (e.getSource() == horiz)
		{
			sx = e.getValue();
			horiz.setValue(sx);
		}
		else
		if (e.getSource() == vert)
		{
			sy = e.getValue();
			vert.setValue(sy);
		}
		repaint();
	}

	private final int SLOW_READ  = 2000;
	private final int SLOW_WRITE = 10000;

	public void read(File file)
	{
		Progress readProgress = null;
		boolean show_progress;
		long max,part,all,tens;
		String line;

		max = file.length();
		tens = max / 10;
		part = all = 0;
		show_progress = max > SLOW_READ;

		if (show_progress)
		{
			if (readProgress == null)
				readProgress = new Progress(parent,"Reading file ... ");

			readProgress.show();
		}

		clear();

		try
		{
			FileReader fr = new FileReader(file);
			BufferedReader br = new BufferedReader(fr);

			while ( (line=br.readLine()) != null)
			{
				if (show_progress)
				{
					part += line.length();
					if (part > tens)
					{
						all += part;
						readProgress.update( (int)(100 * all / max) );
						part = 0;
					}
				}
				appendLine(line);
			}

			br.close();
			fr.close();
		}

		catch (IOException e)
		{
			notify("Error - could not read file", NOACTION);
		}
		if (readProgress != null)
			readProgress.dispose();
	}

	public void write(File file)
	{
		Progress writeProgress = null;
		boolean show_progress;
		long max,part,all,tens;
		String line;

		max = lines.size();
		tens = max / 10;
		part = all = 0;
		show_progress = max > SLOW_WRITE;

		if (show_progress)
		{
			if (writeProgress == null)
				writeProgress = new Progress(parent,"Writing file ... ");

			writeProgress.show();
		}

		try
		{
			FileWriter fw = new FileWriter(file);
			BufferedWriter bw=new BufferedWriter(fw);

			for (int i=0;i<max;i++)
			{
				line = (String)lines.elementAt(i);
				bw.write(line,0,line.length());
				bw.newLine();
				if (show_progress)
				{
					part++;
					if (part > tens)
					{
						all += part;
						writeProgress.update( (int)(100 * all / max) );
						part = 0;
					}
				}
			}

			bw.close();
			fw.close();

			dirty = false;
		}

		catch (IOException e)
		{
			notify("Error - could not write file", NOACTION);
		}
		if (writeProgress != null)
			writeProgress.dispose();
	}

	public void update(Graphics g)
	{ 
		paint(g);
	}

	public void paint(Graphics g)
	{ 
		int i;
		Rectangle clip;
		Rectangle curs = new Rectangle();
		
		pause_cursor();

		g.setPaintMode();
		updateFonts(g);

		i = lines.size();
		if (oldlines != i)
		{
			oldlines = i;
			redoControls(horiz.getValue(), vert.getValue(),true);
		}
		
		if (dimension == null)
			dimension = getSize();

		clip = g.getClipBounds();

		i = clip.y / fontHeight;
		
		if ( clip.height / fontHeight == 1)
			drawLine(g, sy + i);
		else
			for (i = 0; i < ny; i++)
				drawLine(g,sy+i);

		release_cursor();
	}
	
	private synchronized void drawLine(Graphics g, int i)
	{
		int x,y,m;
		int bx,ex;
		String s;
		Rectangle curs = null;
				
		g.setPaintMode();
		g.setColor(Color.black);

		x = EDGE-sx;
		y = (i - sy + 1) * fontHeight;
		m = fontHeight;// + fontDescent;

		if (hactive && (i > oline) && (i < oeline))
			g.fillRect(0,y-fontHeight,dimension.width,m);
		else
			g.clearRect(0,y-fontHeight,dimension.width,m);
		
		if (i >= lines.size())
			return;
			
		s = detabbed((String)lines.elementAt(i));
		
		if (hactive && (i > oline) && (i < oeline))
		{
			g.setColor(Color.white);
			g.drawString(s, x, y-fontDescent);
			g.setColor(Color.black);
		}
		else
			g.drawString(s, x, y-fontDescent);
	
		if (hactive && ((i == oline) || (i == oeline)))
		{
			g.setXORMode(Color.white);
			bx = 0;
			ex = 5000;
			if (i == oline)
			{
				bx = opix - sx;
			}
			if (i == oeline)
			{
				ex = oepix - sx;
				if (oline == oeline)
					ex -= bx;
			}
			g.fillRect(bx,y-fontHeight,ex,m);
		}
	}

	private void notify(String s, int action)
	{
		Note note = new Note(parent,s,action);
		note.setVisible(true);
	}
}
