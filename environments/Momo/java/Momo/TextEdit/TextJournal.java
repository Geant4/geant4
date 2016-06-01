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
import java.util.*;

class	TextJournal
{
	private TextCanvas			canvas;		// the parent
	private Vector				lines;		// the document
	private TextJournalItem		journal;	// the journal
	private TextJournalItem		undo_list;	// where we are in undo
	private TextJournalItem		redo_list;	// where we are in redo

	private	boolean had_undo;
	private	boolean had_redo;

	private String	lineSeparator;

	private final int REPLACE_LINE	= 1;
	private final int SPLIT_LINE	= 2;
	private final int JOIN_LINE		= 3;
	private final int INSERT		= 4;
	private final int DELETE		= 5;
	
	public TextJournal(TextCanvas c, Vector v, String ls)
	{
		canvas = c;
		lines = v;
		lineSeparator = ls;
		undo_list = redo_list = null;
		had_undo = had_redo = false;
	}

	public void clear()
	{
		undo_list = redo_list = null;			
		had_undo = had_redo = false;
	}

	private void updateMenus()
	{
		boolean got_undo;
		boolean got_redo;

		got_undo = undo_list != null;
		got_redo = redo_list != null;

		if ((had_undo != got_undo) || (had_redo != got_redo))
			canvas.updateUndoItems(got_undo, got_redo);

		had_undo = got_undo;
		had_redo = got_redo;
	}

	public TextPosition undo()
	{
		TextJournalItem temp;
		TextPosition tp = null;

		if (undo_list != null)
		{
			temp = undo_list;
			undo_list = temp.next;

			switch (temp.action)
			{
			case REPLACE_LINE:
				redo_line(temp);
				break;
			case SPLIT_LINE:
				split_or_join(temp,true);
				break;
			case JOIN_LINE:
				split_or_join(temp,false);
				break;
			case INSERT:
				copy_or_cut(temp,true);
				break;
			case DELETE:
				insert(temp);
				break;
			default:;
			}

			tp = new TextPosition(temp.line,temp.column);

			temp.next = redo_list;
			redo_list = temp;
		}

		updateMenus();
		return tp;
	}

	public TextPosition redo()
	{
		TextJournalItem temp;
		TextPosition tp = null;
		int line, column;

		if (redo_list != null)
		{
			temp = redo_list;
			redo_list = temp.next;

			line = temp.line;
			column = temp.column;

			switch (temp.action)
			{
			case REPLACE_LINE:
				column = temp.text.length();
				redo_line(temp);
				break;
			case SPLIT_LINE:
				split_or_join(temp,false);
				column = 0;
				break;
			case JOIN_LINE:
				split_or_join(temp,true);
				break;
			case INSERT:
				insert(temp);
				line = temp.eline;
				column = temp.ecolumn;
				break;
			case DELETE:
				copy_or_cut(temp,true);
				break;
			default:;
			}

			tp = new TextPosition(line,column);

			temp.next = undo_list;
			undo_list = temp;
		}

		updateMenus();
		return tp;
	}

	public void insert_char(int line, int column, char c)
	{
		String  s;

		s = (String)lines.elementAt(line);

		remember_line(line,column,s);

		s = s.substring(0, column) + c + s.substring(column, s.length());
		lines.setElementAt(s,line);
	}

	public void delete_prev_char(int line, int column)
	{
		String s;

		s = (String)lines.elementAt(line);

		remember_line(line,column,s);

		s = s.substring(0, column-1) + s.substring(column, s.length());
		lines.setElementAt(s,line);
	}

	public void delete_next_char(int line, int column)
	{
		String s;

		s = (String)lines.elementAt(line);

		remember_line(line,column,s);

		s = s.substring(0, column) + s.substring(column+1, s.length());
		lines.setElementAt(s,line);
	}

	private void remember_line(int line, int column, String s)
	{
		boolean new_line = true;

		if (undo_list != null)
			new_line = (undo_list.action != REPLACE_LINE) || (undo_list.line != line);

		if (new_line)
		{
			TextJournalItem new_journal = new TextJournalItem();
			new_journal.action = REPLACE_LINE;
			new_journal.line = line;
			new_journal.column = column;
			new_journal.text = new String(s);
			new_journal.next = undo_list;
			undo_list = new_journal;
		}

		redo_list = null;
		updateMenus();
	}

	private void redo_line(TextJournalItem i)
	{
		String s;
		s = (String)lines.elementAt(i.line);
		lines.setElementAt(i.text,i.line);
		i.text = s;
	}

	public void split_line(int line, int column)
	{
		TextJournalItem new_journal = new TextJournalItem();
		new_journal.action = SPLIT_LINE;
		new_journal.line = line;
		new_journal.column = column;

		split_or_join(new_journal,false);

		new_journal.next = undo_list;
		undo_list = new_journal;
		redo_list = null;
		updateMenus();
	}

	public void join_line(int line, int column)
	{
		TextJournalItem new_journal = new TextJournalItem();
		new_journal.action = JOIN_LINE;
		new_journal.line = line;
		new_journal.column = column;

		split_or_join(new_journal,true);

		new_journal.next = undo_list;
		undo_list = new_journal;
		redo_list = null;
		updateMenus();
	}

	public void split_or_join(TextJournalItem i, boolean join)
	{
		String s;
		int line = i.line;
		int column = i.column;

		if (join)
		{
			s = (String)lines.elementAt(line);
			s = s.concat((String)lines.elementAt(line+1));
			lines.setElementAt(s,line);
			lines.removeElementAt(line+1);
		}
		else
		{
			s = (String)lines.elementAt(line);
			lines.setElementAt(s.substring(0,column),line);
			lines.insertElementAt(s.substring(column,s.length()),++line);
		}
	}

	public TextPosition insert_section(int line, int column, String s)
	{
		int charCt;
		int saveCt;
		int charMax;
		char c,c2;
		String s2 = null;
		TextPosition tp;

		// start creating journal item

		TextJournalItem new_journal = new TextJournalItem();
		new_journal.action = INSERT;
		new_journal.text = new String(s);
		new_journal.line = line;		// starting line and column
		new_journal.column = column;

		// insert the text

		tp = insert(new_journal);
	
		// finish creating journal entry

		new_journal.eline = tp.line;		// updated line and column
		new_journal.ecolumn = tp.column;
		new_journal.next = undo_list;
		undo_list = new_journal;
		redo_list = null;
		updateMenus();

		return tp;
	}

	public TextPosition insert(TextJournalItem i)
	{
		int line = i.line;
		int column = i.column;
		String s = i.text;

		int charCt;
		int saveCt;
		int charMax;
		char c,c2;
		String s2 = null;

		// insert the text

		charMax = s.length();
		charCt = saveCt = 0;
		
		while (charCt < charMax)
		{
			c = s.charAt(charCt);
			charCt++;
			if ((c == '\r') || (c == '\n'))
			{
				s2 = (String)lines.elementAt(line);
				lines.setElementAt(s2.substring(0,column) + s.substring(saveCt,charCt-1),line);
				lines.insertElementAt(s2.substring(column,s2.length()),++line);
				column = 0;
				if (charCt < charMax)
				{
					c2 = s.charAt(charCt);
					if (((c == '\r') && (c2 == '\n')) || ((c2 == '\r') && (c == '\n')))
						charCt++;
				}
				saveCt = charCt;
			}
		}
		
		if (saveCt < charCt)
		{
			s2 = (String)lines.elementAt(line);
			s = s.substring(saveCt,charCt);
			if (column == 0)
				s2 = s + s2;
			else
				s2 = s2.substring(0, column) + s + s2.substring(column, s2.length());
			column += s.length();
			lines.setElementAt(s2,line);
		}

		return new TextPosition(line,column);
	}

	public String delete_section(int line, int column, int eline, int ecolumn, boolean cut)
	{
		String text;

		TextJournalItem new_journal = new TextJournalItem();
		new_journal.action = DELETE;
		new_journal.line = line;
		new_journal.column = column;
		new_journal.eline = eline;
		new_journal.ecolumn = ecolumn;

		text = copy_or_cut(new_journal,cut);

		if (cut)
		{
			new_journal.text = new String(text);
			new_journal.next = undo_list;
			undo_list = new_journal;
			redo_list = null;
			updateMenus();
		}

		return text;
	}

	private String copy_or_cut(TextJournalItem i, boolean cut)
	{
		String s,s2;
		String text = null;

		int line = i.line;
		int column = i.column;
		int eline = i.eline;
		int ecolumn = i.ecolumn;
	
		if (line == eline)
		{
			s = (String)lines.elementAt(line);
			text = s.substring(column,ecolumn);
			if (cut)
			{
				s = s.substring(0, column) + s.substring(ecolumn, s.length());
				lines.setElementAt(s,line);
			}
		}
		else
		{
			s = (String)lines.elementAt(line);
			text = s.substring(column,s.length());
			s = s.substring(0, column);

			int diff = eline - line;
			int k = line + 1;

			for (int j = 1; j < diff; j++)
			{
				s2 = (String)lines.elementAt(k);
				text = text + lineSeparator + s2;
				if (cut)
					lines.removeElementAt(k);
				else
					k++;
			}

			if (k != lines.size())
    		{
				s2 = (String)lines.elementAt(k);
				text = text + lineSeparator + s2.substring(0,ecolumn);
				s = s + s2.substring(ecolumn,s2.length());
				if (cut)
					lines.removeElementAt(k);
			}
			
			if (cut)
				lines.setElementAt(s,line);
		}

		return text;
	}

}	

class TextJournalItem
{
	public TextJournalItem	next;
	public String			text;
	public int				action;
	public int				line;
	public int				column;
	public int				eline;
	public int				ecolumn;
}
