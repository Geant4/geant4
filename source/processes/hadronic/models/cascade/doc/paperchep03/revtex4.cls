%%
%% This is file `revtex4.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% revtex4.dtx  (with options: `class-pre')
%% ltxutil.dtx  (with options: `ltxutil-krn')
%% ltxgrid.dtx  (with options: `ltxgrid-krn')
%% textcase.dtx  (with options: `package')
%% revtex4.dtx  (with options: `class-pst')
%% 
%% This file is part of the APS files in the REVTeX 4 distribution.
%% For the version number, search on the string 
%% Original version by David Carlisle
%% Modified by Arthur Ogawa (mailto:ogawa@teleport.com)
%% 
%% Copyright (c) 1999 The American Physical Society.
%% http://publish.aps.org/revtex4/
%% mailto:revtex4@aps.org
%% 
%% See the REVTeX 4 README file for restrictions and more information.
%% 
\NeedsTeXFormat{LaTeX2e}[1996/12/01]
\ProvidesClass{revtex4}
          [2001/08/03 v4.0 (http://publish.aps.org/revtex4/ for documentation)]%FileInfo
\def\class@name{revtex4}%
\GenericInfo{}{\space
 Copyright (c) 2001 The American Physical Society.^^J
 mailto:revtex@aps.org^^J
 Licensed under the LPPL:^^Jhttp://www.ctan.org/tex-archive/macros/latex/base/lppl.txt^^J
 Arthur Ogawa <ogawa@teleport.com>^^J
 Based on work by David Carlisle <david@dcarlisle.demon.co.uk>%
 \@gobble
}%
\if@compatibility
  \edef\reserved@a{\errhelp{%
   Change your \string\documentstyle\space statement to
   \string\documentclass\space and rerun.
  }}\reserved@a
  \errmessage{You cannot run \class@name\space in compatability mode}%
  \expandafter\@@end
\fi
\let\DeclareOption@latex\DeclareOption
\let\DeclareOption\@gobbletwo
\let\ProcessOptions@latex\ProcessOptions
\let\ProcessOptions\@empty
\let\ProvidesPackage@latex\ProvidesPackage
\def\ProvidesPackage#1[#2]{%
 \class@info{Incorporating package #1 [#2].}%
}%
\typeout{%
 ltxutil: portions licensed from W. E. Baxter (web@superscript.com)%
}%
\def\class@err#1{\ClassError{\class@name}{#1}\@eha}%
\def\class@warn#1{\ClassWarningNoLine{\class@name}{#1}}%
\def\class@info#1{\ClassInfo{\class@name}{#1}}%
\def\obsolete@command#1{%
 \class@warn@end{Command \string#1\space is obsolete.^^JPlease remove from your document}%
 \global\let#1\@empty
 #1%
}%
\def\replace@command#1#2{%
 \class@warn@end{Command \string#1\space is obsolete;^^JUse \string#2\space instead}%
 \global\let#1#2%
 #1%
}%
\def\replace@environment#1#2{%
 \class@warn@end{Environment #1 is obsolete;^^JUse #2 instead}%
 \glet@environment{#1}{#2}%
 \@nameuse{#1}%
}%
\def\incompatible@package#1{%
 \@ifpackageloaded{#1}{%
  \def\@tempa{I cannot continue. You must remove the \string\usepackage\ statement that caused that package to be loaded.}%
  \ClassError{\class@name}{The #1 package cannot be used with \class@name}%
  \@tempa\stop
 }{%
  \class@info{#1 was not loaded (OK!)}%
 }%
}%
\def\class@warn@end#1{%
 \gappdef\class@enddocumenthook{\class@warn{#1}}%
}%
\AtEndOfClass{%
 \@ifxundefined\class@name{\def\class@name{Generic Class}}{}%
}%
\def\t@{to}%
\dimendef\dimen@iii\thr@@
\def\halignt@{\halign\t@}%
\chardef\f@ur=4\relax
\chardef\cat@letter=11\relax
\chardef\other=12\relax
\def\let@environment#1#2{%
 \expandafter\let
 \csname#1\expandafter\endcsname\csname#2\endcsname
 \expandafter\let
 \csname end#1\expandafter\endcsname\csname end#2\endcsname
}%
\def\glet@environment#1#2{%
 \global\expandafter\let
 \csname#1\expandafter\endcsname\csname#2\endcsname
 \global\expandafter\let
 \csname end#1\expandafter\endcsname\csname end#2\endcsname
}%
\newcommand\tracingplain{%
 \tracingonline\z@\tracingcommands\z@\tracingstats\z@
 \tracingpages\z@\tracingoutput\z@\tracinglostchars\@ne
 \tracingmacros\z@\tracingparagraphs\z@\tracingrestores\z@
 \showboxbreadth5\showboxdepth3\relax %\errorstopmode
 }%
\newcommand\traceoutput{%
 \appdef\@resetactivechars{\showoutput}%
}%
\newcommand\say[1]{\typeout{<\noexpand#1=\meaning#1>}}%
\newcommand\saythe[1]{\typeout{<\noexpand#1=\the#1>}}%
\def\fullinterlineskip{\prevdepth\z@}%
\countdef\count@i\@ne
\countdef\count@ii\tw@
\long\def\prepdef#1#2{%
 \@ifxundefined#1{\toks@{}}{\toks@\expandafter{#1}}%
 \toks@ii{#2}%
 \edef#1{\the\toks@ii\the\toks@}%
}%
\long\def\appdef#1#2{%
 \@ifxundefined#1{\toks@{}}{\toks@\expandafter{#1}}%
 \toks@ii{#2}%
 \edef#1{\the\toks@\the\toks@ii}%
}%
\long\def\gappdef#1#2{%
 \@ifxundefined#1{\toks@{}}{\toks@\expandafter{#1}}%
 \toks@ii{#2}%
 \global\edef#1{\the\toks@\the\toks@ii}%
}%
\long\def\appdef@val#1#2{%
 \appdef#1{{#2}}%
}%
\long\def\appdef@e#1#2{%
 \expandafter\appdef
 \expandafter#1%
 \expandafter{#2}%
}%
\long\def\appdef@eval#1#2{%
 \expandafter\appdef@val
 \expandafter#1%
 \expandafter{#2}%
}%
\toksdef\toks@ii=\tw@
\long\def\@ifxundefined#1{\@ifx{\undefined#1}}%
\long\def\@ifnotrelax#1#2#3{\@ifx{\relax#1}{#3}{#2}}%
\long\def\@argswap#1#2{#2#1}%
\long\def\@argswap@val#1#2{#2{#1}}%
\def\@ifxundefined@cs#1{\expandafter\@ifx\expandafter{\csname#1\endcsname\relax}}%
\def\@boolean#1#2{%
  \long\def#1{%
    #2% \if<something>
      \expandafter\true@sw
    \else
      \expandafter\false@sw
    \fi
  }%
}%
\def\@boole@def#1#{\@boolean{#1}}% Implicit #2
\def\@booleantrue#1{\let#1\true@sw}%
\def\@booleanfalse#1{\let#1\false@sw}%
\@boole@def\@ifx#1{\ifx#1}%
\@boole@def\@ifx@empty#1{\ifx\@empty#1}%
\@boole@def\@if@empty#1{\if!#1!}%
\def\@if@sw#1#2{#1\expandafter\true@sw\else\expandafter\false@sw#2}%
\@boole@def\@ifdim#1{\ifdim#1}%
\@boole@def\@ifeof#1{\ifeof#1}%
\@boole@def\@ifhbox#1{\ifhbox#1}%
\@boole@def\@ifhmode{\ifhmode}%
\@boole@def\@ifinner{\ifinner}%
\@boole@def\@ifmmode{\ifmmode}%
\@boole@def\@ifnum#1{\ifnum#1}%
\@boole@def\@ifodd#1{\ifodd#1}%
\@boole@def\@ifvbox#1{\ifvbox#1}%
\@boole@def\@ifvmode{\ifvmode}%
\@boole@def\@ifvoid#1{\ifvoid#1}%
\long\def\true@sw#1#2{#1}%
\long\def\false@sw#1#2{#2}%
\long\def\loopuntil#1{#1{}{\loopuntil{#1}}}%
\long\def\loopwhile#1{#1{\loopwhile{#1}}{}}%
\def\@provide#1{%
 \@ifx{\undefined#1}{\true@sw}{\@ifx{\relax#1}{\true@sw}{\false@sw}}%
 {\def#1}{\def\j@nk}%
}%
\prepdef\document{%
 \endgroup
 \init@documenthook
 \set@typesize@hook
 \normalsize
 \set@pica@hook
 \true@sw{}%
}%
\def\init@documenthook{}%
\AtBeginDocument{%
 \class@documenthook
}%
\AtEndDocument{%
 \class@enddocumenthook
}%
\def\class@documenthook{}%
\def\class@enddocumenthook{}%
\def\set@typesize@hook{}%
\def\set@pica@hook{}%
\def\enddocument{%
 \@enddocumenthook
 \@checkend{document}%
 \clear@document
 \check@aux
 \deadcycles\z@
 \@@end
}%
\def\clear@document{\clearpage}%
\def\check@aux{\do@check@aux}%
\def\do@check@aux{%
 \@if@sw\if@filesw\fi{%
  \immediate\closeout\@mainaux
  \let\@setckpt\@gobbletwo
  \let\@newl@bel\@testdef
  \@tempswafalse
  \makeatletter
  \input\jobname.aux\relax
 }{}%
 \@dofilelist
 \@ifdim{\font@submax >\fontsubfuzz\relax}{%
  \@font@warning{%
   Size substitutions with differences\MessageBreak
   up to \font@submax\space have occured.\@gobbletwo
  }%
 }{}%
 \@defaultsubs
 \@refundefined
 \@if@sw\if@filesw\fi{%
  \@ifx{\@multiplelabels\relax}{%
   \@if@sw\if@tempswa\fi{%
    \@latex@warning@no@line{%
     Label(s) may have changed.
     Rerun to get cross-references right
    }%
   }{}%
  }{%
    \@multiplelabels
  }%
 }{}%
}%
\def\flushing{%
  \let\\\@normalcr
  \leftskip\z@skip
  \rightskip\z@skip
  \@rightskip\z@skip
  \parfillskip\@flushglue
}%
\def\eqnarray@LaTeX{%
   \stepcounter{equation}%
   \def\@currentlabel{\p@equation\theequation}%
   \global\@eqnswtrue
   \m@th
   \global\@eqcnt\z@
   \tabskip\@centering
   \let\\\@eqncr
   $$\everycr{}\halign\t@\displaywidth\bgroup
       \hskip\@centering$\displaystyle\tabskip\z@skip{##}$\@eqnsel
      &\global\@eqcnt\@ne\hskip \tw@\arraycolsep \hfil${##}$\hfil
      &\global\@eqcnt\tw@ \hskip \tw@\arraycolsep
         $\displaystyle{##}$\hfil\tabskip\@centering
      &\global\@eqcnt\thr@@ \hb@xt@\z@\bgroup\hss##\egroup
         \tabskip\z@skip
      \cr
}
\long\def\eqnarray@fleqn@fixed{%
 \stepcounter{equation}\def\@currentlabel{\p@equation\theequation}%
 \global\@eqnswtrue\m@th\global\@eqcnt\z@
 \tabskip\mathindent
 \let\\=\@eqncr
 \setlength\abovedisplayskip{\topsep}%
 \ifvmode\addtolength\abovedisplayskip{\partopsep}\fi
 \addtolength\abovedisplayskip{\parskip}%
 \setlength\belowdisplayskip{\abovedisplayskip}%
 \setlength\belowdisplayshortskip{\abovedisplayskip}%
 \setlength\abovedisplayshortskip{\abovedisplayskip}%
 $$%
 \everycr{}%
 \halignt@\linewidth\bgroup
  \hskip\@centering$\displaystyle\tabskip\z@skip{##}$\@eqnsel
  &\global\@eqcnt\@ne
   \hskip\tw@\eqncolsep
   \hfil${{}##{}}$\hfil
  &\global\@eqcnt\tw@
   \hskip\tw@\eqncolsep
   $\displaystyle{##}$\hfil\tabskip\@centering
  &\global\@eqcnt\thr@@\hb@xt@\z@\bgroup\hss##\egroup
   \tabskip\z@skip
  \cr
}%
\@ifx{\eqnarray\eqnarray@LaTeX}{%
 \class@info{Repairing broken LaTeX eqnarray}%
 \let\eqnarray\eqnarray@fleqn@fixed
 \newlength\eqncolsep
 \setlength\eqncolsep\z@
 \let\eqnarray@LaTeX\relax
 \let\eqnarray@fleqn@fixed\relax
}{}%
\def\mathindent{\@centering}%
\def\set@eqnarray@skips{}%
\def\footnote{%
 \@ifnextchar[\@xfootnote{\@yfootnote\@footnotetext}%
}%
\def\footnotemark{%
 \@ifnextchar[\@xfootnotemark{\@yfootnote}%
}%
\def\@xfootnote[#1]{%
 \@xfootnotemark[#1]%
 \@footnotetext
}%
\def\@xfootnotemark@ltx[#1]{%
 \begingroup
   \csname c@\@mpfn\endcsname #1\relax
   \unrestored@protected@xdef\@thefnmark{\thempfn}%
 \endgroup
 \H@@footnotemark
}%
\def\@yfootnote{%
 \stepcounter\@mpfn
 \protected@xdef\@thefnmark{\thempfn}%
 \H@@footnotemark
}%
\appdef\class@documenthook{%
 \@ifxundefined\H@@footnotemark{%
  \let\H@@footnotemark\@footnotemark
 }{}%
 \let\@xfootnotemark\@xfootnotemark@ltx
}%
\long\def\@footnotetext{%
 \insert\footins\bgroup
  \make@footnotetext
}%
\long\def\@mpfootnotetext{%
 \minipagefootnote@pick
  \make@footnotetext
}%
\def\make@footnotetext#1{%
  \reset@font\footnotesize
  \interlinepenalty\interfootnotelinepenalty
  \splittopskip\footnotesep
  \splitmaxdepth\dp\strutbox
  \set@footnotewidth
  \@parboxrestore
  \protected@edef\@currentlabel{%
   \csname p@footnote\endcsname\@thefnmark
  }%
  \color@begingroup
   \@makefntext{%
    \rule\z@\footnotesep\ignorespaces#1\@finalstrut\strutbox
   }%
  \color@endgroup
 \minipagefootnote@drop
}%
\def\set@footnotewidth{%
 \hsize\columnwidth
 \linewidth\hsize
}%
\def\set@footnotewidth@ii{%
 \hsize\textwidth
 \advance\hsize\columnsep
 \divide\hsize\tw@
 \advance\hsize-\columnsep
 \linewidth\hsize
}%
\def\addtocontents#1#2{%
 \protected@write\@auxout{%
  \let \label \@gobble \let \index \@gobble \let \glossary \@gobble
  \def\({\string\(}%
  \def\){\string\)}%
  \def\\{\string\\}%
 }{\string \@writefile {#1}{#2}}%
}%
\def\addcontentsline#1#2#3{%
 \addtocontents{#1}{%
  \protect\contentsline{#2}{#3}{\thepage}{}%
 }%
}%
\def\contentsline#1#2#3#4{%
 \csname l@#1\endcsname{#2}{#3}%
}%
\def\label#1{%
 \@bsphack
  \protected@write\@auxout{}{%
   \string\newlabel{#1}{{\@currentlabel}{\thepage}{}{}{}}%
  }%
 \@esphack
}%
\appdef\class@documenthook{%
 \prepdef\caption{\minipagefootnote@here}%
}%
\def\minipagefootnote@init{%
 \setbox\@mpfootins\box\voidb@x
}%
\def\minipagefootnote@pick{%
 \global\setbox\@mpfootins\vbox\bgroup
  \unvbox\@mpfootins
}%
\def\minipagefootnote@drop{%
 \egroup
}%
\def\minipagefootnote@here{%
    \par
    \@ifvoid\@mpfootins{}{%
      \vskip\skip\@mpfootins
      \fullinterlineskip
      \@ifinner{%
       \vtop{\unvcopy\@mpfootins}%
       {\setbox\z@\lastbox}%
      }{}%
      \unvbox\@mpfootins
    }%
}%
\def\minipagefootnote@foot{%
 \@ifvoid\@mpfootins{}{%
  \insert\footins\bgroup\unvbox\@mpfootins\egroup
 }%
}%
\def\endminipage{%
    \par
    \unskip
    \minipagefootnote@here
    \@minipagefalse   %% added 24 May 89
  \color@endgroup
  \egroup
  \expandafter\@iiiparbox\@mpargs{\unvbox\@tempboxa}%
}%
\let\@xfloat@LaTeX\@xfloat
\def\@xfloat#1[#2]{%
  \@xfloat@prep
  \@nameuse{fp@proc@#2}%
  \@ifxundefined\floats@sw{\global\@booleantrue\floats@sw}{}%
  \floats@sw{\@xfloat@LaTeX{#1}[#2]}{\@xfloat@anchored{#1}[]}%
}%
\def\@xfloat@prep{%
  \let\footnote\footnote@latex
  \def\@mpfn{mpfootnote}%
  \def\thempfn{\thempfootnote}%
  \c@mpfootnote\z@
  \let\@footnotetext\@mpfootnotetext
  \let\H@@footnotetext\@mpfootnotetext
  \let\@makefntext\@mpmakefntext
}%
\appdef\class@documenthook{%
 \let\footnote@latex\footnote
}%
\def\@xfloat@anchored#1[#2]{%
 \def\@captype{#1}%
 \begin@float@pagebreak
  \let\end@float\end@float@anchored
  \let\end@dblfloat\end@float@anchored
        \hsize\columnwidth
        \@parboxrestore
        \@floatboxreset
  \minipagefootnote@init
}%
\def\end@float@anchored{%
  \minipagefootnote@here
  \par\vskip\z@skip      %% \par\vskip\z@ added 15 Dec 87
 \par
 \end@float@pagebreak
}%
\def\begin@float@pagebreak{\par\addvspace\intextsep}%
\def\end@float@pagebreak{\par\addvspace\intextsep}%
\def\@mpmakefntext#1{%
 \parindent=1em
 \noindent
 \hb@xt@1em{\hss\@makefnmark}%
 #1%
}%
\def\do@if@floats#1#2{%
 \@ifxundefined\floats@sw{\global\@booleantrue\floats@sw}{}%
 \floats@sw{}{%
  \expandafter\newwrite
              \csname#1write\endcsname
  \expandafter\def
              \csname#1@stream\endcsname{\jobname#2}%
  \expandafter\immediate
  \expandafter\openout
              \csname#1write\endcsname
              \csname#1@stream\endcsname\relax
  \@ifxundefined\@float@LaTeX{%
   \let\@float@LaTeX\@float
   \let\@dblfloat@LaTeX\@dblfloat
   \let\@float\write@float
   \let\@dblfloat\write@floats
  }{}%
  \let@environment{#1@float}{#1}%
  \let@environment{#1@floats}{#1*}%
  \@ifxundefined@cs{#1@write}{}{%
   \let@environment{#1}{#1@write}%
  }%
 }%
}%
\def\triggerpar{\leavevmode\@@par}%
\def\oneapage{\def\begin@float@pagebreak{\newpage}\def\end@float@pagebreak{\newpage}}%
\def\print@float#1#2{%
 \@ifxundefined@cs{#1write}{}{%
  \begingroup
   \@booleanfalse\floats@sw
   #2%
   \raggedbottom
   \def\array@default{v}% floats must
   \let\@float\@float@LaTeX
   \let\@dblfloat\@dblfloat@LaTeX
   \let\trigger@float@par\triggerpar
   \let@environment{#1}{#1@float}%
   \let@environment{#1*}{#1@floats}%
   \expandafter\prepdef\csname#1\endcsname{\trigger@float@par}%
   \expandafter\prepdef\csname#1*\endcsname{\trigger@float@par}%
   \@namedef{fps@#1}{h!}%
   \expandafter\immediate
   \expandafter\closeout
               \csname#1write\endcsname
   \everypar{%
    \global\let\trigger@float@par\relax
    \global\everypar{}\setbox\z@\lastbox
    \@ifxundefined@cs{#1sname}{}{%
     \begin@float@pagebreak
     \expandafter\section
     \expandafter*%
     \expandafter{%
                  \csname#1sname\endcsname
                 }%
    }%
   }%
   \input{\csname#1@stream\endcsname}%
  \endgroup
  \global\expandafter\let\csname#1write\endcsname\relax
 }%
}%
\def\write@float#1{\write@@float{#1}{#1}}%
\def\endwrite@float{\@Esphack}%
\def\write@floats#1{\write@@float{#1*}{#1}}%
\def\endwrite@floats{\@Esphack}%
\def\write@@float#1#2{%
  \ifhmode
     \@bsphack
  \fi
  \chardef\@tempc\csname#2write\endcsname
  \toks@{\begin{#1}}%
  \def\@tempb{#1}%
  \expandafter\let\csname end#1\endcsname\endwrite@float
  \catcode`\^^M\active
  \@makeother\{\@makeother\}\@makeother\%
  \write@floatline
}%
\begingroup
 \catcode`\[\the\catcode`\{\catcode`\]\the\catcode`\}\@makeother\{\@makeother\}%
 \gdef\float@end@tag#1\end{#2}#3\@nul[%
  \def\@tempa[#2]%
  \@ifx[\@tempa\@tempb][\end[#2]][\write@floatline]%
 ]%
 \obeylines%
 \gdef\write@floatline#1^^M[%
  \begingroup%
   \newlinechar`\^^M%
   \toks@\expandafter[\the\toks@#1]\immediate\write\@tempc[\the\toks@]%
  \endgroup%
  \toks@[]%
  \float@end@tag#1\end{}\@nul%
 ]%
\endgroup
\def\@alph#1{\ifcase#1\or a\or b\or c\or d\else\@ialph{#1}\fi}
\def\@ialph#1{\ifcase#1\or \or \or \or \or e\or f\or g\or h\or i\or j\or
  k\or l\or m\or n\or o\or p\or q\or r\or s\or t\or u\or v\or w\or x\or
  y\or z\or aa\or bb\or cc\or dd\or ee\or ff\or gg\or hh\or ii\or jj\or
  kk\or ll\or mm\or nn\or oo\or pp\or qq\or rr\or ss\or tt\or uu\or
  vv\or ww\or xx\or yy\or zz\else\@ctrerr\fi}
\def\@startsection#1#2#3#4#5#6{%
 \@startsection@hook
 \if@noskipsec \leavevmode \fi
 \par
 \@tempskipa #4\relax
 \@afterindenttrue
 \ifdim \@tempskipa <\z@
  \@tempskipa -\@tempskipa \@afterindentfalse
 \fi
 \if@nobreak
  \everypar{}%
 \else
  \addpenalty\@secpenalty\addvspace\@tempskipa
 \fi
 \@ifstar
  {\@dblarg{\@ssect@ltx{#1}{#2}{#3}{#4}{#5}{#6}}}%
  {\@dblarg{\@sect@ltx {#1}{#2}{#3}{#4}{#5}{#6}}}%
}%
\def\@startsection@hook{}%
\class@info
  {Repairing broken LateX \string\@sect}%
\def\@sect@ltx#1#2#3#4#5#6[#7]#8{%
  \@ifnum{#2>\c@secnumdepth}{%
    \def\H@svsec{\phantomsection}%
    \let\@svsec\@empty
  }{%
    \H@refstepcounter{#1}%
    \def\H@svsec{%
     \phantomsection
    }%
    \protected@edef\@svsec{{#1}}%
    \@ifundefined{@#1cntformat}{%
     \prepdef\@svsec\@seccntformat
    }{%
     \expandafter\prepdef
     \expandafter\@svsec
                 \csname @#1cntformat\endcsname
    }%
  }%
  \@tempskipa #5\relax
  \@ifdim{\@tempskipa>\z@}{%
    \begingroup
      \interlinepenalty \@M
      #6{%
       \@ifundefined{@hangfrom@#1}{\@hang@from}{\csname @hangfrom@#1\endcsname}%
       {\hskip#3\relax\H@svsec}{\@svsec}{#8}%
      }%
      \@@par
    \endgroup
    \@ifundefined{#1mark}{\@gobble}{\csname #1mark\endcsname}{#7}%
    \addcontentsline{toc}{#1}{%
      \@ifnum{#2>\c@secnumdepth}{%
       \protect\numberline{}%
      }{%
       \protect\numberline{\csname the#1\endcsname}%
      }%
      #8}%
  }{%
    \def\@svsechd{%
      #6{%
       \@ifundefined{@runin@to@#1}{\@runin@to}{\csname @runin@to@#1\endcsname}%
       {\hskip#3\relax\H@svsec}{\@svsec}{#8}%
      }%
      \@ifundefined{#1mark}{\@gobble}{\csname #1mark\endcsname}{#7}%
      \addcontentsline{toc}{#1}{%
        \@ifnum{#2>\c@secnumdepth}{%
         \protect\numberline{}%
        }{%
         \protect\numberline{\csname the#1\endcsname}%
        }%
        #8}%
    }%
  }%
  \@xsect{#5}%
}%
\def\@hang@from#1#2#3{\@hangfrom{#1#2}#3}%
\def\@runin@to #1#2#3{#1#2#3}%
\def\@ssect@ltx#1#2#3#4#5#6[#7]#8{%
  \def\H@svsec{\phantomsection}%
  \@tempskipa #5\relax
  \@ifdim{\@tempskipa>\z@}{%
    \begingroup
      \interlinepenalty \@M
      #6{%
       \@ifundefined{@hangfroms@#1}{\@hang@froms}{\csname @hangfroms@#1\endcsname}%
       {\hskip#3\relax\H@svsec}{#8}%
      }%
      \@@par
    \endgroup
    \@ifundefined{#1smark}{\@gobble}{\csname #1smark\endcsname}{#7}%
    \addcontentsline{toc}{#1}{\protect\numberline{}#8}%
  }{%
    \def\@svsechd{%
      #6{%
       \@ifundefined{@runin@tos@#1}{\@runin@tos}{\csname @runin@tos@#1\endcsname}%
       {\hskip#3\relax\H@svsec}{#8}%
      }%
      \@ifundefined{#1smark}{\@gobble}{\csname #1smark\endcsname}{#7}%
      \addcontentsline{toc}{#1}{\protect\numberline{}#8}%
    }%
  }%
  \@xsect{#5}%
}%
\def\@hang@froms#1#2{#1#2}%
\def\@runin@tos #1#2{#1#2}%
\appdef\init@documenthook{%
 \providecommand\phantomsection{}%
 \providecommand\hyper@anchor[1]{}%
 \providecommand\hyper@last{}%
 \providecommand\Hy@raisedlink[1]{#1}%
 \providecommand\hyper@anchorstart[1]{}%
 \providecommand\hyper@anchorend{}%
 \providecommand\hyper@linkstart[2]{}%
 \providecommand\hyper@linkend{}%
}%
\let\H@refstepcounter\refstepcounter
\def\sec@upcase#1{\relax{#1}}%
\appdef\class@documenthook{%
 \@ifpackageloaded{array}{\switch@array}{\switch@tabular}%
 \prepdef\endtabular{\endtabular@hook}%
 \@provide\endtabular@hook{}%
 \prepdef\endarray{\endarray@hook}%
 \@provide\endarray@hook{}%
 \providecommand\array@hook{}%
 \prepdef\@tabular{\tabular@hook}%
 \@provide\tabular@hook{}%
}%
\def\switch@tabular{%
 \let\@array@sw\@array@sw@array
 \@ifx{\@array\@array@LaTeX}{%
  \@ifx{\multicolumn\multicolumn@LaTeX}{%
   \@ifx{\@tabular\@tabular@LaTeX}{%
    \@ifx{\@tabarray\@tabarray@LaTeX}{%
     \@ifx{\array\array@LaTeX}{%
      \@ifx{\endarray\endarray@LaTeX}{%
       \@ifx{\endtabular\endtabular@LaTeX}{%
        \@ifx{\@mkpream\@mkpream@LaTeX}{%
         \@ifx{\@addamp\@addamp@LaTeX}{%
          \@ifx{\@arrayacol\@arrayacol@LaTeX}{%
           \@ifx{\@tabacol\@tabacol@LaTeX}{%
            \@ifx{\@arrayclassz\@arrayclassz@LaTeX}{%
             \@ifx{\@tabclassiv\@tabclassiv@LaTeX}{%
              \@ifx{\@arrayclassiv\@arrayclassiv@LaTeX}{%
               \@ifx{\@tabclassz\@tabclassz@LaTeX}{%
                \@ifx{\@classv\@classv@LaTeX}{%
                 \@ifx{\hline\hline@LaTeX}{%
                  \@ifx{\@tabularcr\@tabularcr@LaTeX}{%
                   \@ifx{\@xtabularcr\@xtabularcr@LaTeX}{%
                    \@ifx{\@xargarraycr\@xargarraycr@LaTeX}{%
                     \@ifx{\@yargarraycr\@yargarraycr@LaTeX}{%
                      \true@sw
                     }{%
                      \false@sw
                     }%
                    }{%
                     \false@sw
                    }%
                   }{%
                    \false@sw
                   }%
                  }{%
                   \false@sw
                  }%
                 }{%
                  \false@sw
                 }%
                }{%
                 \false@sw
                }%
               }{%
                \false@sw
               }%
              }{%
               \false@sw
              }%
             }{%
              \false@sw
             }%
            }{%
             \false@sw
            }%
           }{%
            \false@sw
           }%
          }{%
           \false@sw
          }%
         }{%
          \false@sw
         }%
        }{%
         \false@sw
        }%
       }{%
        \false@sw
       }%
      }{%
       \false@sw
      }%
     }{%
      \false@sw
     }%
    }{%
     \false@sw
    }%
   }{%
    \false@sw
   }%
  }{%
   \false@sw
  }%
 }{%
  \false@sw
 }%
 {%
  \class@info{Patching LaTeX tabular.}%
 }{%
  \class@info{Unrecognized LaTeX tabular. Please update this document class! (Proceeding with fingers crossed.)}%
 }%
 \let\@array\@array@ltx
 \let\multicolumn\multicolumn@ltx
 \let\@tabular\@tabular@ltx
 \let\@tabarray\@tabarray@ltx
 \let\array\array@ltx
 \let\endarray\endarray@ltx
 \let\endtabular\endtabular@ltx
 \let\@mkpream\@mkpream@ltx
 \let\@addamp\@addamp@ltx
 \let\@arrayacol\@arrayacol@ltx
 \let\@tabacol\@tabacol@ltx
 \let\@arrayclassz\@arrayclassz@ltx
 \let\@tabclassiv\@tabclassiv@ltx
 \let\@arrayclassiv\@arrayclassiv@ltx
 \let\@tabclassz\@tabclassz@ltx
 \let\@classv\@classv@ltx
 \let\hline\hline@ltx
 \let\@tabularcr\@tabularcr@ltx
 \let\@xtabularcr\@xtabularcr@ltx
 \let\@xargarraycr\@xargarraycr@ltx
 \let\@yargarraycr\@yargarraycr@ltx
}%
\def\switch@array{%
 \let\@array@sw\@array@sw@LaTeX
 \@ifx{\@array\@array@array}{%
  \@ifx{\@tabular\@tabular@array}{%
   \@ifx{\@tabarray\@tabarray@array}{%
    \@ifx{\array\array@array}{%
     \@ifx{\endarray\endarray@array}{%
      \@ifx{\endtabular\endtabular@array}{%
       \@ifx{\@mkpream\@mkpream@array}{%
        \@ifx{\@classx\@classx@array}{%
         \@ifx{\insert@column\insert@column@array}{%
          \@ifx{\@arraycr\@arraycr@array}{%
           \@ifx{\@xarraycr\@xarraycr@array}{%
            \@ifx{\@xargarraycr\@xargarraycr@array}{%
             \@ifx{\@yargarraycr\@yargarraycr@array}{%
              \true@sw
         }{%
          \false@sw
         }%
         }{%
          \false@sw
         }%
         }{%
          \false@sw
         }%
         }{%
          \false@sw
         }%
         }{%
          \false@sw
         }%
        }{%
         \false@sw
        }%
       }{%
        \false@sw
       }%
      }{%
       \false@sw
      }%
     }{%
      \false@sw
     }%
    }{%
     \false@sw
    }%
   }{%
    \false@sw
   }%
  }{%
   \false@sw
  }%
 }{%
  \false@sw
 }{%
  \class@info{Patching array package.}%
 }{%
  \class@info{Unrecognized array package. Please update this document class! (Proceeding with fingers crossed.)}%
 }%
 \let\@array    \@array@array@new
 \let\@@array   \@array % Cosi fan tutti
 \let\@tabular  \@tabular@array@new
 \let\@tabarray \@tabarray@array@new
 \let\array     \array@array@new
 \let\endarray  \endarray@array@new
 \let\endtabular\endtabular@array@new
 \let\@mkpream  \@mkpream@array@new
 \let\@classx   \@classx@array@new
 \let\@arrayacol\@arrayacol@ltx
 \let\@tabacol  \@tabacol@ltx
 \let\insert@column\insert@column@array@new
 \expandafter\let\csname endtabular*\endcsname\endtabular % Cosi fan tutti
 \let\@arraycr  \@arraycr@new
 \let\@xarraycr \@xarraycr@new
 \let\@xargarraycr\@xargarraycr@new
 \let\@yargarraycr\@yargarraycr@new
}%
\def\@array@sw@LaTeX{\@ifx{\\\@tabularcr}}%
\def\@array@sw@array{\@ifx{\d@llarbegin\begingroup}}%
\def\@tabular@LaTeX{%
 \leavevmode
 \hbox\bgroup$%
  \let\@acol\@tabacol
  \let\@classz\@tabclassz
  \let\@classiv\@tabclassiv
  \let\\\@tabularcr
  \@tabarray
}%
\def\@tabular@ltx{%
  \let\@acoll\@tabacoll
  \let\@acolr\@tabacolr
  \let\@acol\@tabacol
  \let\@classz\@tabclassz
  \let\@classiv\@tabclassiv
  \let\\\@tabularcr
  \@tabarray
}%
\def\@tabular@array{%
 \leavevmode
 \hbox\bgroup$%
  \col@sep\tabcolsep
  \let\d@llarbegin\begingroup
  \let\d@llarend\endgroup
  \@tabarray
}%
\def\@tabular@array@new{%
  \let\@acoll\@tabacoll
  \let\@acolr\@tabacolr
  \let\@acol\@tabacol
  \let\col@sep\@undefined
  \let\d@llarbegin\begingroup
  \let\d@llarend\endgroup
  \@tabarray
}%
\def\@tabarray@LaTeX{%
 \m@th\@ifnextchar[\@array{\@array[c]}%
}%
\def\@tabarray@ltx{%
 \m@th\@ifnextchar[\@array{\expandafter\@array\expandafter[\array@default]}%
}%
\def\@tabarray@array{%
 \@ifnextchar[{\@@array}{\@@array[c]}%
}%
\def\@tabarray@array@new{%
 \@ifnextchar[{\@@array}{\expandafter\@@array\expandafter[\array@default]}%
}%
\newcount\intertabularlinepenalty
\intertabularlinepenalty=100
\newcount\@tbpen
\appdef\samepage{\intertabularlinepenalty\@M}%
\def\@tabularcr@LaTeX{{\ifnum 0=`}\fi \@ifstar \@xtabularcr \@xtabularcr}%
\def\@tabularcr@ltx{{\ifnum 0=`}\fi \@ifstar {\global \@tbpen \@M \@xtabularcr }{\global \@tbpen \intertabularlinepenalty \@xtabularcr }}%
\def\@xtabularcr@LaTeX{\@ifnextchar [\@argtabularcr {\ifnum 0=`{\fi }\cr }}%
\def\@xtabularcr@ltx{\@ifnextchar [\@argtabularcr {\ifnum 0=`{\fi }\cr \noalign {\penalty \@tbpen }}}%
\def\@xargarraycr@LaTeX#1{\@tempdima #1\advance \@tempdima \dp \@arstrutbox \vrule \@height \z@ \@depth \@tempdima \@width \z@ \cr}%
\def\@xargarraycr@ltx#1{\@tempdima #1\advance \@tempdima \dp \@arstrutbox \vrule \@height \z@ \@depth \@tempdima \@width \z@ \cr \noalign {\penalty \@tbpen }}%
\def\@yargarraycr@LaTeX#1{\cr \noalign {\vskip #1}}%
\def\@yargarraycr@ltx#1{\cr \noalign {\penalty \@tbpen \vskip #1}}%
\def\@arraycr@array{%
 \relax
 \iffalse{\fi\ifnum 0=`}\fi
 \@ifstar \@xarraycr \@xarraycr
}%
\def\@arraycr@new{%
 \relax
 \iffalse{\fi\ifnum 0=`}\fi
 \@ifstar {\global \@tbpen \@M \@xarraycr }{\global \@tbpen \intertabularlinepenalty \@xarraycr }%
}%
\def\@xarraycr@array{%
 \@ifnextchar [%]
 \@argarraycr {\ifnum 0=`{}\fi\cr}%
}%
\def\@xarraycr@new{%
 \@ifnextchar [%]
 \@argarraycr {\ifnum 0=`{}\fi\cr \noalign {\penalty \@tbpen }}%
}%
\def\@xargarraycr@array#1{%
 \unskip
 \@tempdima #1\advance\@tempdima \dp\@arstrutbox
 \vrule \@depth\@tempdima \@width\z@
 \cr
}%
\def\@xargarraycr@new#1{%
 \unskip
 \@tempdima #1\advance\@tempdima \dp\@arstrutbox
 \vrule \@depth\@tempdima \@width\z@
 \cr
 \noalign {\penalty \@tbpen }%
}%
\def\@yargarraycr@array#1{%
 \cr
 \noalign{\vskip #1}%
}%
\def\@yargarraycr@new#1{%
 \cr
 \noalign{\penalty \@tbpen \vskip #1}%
}%
\def\array@LaTeX{%
 \let\@acol\@arrayacol
 \let\@classz\@arrayclassz
 \let\@classiv\@arrayclassiv
 \let\\\@arraycr
 \let\@halignto\@empty
 \@tabarray
}%
\def\array@ltx{%
 \@ifmmode{}{\@badmath$}%
 \let\@acoll\@arrayacol
 \let\@acolr\@arrayacol
 \let\@acol\@arrayacol
 \let\@classz\@arrayclassz
 \let\@classiv\@arrayclassiv
 \let\\\@arraycr
 \let\@halignto\@empty
 \@tabarray
}%
\def\array@array{%
 \col@sep\arraycolsep
 \def\d@llarbegin{$}\let\d@llarend\d@llarbegin\gdef\@halignto{}%
 \@tabarray
}
\def\array@array@new{%
 \@ifmmode{}{\@badmath$}%
 \let\@acoll\@arrayacol
 \let\@acolr\@arrayacol
 \let\@acol\@arrayacol
\let\col@sep\@undefined
 \def\d@llarbegin{$}%
 \let\d@llarend\d@llarbegin
 \gdef\@halignto{}%
 \@tabarray
}%
\def\@array@LaTeX[#1]#2{%
  \if #1t\vtop \else \if#1b\vbox \else \vcenter \fi\fi
  \bgroup
  \setbox\@arstrutbox\hbox{%
    \vrule \@height\arraystretch\ht\strutbox
           \@depth\arraystretch \dp\strutbox
           \@width\z@}%
  \@mkpream{#2}%
  \edef\@preamble{%
    \ialign \noexpand\@halignto
      \bgroup \@arstrut \@preamble \tabskip\z@skip \cr}%
  \let\@startpbox\@@startpbox \let\@endpbox\@@endpbox
  \let\tabularnewline\\%
    \let\par\@empty
    \let\@sharp##%
    \set@typeset@protect
    \lineskip\z@skip\baselineskip\z@skip
    \ifhmode \@preamerr\z@ \@@par\fi
    \@preamble
}%
\def\@array@ltx[#1]#2{%
 \@nameuse{@array@align@#1}%
  \set@arstrutbox
  \@mkpream{#2}%
  \prepdef\@preamble{%
    \tabskip\tabmid@skip
    \@arstrut
  }%
  \appdef\@preamble{%
    \tabskip\tabright@skip
    \cr
    \array@row@pre
  }%
  \let\tabularnewline\\%
  \let\par\@empty
  \let\@sharp##%
  \set@typeset@protect
  \lineskip\z@skip\baselineskip\z@skip
  \tabskip\tableft@skip\relax
  \ifhmode \@preamerr\z@ \@@par\fi
  \everycr{}%
  \expandafter\halign\expandafter\@halignto\expandafter\bgroup\@preamble
}%
\def\set@arstrutbox{%
  \setbox\@arstrutbox\hbox{%
    \vrule \@height\arraystretch\ht\strutbox
           \@depth\arraystretch \dp\strutbox
           \@width\z@
  }%
}%
\def\@array@array[#1]#2{%
  \@tempdima \ht \strutbox
  \advance \@tempdima by\extrarowheight
  \setbox \@arstrutbox \hbox{\vrule
             \@height \arraystretch \@tempdima
             \@depth \arraystretch \dp \strutbox
             \@width \z@}%
  \begingroup
  \@mkpream{#2}%
  \xdef\@preamble{\noexpand \ialign \@halignto
                  \bgroup \@arstrut \@preamble
                          \tabskip \z@ \cr}%
  \endgroup
  \@arrayleft
  \if #1t\vtop \else \if#1b\vbox \else \vcenter \fi \fi
  \bgroup
  \let \@sharp ##\let \protect \relax
  \lineskip \z@
  \baselineskip \z@
  \m@th
  \let\\\@arraycr \let\tabularnewline\\\let\par\@empty \@preamble
}%
\def\@array@array@new[#1]#2{%
  \@tempdima\ht\strutbox
  \advance\@tempdima by\extrarowheight
  \setbox\@arstrutbox\hbox{%
   \vrule \@height\arraystretch\@tempdima
          \@depth \arraystretch\dp\strutbox
          \@width \z@
  }%
  \begingroup
   \@mkpream{#2}%
   \xdef\@preamble{\@preamble}%
  \endgroup
  \prepdef\@preamble{%
   \tabskip\tabmid@skip
    \@arstrut
  }%
  \appdef\@preamble{%
   \tabskip\tabright@skip
   \cr
   \array@row@pre
  }%
  \@arrayleft
  \@nameuse{@array@align@#1}%
  \m@th
  \let\\\@arraycr
  \let\tabularnewline\\%
  \let\par\@empty
  \let\@sharp##%
  \set@typeset@protect
  \lineskip\z@\baselineskip\z@
  \tabskip\tableft@skip
  \everycr{}%
  \expandafter\halign\expandafter\@halignto\expandafter\bgroup\@preamble
}%
\def\endarray@LaTeX{%
 \crcr\egroup\egroup
}%
\def\endarray@ltx{%
 \crcr\array@row@pst\egroup\egroup
}%
\def\endarray@array{%
 \crcr \egroup \egroup \@arrayright \gdef\@preamble{}%
}%
\def\endarray@array@new{%
 \crcr\array@row@pst\egroup\egroup % Same as \endarray@ltx
 \@arrayright
 \global\let\@preamble\@empty
}%
\def\endtabular@LaTeX{%
 \crcr\egroup\egroup $\egroup
}%
\def\endtabular@ltx{%
 \endarray
}%
\def\endtabular@array{%
 \endarray $\egroup
}%
\def\endtabular@array@new{%
 \endarray
}%
\@namedef{endtabular*}{\endtabular}%
\long\def\multicolumn@LaTeX#1#2#3{%
 \multispan{#1}\begingroup
  \@mkpream{#2}%
  \def\@sharp{#3}\set@typeset@protect
  \let\@startpbox\@@startpbox\let\@endpbox\@@endpbox
  \@arstrut \@preamble\hbox{}\endgroup\ignorespaces
}%
\long\def\multicolumn@ltx#1#2#3{%
 \multispan{#1}%
 \begingroup
  \@mkpream{#2}%
  \def\@sharp{#3}%
  \set@typeset@protect
 %\let\@startpbox\@@startpbox\let\@endpbox\@@endpbox
  \@arstrut
  \@preamble
  \hbox{}%
 \endgroup
 \ignorespaces
}%
\def\@array@align@t{\leavevmode\vtop\bgroup}%
\def\@array@align@b{\leavevmode\vbox\bgroup}%
\def\@array@align@c{\leavevmode\@ifmmode{\vcenter\bgroup}{$\vcenter\bgroup\aftergroup$\aftergroup\relax}}%
\def\@array@align@v{%
 \@ifmmode{%
  \@badmath
  \vcenter\bgroup
 }{%
  \@ifinner{%
   $\vcenter\bgroup\aftergroup$
  }{%
   \@@par\bgroup
  }%
 }%
}%
\def\array@default{c}%
\def\array@row@rst{%
 \let\@array@align@v\@array@align@c
}%
\def\array@row@pre{}%
\def\array@row@pst{}%
\newcommand\toprule{\tab@rule{\column@font}{\column@fil}{\frstrut}}%
\newcommand\colrule{\unskip\lrstrut\\\tab@rule{\body@font}{}{\frstrut}}%
\newcommand\botrule{\unskip\lrstrut\\\noalign{\hline@rule}{}}%
\def\hline@LaTeX{%
 \noalign{\ifnum0=`}\fi\hrule \@height \arrayrulewidth \futurelet
   \reserved@a\@xhline
}%
\def\hline@ltx{%
 \noalign{%
  \ifnum0=`}\fi
  \hline@rule
  \futurelet\reserved@a\@xhline
 % \noalign ended in \@xhline
}%
\def\@xhline@unneeded{%
 \say\reserved@a
 \ifx\reserved@a\hline
  \vskip\doublerulesep
  \vskip-\arrayrulewidth
 \fi
 \ifnum0=`{\fi}%
}%
\def\tab@rule#1#2#3{%
 \crcr
 \noalign{%
  \hline@rule
  \gdef\@arstrut@hook{%
   \global\let\@arstrut@hook\@empty
   #3%
  }%
  \gdef\cell@font{#1}%
  \gdef\cell@fil{#2}%
 }%
}%
\def\column@font{}%
\def\column@fil{}%
\def\body@font{}%
\def\cell@font{}%
\def\frstrut{}%
\def\lrstrut{}%
\def\@arstrut@hline{%
 \relax
 \@ifmmode{\copy}{\unhcopy}\@arstrutbox@hline
 \@arstrut@hook
}%
\let\@arstrut@org\@arstrut
\def\@arstrut@hook{%
 \global\let\@arstrut\@arstrut@org
}%
\newbox\@arstrutbox@hline
\appdef\set@arstrutbox{%
  \setbox\@arstrutbox@hline\hbox{%
    \setbox\z@\hbox{$0^{0}_{}$}%
    \dimen@\ht\z@\advance\dimen@\@arstrut@hline@clnc
    \@ifdim{\dimen@<\arraystretch\ht\strutbox}{\dimen@=\arraystretch\ht\strutbox}{}%
    \vrule \@height\dimen@
           \@depth\arraystretch \dp\strutbox
           \@width\z@
  }%
}%
\def\hline@rule{%
 \hrule \@height \arrayrulewidth
 \global\let\@arstrut\@arstrut@hline
}%
\def\@arstrut@hline@clnc{2\p@}% % Klootch: magic number
\def\tableft@skip{\z@skip}%
\def\tabmid@skip{\z@skip}%\@flushglue
\def\tabright@skip{\z@skip}%
\def\tableftsep{\tabcolsep}%
\def\tabmidsep{\tabcolsep}%
\def\tabrightsep{\tabcolsep}%
\def\cell@fil{}%
\def\pbox@hook{}%
\appdef\@arstrut{\@arstrut@hook}%
\let\@arstrut@hook\@empty
\def\@addtopreamble{\appdef\@preamble}%
\def\@mkpream@LaTeX#1{%
  \@firstamptrue\@lastchclass6
  \let\@preamble\@empty
  \let\protect\@unexpandable@protect
  \let\@sharp\relax
  \let\@startpbox\relax\let\@endpbox\relax
  \@expast{#1}%
  \expandafter\@tfor \expandafter
    \@nextchar \expandafter:\expandafter=\reserved@a\do
       {\@testpach\@nextchar
    \ifcase \@chclass \@classz \or \@classi \or \@classii \or \@classiii
      \or \@classiv \or\@classv \fi\@lastchclass\@chclass}%
  \ifcase \@lastchclass \@acol
      \or \or \@preamerr \@ne\or \@preamerr \tw@\or \or \@acol \fi
}%
\def\@mkpream@ltx#1{%
 \@firstamptrue
 \@lastchclass6
 \let\@preamble\@empty
 \let\protect\@unexpandable@protect
 \let\@sharp\relax
 \@expast{#1}%
 \expandafter\@tfor\expandafter\@nextchar\expandafter:\expandafter=\reserved@a
 \do{%
  \expandafter\@testpach\expandafter{\@nextchar}%
  \ifcase\@chclass
   \@classz
  \or
   \@classi
  \or
   \@classii
  \or
   \@classiii
  \or
   \@classiv
  \or
   \@classv
  \fi
  \@lastchclass\@chclass
 }%
 \ifcase\@lastchclass
  \@acolr % right-hand column
 \or
 \or
  \@preamerr\@ne
 \or
  \@preamerr\tw@
 \or
 \or
  \@acolr % right-hand column
 \fi
}%
\def\insert@column@array{%
   \the@toks \the \@tempcnta
   \ignorespaces \@sharp \unskip
   \the@toks \the \count@ \relax
}%
\def\insert@column@array@new{%
 \the@toks\the\@tempcnta
 \array@row@rst\cell@font
 \ignorespaces\@sharp\unskip
 \the@toks\the\count@
 \relax
}%
\def\@mkpream@relax{%
 \let\tableftsep\relax
 \let\tabmidsep\relax
 \let\tabrightsep\relax
 \let\array@row@rst\relax
 \let\cell@font\relax
 \let\@startpbox\relax
}%
\def\@mkpream@array#1{%
   \gdef\@preamble{}\@lastchclass 4 \@firstamptrue
   \let\@sharp\relax \let\@startpbox\relax \let\@endpbox\relax
   \@temptokena{#1}\@tempswatrue
   \@whilesw\if@tempswa\fi{\@tempswafalse\the\NC@list}%
   \count@\m@ne
   \let\the@toks\relax
   \prepnext@tok
   \expandafter \@tfor \expandafter \@nextchar
    \expandafter :\expandafter =\the\@temptokena \do
   {\@testpach
   \ifcase \@chclass \@classz \or \@classi \or \@classii
     \or \save@decl \or \or \@classv \or \@classvi
     \or \@classvii \or \@classviii
     \or \@classx
     \or \@classx \fi
   \@lastchclass\@chclass}%
   \ifcase\@lastchclass
   \@acol \or
   \or
   \@acol \or
   \@preamerr \thr@@ \or
   \@preamerr \tw@ \@addtopreamble\@sharp \or
   \or
   \else  \@preamerr \@ne \fi
   \def\the@toks{\the\toks}%
}%
\def\@mkpream@array@new#1{%
 \gdef\@preamble{}%
 \@lastchclass\f@ur
 \@firstamptrue
 \let\@sharp\relax
 \@mkpream@relax
 \@temptokena{#1}\@tempswatrue
 \@whilesw\if@tempswa\fi{\@tempswafalse\the\NC@list}%
 \count@\m@ne
 \let\the@toks\relax
 \prepnext@tok
 \expandafter\@tfor\expandafter\@nextchar\expandafter:\expandafter=\the\@temptokena
 \do{%
  \@testpach
  \ifcase\@chclass
   \@classz
  \or
   \@classi
  \or
   \@classii
  \or
   \save@decl
  \or
  \or
   \@classv
  \or
   \@classvi
  \or
   \@classvii
  \or
   \@classviii
  \or
   \@classx
  \or
   \@classx
  \fi
  \@lastchclass\@chclass
 }%
 \ifcase\@lastchclass
  \@acolr % right-hand column
 \or
 \or
  \@acolr % right-hand column
 \or
  \@preamerr\thr@@
 \or
  \@preamerr\tw@\@addtopreamble\@sharp
 \or
 \or
 \else
  \@preamerr\@ne
 \fi
 \def\the@toks{\the\toks}%
}%
\def\@addamp@LaTeX{%
  \if@firstamp\@firstampfalse\else\edef\@preamble{\@preamble &}\fi
}%
\def\@addamp@ltx{%
 \if@firstamp\@firstampfalse\else\@addtopreamble{&}\fi
}%
\def\@arrayacol@LaTeX{%
 \edef\@preamble{\@preamble \hskip \arraycolsep}%
}%
\def\@arrayacol@ltx{%
 \@addtopreamble{\hskip\arraycolsep}%
}%
\def\@tabacoll{%
 \@addtopreamble{\hskip\tableftsep\relax}%
}%
\def\@tabacol@LaTeX{%
 \edef\@preamble{\@preamble \hskip \tabcolsep}%
}%
\def\@tabacol@ltx{%
 \@addtopreamble{\hskip\tabmidsep\relax}%
}%
\def\@tabacolr{%
 \@addtopreamble{\hskip\tabrightsep\relax}%
}%
\def\@arrayclassz@LaTeX{%
 \ifcase \@lastchclass \@acolampacol \or \@ampacol \or
   \or \or \@addamp \or
   \@acolampacol \or \@firstampfalse \@acol \fi
 \edef\@preamble{\@preamble
  \ifcase \@chnum
     \hfil$\relax\@sharp$\hfil \or $\relax\@sharp$\hfil
    \or \hfil$\relax\@sharp$\fi}%
}%
\def\@arrayclassz@ltx{%
 \ifcase\@lastchclass
  \@acolampacol
 \or
  \@ampacol
 \or
 \or
 \or
  \@addamp
 \or
  \@acolampacol
 \or
  \@firstampfalse\@acoll
 \fi
 \ifcase\@chnum
  \@addtopreamble{%
   \hfil\array@row@rst$\relax\@sharp$\hfil
  }%
 \or
  \@addtopreamble{%
   \array@row@rst$\relax\@sharp$\hfil
  }%
 \or
  \@addtopreamble{%
   \hfil\array@row@rst$\relax\@sharp$%
  }%
 \fi
}%
\def\@tabclassz@LaTeX{%
  \ifcase\@lastchclass
    \@acolampacol
  \or
    \@ampacol
  \or
  \or
  \or
    \@addamp
  \or
    \@acolampacol
  \or
    \@firstampfalse\@acol
  \fi
  \edef\@preamble{%
    \@preamble{%
      \ifcase\@chnum
        \hfil\ignorespaces\@sharp\unskip\hfil
      \or
        \hskip1sp\ignorespaces\@sharp\unskip\hfil
      \or
        \hfil\hskip1sp\ignorespaces\@sharp\unskip
      \fi}}%
}%
\def\@tabclassz@ltx{%
 \ifcase\@lastchclass
  \@acolampacol
 \or
  \@ampacol
 \or
 \or
 \or
  \@addamp
 \or
  \@acolampacol
 \or
  \@firstampfalse\@acoll
 \fi
 \ifcase\@chnum
  \@addtopreamble{%
   {\hfil\array@row@rst\cell@font\ignorespaces\@sharp\unskip\hfil}%
  }%
 \or
  \@addtopreamble{%
   {\cell@fil\hskip1sp\array@row@rst\cell@font\ignorespaces\@sharp\unskip\hfil}%
  }%
 \or
  \@addtopreamble{%
   {\hfil\hskip1sp\array@row@rst\cell@font\ignorespaces\@sharp\unskip\cell@fil}%
  }%
 \fi
}%
\def\@tabclassiv@LaTeX{%
 \@addtopreamble\@nextchar
}%
\def\@tabclassiv@ltx{%
 \expandafter\@addtopreamble\expandafter{\@nextchar}%
}%
\def\@arrayclassiv@LaTeX{%
 \@addtopreamble{$\@nextchar$}%
}%
\def\@arrayclassiv@ltx{%
 \expandafter\@addtopreamble\expandafter{\expandafter$\@nextchar$}%
}%
\def\@classv@LaTeX{%
 \@addtopreamble{\@startpbox{\@nextchar}\ignorespaces
 \@sharp\@endpbox}%
}%
\def\@classv@ltx{%
 \expandafter\@addtopreamble
 \expandafter{%
 \expandafter \@startpbox
 \expandafter {\@nextchar}%
 \pbox@hook\array@row@rst\cell@font\ignorespaces\@sharp\@endpbox
 }%
}%
\def\@classx@array{%
  \ifcase \@lastchclass
  \@acolampacol \or
  \@addamp \@acol \or
  \@acolampacol \or
  \or
  \@acol \@firstampfalse \or
  \@addamp
  \fi
}%
\def\@classx@array@new{%
 \ifcase \@lastchclass
  \@acolampacol
 \or
  \@addamp \@acol
 \or
  \@acolampacol
 \or
 \or
  \@firstampfalse\@acoll
 \or
  \@addamp
 \fi
}%
\def\@xbitor@LaTeX #1{\@tempcntb \count#1
   \ifnum \@tempcnta =\z@
   \else
     \divide\@tempcntb\@tempcnta
     \ifodd\@tempcntb \@testtrue\fi
   \fi}%
\def\@xbitor@ltx#1{%
 \@tempcntb\count#1%
 \@ifnum{\@tempcnta=\z@}{}{%
  \divide\@tempcntb\@tempcnta
  \@ifodd\@tempcntb{\@testtrue}{}%
 }%
}%
\@ifx{\@xbitor\@xbitor@LaTeX}{%
  \class@info{Repairing broken LaTeX \string\@xbitor}%
}{%
  \class@info{Unrecognized LaTeX \string\@xbitor. Please update this document class! (Proceeding with fingers crossed.)}%
}%
\let\@xbitor\@xbitor@ltx
\newcommand*\@gobble@opt@one[2][]{}%
\def\@starttoc#1{%
  \begingroup
    \toc@pre
    \makeatletter
    \@input{\jobname.#1}%
    \if@filesw
      \expandafter\newwrite\csname tf@#1\endcsname
      \immediate\openout \csname tf@#1\endcsname \jobname.#1\relax
    \fi
    \@nobreakfalse
    \toc@post
  \endgroup
}%
\def\toc@pre{}%
\def\toc@post{}%
\def\toc@@font{}%{\footnotesize\rmfamily}%
\def\@dotsep{\z@}%{5.5pt}%
\let\tocdim@section       \leftmargini
\let\tocdim@subsection    \leftmarginii
\let\tocdim@subsubsection \leftmarginiii
\let\tocdim@paragraph     \leftmarginiv
\let\tocdim@appendix      \leftmarginv
\let\tocdim@pagenum       \leftmarginvi
\def\toc@pre@auto{%
  \toc@@font
  \@tempdima\z@
  \toc@setindent\@tempdima{section}%
  \toc@setindent\@tempdima{subsection}%
  \toc@setindent\@tempdima{subsubsection}%
  \toc@setindent\@tempdima{paragraph}%
  \toc@letdimen{appendix}%
  \toc@letdimen{pagenum}%
}%
\def\toc@post@auto{%
  \if@filesw
   \begingroup
    \toc@writedimen{section}%
    \toc@writedimen{subsection}%
    \toc@writedimen{subsubsection}%
    \toc@writedimen{paragraph}%
    \toc@writedimen{appendix}%
    \toc@writedimen{pagenum}%
   \endgroup
  \fi
}%
\def\toc@setindent#1#2{%
 \csname tocdim@#2\endcsname\tocdim@min\relax
 \@ifundefined{tocmax@#2}{\@namedef{tocmax@#2}{\z@}}{}%
 \advance#1\@nameuse{tocmax@#2}\relax
 \expandafter\edef\csname tocleft@#2\endcsname{\the#1}%
}%
\def\toc@letdimen#1{%
 \csname tocdim@#1\endcsname\tocdim@min\relax
 \@ifundefined{tocmax@#1}{\@namedef{tocmax@#1}{\z@}}{}%
 \expandafter\let\csname tocleft@#1\expandafter\endcsname\csname tocmax@#1\endcsname
}%
\def\toc@writedimen#1{%
 \immediate\write\@auxout{%
  \gdef\expandafter\string\csname tocmax@#1\endcsname{%
   \expandafter\the\csname tocdim@#1\endcsname
  }%
 }%
}%
\def\l@@sections#1#2#3#4{%
 % #1 - superior section
 % #2 - this section
 % #3 - content, including possible \numberline
 % #4 - page number
 \begingroup
  \everypar{}%
  \set@tocdim@pagenum{#4}%
  \global\@tempdima\csname tocdim@#2\endcsname
  \leftskip\csname tocleft@#2\endcsname\relax
  \dimen@\csname tocleft@#1\endcsname\relax
  \parindent-\leftskip\advance\parindent\dimen@
  \rightskip\tocleft@pagenum plus 1fil\relax
  \skip@\parfillskip\parfillskip\z@
  \let\numberline\numberline@@sections
  \@nameuse{l@f@#2}%
  \ignorespaces#3\unskip\nobreak\hskip\skip@
  \hb@xt@\rightskip{\hfil\unhbox\@tempboxa}\hskip-\rightskip\hskip\z@skip
  \par
  \expandafter\aftergroup\csname tocdim@#2\endcsname\expandafter
 \endgroup\the\@tempdima\relax
}%
\def\set@tocdim@pagenum#1{%
 \setbox\@tempboxa\hbox{\ignorespaces#1}%
 \@ifdim{\tocdim@pagenum<\wd\z@}{\global\tocdim@pagenum\wd\z@}{}%
}%
\def\numberline@@sections#1{%
 \leavevmode\hb@xt@-\parindent{%
  \hfil
  \@if@empty{#1}{}{%
   \setbox\z@\hbox{#1.\kern\@dotsep}%
   \@ifdim{\@tempdima<\wd\z@}{\global\@tempdima\wd\z@}{}%
   \unhbox\z@
  }%
 }%
 \ignorespaces
}%
\def\tocdim@min{\z@}%
\def\list#1#2{%
  \ifnum \@listdepth >5\relax
    \@toodeep
  \else
    \global\advance\@listdepth\@ne
  \fi
  \rightmargin\z@
  \listparindent\z@
  \itemindent\z@
  \csname @list\romannumeral\the\@listdepth\endcsname
  \def\@itemlabel{#1}%
  \let\makelabel\@mklab
  \@nmbrlistfalse
  #2\relax
  \@trivlist
  \parskip\parsep
  \set@listindent
  \ignorespaces
}%
\def\set@listindent@parshape{%
 \parindent\listparindent
 \advance\@totalleftmargin\leftmargin
 \advance\linewidth-\rightmargin
 \advance\linewidth-\leftmargin
 \parshape\@ne\@totalleftmargin\linewidth
}%
\def\set@listindent@{%
 \parindent\listparindent
 \advance\@totalleftmargin\leftmargin
 \advance\rightskip\rightmargin
 \advance\leftskip\@totalleftmargin
}%
\let\set@listindent\set@listindent@parshape
\typeout{%
 ltxgrid: portions licensed from W. E. Baxter (web@superscript.com)%
}%
\newcounter{linecount}
\def\lineloop#1{%
 \loop
 \ifnum\c@linecount<#1\relax
  \global\advance\c@linecount\@ne
  \par
  \hb@xt@\hsize{%
   \ifnum\c@linecount<100 0\fi\ifnum\c@linecount<10 0\fi\number\c@linecount
   \vrule depth2.5\p@
   \leaders\hrule\hfil
  }%
  \penalty\interlinepenalty
 \repeat
}%
\let\@@mark\mark
\let\@@topmark\topmark
\let\@@firstmark\firstmark
\let\@@botmark\botmark
\let\@@splitfirstmark\splitfirstmark
\let\@@splitbotmark\splitbotmark
\def\@themark{{}{}{}{}}%
\def\nul@mark{{}{}{}{}\@@nul}%
\def\set@mark@netw@#1#2#3#4#5#6#7{\gdef#1{{#6}{#7}{#4}{#5}}\do@mark}%
\def\set@marktw@#1#2#3#4#5#6{\gdef#1{{#2}{#6}{#4}{#5}}\do@mark}%
\def\set@markthr@@#1#2#3#4#5#6{\gdef#1{{#2}{#3}{#6}{#5}}\do@mark}%
\def\get@mark@@ne#1#2#3#4#5\@@nul{#1}%
\def\get@mark@tw@#1#2#3#4#5\@@nul{#2}%
\def\get@mark@thr@@#1#2#3#4#5\@@nul{#3}%
\def\get@mark@f@ur#1#2#3#4#5\@@nul{#4}%
\def\mark@netw@{\expandafter\set@mark@netw@\expandafter\@themark\@themark}%
\def\marktw@{\expandafter\set@marktw@\expandafter\@themark\@themark}%
\def\markthr@@{\expandafter\set@markthr@@\expandafter\@themark\@themark}%
\def\do@mark{\do@@mark\@themark\nobreak@mark}%
\def\do@@mark#1{%
 \begingroup
  \let@mark
  \@@mark{#1}%
 \endgroup
}%
\def\let@mark{%
 \let\protect\@unexpandable@protect
 \let\label\relax
 \let\index\relax
 \let\glossary\relax
}%
\def\nobreak@mark{%
 \@if@sw\if@nobreak\fi{\@ifvmode{\nobreak}{}}{}%
}%
\def\mark@envir{\markthr@@}%
\def\bot@envir{%
 \expandafter\expandafter
 \expandafter\get@mark@thr@@
 \expandafter\@@botmark
             \nul@mark
}%
\def\markboth{\mark@netw@}%
\def\markright{\marktw@}%
\def\leftmark{%
 \expandafter\expandafter
 \expandafter\get@mark@@ne
 \expandafter\saved@@botmark
             \nul@mark
}%
\def\rightmark{%
 \expandafter\expandafter
 \expandafter\get@mark@tw@
 \expandafter\saved@@firstmark
             \nul@mark
}%
\let\primitive@output\output
\long\def\@tempa#1\@@nil{#1}%
\toks@
\expandafter\expandafter
\expandafter{%
\expandafter \@tempa
             \the\output
             \@@nil
             }%
\newtoks\output
\output\expandafter{\the\toks@}%
\primitive@output{\dispatch@output}%
\def\dispatch@output{%
 \let\par\@@par
 \expandafter\let\expandafter\@tempa\csname output@\the\outputpenalty\endcsname
 \outputdebug@sw{%
  \saythe\badness
  \saythe\outputpenalty
  \saythe\holdinginserts
  \say\thepagegrid
  \saythe\pagegrid@col
  \saythe\pagegrid@cur
 %\say\bot@envir
  \saythe\insertpenalties
 %\say\@@topmark
 %\say\saved@@topmark
 %\say\@@firstmark
 %\say\saved@@firstmark
  \say\@@botmark
 %\say\saved@@botmark
  \saythe\pagegoal
  \saythe\pagetotal
  \saythe{\badness\@cclv}%
  \expandafter\@ifx\expandafter{\csname output@-\the\execute@message@pen\endcsname\@tempa}{%
   \say\@message@saved
  }{%
   \expandafter\say\csname output@\the\outputpenalty\endcsname
  }%
  \say\@toplist
  \say\@botlist
  \say\@dbltoplist
  \say\@deferlist
  {\tracingall\scrollmode
   \showbox\@cclv
   \showbox\@cclv@saved
   \showbox\pagesofar
   \showbox\footbox
   \showbox\footins@saved
   \showbox\footins
   \showlists
  }%
 }{}%
 \@ifnotrelax\@tempa{\@tempa}{\the\output}%
}%
\@ifxundefined{\outputdebug@sw}{%
 \@booleanfalse\outputdebug@sw
}{}%
\output={\toggle@insert\output@holding\output@moving}%
\def\output@holding{%
\csname output@init@\bot@envir\endcsname
\@if@exceed@pagegoal{\unvcopy\@cclv}{%
  \setbox\z@\vbox{\unvcopy\@cclv}%
\outputdebug@sw{{\tracingall\scrollmode\showbox\z@}}{}%
\dimen@\ht\@cclv\advance\dimen@-\ht\z@
\dead@cycle@repair\dimen@
}{%
\dead@cycle
}%
}%
\def\@if@exceed@pagegoal#1{%
 \begingroup
  \setbox\z@\vbox{#1}%
 \dimen@\ht\z@\advance\dimen@\dp\z@
  \outputdebug@sw{\saythe\dimen@}{}%
  \@ifdim{\dimen@>\pagegoal}{%
  \setbox\z@\vbox{\@@mark{}\unvbox\z@}%
   \splittopskip\topskip
   \splitmaxdepth\maxdepth
   \vbadness\@M
   \vfuzz\maxdimen
   \setbox\tw@\vsplit\z@ to\pagegoal
  \outputdebug@sw{{\tracingall\scrollmode\showbox\tw@\showbox\z@}}{}%
  \setbox\tw@\vbox{\unvbox\tw@}%
  \@ifdim{\ht\tw@=\z@}{%
   \ltxgrid@info{Found overly large chunk while preparing to move insertions. Attempting repairs}%
    \aftergroup\true@sw
  }{%
    \aftergroup\false@sw
  }%
 }{%
 \aftergroup\false@sw
 }%
\endgroup
}%
%% \item
%%  Put down the same interrupts as for the non-trivial case above.
\def\output@moving{%
 \set@top@firstmark
 \@ifnum{\outputpenalty=\do@newpage@pen}{%
  \setbox\@cclv\vbox{%
   \unvbox\@cclv
   \setbox\z@\lastbox
   \@ifdim{\ht\z@=\ht\@protection@box}{\box\lastbox}{\unskip}%
  }%
 }{}%
 \@cclv@nontrivial@sw{%
  \csname output@prep@\bot@envir \endcsname
  \@makecol\csname output@column@\thepagegrid\endcsname
  \protect@penalty\do@startcolumn@pen
  \clearpage@sw{%
   \protect@penalty\do@endpage@pen
  }{}%
  \csname output@post@\bot@envir \endcsname
 }{%
  {\setbox\z@\box\@cclv}%
 }%
 \set@colroom
 \global\@mparbottom\z@
 \global\@textfloatsheight\z@ %FIXME: this legacy LaTeX variable is set, but never queried!
}%
\def\@cclv@nontrivial@sw{%
\@ifx@empty\@toplist{%
\@ifx@empty\@botlist{%
\@ifvoid\footins{%
 \@ifvoid\@cclv{%
  \false@sw
 }{%
\setbox\z@\vbox{\unvcopy\@cclv}%
\@ifdim{\ht\z@=\topskip}{%
\setbox\z@\vbox{%
\unvbox\z@
\setbox\z@\lastbox\dimen@\lastskip\unskip
\@ifdim{\ht\z@=\ht\@protection@box}{%
\advance\dimen@\ht\z@
\@ifdim{\dimen@=\topskip}{%
\aftergroup\true@sw
}{%
\aftergroup\false@sw
}%
}{%
\aftergroup\false@sw
}%
}%
{%
\false@sw % Normal for \clearpage
}{%
\true@sw
}%
}{%
\@ifdim{\ht\z@=\z@}{%
\ltxgrid@info{Found trivial column. Discarding it}%
\outputdebug@sw{{\tracingall\scrollmode\showbox\@cclv}}{}%
\false@sw
}{%
\true@sw
}%
}%
 }%
}{%
\true@sw
}%
}{%
\true@sw
}%
}{%
\true@sw
}%
}%
\def\protect@penalty#1{\protection@box\penalty-#1\relax}%
\newbox\@protection@box
\setbox\@protection@box\vbox to1986sp{\vfil}%
\def\protection@box{\nointerlineskip\copy\@protection@box}%
\def\dead@cycle@repair#1{%
\expandafter\do@@mark
\expandafter{%
\@@botmark
}%
\unvbox\@cclv
\nointerlineskip
\vbox to#1{\vss}%
\@ifnum{\outputpenalty<\@M}{\penalty\outputpenalty}{}%
}%
\def\dead@cycle@repair@protected#1{%
\expandafter\do@@mark
\expandafter{%
\@@botmark
}%
\begingroup
  \unvbox\@cclv
  \setbox\z@\lastbox % Remove protection box
 \nointerlineskip
  \advance#1-\ht\@protection@box
 \vbox to#1{\vss}%
  \protection@box % Reinsert protection box
 \@ifnum{\outputpenalty<\@M}{\penalty\outputpenalty}{}%
\endgroup
}%
\def\dead@cycle{%
 \expandafter\do@@mark
 \expandafter{%
              \@@botmark
             }%
 \unvbox\@cclv
 \@ifnum{\outputpenalty<\@M}{\penalty\outputpenalty}{}%
}%
\def\output@init@document{%
 \@ifvoid\footbox{}{%
  \global\advance\vsize-\ht\footbox
  \global\advance\vsize-\dp\footbox
 }%
}%
\def\output@prep@document{%
 \@ifvoid\footbox{}{%
  \setbox\footins\vbox{\unvbox\footbox\unvbox\footins}%
 }%
}%
\def\output@post@document{}%
\let\@opcol\@undefined
\def\@makecol{%
 \setbox\@outputbox\vbox{%
  \boxmaxdepth\@maxdepth
  \@tempdima\dp\@cclv
  \unvbox\@cclv
  \vskip-\@tempdima
 }%
 \xdef\@freelist{\@freelist\@midlist}\global\let\@midlist\@empty
 \@combinefloats
 \@combineinserts\@outputbox\footins
  \set@adj@colht\dimen@
  \count@\vbadness
  \vbadness\@M
  \setbox\@outputbox\vbox to\dimen@{%
   \@texttop
   \dimen@\dp\@outputbox
   \unvbox\@outputbox
   \vskip-\dimen@
   \@textbottom
  }%
  \vbadness\count@
 \global\maxdepth\@maxdepth
}%
\let\@makespecialcolbox\@undefined
\def\@combineinserts#1#2{%
 \setbox#1\vbox{%
  \unvbox#1%
  \vbox{%
   \@ifvoid#2{}{%
    \vskip\skip\footins
    \color@begingroup
    \normalcolor
    \footnoterule
    \nointerlineskip
    \box#2%
    \color@endgroup
   }{}%
  }%
 }%
}%
\appdef\@floatplacement{%
 \global\@fpmin\@fpmin
}%
\mathchardef\pagebreak@pen=\@M
\expandafter\let\csname output@-\the\pagebreak@pen\endcsname\relax
\mathchardef\do@startcolumn@pen=10005
\@namedef{output@-\the\do@startcolumn@pen}{\do@startcolumn}%
\def\do@startcolumn{%
 \setbox\@cclv\vbox{\unvbox\@cclv\setbox\z@\lastbox\unskip}%
 \clearpage@sw{\@clearfloatplacement}{\@floatplacement}%
 \set@colroom
 \@booleanfalse\pfloat@avail@sw
 \begingroup
  \@colht\@colroom
  \@booleanfalse\float@avail@sw
  \@tryfcolumn\test@colfloat
  \float@avail@sw{\aftergroup\@booleantrue\aftergroup\pfloat@avail@sw}{}%
 \endgroup
 \fcolmade@sw{%
  \setbox\@cclv\vbox{\unvbox\@outputbox\unvbox\@cclv}%
  \outputpenalty-\pagebreak@pen % ask for a return visit, this time with insertions and all.
  \dead@cycle
 }{%
  \begingroup
   \let\@elt\@scolelt
   \let\reserved@b\@deferlist\global\let\@deferlist\@empty\reserved@b
  \endgroup
  \clearpage@sw{%
   \outputpenalty\@M
  }{%
   \outputpenalty\do@newpage@pen
  }%
  \dead@cycle
 }%
 \check@deferlist@stuck\do@startcolumn
 \set@vsize
}%
\def\@scolelt#1{\def\@currbox{#1}\@addtonextcol}%
\def\test@colfloat#1{%
 \csname @floatselect@sw@\thepagegrid\endcsname#1{}{\@testtrue}%
 \@if@sw\if@test\fi{}{\aftergroup\@booleantrue\aftergroup\float@avail@sw}%
}%
\def\@addtonextcol{%
 \begingroup
  \@insertfalse
  \@setfloattypecounts
  \csname @floatselect@sw@\thepagegrid\endcsname\@currbox{%
   \@ifnum{\@fpstype=8 }{}{%
     \@ifnum{\@fpstype=24 }{}{%
       \@flsettextmin
       \@reqcolroom \ht\@currbox
       \advance \@reqcolroom \@textmin
       \advance \@reqcolroom \vsize % take into account split insertions
       \advance \@reqcolroom -\pagegoal
       \@ifdim{\@colroom>\@reqcolroom}{%
         \@flsetnum \@colnum
         \@ifnum{\@colnum>\z@}{%
            \@bitor\@currtype\@deferlist
            \@if@sw\if@test\fi{}{%
              \@addtotoporbot
            }%
         }{}%
       }{}%
     }%
   }%
  }{}%
  \@if@sw\if@insert\fi{}{%
    \@cons\@deferlist\@currbox
  }%
 \endgroup
}%
\mathchardef\do@startpage@pen=10006
\@namedef{output@-\the\do@startpage@pen}{\do@startpage}%
\def\do@startpage{%
 \setbox\@cclv\vbox{\unvbox\@cclv\setbox\z@\lastbox\unskip}%
 \clearpage@sw{\@clearfloatplacement}{\@dblfloatplacement}%
 \set@colht
 \@booleanfalse\pfloat@avail@sw
 \begingroup
  \@booleanfalse\float@avail@sw
  \@tryfcolumn\test@dblfloat
  \float@avail@sw{\aftergroup\@booleantrue\aftergroup\pfloat@avail@sw}{}%
 \endgroup
 \fcolmade@sw{%
  \global\setbox\pagesofar\vbox{\unvbox\pagesofar\unvbox\@outputbox}%
  \@combinepage
  \@combinedblfloats
  \@outputpage
  \global\pagegrid@cur\@ne
  \protect@penalty\do@startpage@pen
 }{%
  \begingroup
   \@booleanfalse\float@avail@sw
   \let\@elt\@sdblcolelt
   \let\reserved@b\@deferlist\global\let\@deferlist\@empty\reserved@b
  \endgroup
  \@ifdim{\@colht=\textheight}{% No luck...
   \pfloat@avail@sw{% ...but a float *was* available!
    \forcefloats@sw{%
     \ltxgrid@warn{Forced dequeueing of floats stalled}%
    }{%
     \ltxgrid@warn{Dequeueing of floats stalled}%
    }%
   }{}%
  }{}%
  \outputpenalty\@M
  \dead@cycle
 }%
 \check@deferlist@stuck\do@startpage
 \set@colht
}%
\def\@sdblcolelt#1{\def\@currbox{#1}\@addtodblcol}%
\def\test@dblfloat#1{%
 \@if@notdblfloat{#1}{\@testtrue}{}%
 \@if@sw\if@test\fi{}{\aftergroup\@booleantrue\aftergroup\float@avail@sw}%
}%
\def\@if@notdblfloat#1{\@ifdim{\wd#1<\textwidth}}%
\@booleanfalse\forcefloats@sw
\def\@addtodblcol{%
 \begingroup
  \@if@notdblfloat{\@currbox}{%
   \false@sw
  }{%
   \@setfloattypecounts
   \@getfpsbit \tw@
   \@bitor \@currtype \@deferlist
   \@if@sw\if@test\fi{%
    \false@sw
   }{%
    \@ifodd\@tempcnta{%
     \aftergroup\@booleantrue\aftergroup\float@avail@sw
     \@flsetnum \@dbltopnum
     \@ifnum{\@dbltopnum>\z@}{%
       \@ifdim{\@dbltoproom>\ht\@currbox}{%
        \true@sw
       }{%
        \@ifnum{\@fpstype<\sixt@@n}{%
         \begingroup
          \advance \@dbltoproom \@textmin
          \@ifdim{\@dbltoproom>\ht\@currbox}{%
           \endgroup\true@sw
          }{%
           \endgroup\false@sw
          }%
        }{%
         \false@sw
        }%
       }%
     }{%
      \false@sw
     }%
    }{%
     \false@sw
    }%
   }%
  }%
  {%
   \@tempdima -\ht\@currbox
   \advance\@tempdima
    -\@ifx{\@dbltoplist\@empty}{\dbltextfloatsep}{\dblfloatsep}%
   \global \advance \@dbltoproom \@tempdima
   \global \advance \@colht \@tempdima
   \global \advance \@dbltopnum \m@ne
   \@cons \@dbltoplist \@currbox
  }{%
   \@cons \@deferlist \@currbox
  }%
 \endgroup
}%
\def\@tryfcolumn#1{%
  \global\@booleanfalse\fcolmade@sw
  \@ifx@empty\@deferlist{}{%
    \global\let\@trylist\@deferlist
    \global\let\@failedlist\@empty
    \begingroup
      \dimen@\vsize\advance\dimen@-\pagegoal\@ifdim{\dimen@>\z@}{%
       \advance\@fpmin-\dimen@
      }{}%
      \def\@elt{\@xtryfc#1}\@trylist
    \endgroup
    \fcolmade@sw{%
      \global\setbox\@outputbox\vbox{\vskip \@fptop}%
      \let \@elt \@wtryfc \@flsucceed
      \global\setbox\@outputbox\vbox{\unvbox\@outputbox
        \unskip \vskip \@fpbot
      }%
      \let \@elt \relax
      \xdef\@deferlist{\@failedlist\@flfail}%
      \xdef\@freelist{\@freelist\@flsucceed}%
    }{}%
  }%
}%
\def\@wtryfc #1{%
  \global\setbox\@outputbox\vbox{\unvbox\@outputbox
    \box #1\vskip\@fpsep
  }%
}%
\def\@xtryfc#1#2{%
  \@next\reserved@a\@trylist{}{}% trim \@trylist. Ugly!
  \@currtype \count #2%
  \divide\@currtype\@xxxii\multiply\@currtype\@xxxii
  \@bitor \@currtype \@failedlist
  \@testfp #2%
  #1#2%
  \@ifdim{\ht #2>\@colht   }{\@testtrue}{}%
  \@if@sw\if@test\fi{%
   \@cons\@failedlist #2%
  }{%
   \begingroup
     \gdef\@flsucceed{\@elt #2}%
     \global\let\@flfail\@empty
     \@tempdima\ht #2%
     \def \@elt {\@ztryfc#1}\@trylist
     \@ifdim{\@tempdima >\@fpmin}{%
       \global\@booleantrue\fcolmade@sw
     }{%
       \@cons\@failedlist #2%
     }%
   \endgroup
   \fcolmade@sw{%
     \let \@elt \@gobble
   }{}%
  }%
}%
\def\@ztryfc #1#2{%
  \@tempcnta \count#2%
  \divide\@tempcnta\@xxxii\multiply\@tempcnta\@xxxii
  \@bitor \@tempcnta {\@failedlist \@flfail}%
  \@testfp #2%
  #1#2%
  \@tempdimb\@tempdima
  \advance\@tempdimb \ht#2\advance\@tempdimb\@fpsep
  \@ifdim{\@tempdimb >\@colht}{%
    \@testtrue
  }{}%
  \@if@sw\if@test\fi{%
    \@cons\@flfail #2%
  }{%
    \@cons\@flsucceed #2%
    \@tempdima\@tempdimb
  }%
}%
\def\newpage@prep{%
  \if@noskipsec
    \ifx \@nodocument\relax
      \leavevmode
      \global \@noskipsecfalse
    \fi
  \fi
  \if@inlabel
    \leavevmode
    \global \@inlabelfalse
  \fi
  \if@nobreak \@nobreakfalse \everypar{}\fi
  \par
}%
\def \newpage {%
 \newpage@prep
 \do@output@MVL{%
  \vfil
  \penalty-\pagebreak@pen
 }%
}%
\def\clearpage{%
 \newpage@prep
 \do@output@MVL{%
  \vfil
  \penalty-\pagebreak@pen
  \global\@booleantrue\clearpage@sw
  \protect@penalty\do@startcolumn@pen
  \protect@penalty\do@endpage@pen
 }%
 \do@output@MVL{%
  \global\@booleanfalse\clearpage@sw
 }%
}%
\def\cleardoublepage{%
 \clearpage
 \@if@sw\if@twoside\fi{%
  \@ifodd\c@page{}{%
   \null\clearpage
  }%
 }{}%
}%
\@booleanfalse\clearpage@sw
\mathchardef\do@endpage@pen=10007
\@namedef{output@-\the\do@endpage@pen}{%
 \csname end@column@\thepagegrid\endcsname
}%
\mathchardef\do@newpage@pen=10001
\expandafter\let\csname output@-\the\do@newpage@pen\endcsname\relax
\def\@clearfloatplacement{%
 \global\@topnum     \maxdimen % \c@topnumber
 \global\@toproom    \maxdimen % \topfraction\@colht
 \global\@botnum     \maxdimen % \c@bottomnumber
 \global\@botroom    \maxdimen % \bottomfraction\@colht
 \global\@colnum     \maxdimen % \c@totalnumber
 \global\@dbltopnum  \maxdimen % \c@dbltopnumber
 \global\@dbltoproom \maxdimen % \dbltopfraction\@colht
 \global\@textmin    \z@       % \@colht\advance \@textmin -\@dbltoproom
 \global\@fpmin      \z@       % \dblfloatpagefraction\textheight
 \let\@testfp\@gobble
 \appdef\@setfloattypecounts{\@fpstype16\advance\@fpstype\m@ne}%
}%
\let\@doclearpage\@undefined
\let\@makefcolumn\@undefined
\def\clr@top@firstmark{%
 \global\let\saved@@topmark\@undefined
 \global\let\saved@@firstmark\@empty
 \global\let\saved@@botmark\@empty
}%
\clr@top@firstmark
\def\set@top@firstmark{%
 \@ifxundefined\saved@@topmark{\expandafter\gdef\expandafter\saved@@topmark\expandafter{\@@topmark}}{}%
 \@if@empty\saved@@firstmark{\expandafter\gdef\expandafter\saved@@firstmark\expandafter{\@@firstmark}}{}%
 \@if@empty\@@botmark{}{\expandafter\gdef\expandafter\saved@@botmark\expandafter{\@@botmark}}%
}%
\appdef\@outputpage{%
 \clr@top@firstmark
}%
\def\@float#1{%
 \@ifnextchar[{%}]{%Brace-matching klootch
  \@yfloat\width@float{#1}%
 }{%
  \@ifxundefined@cs{fps@#1}{%
   \edef\reserved@a{\noexpand\@yfloat\noexpand\width@float{#1}[\csname fps@\endcsname]}\reserved@a
  }{%
   \edef\reserved@a{\noexpand\@yfloat\noexpand\width@float{#1}[\csname fps@#1\endcsname]}\reserved@a
  }%
 }%
}%
\def\@dblfloat#1{%
 \@ifnum{\pagegrid@col=\@ne}{%
  \@float{#1}%
 }{%
  \@ifnextchar[{%}]{%Brace-matching klootch
   \@yfloat\widthd@float{#1}%
  }{%
   \@ifxundefined@cs{fpsd@#1}{%
    \edef\reserved@a{\noexpand\@yfloat\noexpand\widthd@float{#1}[\csname fpsd@\endcsname]}\reserved@a
   }{%
    \edef\reserved@a{\noexpand\@yfloat\noexpand\widthd@float{#1}[\csname fpsd@#1\endcsname]}\reserved@a
   }%
  }%
 }%
}%
\def\@yfloat#1#2[#3]{%
 \@xfloat{#2}[#3]%
 \hsize#1\linewidth\hsize
 \minipagefootnote@init
}%
\def\fps@{tbp}%
\def\fpsd@{tp}%
\def\width@float{\columnwidth}%
\def\widthd@float{\textwidth}%
\def\end@float{%
 \end@@float{%
  \check@currbox@count
 }%
}%
\def\end@dblfloat{%
 \@ifnum{\pagegrid@col=\@ne}{%
  \end@float
 }{%
  \end@@float{%
   \@boxfpsbit\@currbox{1}\@ifodd\@tempcnta{\global\advance\count\@currbox\m@ne}{}%
   \@boxfpsbit\@currbox{4}\@ifodd\@tempcnta{\global\advance\count\@currbox-4\relax}{}%
   \global\wd\@currbox\textwidth % Klootch
   \check@currbox@count
  }%
 }%
}%
\def\end@@float#1{%
 \minipagefootnote@here
 \@endfloatbox
 #1%
 \@ifnum{\@floatpenalty <\z@}{%
  \@largefloatcheck
  \@cons\@currlist\@currbox
  \@ifnum{\@floatpenalty <-\@Mii}{%
   \do@output@cclv{\@add@float}%
  }{%
   \vadjust{\do@output@cclv{\@add@float}}%
   \@Esphack
  }%
 }{}%
}%
\def\check@currbox@count{%
 \@ifnum{\count\@currbox>\z@}{%
\count@\count\@currbox\divide\count@\sixt@@n\multiply\count@\sixt@@n
\@tempcnta\count\@currbox\advance\@tempcnta-\count@
\@ifnum{\@tempcnta=\z@}{%
   \ltxgrid@warn{Float cannot be placed}%
}{}%
 }{%
  % Is a \marginpar
 }%
}%
\providecommand\minipagefootnote@init{}%
\providecommand\minipagefootnote@here{}%
\let\@specialoutput\@undefined
\def\@add@float{%
 \@pageht\ht\@cclv\@pagedp\dp\@cclv
 \unvbox\@cclv
 \@next\@currbox\@currlist{%
  \csname @floatselect@sw@\thepagegrid\endcsname\@currbox{%
   \@ifnum{\count\@currbox>\z@}{%
    \advance \@pageht \@pagedp
    \advance \@pageht \vsize \advance \@pageht -\pagegoal % do not assume \holdinginserts is cleared!
    \@addtocurcol % Commit an h float
  }{%
   \@addmarginpar
  }%
  }{%
   \@resethfps
   \@cons\@deferlist\@currbox
  }%
 }{\@latexbug}%
 \@ifnum{\outputpenalty<\z@}{%
  \@if@sw\if@nobreak\fi{%
   \nobreak
  }{%
   \addpenalty \interlinepenalty
  }%
 }{}%
 \set@vsize
}%
\let\@reinserts\@undefined
\def \@addtocurcol {%
   \@insertfalse
   \@setfloattypecounts
   \ifnum \@fpstype=8
   \else
     \ifnum \@fpstype=24
     \else
       \@flsettextmin
       \advance \@textmin \@textfloatsheight
       \@reqcolroom \@pageht
       \ifdim \@textmin>\@reqcolroom
         \@reqcolroom \@textmin
       \fi
       \advance \@reqcolroom \ht\@currbox
       \ifdim \@colroom>\@reqcolroom
         \@flsetnum \@colnum
         \ifnum \@colnum>\z@
           \@bitor\@currtype\@deferlist
           \if@test
           \else
             \@bitor\@currtype\@botlist
             \if@test
               \@addtobot
             \else
               \ifodd \count\@currbox
                 \advance \@reqcolroom \intextsep
                 \ifdim \@colroom>\@reqcolroom
                   \global \advance \@colnum \m@ne
                   \global \advance \@textfloatsheight \ht\@currbox
                   \global \advance \@textfloatsheight 2\intextsep
                   \@cons \@midlist \@currbox
                   \if@nobreak
                     \nobreak
                     \@nobreakfalse
                     \everypar{}%
                   \else
                     \addpenalty \interlinepenalty
                   \fi
                   \vskip \intextsep
                   \unvbox\@currbox %AO
                   \penalty\interlinepenalty
                   \vskip\intextsep
                   \ifnum\outputpenalty <-\@Mii \vskip -\parskip\fi
                   \outputpenalty \z@
                   \@inserttrue
                 \fi
               \fi
               \if@insert
               \else
                 \@addtotoporbot
               \fi
             \fi
           \fi
         \fi
       \fi
     \fi
   \fi
   \if@insert
   \else
     \@resethfps
     \@cons\@deferlist\@currbox
   \fi
}%
\@twocolumnfalse
\let\@twocolumntrue\@twocolumnfalse
\def\@addmarginpar{%
 \@next\@marbox\@currlist{%
  \@cons\@freelist\@marbox\@cons\@freelist\@currbox
 }\@latexbug
 \setbox\@marbox\hb@xt@\columnwidth{%
  \csname @addmarginpar@\thepagegrid\endcsname{%
   \hskip-\marginparsep\hskip-\marginparwidth
   \box\@currbox
  }{%
   \hskip\columnwidth\hskip\marginparsep
   \box\@marbox
  }%
  \hss
 }%
 \setbox\z@\box\@currbox
    \@tempdima\@mparbottom
    \advance\@tempdima -\@pageht
    \advance\@tempdima\ht\@marbox
 \@ifdim{\@tempdima >\z@}{%
   \@latex@warning@no@line {Marginpar on page \thepage\space moved}%
 }{%
   \@tempdima\z@
 }%
    \global\@mparbottom\@pageht
    \global\advance\@mparbottom\@tempdima
    \global\advance\@mparbottom\dp\@marbox
    \global\advance\@mparbottom\marginparpush
    \advance\@tempdima -\ht\@marbox
    \global\setbox \@marbox
                   \vbox {\vskip \@tempdima
                          \box \@marbox}%
    \global \ht\@marbox \z@
    \global \dp\@marbox \z@
    \kern -\@pagedp
    \nointerlineskip
  \box\@marbox
    \nointerlineskip
    \hbox{\vrule \@height\z@ \@width\z@ \@depth\@pagedp}%
}%
\newenvironment{turnpage}{%
 \def\width@float{\textheight}%
 \def\widthd@float{\textheight}%
 \appdef\@endfloatbox{%
  \@ifxundefined\@currbox{%
   \ltxgrid@warn{Cannot rotate! Not a float}%
  }{%
   \setbox\@currbox\vbox to\textwidth{\vfil\unvbox\@currbox\vfil}%
   \global\setbox\@currbox\vbox{\rotatebox{90}{\box\@currbox}}%
  }%
 }%
}{%
}%
\def\rotatebox@dummy#1#2{%
 \ltxgrid@warn{You must load the graphics or graphicx package in order to use the turnpage environment}%
 #2%
}%
\AtBeginDocument{%
 \@ifxundefined\rotatebox{\let\rotatebox\rotatebox@dummy}{}%
}%
\@namedef{output@-1073741824}{%"40000000
 \deadcycles\z@
 \setbox\z@\box\@cclv
}%
\mathchardef\save@column@pen=10016
\@namedef{output@-\the\save@column@pen}{\save@column}%
\let \@cclv@saved \@holdpg
\let \@holdpg \@undefined
\def\save@column{%
 \@ifvoid\@cclv@saved{%
\set@top@firstmark
\global\@topmark@saved\expandafter{\@@topmark}%
 }{}%
\global\setbox\@cclv@saved\vbox{%
  \@ifvoid\@cclv@saved{}{%
 \unvbox\@cclv@saved
 \marry@baselines
}%
\unvbox\@cclv
\lose@breaks
\setbox\z@\lastbox
}%
}%
\newtoks\@topmark@saved
\def\prep@cclv{%
 \setbox\z@\box\@cclv
 \setbox\@cclv\box\@cclv@saved
 \vbadness\@M
}%
\mathchardef\save@column@insert@pen=10017
\@namedef{output@-\the\save@column@insert@pen}{\toggle@insert\savecolumn@holding\savecolumn@moving}%
\def\savecolumn@holding{%
\@if@exceed@pagegoal{\unvcopy\@cclv\setbox\z@\lastbox}{%
  \setbox\z@\vbox{\unvcopy\@cclv\setbox\z@\lastbox}%
  \outputdebug@sw{{\tracingall\scrollmode\showbox\z@}}{}%
\dimen@\ht\@cclv\advance\dimen@-\ht\z@
\dead@cycle@repair@protected\dimen@
}{%
\dead@cycle
}%
}%
\def\savecolumn@moving{%
 \@cclv@nontrivial@sw{%
  \save@column
 }{%
  {\setbox\z@\box\@cclv}%
 }%
\@ifvoid\footins{}{%
\outputdebug@sw{{\tracingall\scrollmode\showbox\footins}}{}%
\global\setbox\footins@saved\vbox{\unvbox\footins@saved\marry@baselines\unvbox\footins}%
\protect@penalty\save@column@insert@pen
}%
}%
\newbox\footins@saved
\mathchardef\save@message@pen=10018
\@namedef{output@-\the\save@message@pen}{\save@message}%
\def\save@message{%
 \setbox\z@\box\@cclv %FIXME: what if \box\@cclv is not empty?
 \toks@\expandafter{\@@firstmark}%
 \expandafter\gdef\expandafter\@message@saved\expandafter{\the\toks@}%
 \expandafter\do@@mark\expandafter{\the\@topmark@saved}%
}%
\gdef\@message@saved{}%
\mathchardef\execute@message@pen=10019
\@namedef{output@-\the\execute@message@pen}{\@message@saved}%
\def\execute@message{%
 \@execute@message\save@column@pen %Implicit #2
}%
\def\execute@message@insert#1{%
 \@execute@message\save@column@insert@pen{\setbox\footins\box\footins@saved#1}%
}%
\long\def\@execute@message#1#2{%
 \begingroup
  \dimen@\prevdepth\@ifdim{\dimen@<\z@}{\dimen@\z@}{}%
  \setbox\z@\vbox{%
   \protect@penalty#1%
   \protection@box
   \toks@{\prep@cclv#2}%
   \@@mark{\the\toks@}%
   \penalty-\save@message@pen
   \setbox\z@\null\dp\z@\dimen@\ht\z@-\dimen@
   \nointerlineskip\box\z@
   \penalty-\execute@message@pen
  }\unvbox\z@
 \endgroup
}%
\def\do@output@cclv{\execute@message}%
\def\do@output@MVL#1{%
 \@ifvmode{%
  \begingroup\execute@message{\unvbox\@cclv#1}\endgroup
 }{%
  \@ifhmode{%
   \vadjust{\execute@message{\unvbox\@cclv#1}}%
  }{%
   \@latexerr{\string\do@output@MVL\space cannot be executed in this mode!}\@eha
  }%
 }%
}%
\def\lose@breaks{%
 \loopwhile{%
  \count@\lastpenalty
  \@ifnum{\count@=\@M}{% 10000 is a TeX magic number!
   \unpenalty\true@sw
  }{%
   \false@sw
  }%
 }%
}%
\def\removestuff{\do@output@MVL{\unskip\unpenalty}}%
\def\removephantombox{%
 \vadjust{%
  \execute@message{%
   \unvbox\@cclv
   \setbox\z@\lastbox
   \unskip
   \unskip
   \unpenalty
   \penalty\predisplaypenalty
   \vskip\abovedisplayskip
  }%
 }%
}%
\def\addstuff#1#2{\edef\@tempa{\noexpand\do@output@MVL{\noexpand\@addstuff{#1}{#2}}}\@tempa}%
\def\@addstuff#1#2{%
 \skip@\lastskip\unskip
 \count@\lastpenalty\unpenalty
 \@if@empty{#1}{}{\penalty#1\relax}%
 \@ifnum{\count@=\z@}{}{\penalty\count@}%
 \vskip\skip@
 \@if@empty{#2}{}{\vskip#2\relax}%
}%
\def\replacestuff#1#2{\edef\@tempa{\noexpand\do@output@MVL{\noexpand\@replacestuff{#1}{#2}}}\@tempa}%
\def\@replacestuff#1#2{%
 \skip@\lastskip\unskip
 \count@\lastpenalty\unpenalty
 \@if@empty{#1}{}{%
 \@ifnum{\count@>\@M}{}{%
   \@ifnum{\count@=\z@}{\count@=#1\relax}{%
    \@ifnum{\count@<#1\relax}{}{%
     \count@=#1\relax
    }%
   }%
 }%
 }%
 \@ifnum{\count@=\z@}{}{\penalty\count@}%
 \@if@empty{#2}{}{%
  \@tempskipa#2\relax
  \@ifdim{\z@>\@tempskipa}{%
   \advance\skip@-\@tempskipa
  }{%
   \@ifdim{\skip@>\@tempskipa}{}{%
    \skip@\@tempskipa
   }%
  }%
 }%
 \vskip\skip@
}%
\def\move@insertions{\global\holdinginserts\z@}%
\def\hold@insertions{\global\holdinginserts\@ne}%
\hold@insertions
\def\move@insert@sw{\@ifnum{\holdinginserts=\z@}}%
\def\toggle@insert#1#2{%
 \@ifnum{\holdinginserts=\z@}{\hold@insertions#2}{\move@insertions#1}%
}%
\def\do@columngrid#1#2{%
 \par
 \expandafter\let\expandafter\@tempa\csname open@column@#1\endcsname
 \@ifx{\relax\@tempa}{%
  \ltxgrid@warn{Unknown page grid #1. No action taken}%
 }{%
  \do@output@MVL{\start@column{#1}{#2}}%
 }%
}%
\def\start@column#1#2{%
 \def\@tempa{#1}\@ifx{\@tempa\thepagegrid}{%
  \ltxgrid@info{Already in page grid \thepagegrid. No action taken}%
 }{%
  \expandafter\execute@message@insert
  \expandafter{%
               \csname shut@column@\thepagegrid\expandafter\endcsname
               \csname open@column@#1\endcsname{#2}%
               \set@vsize
             }%
 }%
}%
\def\thepagegrid{one}%
\newbox\pagesofar
\newbox\footbox
\newcommand\onecolumngrid{\do@columngrid{one}{\@ne}}%
\let\onecolumn\@undefined
\def\open@column@one#1{%
 \unvbox\pagesofar
 \gdef\thepagegrid{one}%
 \global\pagegrid@col#1%
 \global\pagegrid@cur\@ne
 \set@colht
 \set@column@hsize\pagegrid@col
}%
\def\shut@column@one{%
 \@makecol
 \global\setbox\pagesofar\vbox{\unvbox\@outputbox\recover@footins}%
 \set@colht
}%
\def\float@column@one{%
 \@makecol
 \@outputpage
}%
\def\end@column@one{%
 \unvbox\@cclv\setbox\z@\lastbox
 \protect@penalty\do@newpage@pen
}%
\def\output@column@one{%
 \@outputpage
}%
\def\@addmarginpar@one{%
 \@if@sw\if@mparswitch\fi{%
  \@ifodd\c@page{\false@sw}{\true@sw}%
 }{\false@sw}{%
  \@if@sw\if@reversemargin\fi{\false@sw}{\true@sw}%
 }{%
  \@if@sw\if@reversemargin\fi{\true@sw}{\false@sw}%
 }%
}%
\def\@floatselect@sw@one#1{\true@sw}%
\def\onecolumngrid@push{%
 \do@output@MVL{%
  \@ifnum{\pagegrid@col=\@ne}{%
   \global\let\restorecolumngrid\@empty
  }{%
   \xdef\restorecolumngrid{%
    \noexpand\start@column{\thepagegrid}{\the\pagegrid@col}%
   }%
   \start@column{one}{\@ne}%
  }%
 }%
}%
\def\onecolumngrid@pop{%
 \do@output@MVL{\restorecolumngrid}%
}%
\newcommand\twocolumngrid{\do@columngrid{mlt}{\tw@}}%
\let\twocolumn\@undefined
\let\@topnewpage\@undefined
\def\open@column@mlt#1{%
 \gdef\thepagegrid{mlt}%
 \global\pagegrid@col#1%
 \global\pagegrid@cur\@ne
 \set@column@hsize\pagegrid@col
 \set@colht
}%
\def\shut@column@mlt{%
 \@cclv@nontrivial@sw{%
\@makecol
\@ifnum{\pagegrid@cur<\pagegrid@col}{%
\expandafter\global\expandafter\setbox\csname col@\the\pagegrid@cur\endcsname\box\@outputbox
\global\advance\pagegrid@cur\@ne
}{}%
 }{%
  {\setbox\z@\box\@cclv}%
 }%
\@ifnum{\pagegrid@cur>\@ne}{%
\csname balance@\the\pagegrid@col\endcsname
\grid@column{}%
\@combinepage
\@combinedblfloats
\global\setbox\pagesofar\box\@outputbox
 }{}%
 \set@colht
}%
\def\float@column@mlt{%
  \@combinepage
  \@combinedblfloats
  \@outputpage
  \global\pagegrid@cur\@ne
  \protect@penalty\do@startpage@pen
}%
\def\end@column@mlt{%
 \@ifx@empty\@toplist{%
  \@ifx@empty\@botlist{%
   \@ifx@empty\@dbltoplist{%
    \@ifx@empty\@deferlist{%
     \@ifnum{\pagegrid@cur=\@ne}{%
      \false@sw
     }{%
      \true@sw
     }%
    }{%
     \true@sw
    }%
   }{%
    \true@sw
   }%
  }{%
   \true@sw
  }%
 }{%
  \true@sw
 }%
 % true = kick out a column and try again
 {%
  \@cclv@nontrivial@sw{%
   \unvbox\@cclv\setbox\z@\lastbox
  }{%
   \unvbox\@cclv\setbox\z@\lastbox\unskip\null
  }%
  \protect@penalty\do@newpage@pen
  \protect@penalty\do@endpage@pen
 }{%
  \unvbox\@cclv\setbox\z@\lastbox
 }%
}%
\def\output@column@mlt{%
 \@ifnum{\pagegrid@cur<\pagegrid@col}{%
  \expandafter\global\expandafter\setbox\csname col@\the\pagegrid@cur\endcsname\box\@outputbox
  \global\advance\pagegrid@cur\@ne
 }{%
  \set@adj@colht\dimen@
  \grid@column{}%{\dimen@}%
  \@combinepage
  \@combinedblfloats
  \@outputpage
  \global\pagegrid@cur\@ne
  \protect@penalty\do@startpage@pen
 }%
}%
\let\@outputdblcol\@undefined
\def\@floatselect@sw@mlt#1{\@if@notdblfloat{#1}}%
\def\@addmarginpar@mlt{% emits a boolean
 \@ifnum{\pagegrid@cur=\@ne}%
}%
\let\pagegrid@cur\col@number
\let\col@number\@undefined
\newcount\pagegrid@col
\pagegrid@cur\@ne
\expandafter\let\csname col@\the\pagegrid@cur\endcsname\@leftcolumn
\let\@leftcolumn\@undefined
\pagegrid@col\tw@
\def\pagegrid@init{%
 \advance\pagegrid@cur\@ne
 \@ifnum{\pagegrid@cur<\pagegrid@col}{%
  \csname newbox\expandafter\endcsname\csname col@\the\pagegrid@cur\endcsname
  \pagegrid@init
 }{%
 }%
}%
\appdef\class@documenthook{%
 \pagegrid@init
}%
\def\grid@column#1{%
 \global\setbox\@outputbox\vbox{%
  \hb@xt@\textwidth{%
   \vrule\@height\z@\@width\z@\@if@empty{#1}{}{\@depth#1}%
   \pagegrid@cur\@ne
   \append@column
   \box@column\@outputbox
  }%
  \vskip\z@skip % FIXME: page depth!
 }%
}%
\def\append@column{%
 \@ifnum{\pagegrid@cur<\pagegrid@col}{%
  \expandafter\box@column\csname col@\the\pagegrid@cur\endcsname
  \hfil
  \vrule \@width\columnseprule
  \hfil
  \advance\pagegrid@cur\@ne
  \append@column
 }{%
 }%
}%
\def\box@column#1{%
 \raise\topskip
 \hb@xt@\columnwidth{%
  \dimen@\ht#1\@ifdim{\dimen@>\@colht}{\dimen@\@colht}{}%
  \count@\vbadness\vbadness\@M
  \dimen@ii\vfuzz\vfuzz\maxdimen
  \outputdebug@sw{\saythe\@colht\saythe\dimen@}{}%
  \vtop to\dimen@
  {\hrule\@height\z@
   \unvbox#1%
   \raggedcolumn@skip
  }%
  \vfuzz\dimen@ii
  \vbadness\count@
  \hss
 }%
}%
\def\marry@baselines{%
 \vskip\marry@skip\relax
}%
\gdef\marry@skip{\z@skip}%
\def\set@marry@skip{%
\begingroup
  \skip@\baselineskip\advance\skip@-\topskip
  \@ifdim{\skip@>\z@}{%
   \xdef\marry@skip{\the\skip@}%
  }{}%
 \endgroup
}%
\AtBeginDocument{%
 \@ifxundefined\raggedcolumn@sw{\@booleanfalse\raggedcolumn@sw}{}%
}%
\def\raggedcolumn@skip{%
 \vskip\z@\raggedcolumn@sw{\@plus.0001fil\@minus.0001fil}{}\relax
}%
\def\@combinepage{%
 \@ifvoid\pagesofar{}{%
  \setbox\@outputbox\vbox{%
   \unvbox\pagesofar
   \marry@baselines
   \unvbox\@outputbox
  }%
 }%
 \@ifvoid\footbox{}{%
  \setbox\@outputbox\vbox{%
   \unvbox\@outputbox
   \marry@baselines
   \unvbox\footbox
  }%
 }%
}%
\def\@combinedblfloats{%
 \@ifx@empty\@dbltoplist{}{%
  \setbox\@tempboxa\vbox{}%
  \let\@elt\@comdblflelt\@dbltoplist
  \let\@elt\relax\xdef\@freelist{\@freelist\@dbltoplist}%
  \global\let\@dbltoplist\@empty
  \setbox\@outputbox\vbox{%
   %\boxmaxdepth\maxdepth   %% probably not needed, CAR
   \unvbox\@tempboxa\unskip
   \@ifnum{\@dbltopnum>\m@ne}{\dblfigrule}{}%FIXME: how is \@dbltopnum maintained?
   \vskip\dbltextfloatsep
   \unvbox\@outputbox
  }%
 }%
}%
\def\set@column@hsize#1{%
 \pagegrid@col#1%
 \global\columnwidth\textwidth
 \global\advance\columnwidth\columnsep
 \global\divide\columnwidth\pagegrid@col
 \global\advance\columnwidth-\columnsep
 \global\hsize\columnwidth
 \global\linewidth\columnwidth
 \skip@\baselineskip\advance\skip@-\topskip
 \@ifnum{\pagegrid@col>\@ne}{\set@marry@skip}{}%
}%
\def\set@colht{%
 \set@adj@textheight\@colht
 \global\let\enlarge@colroom\@empty
 \set@colroom
}%
\def\set@adj@textheight#1{%
 #1\textheight
 \def\@elt{\adj@page#1}%
 \@booleantrue\firsttime@sw\@dbltoplist
 \let\@elt\relax
 \global#1#1\relax
 \outputdebug@sw{\saythe#1}{}%
}%
\def\set@colroom{%
 \set@adj@colht\@colroom
 \@if@empty\enlarge@colroom{}{%
  \global\advance\@colroom\enlarge@colroom\relax
 }%
 \outputdebug@sw{\saythe\@colroom}{}%
 \@ifdim{\@colroom>\topskip}{}{%
  \ltxgrid@info{Not enough room: \string\@colroom=\the\@colroom; increasing to \the\topskip}%
  \@colroom\topskip
 }%
 \global\@colroom\@colroom
 \set@vsize
}%
\def\set@vsize{%
 \global\vsize\@colroom
 \outputdebug@sw{\saythe\vsize}{}%
}%
\def\set@adj@colht#1{%
 #1\@colht
 \@ifvoid\pagesofar{}{%
  \advance#1-\ht\pagesofar\advance#1-\dp\pagesofar
 }%
 \@ifvoid\footbox{}{%
  \advance#1-\ht\footbox\advance#1-\dp\footbox
 }%
 \def\@elt{\adj@column#1}%
 \@booleantrue\firsttime@sw\@toplist
 \@booleantrue\firsttime@sw\@botlist
 \let\@elt\relax
 \outputdebug@sw{\saythe#1}{}%
}%
\def\adj@column#1#2{%
 \advance#1-\ht#2%
 \advance#1-\firsttime@sw{\textfloatsep\@booleanfalse\firsttime@sw}{\floatsep}%
}%
\def\adj@page#1#2{%
 \advance#1-\ht#2%
 \advance#1-\firsttime@sw{\dbltextfloatsep\@booleanfalse\firsttime@sw}{\dblfloatsep}%
}%
\appdef\@outputpage{%
 \set@colht          % FIXME: needed?
 \@floatplacement    % FIXME: needed?
 \@dblfloatplacement % FIXME: needed?
}%
\@namedef{balance@2}{%
 \expandafter\balance@two\csname col@1\endcsname\@outputbox
 % Avoid a bug by preventing a restore when leaving this group
 \global\setbox\csname col@1\endcsname\box\csname col@1\endcsname
 \@ifvoid\footbox{}{%
  \global\setbox\footbox\vbox{%
   \setbox\z@\box\@tempboxa
   \let\recover@footins\relax
   \balance@two\footbox\@tempboxa
   \hb@xt@\textwidth{\box\footbox\hfil\box\@tempboxa}%
  }%
 }%
}%
\def\balance@two#1#2{%
\outputdebug@sw{{\tracingall\scrollmode\showbox#1\showbox#2}}{}%
 \setbox\@ne\vbox{%
  \@ifvoid#1{}{%
   \unvcopy#1\recover@footins
   \@ifvoid#2{}{\marry@baselines}%
  }%
  \@ifvoid#2{}{%
   \unvcopy#2\recover@footins
  }%
 }%
 \dimen@\ht\@ne\divide\dimen@\tw@
 \dimen@i\dimen@
 \vbadness\@M
 \vfuzz\maxdimen
 \loopwhile{%
  \dimen@i=.5\dimen@i
  \outputdebug@sw{\saythe\dimen@\saythe\dimen@i\saythe\dimen@ii}{}%
  \setbox\z@\copy\@ne\setbox\tw@\vsplit\z@ to\dimen@
  \setbox\z@ \vbox{%
   \unvcopy\z@
   \setbox\z@\vbox{\unvbox\z@ \setbox\z@\lastbox\aftergroup\vskip\aftergroup-\expandafter}\the\dp\z@\relax
  }%
  \setbox\tw@\vbox{%
   \unvcopy\tw@
   \setbox\z@\vbox{\unvbox\tw@\setbox\z@\lastbox\aftergroup\vskip\aftergroup-\expandafter}\the\dp\z@\relax
  }%
  \dimen@ii\ht\tw@\advance\dimen@ii-\ht\z@
  \@ifdim{\dimen@i>.5\p@}{%
   \advance\dimen@\@ifdim{\dimen@ii<\z@}{}{-}\dimen@i
   \true@sw
  }{%
   \@ifdim{\dimen@ii<\z@}{%
    \advance\dimen@\tw@\dimen@i
    \true@sw
   }{%
    \false@sw
   }%
  }%
 }%
 \outputdebug@sw{\saythe\dimen@\saythe\dimen@i\saythe\dimen@ii}{}%
\@ifdim{\ht\z@=\z@}{%
\@ifdim{\ht\tw@=\z@}{%
\true@sw
}{%
\false@sw
}%
}{%
\true@sw
}%
{%
}{%
\ltxgrid@info{Unsatifactorily balanced columns: giving up}%
\setbox\tw@\box#1%
\setbox\z@ \box#2%
}%
 \setbox\tw@\vbox{\unvbox\tw@\vskip\z@skip}%
 \setbox\z@ \vbox{\unvbox\z@ \vskip\z@skip}%
 \set@colroom
\dimen@\ht\z@\@ifdim{\dimen@<\ht\tw@}{\dimen@\ht\tw@}{}%
\@ifdim{\dimen@>\@colroom}{\dimen@\@colroom}{}%
 \outputdebug@sw{\saythe{\ht\z@}\saythe{\ht\tw@}\saythe\@colroom\saythe\dimen@}{}%
\setbox#1\vbox to\dimen@{\unvbox\tw@\unskip\raggedcolumn@skip}%
\setbox#2\vbox to\dimen@{\unvbox\z@ \unskip\raggedcolumn@skip}%
\outputdebug@sw{{\tracingall\scrollmode\showbox#1\showbox#2}}{}%
}%
\def\recover@footins{%
 \skip\z@ \lastskip\unskip
 \skip\@ne\lastskip\unskip
 \setbox\z@\lastbox
 \@ifvbox\z@{%
\setbox\z@\vbox{%
\unvbox\z@
\setbox\z@\lastbox
\@ifvoid\z@{}{%
\global\setbox\footbox\vbox{%
 \unvbox\footbox
     \@ifvbox\z@{%
   \unvbox\z@
  }{%
   \box\z@
  }%
}%
}%
}%
 }{}%
 \outputdebug@sw{{\tracingall\scrollmode\showbox\footbox}}{}%
}%
\prepdef\@begindocumenthook{%
 \open@column@one\@ne
 \set@colht
 \@floatplacement
 \@dblfloatplacement
}%
\def\longtable@longtable{%
 \par
 \ifx\multicols\@undefined\else\ifnum\col@number>\@ne\@twocolumntrue\fi\fi
 \if@twocolumn\LT@err{longtable not in 1-column mode}\@ehc\fi
 \begingroup
 \@ifnextchar[\LT@array{\LT@array[x]}%
}%
\def\longtable@new{%
 \par
  \@ifnextchar[\LT@array{\LT@array[x]}%
}%
\def\endlongtable@longtable{%
  \crcr
  \noalign{%
    \let\LT@entry\LT@entry@chop
    \xdef\LT@save@row{\LT@save@row}}%
  \LT@echunk
  \LT@start
  \unvbox\z@
  \LT@get@widths
  \if@filesw
    {\let\LT@entry\LT@entry@write\immediate\write\@auxout{%
      \gdef\expandafter\noexpand
        \csname LT@\romannumeral\c@LT@tables\endcsname
          {\LT@save@row}}}%
  \fi
  \ifx\LT@save@row\LT@@save@row
  \else
    \LT@warn{Column \@width s have changed\MessageBreak
             in table \thetable}%
    \LT@final@warn
  \fi
  \endgraf\penalty -\LT@end@pen
  \endgroup
  \global\@mparbottom\z@
  \pagegoal\vsize
  \endgraf\penalty\z@\addvspace\LTpost
  \ifvoid\footins\else\insert\footins{}\fi
}%
\def\endlongtable@new{%
  \crcr
  \noalign{%
   \let\LT@entry\LT@entry@chop
   \xdef\LT@save@row{\LT@save@row}%
  }%
  \LT@echunk
  \LT@start
  \unvbox\z@
  \LT@get@widths
  \@if@sw\if@filesw\fi{%
   {%
    \let\LT@entry\LT@entry@write
    \immediate\write\@auxout{%
     \gdef\expandafter\noexpand\csname LT@\romannumeral\c@LT@tables\endcsname
     {\LT@save@row}%
    }%
   }%
  }{}%
  \@ifx\LT@save@row\LT@@save@row{}{%
   \LT@warn{%
    Column \@width s have changed\MessageBreak in table \thetable
   }\LT@final@warn
  }%
  \endgraf
  \nobreak
  \box\@ifvoid\LT@lastfoot{\LT@foot}{\LT@lastfoot}%
 \global\@mparbottom\z@
 \endgraf
 \LT@post
}%
\def\LT@start@longtable{%
 \let\LT@start\endgraf
 \endgraf
 \penalty\z@
 \vskip\LTpre
 \dimen@\pagetotal
 \advance\dimen@ \ht\ifvoid\LT@firsthead\LT@head\else\LT@firsthead\fi
 \advance\dimen@ \dp\ifvoid\LT@firsthead\LT@head\else\LT@firsthead\fi
 \advance\dimen@ \ht\LT@foot
 \dimen@ii\vfuzz\vfuzz\maxdimen
   \setbox\tw@\copy\z@
   \setbox\tw@\vsplit\tw@ to \ht\@arstrutbox
   \setbox\tw@\vbox{\unvbox\tw@}%
 \vfuzz\dimen@ii
 \advance\dimen@ \ht
       \ifdim\ht\@arstrutbox>\ht\tw@\@arstrutbox\else\tw@\fi
 \advance\dimen@\dp
       \ifdim\dp\@arstrutbox>\dp\tw@\@arstrutbox\else\tw@\fi
 \advance\dimen@ -\pagegoal
 \ifdim \dimen@>\z@\vfil\break\fi
     \global\@colroom\@colht
 \ifvoid\LT@foot\else
   \advance\vsize-\ht\LT@foot
   \global\advance\@colroom-\ht\LT@foot
   \dimen@\pagegoal\advance\dimen@-\ht\LT@foot\pagegoal\dimen@
   \maxdepth\z@
 \fi
 \ifvoid\LT@firsthead\copy\LT@head\else\box\LT@firsthead\fi
 \output{\LT@output}%
}%
\def\LT@start@new{%
 \let\LT@start\endgraf
 \endgraf
 \markthr@@{}%
 \LT@pre
 \@ifvoid\LT@firsthead{\LT@top}{\box\LT@firsthead\nobreak}%
 \mark@envir{longtable}%
}%
\def\LT@end@hd@ft@longtable#1{%
 \LT@echunk
 \ifx\LT@start\endgraf
  \LT@err{Longtable head or foot not at start of table}{Increase LTchunksize}%
 \fi
 \setbox#1\box\z@
 \LT@get@widths\LT@bchunk
}%
\def\LT@end@hd@ft@new#1{%
 \LT@echunk
 \@ifx{\LT@start\endgraf}{%
  \LT@err{Longtable head or foot not at start of table}{Increase LTchunksize}%
 }%
 \global\setbox#1\box\z@
 \LT@get@widths
 \LT@bchunk
}%
\def\LT@array@longtable[#1]#2{%
  \refstepcounter{table}\stepcounter{LT@tables}%
  \if l#1%
    \LTleft\z@ \LTright\fill
  \else\if r#1%
    \LTleft\fill \LTright\z@
  \else\if c#1%
    \LTleft\fill \LTright\fill
  \fi\fi\fi
  \let\LT@mcol\multicolumn
  \let\LT@@tabarray\@tabarray
  \let\LT@@hl\hline
  \def\@tabarray{%
    \let\hline\LT@@hl
    \LT@@tabarray}%
  \let\\\LT@tabularcr\let\tabularnewline\\%
  \def\newpage{\noalign{\break}}%
  \def\pagebreak{\noalign{\ifnum`}=0\fi\@testopt{\LT@no@pgbk-}4}%
  \def\nopagebreak{\noalign{\ifnum`}=0\fi\@testopt\LT@no@pgbk4}%
  \let\hline\LT@hline \let\kill\LT@kill\let\caption\LT@caption
  \@tempdima\ht\strutbox
  \let\@endpbox\LT@endpbox
  \ifx\extrarowheight\@undefined
    \let\@acol\@tabacol
    \let\@classz\@tabclassz \let\@classiv\@tabclassiv
    \def\@startpbox{\vtop\LT@startpbox}%
    \let\@@startpbox\@startpbox
    \let\@@endpbox\@endpbox
    \let\LT@LL@FM@cr\@tabularcr
  \else
    \advance\@tempdima\extrarowheight
    \col@sep\tabcolsep
    \let\@startpbox\LT@startpbox\let\LT@LL@FM@cr\@arraycr
  \fi
  \setbox\@arstrutbox\hbox{\vrule
    \@height \arraystretch \@tempdima
    \@depth \arraystretch \dp \strutbox
    \@width \z@}%
  \let\@sharp##\let\protect\relax
   \begingroup
    \@mkpream{#2}%
    \xdef\LT@bchunk{%
       \global\advance\c@LT@chunks\@ne
       \global\LT@rows\z@\setbox\z@\vbox\bgroup
       \LT@setprevdepth
       \tabskip\LTleft\halign to\hsize\bgroup
      \tabskip\z@ \@arstrut \@preamble \tabskip\LTright \cr}%
  \endgroup
  \expandafter\LT@nofcols\LT@bchunk&\LT@nofcols
  \LT@make@row
  \m@th\let\par\@empty
  \everycr{}\lineskip\z@\baselineskip\z@
  \LT@bchunk
}%
\def\LT@LR@l{\LTleft\z@   \LTright\fill}%
\def\LT@LR@r{\LTleft\fill \LTright\z@  }%
\def\LT@LR@c{\LTleft\fill \LTright\fill}%
\def\LT@array@new[#1]#2{%
 \refstepcounter{table}\stepcounter{LT@tables}%
 \table@hook
 \LTleft\fill \LTright\fill
 \csname LT@LR@#1\endcsname
 \let\LT@mcol\multicolumn
 \let\LT@@hl\hline
 \prepdef\@tabarray{\let\hline\LT@@hl}%
 \let\\\LT@tabularcr
 \let\tabularnewline\\%
 \def\newpage{\noalign{\break}}%
 \def\pagebreak{\noalign{\ifnum`}=0\fi\@testopt{\LT@no@pgbk-}4}%
 \def\nopagebreak{\noalign{\ifnum`}=0\fi\@testopt\LT@no@pgbk4}%
 \let\hline\LT@hline
 \let\kill\LT@kill
 \let\caption\LT@caption
 \@tempdima\ht\strutbox
 \let\@endpbox\LT@endpbox
 \@ifxundefined\extrarowheight{%
  \let\@acol\@tabacol
  \let\@classz\@tabclassz
  \let\@classiv\@tabclassiv
  \def\@startpbox{\vtop\LT@startpbox}%
  \let\@@startpbox\@startpbox
  \let\@@endpbox\@endpbox
  \let\LT@LL@FM@cr\@tabularcr
 }{%
  \advance\@tempdima\extrarowheight
  \col@sep\tabcolsep
  \let\@startpbox\LT@startpbox
  \let\LT@LL@FM@cr\@arraycr
 }%
 \let\@acoll\@tabacoll
 \let\@acolr\@tabacolr
 \let\@acol\@tabacol
 \setbox\@arstrutbox\hbox{%
  \vrule
  \@height \arraystretch \@tempdima
  \@depth \arraystretch \dp \strutbox
  \@width \z@
 }%
 \let\@sharp##%
 \let\protect\relax
 \begingroup
  \@mkpream{#2}%
  \@mkpream@relax
  \edef\@preamble{\@preamble}%
  \prepdef\@preamble{%
   \global\advance\c@LT@chunks\@ne
   \global\LT@rows\z@
   \setbox\z@\vbox\bgroup
    \LT@setprevdepth
    \tabskip\LTleft
    \halign to\hsize\bgroup
     \tabskip\z@
     \@arstrut
  }%
  \appdef\@preamble{%
     \tabskip\LTright
     \cr
  }%
  \global\let\LT@bchunk\@preamble
 \endgroup
 \expandafter\LT@nofcols\LT@bchunk&\LT@nofcols
 \LT@make@row
 \m@th
 \let\par\@empty
 \everycr{}%
 \lineskip\z@
 \baselineskip\z@
 \LT@bchunk
}%
\appdef\table@hook{}%
%% Note that it is not enough to define the environment itself; we also have to create the corresponding
%% \cmd\output\ routine procedures, which provide for continued footers and headers
%% (the very feature of \env{longtable} requiring support in the output routine).
%% This same consideration would arise in defining any syntactic extension to \env{longtable}, because
%% the environment name itself is exposed in the output routine.
\def\switch@longtable{%
 \@ifpackageloaded{longtable}{%
  \@ifx{\longtable\longtable@longtable}{%
   \@ifx{\endlongtable\endlongtable@longtable}{%
    \@ifx{\LT@start\LT@start@longtable}{%
     \@ifx{\LT@end@hd@ft\LT@end@hd@ft@longtable}{%
      \@ifx{\LT@array\LT@array@longtable}{%
       \true@sw
      }{\false@sw}%
     }{\false@sw}%
    }{\false@sw}%
   }{\false@sw}%
  }{\false@sw}%
  {%
   \class@info{Patching longtable package}%
  }{%
   \class@info{Patching unrecognized longtable package. (Proceeding with fingers crossed)}%
  }%
  \let\longtable\longtable@new
  \let\endlongtable\endlongtable@new
  \let\LT@start\LT@start@new
  \let\LT@end@hd@ft\LT@end@hd@ft@new
  \let\LT@array\LT@array@new
  \newenvironment{longtable*}{%
   \onecolumngrid@push
   \longtable
  }{%
   \endlongtable
   \onecolumngrid@pop
  }%
 }{}%
}%
\def\LT@pre{\penalty\z@\vskip\LTpre}%
\def\LT@bot{\nobreak\copy\LT@foot\vfil}%
\def\LT@top{\copy\LT@head\nobreak}%
\def\LT@post{\penalty\z@\addvspace\LTpost\mark@envir{\curr@envir}}%
\def\LT@adj{%
 \setbox\z@\vbox{\null}\dimen@-\ht\z@
 \setbox\z@\vbox{\unvbox\z@\LT@bot}\advance\dimen@\ht\z@
 \global\advance\vsize-\dimen@
}%
\def\output@init@longtable{\LT@adj}%
\def\output@prep@longtable{\setbox\@cclv\vbox{\unvbox\@cclv\LT@bot}}%
\def\output@post@longtable{\LT@top}%
\let\output@init@theindex\@empty
\let\output@prep@theindex\@empty
\def\output@post@theindex{%
 \@ifodd\c@page{}{%
  \@ifnum{\pagegrid@cur=\@ne}{% we have the leftmost column of a verso page
   % insert the current top-level continued head
  }%
 }%
}%
\def\check@aux{\do@output@MVL{\do@check@aux}}%
\def\check@deferlist@stuck#1{%
 \@ifx{\@deferlist@postshipout\@empty}{}{%
  \@ifx{\@deferlist@postshipout\@deferlist}{%
   \@fltstk
   \clearpage@sw{%
    \ltxgrid@warn{Deferred float stuck during \string\clearpage\space processing}%
   }{%
    \force@deferlist@stuck#1%
   }%
  }{%
   %Successfully committed float(s)
  }%
  \global\let\@deferlist@postshipout\@empty
 }%
}%
\def\@fltstk{%
 \@latex@warning{A float is stuck (cannot be placed without \string\clearpage)}%
}%
\appdef\@outputpage{%
 \global\let\@deferlist@postshipout\@deferlist
}%
\def\@next#1#2{%
 \@ifx{#2\@empty}{\false@sw}{%
  \expandafter\@xnext#2\@@#1#2%
  \true@sw
 }%
}%
\def\@xnext\@elt#1#2\@@#3#4{%
 \def#3{#1}%
 \gdef#4{#2}%
 \def\@tempa{#4}\def\@tempb{\@freelist}%
 \@ifx{\@tempa\@tempb}{%
  \@ifx{#4\@empty}{%
   \force@deferlist@empty%{Float register pool exhausted}%
  }{}%
 }{}%
}%
\def\force@deferlist@stuck#1{%
\force@deferlist@sw{%
  \@booleantrue\clearpage@sw
  \@booleantrue\forcefloats@sw
  #1%
}{%
}%
}%
\def\force@deferlist@empty{%
 \force@deferlist@sw{%
  \penalty-\pagebreak@pen
  \protect@penalty\do@forcecolumn@pen
 }{%
 }%
}%
\@booleanfalse\force@deferlist@sw
\mathchardef\do@forcecolumn@pen=10009
\@namedef{output@-\the\do@forcecolumn@pen}{\do@forcecolumn}%
\def\do@forcecolumn{%
 \@booleantrue\clearpage@sw
 \@booleantrue\forcefloats@sw
 \do@startcolumn
}%
\def\enlargethispage{%
\@ifstar{%
\@enlargethispage{}%
}{%
\@enlargethispage{}%
}%
}%
\def\@enlargethispage#1#2{%
 \begingroup
  \dimen@#2\relax
  \edef\@tempa{#1}%
  \edef\@tempa{\noexpand\@@enlargethispage{\@tempa}{\the\dimen@}}%
  \expandafter\do@output@MVL\expandafter{\@tempa}%
 \endgroup
}%
\def\@@enlargethispage#1#2{%
 \def\@tempa{one}%
 \@ifx{\thepagegrid\@tempa}{%
  \true@sw
 }{%
  \def\@tempa{mlt}%
  \@ifx{\thepagegrid\@tempa}{%
   \@ifnum{\pagegrid@cur=\@ne}{% OK to adjust this page
    \gdef\enlarge@colroom{#2}%
    \true@sw
   }{% Can only adjust this column; give up
    \ltxgrid@warn{Too late to enlarge this page; move the command to the first column.}%
    \false@sw
   }%
  }{% Unknown page grid
   \ltxgrid@warn{Unable to enlarge a page of this kind.}%
   \false@sw
  }%
 }%
 {%
  \class@info{Enlarging page \thepage\space by #2}%
  \global\advance\@colroom#2\relax
  \set@vsize
 }{%
  % Could not adjust this page
 }%
}%
\let\enlarge@colroom\@empty
\let\@kludgeins\@undefined
\@booleantrue\textheight@sw
\prepdef\@outputpage{%
 \textheight@sw{%
  \count@\vbadness\vbadness\@M
  \dimen@\vfuzz\vfuzz\maxdimen
  \setbox\@outputbox\vbox to\textheight{\unvbox\@outputbox}%
  \vfuzz\dimen@
  \vbadness\count@
 }{}%
}%
\def\ltxgrid@info{%
 \ltxgrid@info@sw{\class@info}{\@gobble}%
}%
\@booleanfalse\ltxgrid@info@sw
\def\ltxgrid@warn{%
 \ltxgrid@warn@sw{\class@warn}{\@gobble}%
}%
\@booleantrue\ltxgrid@warn@sw
%% Source File: textcase.dtx
%% Copyright 1997 1998 David Carlisle
%%
%% This file may be distributed under the terms of the LPPL.
%% See 00readme.txt for details.
\NeedsTeXFormat{LaTeX2e}
\ProvidesPackage{textcase}
       [1998/11/12 v0.06 Text only upper/lower case changing (DPC)]
\def\@uclcnotmath#1#2#3#4{\begingroup
      #1%
      \def\({$}\let\)\(%
      \def\NoCaseChange##1{\noexpand\NoCaseChange{\noexpand##1}}%
      \@nonchangecase\label
      \@nonchangecase\ref
      \@nonchangecase\ensuremath
      \def\cite##1##{\toks@{\noexpand\cite##1}\@citex}%
      \def\@citex##1{\NoCaseChange{\the\toks@{##1}}}%
      \def\reserved@a##1##2{\let#2\reserved@a}%
      \expandafter\reserved@a\@uclclist\reserved@b{\reserved@b\@gobble}%
      \protected@edef\reserved@a{\endgroup
          \noexpand\@skipmath#3#4$\valign$}%
      \reserved@a}
\def\@nonchangecase#1{\def#1##1{\NoCaseChange{#1{##1}}}}
\let\NoCaseChange\@firstofone
\def\@skipmath#1#2$#3${%
  \@skip@nonchangecase#1#2\NoCaseChange\valign
  \ifx\valign#3%
  \else
    $#3$%
    \expandafter\@skipmath\expandafter#1%
  \fi}
\def\@skip@nonchangecase#1#2\NoCaseChange#3{%
  #1{#2}%
  \ifx\valign#3%
  \else
    #3%
    \expandafter\@skip@nonchangecase\expandafter#1%
  \fi}
\DeclareRobustCommand\MakeTextUppercase{%
  \@uclcnotmath{\def\i{I}\def\j{J}}{##1##2}\uppercase}
\protected@edef\MakeTextUppercase#1{\MakeTextUppercase{#1}}
\DeclareRobustCommand\MakeTextLowercase{%
  \@uclcnotmath{}{##2##1}\lowercase}
\protected@edef\MakeTextLowercase#1{\MakeTextLowercase{#1}}
\DeclareOption{overload}{%
  \expandafter\let\csname MakeUppercase \expandafter\endcsname
                  \csname MakeTextUppercase \endcsname
  \expandafter\let\csname MakeLowercase \expandafter\endcsname
                  \csname MakeTextLowercase \endcsname}
\ProcessOptions
\let\ProvidesPackage\ProvidesPackage@latex
\let\ProcessOptions\ProcessOptions@latex
\let\DeclareOption\DeclareOption@latex
\expandafter
\let\csname MakeUppercase \expandafter\endcsname
    \csname MakeTextUppercase \endcsname
\expandafter
\let\csname MakeLowercase \expandafter\endcsname
    \csname MakeTextLowercase \endcsname
\appdef\class@documenthook{%
 \switch@longtable
}%
\appdef\class@documenthook{%
 \expandafter\@ifnotrelax\csname boldmath \endcsname{}{%
  \expandafter\let\csname boldmath \endcsname\boldmath
  \expandafter\def\expandafter\boldmath\expandafter{\expandafter\protect\csname boldmath \endcsname}%
 }%
}%
\appdef\class@documenthook{%
 \class@inithook
}%
\def\class@inithook{}%
\appdef\class@inithook{%
 \@ifxundefined\twocolumn@sw{\@booleanfalse\twocolumn@sw}{}%
}%
\appdef\class@inithook{%
 \@ifxundefined\footinbib@sw{\@booleanfalse\footinbib@sw}{}%
}%
\appdef\class@inithook{%
 \@ifxundefined\preprintsty@sw{\@booleanfalse\preprintsty@sw}{}%
}%
\appdef\class@inithook{%
 \@ifxundefined\preprint@sw{\@booleanfalse\preprint@sw}{}%
 \preprintsty@sw{\@booleantrue\preprint@sw}{}%
}%
\appdef\class@inithook{%
 \@ifxundefined\galley@sw{\@booleanfalse\galley@sw}{}%
}%
\appdef\class@inithook{%
 \@ifxundefined\showPACS@sw{\@booleanfalse\showPACS@sw}{}%
 \@ifxundefined\showKEYS@sw{\@booleanfalse\showKEYS@sw}{}%
}%
\appdef\class@inithook{%
 \@ifxundefined\twoside@sw{\@booleanfalse\twoside@sw}{}%
 \twoside@sw{\@twosidetrue}{\@twosidefalse}%
}%
\appdef\class@inithook{%
 \@ifxundefined\draft@sw{\@booleanfalse\draft@sw}{}%
 \draft@sw{\overfullrule 5\p@}{\overfullrule\z@}%
}%
\appdef\class@inithook{%
 \@ifxundefined\tightenlines@sw{\@booleanfalse\tightenlines@sw}{}%
}%
\appdef\class@inithook{%
 \@ifxundefined\lengthcheck@sw{\@booleanfalse\lengthcheck@sw}{}%
}%
\appdef\class@inithook{%
 \@ifxundefined\eqsecnum@sw{\@booleanfalse\eqsecnum@sw}{}%
 \eqsecnum@sw{%
  \@addtoreset{equation}{section}%
  \def\theequation@prefix{\arabic{section}.}%
 }{}%
}%
\appdef\class@inithook{%
 \@ifxundefined\byrevtex@sw{\@booleanfalse\byrevtex@sw}{}%
}%
\DeclareOption{checkin}{%
  \@booleantrue\dateinRH@sw
  \@booleantrue\preprintsty@sw
  \def\@pointsize{12}%
  \@booleantrue\showPACS@sw
  \@booleantrue\showKEYS@sw
  \def\fp@proc@h{\allow@breaking@tables}%
  \def\fp@proc@H{\allow@breaking@tables}%
}%
\@booleanfalse\dateinRH@sw
\def\checkindate{\dateinRH@sw{{\tiny(\today)}}{}}%
\def\allow@breaking@tables{%
 \def\array@default{v}% tabular can break over pages
 \@booleanfalse\floats@sw % table can break over pages
}%
\DeclareOption{preprint}{%
  \@booleantrue\preprintsty@sw
  \def\@pointsize{12}%
}%
\DeclareOption{manuscript}{%
 \class@warn{Document class option manuscript is obsolete; use preprint instead}%
 \ExecuteOptions{preprint}%
}%
\DeclareOption{showpacs}{%
  \@booleantrue\showPACS@sw
}%
\DeclareOption{noshowpacs}{%
  \@booleanfalse\showPACS@sw
}%
\DeclareOption{showkeys}{%
  \@booleantrue\showKEYS@sw
}%
\DeclareOption{noshowkeys}{%
  \@booleanfalse\showKEYS@sw
}%
\DeclareOption{balancelastpage}{%
  \@booleantrue\balancelastpage@sw
}%
\DeclareOption{nobalancelastpage}{%
  \@booleanfalse\balancelastpage@sw
}%
\appdef\class@inithook{%
 \@ifxundefined\balancelastpage@sw{%
  \@booleantrue\balancelastpage@sw
 }{}%
}%
\DeclareOption{nopreprintnumbers}{\@booleanfalse\preprint@sw}%
\DeclareOption{preprintnumbers}{\@booleantrue\preprint@sw}%
\DeclareOption{hyperref}{%
 \class@warn{Class option "hyperref" is no longer supported.^^JEmploy \string\usepackage{hyperref} instead}%
}%
\DeclareOption{10pt}{\def\@pointsize{10}}
\DeclareOption{11pt}{\def\@pointsize{11}}
\DeclareOption{12pt}{\def\@pointsize{12}}
\DeclareOption{a4paper}
   {\setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}}
\DeclareOption{a5paper}
   {\setlength\paperheight {210mm}%
    \setlength\paperwidth  {148mm}}
\DeclareOption{b5paper}
   {\setlength\paperheight {250mm}%
    \setlength\paperwidth  {176mm}}
\DeclareOption{letterpaper}
   {\setlength\paperheight {11in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{legalpaper}
   {\setlength\paperheight {14in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{executivepaper}
   {\setlength\paperheight {10.5in}%
    \setlength\paperwidth  {7.25in}}
\DeclareOption{landscape}
   {\setlength\@tempdima   {\paperheight}%
    \setlength\paperheight {\paperwidth}%
    \setlength\paperwidth  {\@tempdima}}
\DeclareOption{bibnotes}{\@booleantrue\bibnotes@sw}%
\DeclareOption{nobibnotes}{\@booleanfalse\bibnotes@sw}%
\appdef\class@enddocumenthook{\auto@bib}%
\appdef\class@inithook{%
 \@ifxundefined\bibnotes@sw{\@booleanfalse\bibnotes@sw}{}%
 \bibnotes@sw{%
  \let\frontmatter@footnote@produce\frontmatter@footnote@produce@endnote
 }{%
  \let\frontmatter@footnote@produce\frontmatter@footnote@produce@footnote
 }%
}%
\DeclareOption{footinbib}{\@booleantrue\footinbib@sw}
\DeclareOption{nofootinbib}{\@booleanfalse\footinbib@sw}
\DeclareOption{altaffilletter}{\@booleantrue\altaffilletter@sw}%
\DeclareOption{altaffilsymbol}{\@booleanfalse\altaffilletter@sw}%
\appdef\class@inithook{%
 \@ifxundefined\altaffilletter@sw{\@booleanfalse\altaffilletter@sw}{}%
}%
\DeclareOption{superbib}{%
 \def\place@bibnumber{\place@bibnumber@sup}%
}%
\DeclareOption{citeautoscript}{%
 \@booleantrue\citeautoscript@sw
}%
\@booleanfalse\citeautoscript@sw
\DeclareOption{twoside}{\@booleantrue\twoside@sw \@mparswitchfalse}
\DeclareOption{oneside}{\@booleanfalse\twoside@sw\@mparswitchtrue}
\AtEndOfClass{%
 \def\clear@document{%
  \close@column
 }%
 \appdef\class@inithook{%
  \@ifxundefined\close@column{%
   \let\close@column\close@column@default
  }{}%
 }%
}%
\DeclareOption{twocolumn}{%
  \@booleantrue\twocolumn@sw
}%
\DeclareOption{onecolumn}{%
  \@booleanfalse\twocolumn@sw
}%
\DeclareOption{galley}{%
  \ExecuteOptions{onecolumn}%
  \@booleantrue\galley@sw
  \appdef\set@pica@hook{%
    \advance\textwidth-\columnsep
    \textwidth.5\textwidth
  }%
}%
\DeclareOption{raggedbottom}{\@booleantrue\raggedcolumn@sw}
\DeclareOption{flushbottom}{\@booleanfalse\raggedcolumn@sw}
\DeclareOption{tightenlines}{\@booleantrue\tightenlines@sw}
\DeclareOption{lengthcheck}{%
 \@booleantrue\lengthcheck@sw
 \@booleantrue\twoside@sw
 \@booleantrue\twocolumn@sw
}%
\DeclareOption{draft}{%
 \@booleantrue\draft@sw
}%
\DeclareOption{final}{%
 \@booleanfalse\draft@sw
}%
\DeclareOption{eqsecnum}{\@booleantrue\eqsecnum@sw}
\DeclareOption{secnumarabic}{\@booleantrue\secnumarabic@sw}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\DeclareOption{floats}{\@booleantrue\floats@sw\@booleanfalse\floatp@sw}
\DeclareOption{endfloats}{\@booleanfalse\floats@sw\@booleanfalse\floatp@sw}
\DeclareOption{endfloats*}{\@booleanfalse\floats@sw\@booleantrue\floatp@sw}
\appdef\class@inithook{%
 \@ifxundefined\floats@sw{\@booleantrue\floats@sw}{}%
 \@ifxundefined\floatp@sw{\@booleantrue\floatp@sw}{}%
}%
\DeclareOption{titlepage}{\@booleantrue\titlepage@sw}
\DeclareOption{notitlepage}{\@booleanfalse\titlepage@sw}
\appdef\defaults@hook{%
 \@ifxundefined\titlepage@sw{\@booleanfalse\titlepage@sw}{}%
}%
\def\change@society#1{%
 \def\@tempa{#1}%
 \@ifxundefined\@society{}{%
  \@ifx{\@tempa\@society}{}{%
   \class@warn{Conflicting society \@tempa<>\@society.}%
  }%
 }%
 \let\@society\@tempa
}%
\DeclareOption{aps}{\change@society{aps}\def\@journal{}}%
\DeclareOption{pra}{\change@society{aps}\def\@journal{pra}}%
\DeclareOption{prb}{\change@society{aps}\def\@journal{prb}}%
\DeclareOption{prc}{\change@society{aps}\def\@journal{prc}}%
\DeclareOption{prd}{\change@society{aps}\def\@journal{prd}}%
\DeclareOption{pre}{\change@society{aps}\def\@journal{pre}}%
\DeclareOption{prl}{\change@society{aps}\def\@journal{prl}}%
\DeclareOption{prstab}{\change@society{aps}\def\@journal{prstab}}%
\DeclareOption{rmp}{\change@society{aps}\def\@journal{rmp}}%
\DeclareOption{osa}{\change@society{osa}\def\@journal{}}%
\DeclareOption{osameet}{\change@society{osa}\def\@journal{osameet}}%
\DeclareOption{opex}{\change@society{osa}\def\@journal{opex}}%
\DeclareOption{tops}{\change@society{osa}\def\@journal{tops}}%
\DeclareOption{josa}{\change@society{osa}\def\@journal{josa}}%
\DeclareOption{amsfonts}{%
  \def\class@amsfonts{\RequirePackage{amsfonts}}%
}%
\DeclareOption{noamsfonts}{%
  \let\class@amsfonts\@empty
}%
\AtEndOfClass{%
 \@ifxundefined\class@amsfonts{}{\class@amsfonts}%
}%
\DeclareOption{amssymb}{%
  \def\class@amssymb{\RequirePackage{amssymb}}%
}%
\DeclareOption{noamssymb}{%
  \let\class@amssymb\@empty
}%
\AtEndOfClass{%
 \@ifxundefined\class@amssymb{}{\class@amssymb}%
}%
\DeclareOption{amsmath}{%
  \def\class@amsmath{\RequirePackage{amsmath}[\ver@amsmath@prefer]}%
}%
\DeclareOption{noamsmath}{%
  \let\class@amsmath\@empty
}%
\AtEndOfClass{%
 \preserve@LaTeX
 \@ifxundefined\class@amsmath{}{\class@amsmath}%
 \appdef\class@enddocumenthook{\test@amsmath@ver}%
}%
\appdef\preserve@LaTeX{%
 \let\@ifl@aded@LaTeX\@ifl@aded
 \let\@ifpackageloaded@LaTeX\@ifpackageloaded
 \let\@pkgextension@LaTeX\@pkgextension
 \let\@ifpackagelater@LaTeX\@ifpackagelater
 \let\@ifl@ter@LaTeX\@ifl@ter
 \let\@ifl@t@r@LaTeX\@ifl@t@r
 \let\@parse@version@LaTeX\@parse@version
}%
\appdef\restore@LaTeX{%
 \let\@ifl@aded\@ifl@aded@LaTeX
 \let\@ifpackageloaded\@ifpackageloaded@LaTeX
 \let\@pkgextension\@pkgextension@LaTeX
 \let\@ifpackagelater\@ifpackagelater@LaTeX
 \let\@ifl@ter\@ifl@ter@LaTeX
 \let\@ifl@t@r\@ifl@t@r@LaTeX
 \let\@parse@version\@parse@version@LaTeX
}%
\def\test@amsmath@ver{%
 \begingroup
  \restore@LaTeX
  \@ifpackageloaded{amsmath}{%
   \@ifpackagelater{amsmath}{\ver@amsmath@prefer}{}{%
     \class@warn{%
      You have loaded amsmath, version "\csname ver@amsmath.sty\endcsname",\MessageBreak
      but this class requires version "\ver@amsmath@prefer", or later.\MessageBreak
      Please update your LaTeX installation.
     }%
   }%
  }{%
  }%
 \endgroup
}%
\def\ver@amsmath@prefer{2000/01/15 v2.05 AMS math features}%
\DeclareOption{groupedaddress}{\clo@groupedaddress}%
\def\clo@groupedaddress{%
 \@booleantrue\groupauthors@sw
 \@booleantrue\@affils@sw
 \@booleanfalse\runinaddress@sw
}%
\appdef\class@inithook{%
 \@ifxundefined\groupauthors@sw{\clo@groupedaddress}{}%
}%
\DeclareOption{unsortedaddress}{\clo@unsortedaddress}%
\def\clo@unsortedaddress{%
 \@booleantrue\groupauthors@sw
 \@booleanfalse\@affils@sw
 \@booleanfalse\runinaddress@sw
}%
\DeclareOption{runinaddress}{\clo@runinaddress}%
\def\clo@runinaddress{%
 \@booleantrue\groupauthors@sw
 \@booleantrue\@affils@sw
 \@booleantrue\runinaddress@sw
}%
\DeclareOption{superscriptaddress}{\clo@superscriptaddress}%
\def\clo@superscriptaddress{%
 \@booleanfalse\groupauthors@sw
 \@booleantrue\@affils@sw
 \@booleanfalse\runinaddress@sw
}%
\DeclareOption{byrevtex}{\@booleantrue\byrevtex@sw}%
\DeclareOption{floatfix}{\@booleantrue\force@deferlist@sw}%
\DeclareOption{nofloatfix}{\@booleanfalse\force@deferlist@sw}%
\appdef\class@inithook{%
 \@ifxundefined\force@deferlist@sw{\@booleanfalse\force@deferlist@sw}{}%
}%
\gdef\@fltovf{%
 \@latex@error{%
  Too many unprocessed floats%
  \force@deferlist@sw{}{; try class option [floatfix]}%
 }\@ehb
}%
\def\@fltstk{%
 \@latex@warning{%
  A float is stuck (cannot be placed)%
  \force@deferlist@sw{}{; try class option [floatfix]}%
 }%
}%
\DeclareOption{ltxgridinfo}{\@booleantrue\ltxgrid@info@sw}%
\DeclareOption{outputdebug}{\@booleantrue\outputdebug@sw}%
\DeclareOption{raggedfooter}{\@booleanfalse\textheight@sw}%
\DeclareOption{noraggedfooter}{\@booleantrue\textheight@sw}%
\DeclareOption{newabstract}{\@booleanfalse\oldabstract@sw}%
\DeclareOption{oldabstract}{\@booleantrue \oldabstract@sw}%
\appdef\class@inithook{%
 \@ifxundefined\oldabstract@sw{\@booleanfalse\oldabstract@sw}{}%
 \oldabstract@sw{%
  \let@environment{abstract}{abstract-old}%
 }{}%
}%
\DeclareOption*{%
 \IfFileExists{\CurrentOption\substyle@ext}{%
  \expandafter\change@society\expandafter{\CurrentOption}%
 }{%
  \OptionNotUsed
 }%
}%
\def\substyle@ext{.rtx}%
\AtEndOfClass{%
 \@ifxundefined\@society{%
  \class@warn{No Society specified, using generic APS style}%
  \def\@society{aps}\def\@journal{}%
 }{}%
 \expandafter\input\expandafter{\@society\substyle@ext}%
}%
\ExecuteOptions{10pt}
\ProcessOptions\relax
  \def\ps@headings{%
      \let\@oddfoot\@empty\let\@evenfoot\@empty
      \def\@evenhead{\thepage\hfil\slshape\leftmark}%
      \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
      \let\@mkboth\markboth
    \def\sectionmark##1{%
      \markboth {\MakeUppercase{%
        \ifnum \c@secnumdepth >\z@
          \thesection\quad
        \fi
        ##1}}{}}%
    \def\subsectionmark##1{%
      \markright {%
        \ifnum \c@secnumdepth >\@ne
          \thesubsection\quad
        \fi
        ##1}}}%
\def\ps@myheadings{%
    \let\@oddfoot\@empty\let\@evenfoot\@empty
    \def\@evenhead{\thepage\hfil\slshape\leftmark}%
    \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
    \let\@mkboth\@gobbletwo
    \let\sectionmark\@gobble
    \let\subsectionmark\@gobble
    }%
\def\ps@article{%
    \@provide\@evenhead{\let\\\heading@cr\thepage\quad\checkindate\hfil{\leftmark}}%
    \@provide\@oddhead{\let\\\heading@cr{\rightmark}\hfil\checkindate\quad\thepage}%
    \@provide\@oddfoot{}%
    \@provide\@evenfoot{}%
    \let\@mkboth\markboth
  \let\sectionmark\@gobble
  \let\subsectionmark\@gobble
}%
\def\ps@article@final{%
    \@provide\@evenhead{\let\\\heading@cr\thepage\quad\checkindate\hfil{\leftmark}}%
    \@provide\@oddhead{\let\\\heading@cr{\rightmark}\hfil\checkindate\quad\thepage}%
    \@provide\@oddfoot{}%
    \@provide\@evenfoot{}%
    \let\@mkboth\markboth
    \def\sectionmark##1{%
      \markboth{%
       \MakeTextUppercase{%
        \@ifnum{\c@secnumdepth >\z@}{\thesection\hskip 1em\relax}{}%
         ##1%
       }%
       }{}%
    }%
    \def\subsectionmark##1{%
      \markright {%
        \@ifnum{\c@secnumdepth >\@ne}{\thesubsection\hskip 1em\relax}{}%
         ##1%
      }%
    }%
}%
\appdef\class@inithook{%
}%
\def\heading@cr{\unskip\space\ignorespaces}%
\def\ps@preprint{%
  \@provide\@oddfoot{\hfil\thepage\quad\checkindate\hfil}%
  \@provide\@evenfoot{\hfil\thepage\quad\checkindate\hfil}%
  \@provide\@oddhead{}%
  \@provide\@evenhead{}%
  \let\@mkboth\@gobbletwo
  \let\sectionmark\@gobble
  \let\subsectionmark\@gobble
}%
\let\@oddhead\@undefined
\let\@evenhead\@undefined
\let\@oddfoot\@undefined
\let\@evenfoot\@undefined
\AtEndOfClass{%
 \appdef\class@documenthook{%
  \@provide\@oddhead{}%
  \@provide\@evenhead{}%
  \@provide\@oddfoot{}%
  \@provide\@evenfoot{}%
 }%
}%
\def\lastpage@putlabel{%
 \if@filesw
  \begingroup
    \advance\c@page\m@ne
    \immediate\write\@auxout{\string\newlabel{LastPage}{{}{\thepage}}}%
  \endgroup
 \fi
}%
\AtEndOfClass{%
 \appdef\clear@document{%
  \clearpage
  \lastpage@putlabel
 }%
}%
\RequirePackage{url}%
\appdef\class@inithook{%
 \incompatible@package{cite}%
 \incompatible@package{multicol}%
}%
\renewcommand\title[2][]{%
 \def\@title{#2}%
 \def\@shorttitle{#1}%
 \let\@AF@join\@title@join
}%
\appdef\class@inithook{%
 \def\@title{\class@warn{No title}}%
 \let\@shorttitle\@empty
 \let\@title@aux\@title@aux@cleared
}%
\def\@title@join{\expandafter\@title@join@\@title@aux}%
\def\@title@join@#1#2{%
 \def\@title@aux{{\@join{\@separator}{#1}{#2}}}%
}%
\def\@title@aux@cleared{{}}%
\newcounter{affil}
\renewcommand\author{% implicit #1
 \@author@def{}% implicit #2
}%
\def\collaboration{% implicit #1
 \groupauthors@sw{%
  \class@warn{%
   \string\collaboration\space is only allowed with superscriptaddress option.
   Use superscriptaddress in place of groupedaddress, unsortedaddress, or runinaddress%
  }%
 }{}%
 \@author@def{\@booleantrue\collaboration@sw}% implicit #2
}%
\def\@author@cleared{{}{}}%
\def\@author@init{%
 \let\@author\@author@cleared
}%
\def\@authorclear@sw{\@ifx{\@author\@author@cleared}}%
\appdef\class@inithook{%
 \@author@init
}%
\def\@author@def#1#2{%
 \move@AU\move@AF\move@AUAF
 \let\@AF@join\@author@join
 #1%
 \def\@author{{#2}{}}%
}%
\def\@author@join@#1#2#3{%
  \def\@author{{#1}{\@join{\@separator}{#2}{#3}}}%
}%
\def\@author@join{\expandafter\@author@join@\@author}%
\def\@author@gobble#1#2{}%
\def\move@AU{%
 \@authorclear@sw{}{%
  \collaboration@sw{%
    \expandafter\appdef
    \expandafter\AU@grp
    \expandafter{%
    \expandafter \CO@opr
                 \@author
                }%
  }{%
    \expandafter\appdef
    \expandafter\AU@grp
    \expandafter{%
    \expandafter \AU@opr
                 \@author
                }%
  }%
 }%
 \@author@init
 \@booleanfalse\collaboration@sw
}%
\def\@author@finish{%
 \move@AU\move@AF
 \@ifx{\AU@grp\@empty}{}{%
  \@ifx{\AF@grp\@empty}{%
   \begingroup
    \let\href\@secondoftwo
    \let\AU@opr\@firstoftwo
    \let\CO@opr\@firstoftwo
    \class@warn{Assuming \string\noaffiliation\space for \AU@grp}%
   \endgroup
   \begingroup\@affiliation{\blankaffiliation}%\endgroup in \@affiliation
   \move@AF
  }{}%
 }%
 \move@AUAF
}%
\def\@join#1#2#3{%
  \@if@empty{#2}{#3}{#2#1#3}%
}%
\def\@separator{;\space}%
\let\surname\@firstofone
\let\firstname\@firstofone
\let\AU@grp\@empty
\renewcommand\and{\class@err{\protect\and\space is not supported}}
\def\cat@comma@active{\catcode`\,\active}%
{\cat@comma@active\gdef,{\active@comma}}%
\def\active@comma{,\penalty-300\relax}%
\newcommand\affiliation{%
 \move@AU\move@AF
 \begingroup
  \cat@comma@active
  \@affiliation
}%
\def\@affiliation#1{%
 \endgroup
 \let\@AF@join\@affil@join
 \@affil@def{#1}%
}%
\newcommand\noaffiliation{%
  \move@AU\move@AF
  \begingroup\@affiliation{\blankaffiliation}%
}%
\def\blankaffiliation{}%
\def\@affil@cleared{{}{}}%
\def\@affil@init{%
 \let\@affil\@affil@cleared
}%
\def\@affilclear@sw{\@ifx{\@affil\@affil@cleared}}%
\appdef\class@inithook{%
 \@affil@init
}%
\def\@affil@def#1{\def\@affil{{#1}{}}}%
\def\@affil@join@#1#2#3{%
  \def\@affil{{#1}{\@join{\@separator}{#2}{#3}}}%
}%
\def\@affil@join{\expandafter\@affil@join@\@affil}%
\def\move@AF{%
 \@affilclear@sw{}{%
  \@booleanfalse\temp@sw
  \let\@tempd\@empty
  \@affils@sw{%
    \expandafter\@firstoftwodef\expandafter\@tempa\@affil
    \def\AFF@opr{\@affil@match\@tempa}%
    \@AFF@list
  }{}\temp@sw
  {%
   \expandafter\@secondoftwodef\expandafter\@tempb\@affil
   \@ifx{\@tempb\@empty}{}{%
    \@ifx{\@tempb\@tempd}{}{%
     \class@warn{%
      Ancillary information for \@tempa\space must not be different!
      Please put it on the first instance%
     }%
    }%
   }%
   \appdef@eval\AF@grp\@tempc
  }%
  {%
    \@if@empty{\expandafter\@firstoftwo\@affil}{%
\edef\@tempc{\z@}%
    }{%
\stepcounter{affil}%
\edef\@tempc{\theaffil}%
    }%
\appdef     \@AFF@list{\AFF@opr}%
\appdef@eval\@AFF@list\@tempc
\appdef@e   \@AFF@list{\@affil}%
\appdef@eval\AF@grp\@tempc
  }%
  \@affil@init
 }%
}%
\def\@firstoftwodef#1#2#3{%
 \def#1{#2}%
}%
\def\@secondoftwodef#1#2#3{%
 \def#1{#3}%
}%
\def\move@AUAF{%
 \@ifx{\AF@grp\@empty}%
 {}{%
    \appdef     \@AAC@list{\AF@opr}%
    \appdef@eval\@AAC@list{\AF@grp}%
    \appdef@e   \@AAC@list{\AU@grp}%
    \@ifx{\@empty\AU@grp}{}{%
     \@ifx{\@empty\AF@grp}{}{%
      \@booleanfalse\temp@sw
      \def\AFS@opr{\x@match\AF@grp}%
      \@AFG@list
      \temp@sw{}{%
        \appdef     \@AFG@list{\AFS@opr}%
        \appdef@eval\@AFG@list{\AF@grp}%
      }%
     }%
    }%
    \let\AU@grp\@empty
    \let\AF@grp\@empty
 }%
}%
\appdef\class@inithook{%
 \let\@AAC@list\@empty
}%
\appdef\class@inithook{%
 \let\@AFG@list\@empty
}%
\appdef\class@inithook{%
 \let\@AFF@list\@empty
}%
\let\AF@grp\@empty
\def\@affil@match#1#2#3#4{%
 \temp@sw{}{%
  \def\@tempifx{#3}%
  \@ifx{\@tempifx#1}{%
    \@booleantrue\temp@sw
    \def\@tempc{#2}%
    \def\@tempd{#4}%
  }{%
  }%
 }%
}%
\appdef\class@inithook{%
 \let\@AF@join\@AF@join@error
}%
\def\@AF@join@error#1{%
 \class@warn{%
  \string\email, \string\homepage, \string\thanks, or \string\altaffiliation\space
  appears in wrong context.
 }%
}%
\def\sanitize@url{%
 \@makeother\%%
 \@makeother\~%
 \@makeother\_%
}%
\newcommand*\email[1][Electronic address: ]{\begingroup\sanitize@url\@email{#1}}%
\def\@email#1#2{%
 \endgroup
 \@AF@join{#1\href{mailto:#2}{#2}}%
}%
\newcommand*\homepage[1][URL:~]{\begingroup\sanitize@url\@homepage{#1}}%
\def\@homepage#1#2{%
 \endgroup
 \@AF@join{#1\url{#2}}%
}%
\appdef\class@documenthook{%
 \providecommand\url{\texttt}%
 \providecommand\href[2]{\url{#2}}%
}%
\def\thanks{% implicit #1
  \@AF@join
}%
\newcommand*\altaffiliation[2][]{%
  \@AF@join{#1#2}%
}%
\def\set@listcomma@list#1{%
  \expandafter\@reset@ac\expandafter#1#1{0}\@reset@ac{%
   \let\@listcomma\relax
  }{%
   \let\@listcomma\@listcomma@comma
  }%
}%
\def\set@listcomma@count#1{%
\@ifnum{#1=\tw@}{%
\let\@listcomma\relax
}{%
\let\@listcomma\@listcomma@comma
}%
}%
\def\@reset@ac#1#2#3\@reset@ac{%
  \def#1{#3}%
  \@tempcnta#2\relax
  \@ifnum{#2=\tw@}%
}%
\def\@listand{\@ifnum{\@tempcnta=\tw@}{\andname\space}{}}%
\def\@listcomma@comma{\@ifnum{\@tempcnta>\@ne}{,}{}}%
\def\@listcomma@comma@UK{\@ifnum{\@tempcnta>\tw@}{,}{}}%
\def\@collaboration@err#1#2{%
}%
\def\doauthor#1#2#3{%
  \ignorespaces#1\unskip\@listcomma
  \def\@tempifx{#3}%
  \@ifx{\@tempifx\@empty}{%
   \def\@tempifx{#2}%
   \@ifx{\@tempifx\@empty}{}{\frontmatter@footnote{#2}}%
  }{%
   #3%
   \def\@tempifx{#2}%
   \@ifx{\@tempifx\@empty}{}{\comma@space\frontmatter@footnote{#2}}%
  }%
  \space \@listand
}%
\def\@thanks{\comma@space\footnote}
\def\x@match#1#2{%
 \temp@sw{}{%
  \def\@tempifx{#2}%
  \@ifx{\@tempifx#1}{%
    \@booleantrue\temp@sw
  }{%
  }%
 }%
}%
\def\y@match#1#2#3{%
 \temp@sw{}{%
  \def\@tempifx{#3}%
  \@ifx{\@tempifx#1}{%
    \@booleantrue\temp@sw
    \def\@tempb{#2}%
  }{%
  }%
 }%
}%
\def\frontmatter@footnote#1{%
 \begingroup
  \@booleanfalse\temp@sw
  \def\@tempa{#1}%
  \let\@tempb\@empty
  \def\@TBN@opr{\y@match\@tempa}%
  \@FMN@list
  \temp@sw{%
   \expandafter\frontmatter@footnotemark
   \expandafter{\@tempb}%
  }{%
   \stepcounter\@mpfn
   \protected@xdef\@tempb{\the\csname c@\@mpfn\endcsname}%
   \expandafter\frontmatter@footnotemark
   \expandafter{\@tempb}%
   \expandafter\g@addto@macro
   \expandafter\@FMN@list
   \expandafter{%
   \expandafter \@TBN@opr
   \expandafter {\@tempb}{#1}%
               }%
  }%
 \endgroup
}%
\appdef\class@inithook{%
 \global\let\@FMN@list\@empty
}%
\def\frontmatter@footnotemark#1{%
 \leavevmode
 \ifhmode\edef\@x@sf{\the\spacefactor}\nobreak\fi
  \begingroup
   \hyper@linkstart {link}{frontmatter@#1}%
    \csname c@\@mpfn\endcsname#1\relax
    \def\@thefnmark{\thempfn}%
    \@makefnmark
   \hyper@linkend
  \endgroup
 \ifhmode\spacefactor\@x@sf\fi
 \relax
}%
\def\keywords#1{%
  \aftermaketitle@chk{\keywords}%
  \gdef\@keywords{#1}%
}%
\def\@keywords{}%
\renewcommand*\date[2][\Dated@name]{\def\@date{#1#2}}%
\def\@date{}%
\newcommand*\received[2][\Received@name]{\def\@received{#1#2}}%
\def\@received{}%
\newcommand*\revised[2][\Revised@name]{\def\@revised{#1#2}}%
\def\@revised{}%
\newcommand*\accepted[2][\Accepted@name]{\def\@accepted{#1#2}}%
\def\@accepted{}%
\newcommand*\published[2][\Published@name]{\def\@published{#1#2}}%
\def\@published{}%
\def\pacs#1{%
  \aftermaketitle@chk{\pacs}%
  \gdef\@pacs{#1}%
}%
\def\@pacs{}%{\class@warn{No PACS numbers}}%
\def\preprint#1{\gappdef\@preprint{\preprint{#1}}}%
\let\@preprint\@empty
\newbox\absbox
\newenvironment{abstract}{%
  \aftermaketitle@chk{\begin{abstract}}%
  \global\setbox\absbox\vbox\bgroup
   \color@begingroup
   \columnwidth\textwidth
   \hsize\columnwidth
   \@parboxrestore
   \def\@mpfn{mpfootnote}\def\thempfn{\thempfootnote}\c@mpfootnote\z@
   \let\@footnotetext\frontmatter@footnotetext
   \minipagefootnote@init
   \let\@listdepth\@mplistdepth \@mplistdepth\z@
   \@minipagerestore
   \@setminipage
    \frontmatter@abstractheading
    \frontmatter@abstractfont
    \let\footnote\mini@note
}{%
    \par
    \unskip
    \minipagefootnote@here
    \@minipagefalse   %% added 24 May 89
    \color@endgroup
  \egroup
}%
\long\def\frontmatter@footnotetext#1{%
  \minipagefootnote@pick
    \reset@font\footnotesize
    \interlinepenalty \interfootnotelinepenalty
    \splittopskip \footnotesep
    \splitmaxdepth \dp \strutbox
    \set@footnotewidth
    \@parboxrestore
    \protected@edef\@currentlabel{\csname p@mpfootnote\endcsname\@thefnmark}%
    \color@begingroup
      \frontmatter@makefntext{%
        \rule\z@\footnotesep\ignorespaces#1\@finalstrut\strutbox
      }%
    \color@endgroup
  \minipagefootnote@drop
}%
\def\frontmatter@abstractheading{%
 \begingroup
  \centering\large
  \abstractname
  \par
 \endgroup
}%
\def\frontmatter@abstractfont{}%
\def\frontmatter@abstractwidth{\textwidth}
\def\frontmatter@abstract@produce{%
  \par
  \preprintsty@sw{%
   \do@output@MVL{%
    \vskip\frontmatter@preabstractspace
    \vskip200\p@\@plus1fil
    \penalty-200\relax
    \vskip-200\p@\@plus-1fil
   }%
  }{%
   \addvspace{\frontmatter@preabstractspace}%
  }%
   \begingroup
    \dimen@\baselineskip
    \setbox\z@\vtop{\unvcopy\absbox}%
    \advance\dimen@-\ht\z@\advance\dimen@-\prevdepth
    \@ifdim{\dimen@>\z@}{\vskip\dimen@}{}%
   \endgroup
   \unvbox\absbox
  \@ifx{\@empty\mini@notes}{}{\mini@notes\par}%
  \addvspace\frontmatter@postabstractspace
}%
\def\frontmatter@preabstractspace{.5\baselineskip}
\def\frontmatter@postabstractspace{.5\baselineskip}
\newenvironment{titlepage}{%
      \twocolumn@sw{\onecolumngrid}{\newpage}%
      \thispagestyle{titlepage}%
      \setcounter{page}\@ne
}{%
     \twocolumn@sw{\twocolumngrid}{\newpage}%
     \twoside@sw{}{%
        \setcounter{page}\@ne
     }%
}%
\def\maketitle{%
  \@author@finish
  \title@column\titleblock@produce
  \suppressfloats[t]%
  \let\and\relax
  \let\affiliation\@gobble
  \let\author\@gobble
  \let\@AAC@list\@empty
  \let\@AFF@list\@empty
  \let\@AFG@list\@empty
  \let\@AF@join\@AF@join@error
  \let\email\@gobble
  \let\@address\@empty
  \let\maketitle\relax
  \let\thanks\@gobble
  \titlepage@sw{%
   \vfil
   \clearpage
  }{}%
}%
\def\titleblock@produce{%
 \begingroup
  \let\footnote\footnote@latex
  \let\@makefnmark\@makefnmark@latex
  \let\@footnotemark\@footnotemark@latex
  \let\thefootnote\frontmatter@thefootnote
  \global\c@footnote\z@
  \let\@makefnmark\frontmatter@makefnmark
  \frontmatter@setup
  \thispagestyle{titlepage}\label{FirstPage}%
  \frontmatter@title@produce
  \groupauthors@sw{%
   \frontmatter@author@produce@group
  }{%
   \frontmatter@author@produce@script
  }%
  \frontmatter@RRAPformat{%
   \expandafter\produce@RRAP\expandafter{\@date}%
   \expandafter\produce@RRAP\expandafter{\@received}%
   \expandafter\produce@RRAP\expandafter{\@revised}%
   \expandafter\produce@RRAP\expandafter{\@accepted}%
   \expandafter\produce@RRAP\expandafter{\@published}%
  }%
  \frontmatter@abstract@produce
  \@pacs@produce\@pacs
  \@keywords@produce\@keywords
  \par
  \frontmatter@finalspace
 \endgroup
}%
\def\frontmatter@title@produce{%
 \begingroup
  \frontmatter@title@above
  \frontmatter@title@format
  \@title
  \@ifx{\@title@aux\@title@aux@cleared}{}{%
   \expandafter\frontmatter@footnote\expandafter{\@title@aux}%
  }%
  \par
%% \class@warn{Unused short title ignored}%
  \frontmatter@title@below
 \endgroup
}%
\appdef\let@mark{\let\\\relax}%
\def\frontmatter@title@above{}%
\def\frontmatter@title@format{}%
\def\frontmatter@title@below{\addvspace{\baselineskip}}%
\def\frontmatter@author@produce@script{%
  \begingroup
    \let\@author@present\@author@present@script
    \let\AU@temp\@empty
    \@tempcnta\z@
    \let\AF@opr \@gobble
    \def\AU@opr{\@author@count\@tempcnta}%
    \def\CO@opr{\@collaboration@count\AU@temp\@tempcnta}%
    \@AAC@list
    \expandafter\CO@opr\@author@cleared
    \begingroup
\frontmatter@authorformat
\let\AF@opr \@affilID@def
\let\AU@opr \@author@present
\def\CO@opr{\@collaboration@present\AU@temp}%
\set@listcomma@list\AU@temp
\@AAC@list
\unskip\unskip
\par
\endgroup
    \begingroup
\frontmatter@above@affiliation@script
\let\AFF@opr \@affil@script
     \@AFF@list
\frontmatter@footnote@produce
\par
\endgroup
  \endgroup
}%
\def\@author@count#1{\advance#1\@ne\@author@gobble}%
\def\@collaboration@present#1#2#3{%
 \par
 \begingroup
  \frontmatter@collaboration@above
  \@author@present{(\ignorespaces#2\unskip)}{#3}%
  \par
 \endgroup
 \set@listcomma@list#1%
}%
\def\frontmatter@collaboration@above{}%
\def\@collaboration@count#1#2{%
 \appdef@eval#1{\the#2}#2\z@
 \@author@gobble
}%
\def\@affilID@def{\def\@affilID@temp}%
\let\@affilID@temp\@empty
\def\affil@script#1#2#3{%
 \def\@tempifx{#1}\@ifx{\@tempifx\@tempa}{%
  \@if@empty{#2}{}{%
   \par
   \begingroup
    \textsuperscript{\normalfont#1}%
    #2%
    \@if@empty{#3}{}{\frontmatter@footnote{#3}}%
    \par
   \endgroup
  }%
 }{}%
}%
\def\@affil@script#1#2#3{%
 \@ifnum{#1=\z@}{}{%
  \par
  \begingroup
   \frontmatter@affiliationfont
   \textsuperscript{#1}%
   #2%
   \@if@empty{#3}{}{\frontmatter@footnote{#3}}%
   \par
  \endgroup
 }%
}%
\def\@author@present@script#1#2{%
 \begingroup
  \gdef\comma@space{\textsuperscript{,\,}}%
 \doauthor{#1}{#2}{%
  \let\@tempa\@empty
  \expandafter\@affil@present@script\@affilID@temp\relax
 }%
 \endgroup
 \advance\@tempcnta\m@ne
}%
\def\@affilcomma#1#2{%
 \@ifx{\z@#1}{%
  \@ifx{\relax#2}{}{%
   \@affilcomma{#2}%
  }%
 }{%
  #1%
  \@ifx{\relax#2}{}{%
   \@ifx{\z@#2}{%
    \@affilcomma
   }{%
    ,\,\@affilcomma{#2}%
   }%
  }%
 }%
}%
\def\@affil@present@script#1{%
 \@ifx{\relax#1}{%
  \@ifx{\@tempa\@empty}{}{%
   \textsuperscript{\expandafter\@affilcomma\@tempa\relax\relax}%
  }%
 }{%
  \@ifx{\z@#1}{}{\appdef\@tempa{{#1}}}%
  \@affil@present@script
 }%
}%
\@provide\@author@parskip{\z@skip}%
\def\frontmatter@author@produce@group{%
  \begingroup
    \let\@author@present\@author@present@group
    \frontmatter@authorformat
    \let\AFS@opr \affils@present@group
    \let\@listcomma\relax
    \@AFG@list
    \frontmatter@footnote@produce
    \par
  \endgroup
  \frontmatter@authorbelow
}%
\@provide\frontmatter@authorbelow{}%
\def\affils@present@group#1{%
 \begingroup
   \def\AF@temp{#1}%               % 0B13
   \@tempcnta\z@
   \let\AU@opr \undefined          % 0B13
   \let\CO@opr \undefined          % 0B13
   \def\AF@opr{\@affilID@count\AF@temp\@tempcnta}% 0B13
   \@AAC@list
   \@ifnum{\@tempcnta=\z@}{}{%
    \begingroup
     \frontmatter@above@affilgroup
     \set@listcomma@count\@tempcnta
     \let\AU@opr \undefined          % 0B13
     \let\CO@opr \undefined          % 0B13
     \def\AF@opr{\@affilID@match\AF@temp}%
     \@AAC@list
    \endgroup
    \begingroup
     \par
     \frontmatter@above@affiliation
     \frontmatter@affiliationfont
     \let\\\frontmatter@addressnewline
     \@tempcnta\z@
     \@tfor\AF@temp:=#1\do{%
      \expandafter\@ifx\expandafter{\z@\AF@temp}{}{%
       \advance\@tempcnta\@ne
      }%
     }%
     \@ifnum{\@tempcnta=\tw@}{%
      \let\@listcomma\relax
     }{}%
     \def\after@address{\@listcomma\ \@listand}%
     \runinaddress@sw{%
     }{%
      \tightenlines@sw{%
      }{%
       \parskip\z@ %space between affilations in grouped mode
      }%
      \appdef\after@address\par
     }%
     \let\AFF@opr \@affil@group
     \do@affil@fromgroup\@AFF@list#1\relax
    \endgroup
   }%
   \par
 \endgroup
}%
\def\@affilID@count#1#2#3{%
  \def\@tempifx{#3}%
  \@ifx{\@tempifx#1}{%
    \def\AU@opr{\@author@count#2}%
  }{%
    \let\AU@opr \@author@gobble
  }%
  \let\CO@opr \@collaboration@err
}%
\def\@affilID@match#1#2{%
 \def\@tempifx{#2}%
 \@ifx{\@tempifx#1}{%
   \let\AU@opr \@author@present
 }{%
   \let\AU@opr \@author@gobble
 }%
  \let\CO@opr \@collaboration@err
}%
\def\do@affil@fromgroup#1#2{%
  \@ifx{\relax#2}{}{%
\count@#2\relax
\@ifnum{\count@=\z@}{}{#1}%
    \do@affil@fromgroup#1%
  }%
}%
\def\@affil@group#1#2#3{%
  \@ifnum{\count@=#1\relax}{%
    \@ifx{#2\blankaffiliation}{}{%
#2%
\@if@empty{#3}{}{%
\frontmatter@footnote{#3}%
}%
\after@address
    }%
    \advance\@tempcnta\m@ne
  }{}%
}%
\def\@author@present@group#1#2{%
  \gdef\comma@space{\gdef\comma@space{\textsuperscript{,\,}}}%
  \doauthor{#1}{#2}{}%
  \advance\@tempcnta\m@ne
}%
\def\@pacs@produce#1{%
 \showPACS@sw{%
  \begingroup
   \frontmatter@PACS@format
   \@pacs@name#1\par
  \endgroup
 }{%
  \@if@empty{#1}{}{%
   \class@warn{\PACS@warn}%
  }%
 }%
}%
\def\PACS@warn{If you want your PACS to appear in your output, use document class option showpacs}%
\def\@keywords@produce#1{%
 \showKEYS@sw{%
  \begingroup
   \frontmatter@keys@format
   \@keys@name#1\par
  \endgroup
 }{%
  \@if@empty{#1}{}{%
   \class@warn{If you want your keywords to appear in your output, use document class option showkeys}%
  }%
 }%
}%
\def\frontmatter@footnote@produce@footnote{%
 \let\@TBN@opr\present@FM@footnote
 \@FMN@list
 \global\let\@FMN@list\@empty
}%
\def\present@FM@footnote#1#2{%
 \begingroup
  \csname c@\@mpfn\endcsname#1\relax
  \def\@thefnmark{\thempfn}%
  \frontmatter@footnotetext{#2}%
 \endgroup
}%
\def\frontmatter@footnote@produce@endnote{%
}%
\def\frontmatter@thefootnote{%
 \altaffilletter@sw{\@alph}{\@fnsymbol}\c@footnote
}%
\def\frontmatter@makefnmark{%
 \@textsuperscript{%
  \normalfont\@thefnmark
 }%
}%
\long\def\frontmatter@makefntext#1{%
 \parindent 1em
 \noindent
 \Hy@raisedlink{\hyper@anchorstart{frontmatter@\the\c@footnote}\hyper@anchorend}%
 \@makefnmark
 #1%
}%
\def\frontmatter@setup{}%
\def\frontmatter@RRAPformat#1{%
 \removelastskip
 \begingroup
  \frontmatter@RRAP@format
  #1\par
 \endgroup
}%
\def\punct@RRAP{; }%
\def\produce@RRAP#1{%
  \@if@empty{#1}{}{%
   \@ifvmode{\leavevmode}{\unskip\punct@RRAP\ignorespaces}%
   #1%
  }%
}%
\def\frontmatter@authorformat{}%
\def\frontmatter@above@affilgroup{}%
\def\frontmatter@above@affiliation{}%
\def\frontmatter@above@affiliation@script{}%
\def\frontmatter@affiliationfont{}%
\def\frontmatter@RRAP@format{}%
\def\frontmatter@PACS@format{}%
\def\frontmatter@keys@format{}%
\def\frontmatter@finalspace{\addvspace{18\p@}}
\def\frontmatter@addressnewline{%
  \@ifhmode{\skip@\lastskip\unskip\unpenalty\break\hskip\skip@}{}%
  % was: \vskip-.5ex
}%
\def\frontmatter@preabstractspace{5.5\p@}
\def\frontmatter@postabstractspace{6.5\p@}
\def\aftermaketitle@chk#1{%
  \@ifx{\maketitle\relax}{%
    \class@err{\protect#1 must be used before \protect\maketitle}%
  }{}%
}%
\def\ps@titlepage{\ps@empty}%
\def\volumeyear#1{\gdef\@volumeyear{#1}}%
\def\@volumeyear{}%
\def\volumenumber#1{\gdef\@volumenumber{#1}}%
\def\@volumenumber{}%
\def\issuenumber#1{\gdef\@issuenumber{#1}}%
\def\@issuenumber{}%
\def\eid#1{\gdef\@eid{#1}}%
\def\@eid{}%
\def\startpage#1{\gdef\@startpage{#1}\c@page#1\relax}%
\def\@startpage{\pageref{FirstPage}}%
\def\endpage#1{\gdef\@endpage{#1}}%
\def\@endpage{\pageref{LastPage}}%
\def\print@toc#1{%
 \begingroup
  \expandafter\section
  \expandafter*%
  \expandafter{%
              \csname#1name\endcsname
              }%
  \let\appendix\appendix@toc
  \@starttoc{#1}%
 \endgroup
}%
\def\appendix@toc{}%
\def\labelenumi{\theenumi.}
\def\theenumi{\arabic{enumi}}
\def\labelenumii{(\theenumii)}
\def\theenumii{\alph{enumii}}
\def\p@enumii{\theenumi}
\def\labelenumiii{\theenumiii.}
\def\theenumiii{\roman{enumiii}}
\def\p@enumiii{\theenumi(\theenumii)}
\def\labelenumiv{\theenumiv.}
\def\theenumiv{\Alph{enumiv}}
\def\p@enumiv{\p@enumiii\theenumiii}
\def\labelitemi{\textbullet}
\def\labelitemii{\normalfont\bfseries\textendash}
\def\labelitemiii{\textasteriskcentered}
\def\labelitemiv{\textperiodcentered}
\pagenumbering{arabic}
\setcounter{topnumber}{2}
\def\topfraction{.9}
\setcounter{bottomnumber}{1}
\def\bottomfraction{.9}
\setcounter{totalnumber}{3}
\def\textfraction{.1}
\def\floatpagefraction{.9}
\setcounter{dbltopnumber}{2}
\def\dbltopfraction{.9}
\def\dblfloatpagefraction{.9}
\newenvironment{verse}{%
  \let\\=\@centercr
  \list{}{%
    \itemsep\z@ \itemindent -1.5em\listparindent \itemindent
    \rightmargin\leftmargin\advance\leftmargin 1.5em}\item[]%
}{%
 \endlist
}%
\newenvironment{quotation}{%
  \list{}{%
    \listparindent 1.5em
    \itemindent\listparindent
    \rightmargin\leftmargin \parsep \z@ \@plus\p@}\item[]%
}{%
 \endlist
}%
\newenvironment{quote}{%
  \list{}{%
  \rightmargin\leftmargin}\item[]%
}{%
 \endlist
}%
\def\descriptionlabel#1{%
  \hspace\labelsep \normalfont\bfseries #1%
}%
\newenvironment{description}{%
 \list{}{%
    \labelwidth\z@ \itemindent-\leftmargin
    \let\makelabel\descriptionlabel
 }%
}{%
 \endlist
}%
\newcounter{part}
\newcounter{section}
\newcounter{subsection}[section]
\newcounter{subsubsection}[subsection]
\newcounter{paragraph}[subsubsection]
\newcounter{subparagraph}[paragraph]
\def\thepart          {\Roman{part}} %
\def\thesection       {\Roman{section}}
\def\p@section        {}
\def\thesubsection    {\Alph{subsection}}
\def\p@subsection     {\thesection\,}
\def\thesubsubsection {\arabic{subsubsection}}
\def\p@subsubsection  {\thesection\,\thesubsection\,}
\def\theparagraph     {\alph{paragraph}}
\def\p@paragraph      {\thesection\,\thesubsection\,\thesubsubsection\,}
\def\thesubparagraph  {\arabic{subparagraph}}
\def\p@subparagraph   {\thesection\,\thesubsection\,\thesubsubsection\,\theparagraph\,}
\newenvironment{acknowledgments}{%
 \acknowledgments@sw{%
  \section*{\acknowledgmentsname}%
 }{%
  \par
 \phantomsection
  \addcontentsline{toc}{section}{\protect\numberline{}\acknowledgmentsname}%
 }%
}{%
 \par
}%
\@booleantrue\acknowledgments@sw
\newenvironment{acknowledgements}{%
 \replace@environment{acknowledgements}{acknowledgments}%
}{%
 \endacknowledgments
}%
\def\part{\par
   \addvspace{4ex}%
   \@afterindentfalse
   \secdef\@part\@spart}%
\def\@part[#1]#2{%
 \@ifnum{\c@secnumdepth >\m@ne}{%
        \refstepcounter{part}%
        \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
 }{%
      \addcontentsline{toc}{part}{#1}%
 }%
 \begingroup
    \parindent \z@ \raggedright
    \interlinepenalty\@M
    \@ifnum{\c@secnumdepth >\m@ne}{%
      \Large \bf \partname~\thepart%
      \par\nobreak
    }{}%
    \huge \bf
    #2%
    \markboth{}{}\par
 \endgroup
   \nobreak
   \vskip 3ex
   \@afterheading
}%
\def\@spart#1{{\parindent \z@ \raggedright
    \interlinepenalty\@M
    \huge \bf
    #1\par}
    \nobreak
    \vskip 3ex
    \@afterheading}
\def\section{%
  \@startsection
    {section}%
    {1}%
    {\z@}%
    {0.8cm \@plus1ex \@minus .2ex}%
    {0.5cm}%
    {\normalfont\small\bfseries}%
}%
\def\subsection{%
  \@startsection
    {subsection}%
    {2}%
    {\z@}%
    {.8cm \@plus1ex \@minus .2ex}%
    {.5cm}%
    {\normalfont\small\bfseries}%
}%
\def\subsubsection{%
  \@startsection
    {subsubsection}%
    {3}%
    {\z@}%
    {.8cm \@plus1ex \@minus .2ex}%
    {.5cm}%
    {\normalfont\small\itshape}%
}%
\def\paragraph{%
  \@startsection
    {paragraph}%
    {4}%
    {\parindent}%
    {\z@}%
    {-1em}%
    {\normalfont\normalsize\itshape}%
}%
\def\subparagraph{%
  \@startsection
    {subparagraph}%
    {5}%
    {\parindent}%
    {3.25ex \@plus1ex \@minus .2ex}%
    {-1em}%
    {\normalfont\normalsize\bfseries}%
}%
\def\theequation{%
 \theequation@prefix\arabic{equation}%
}%
\def\theequation@prefix{}%
\AtEndOfClass{%
 \expandafter\input\expandafter{\@pointsize pt\substyle@ext}%
}%
\setcounter{secnumdepth}{4}
\lineskip 1pt
\normallineskip 1pt
\def\baselinestretch{1}%
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
\arraycolsep 3pt
\tabcolsep 2pt
\arrayrulewidth .4pt
\doublerulesep 2pt
\skip\@mpfootins = 0pt
\fboxsep  = 3.0pt
\fboxrule = 0.4pt
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{10\p@}
\setlength\belowcaptionskip{2\p@}
\long\def\@makecaption#1#2{%
  \par
  \vskip\abovecaptionskip
  \begingroup
   \small\rmfamily
   \sbox\@tempboxa{%
    \let\\\heading@cr
    #1: #2%
   }%
   \@ifdim{\wd\@tempboxa >\hsize}{%
    \begingroup
     \samepage
     \flushing
     \let\footnote\@footnotemark@gobble
     #1: #2\par
    \endgroup
   }{%
     \global \@minipagefalse
     \hb@xt@\hsize{\hfil\unhbox\@tempboxa\hfil}%
   }%
  \endgroup
  \vskip\belowcaptionskip
}%
\def\@footnotemark@gobble{%
 \@footnotemark
 \@ifnextchar[{\@gobble@opt@i}{\@gobble}%
}%
\def\@gobble@opt@i[#1]#2{}%
\def\@mpmakefntext#1{%
 \flushing
 \parindent=1em
 \noindent
 \hb@xt@1em{\hss\@makefnmark}%
 #1%
}%
\newcounter{figure}
\renewcommand \thefigure {\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\appdef\class@documenthook{%
 \do@if@floats{figure}{.fgx}%
}%
\appdef\class@enddocumenthook{%
 \printfigures\relax
}%
\newcommand\printfigures{%
 \@ifstar{\true@sw}{\floatp@sw{\true@sw}{\false@sw}}%
 {%
  \print@float{figure}{\oneapage}%
 }{%
  \print@float{figure}{}%
 }%
}%
\appdef\@xfloat@prep{%
 \appdef\@parboxrestore{\centering}%
 \let\@makefnmark\@makefnmark@latex
 \let\@footnotemark\@footnotemark@latex
}%
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
\newcounter{table}
\renewcommand\thetable{\@Roman\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}
\def\table@hook{\small}%
\def\squeezetable{\def\table@hook{\scriptsize}}%
\appdef\@floatboxreset{\table@hook}%
\appdef\class@documenthook{%
 \floats@sw{}{%
  \let@environment{longtable@float}{longtable}%
  \let@environment{longtable}{longtable@write}%
  \let@environment{longtable*@float}{longtable*}%
  \let@environment{longtable*}{longtable*@write}%
  \let@environment{turnpage@float}{turnpage}%
  \let@environment{turnpage}{turnpage@write}%
 }%
 \do@if@floats{table}{.tbx}%
}%
\appdef\class@enddocumenthook{%
 \printtables\relax
}%
\newenvironment{longtable@write}{%
 \write@@float{longtable}{table}%
}{%
 \endwrite@float
}%
\newenvironment{longtable*@write}{%
 \write@@float{longtable*}{table}%
}{%
 \endwrite@float
}%
\newenvironment{turnpage@write}{%
 \immediate\write\tablewrite{\string\begin{turnpage}}%
}{%
 \immediate\write\tablewrite{\string\end{turnpage}}%
}%
\newcommand\printtables{%
 \begingroup
  \let@environment{longtable}{longtable@float}%
  \let@environment{longtable*}{longtable*@float}%
  \let@environment{turnpage}{turnpage@anchored}%
  \prepdef\longtable{\trigger@float@par}%
  \expandafter\prepdef\csname longtable*\endcsname{\trigger@float@par}%
  \expandafter\prepdef\csname table@floats\endcsname{%
   \onecolumngrid@push
  }%
  \expandafter\appdef\csname endtable@floats\endcsname{%
   \onecolumngrid@pop
  }%
  \@ifstar{\true@sw}{\floatp@sw{\true@sw}{\false@sw}}%
  {%
   \print@float{table}{\oneapage}%
  }{%
   \print@float{table}{}%
  }%
 \endgroup
}%
\newenvironment{turnpage@anchored}{%
 \onecolumngrid@push
 \setbox\z@\vbox to\textwidth\bgroup
  \columnwidth\textheight
}{%
  \vfil
 \egroup
 \rotatebox{90}{\box\z@}%
 \onecolumngrid@pop
}%
\def\endtabular@hook{}%
\appdef\class@documenthook{%
 \@ifpackageloaded{dcolumn}{%
  \newcolumntype{d}{D{.}{.}{-1}}%
 }{}%
}%
\def\toprule{\hline\hline}%
\def\colrule{\hline}%
\def\botrule{\hline\hline}%
\newenvironment{ruledtabular}{%
 \def\array@default{v}%
 \appdef\tabular@hook{\def\@halignto{to\hsize}}%
 \let\tableft@skip@default\tableft@skip
 \let\tableft@skip\tableft@skip@float
 \let\tabmid@skip@default\tabmid@skip
 \let\tabmid@skip\tabmid@skip@float
 \let\tabright@skip@default\tabright@skip
 \let\tabright@skip\tabright@skip@float
 \let\array@row@pre@default\array@row@pre
 \let\array@row@pre\array@row@pre@float
 \let\array@row@pst@default\array@row@pst
 \let\array@row@pst\array@row@pst@float
 \appdef\array@row@rst{%
  \let\array@row@pre\array@row@pre@default
  \let\array@row@pst\array@row@pst@default
  \let\tableft@skip\tableft@skip@default
  \let\tabmid@skip\tabmid@skip@default
  \let\tabright@skip\tabright@skip@default
  \appdef\tabular@hook{\let\@halignto\@empty}%
 }%
}{%
}%
\def\@makefntext#1{%
  \def\baselinestretch{1}%
  \reset@font\footnotesize
  \parindent 1em%
  \noindent
  \hb@xt@1.8em{%
   \Hy@raisedlink{\hyper@anchorstart{footnote@\the\c@footnote}\hyper@anchorend}%
   \hss\@makefnmark
  }%
  #1\par
}%
\newenvironment{thebibliography}{}{}%
\let\@listi\@empty
\RequirePackage[sort&compress]{natbib}[1999/05/28 7.0]%
\let\NAT@thebibliography\thebibliography
\let\NAT@endthebibliography\endthebibliography
\let\bibliographystyle@latex\bibliographystyle
\renewcommand\NAT@citesuper[3]{%
 \ifNAT@swa
  \leavevmode
  \unskip
  \textsuperscript{#1}%
  \if*#3*\else\ (#3)\fi
 \else
  #1%
 \fi
 \endgroup
}%
\def\@makefnmark@cite{\begingroup\NAT@swatrue\@cite{{\@thefnmark}}{}{}}%
\let\bibliographystyle@latex\bibliographystyle
\def\bibliographystyle{\def\@bibstyle}%
\def\bibsection{%
 \@ifx@empty\refname{%
  \par
 }{%
  \let\@hangfroms@section\@hang@froms
  \section*{\refname}%
  \@nobreaktrue
 }%
}%
\let\NATx@bibsetnum\NAT@bibsetnum
\def\NAT@bibsetnum#1{%
 \setlength{\topsep}{\z@}%
 \NATx@bibsetnum{\ref{LastBibItem}}%
}%
\let\NATx@bibsetup\NAT@bibsetup
\def\NAT@bibsetup{%
 \setlength{\labelwidth}{\z@}%
 \setlength{\labelsep}{\z@}%
 \setlength{\itemindent}{\z@}%
 \setlength{\listparindent}{\z@}%
 \setlength{\topsep}{\z@}%
 \setlength{\parsep}{\z@}%
 \NATx@bibsetup
}%
\let\bibpreamble\@empty
\def\newblock{\ }%
\let\NATx@bibnumfmt\bibnumfmt
\def\bibnumfmt{\place@bibnumber}%
\AtEndOfClass{%
 \appdef\class@documenthook{%
  \@ifxundefined\place@bibnumber{%
   \let\place@bibnumber\NATx@bibnumfmt
  }{}%
 }%
}%
\let\NAT@citeyear\citeyear
\def\NAT@onlinecite#1{%
 \begingroup
  \let\@cite\NAT@citenum
  \let\NAT@mbox\mbox
  \let\citeyear\NAT@citeyear
  \def\NAT@space{ }%
  \cite{#1}%
 \endgroup
}%
\appdef\class@inithook{%
 \@ifxundefined\@onlinecite{\def\@onlinecite{\citealp}}{}%
 \@ifxundefined\@textcite{\def\@textcite{\citet}}{}%
}%
\DeclareRobustCommand\onlinecite{\@onlinecite}%
\DeclareRobustCommand\textcite{\@textcite}%
\renewenvironment{thebibliography}[1]{%
 \NAT@thebibliography{#1}%
 \@clubpenalty\clubpenalty
 \let\@TBN@opr\present@bibnote
 \@FMN@list
}{%
 \@endnotesinbib
 \edef\@currentlabel{\arabic{NAT@ctr}}%
 \label{LastBibItem}%
 \NAT@endthebibliography
 \global\let\auto@bib\@empty
}%
\def\present@bibnote#1#2{%
 \item[%
  \textsuperscript{%
   \normalfont
   \Hy@raisedlink{\hyper@anchorstart{frontmatter@#1}\hyper@anchorend}%
   \begingroup
    \c@footnote#1\relax
    \frontmatter@thefootnote
   \endgroup
  }%
 ]#2\par
}%
\def\place@bibnumber@sup#1{\textsuperscript{#1}}%
\def\place@bibnumber@inl#1{[#1]}%
%% \begin{macro}{\@cite}
%% \begin{macro}{\nocite}
%% For the sake of Bib\TeX, we must execute the \cmd\bibliographystyle\ command
%% (putting a \cmd\bibstyle\ command into the auxiliary file)
%% somewhere before the first \cmd\cite\ command (or any other like command),
%% which puts a \cmd\citation\ command into the auxiliary file.
%%
%% Therefore, we headpatch the two procedures in \classname{natbib} that can result in
%% the appearance in the auxiliary file of a \cmd\citation\ command.
%% These two procedures are \cmd\@cite\ and \cmd\nocite.
%%
%% Not!
\def\write@bibliographystyle{%
 \@ifxundefined\@bibstyle{}{%
  \expandafter\bibliographystyle@latex\expandafter{\@bibstyle}%
 }%
 \global\let\write@bibliographystyle\relax
}%
\AtEndDocument{\write@bibliographystyle}%
%% \end{macro}
%% \end{macro}
\def\super@cite[#1][#2]#3{%
 \leavevmode
 \skip@\lastskip
 \unskip
 \begingroup
  \def\@tempa{[#1][#2]{#3}}%
  \super@cite@
}%
\def\super@cite@{%
 \futurelet\@let@token\super@cite@@
}%
\long\def\super@cite@@#1{%
 \@ifx{\@let@token\@sptoken}{%
  \aftergroup\@sptoken\false@sw
 }{%
  \@ifx{\@let@token.}{\true@sw}{%
   \@ifx{\@let@token,}{\true@sw}{%
    \@ifx{\@let@token'}{\true@sw}{%
     \@ifx{\@let@token)}{\true@sw}{%
      \false@sw
     }%
    }%
   }%
  }%
 }%
 {%
  #1%
  \super@cite@
 }{%
  \expandafter\@citex@NAT\@tempa\endgroup
  #1%
 }%
}%
\AtBeginDocument{%
 \citeautoscript@sw{%
  \@ifx{\place@bibnumber\place@bibnumber@sup}{%
   \let\@citex@NAT\@citex\let\@citex\super@cite
  }{}%
 }{}%
}%
\@ifx{\chapter\relax}{\let\chapter\@undefined}{}%
\def\mini@note{\save@note\mini@notes}%Implicit #1
\def\save@note#1#2{%
  \stepcounter{footnote}%
  \protected@xdef\@thefnmark{\the\c@footnote}%
  \@footnotemark
  \expandafter\g@addto@macro
  \expandafter#1%
  \expandafter{%
  \expandafter \@@footnotetext
  \expandafter {\@thefnmark}{#2}%
              }%
}%
\long\def\@@footnotetext#1{\def\@thefnmark{#1}\@footnotetext}%
\let\mini@notes\@empty
\def\endnote{%
 \begingroup
  \aftergroup\@footnotemark
  \aftergroup\@endnotetext
  \@ifnextchar[{%
   \@xendnote
  }{%
   \stepcounter{footnote}%
   \protected@xdef\@tempa{\thefootnote}%
   \expandafter\@xendnote\expandafter[\@tempa]%
 }%
}%
\def\@xendnote[#1]{%
  \begingroup
   \c@footnote#1\relax
   \unrestored@protected@xdef\@thefnmark{endnote\thefootnote}%
  \endgroup
 \endgroup
}%
\def\@endnotemark{%
 \expandafter\cite\expandafter{\@thefnmark}%
}%
\def\endnote@ext{.end}%
\long\def\@endnotetext#1{%
  \@ifxundefined\@endnoteout{%
    \newwrite\@endnoteout
    \gdef\endnote@stream{\jobname\endnote@ext}%
    \immediate\openout\@endnoteout\endnote@stream\relax
  }{}%
  \begingroup
    \let\label\relax \let\index\relax \let\glossary\relax
    \let\cite \relax \let\ref  \relax \let\pageref \relax
    \let\(    \relax \let\)    \relax \let\\    \relax
    \let \protect \@unexpandable@protect
    \newlinechar`\^^M%
    \immediate\write\@endnoteout{\string\@doendnote{\@thefnmark}{#1}}%
  \endgroup
}%
\def\endnotetext{%
  \@ifnextchar[{%
   \@xendnotenext
  }{%
   \protected@xdef\@thefnmark{\thefootnote}%
   \@endnotetext
  }%
}%
\def\@xendnotenext[#1]{%
  \begingroup
     \c@footnote#1\relax
     \unrestored@protected@xdef\@thefnmark{\thefootnote}%
  \endgroup
  \@endnotetext
}%
\def\printendnotes{%
 \@ifx{\@endnoteout\@undefined}{}{%
  \begingroup
   \c@secnumdepth-\maxdimen
   \section{%
    \notesname
   }%
   \immediate\closeout\@endnoteout
   \global\let\@endnoteout\@undefined
   \begin{enumerate}%
   \makeatletter
   \input{\endnote@stream}%
   \end{enumerate}%
  \endgroup
 }%
}%
\appdef\class@enddocumenthook{%
 \printendnotes
}%
\mathchardef\firstnote@num\z@
\appdef\class@documenthook{%
 \global\c@footnote\firstnote@num
 \appdef\maketitle{\global\c@footnote\firstnote@num}%
}%
\def\set@firstnote#1{%
 \@ifnum{\firstnote@num=#1\relax}{}{%
  \class@warn@end{Endnote numbers changed: rerun LaTeX}%
 }%
 \immediate\write\@mainaux{%
   \global\chardef\string\firstnote@num#1\relax
 }%
}%
\def\@endnotesinbib{%
 \set@firstnote\z@
}%
\def\@endnotesinbibliography{%
 \expandafter\set@firstnote\expandafter{\the\c@NAT@ctr}%
 \@ifx{\@endnoteout\@undefined}{}{%
  \immediate\closeout\@endnoteout
  \global\let\@endnoteout\@undefined
  \makeatletter
  \input{\endnote@stream}%
 }%
}%
\appdef\class@documenthook{%
 \let\@makefnmark@latex\@makefnmark
 \let\@footnotemark@latex\@footnotemark
 \footinbib@sw{%
  \let\footnote\endnote
  \let\@footnotemark\@endnotemark
  \let\@makefnmark\@makefnmark@cite
  \let\@endnotesinbib\@endnotesinbibliography
  \appdef\class@enddocumenthook{\auto@bib}%
  \let\printendnotes\relax
 }{}%
}%
\def\auto@bib{%
 \@ifx@empty\@FMN@list{%
  \@ifxundefined\@endnoteout{%
   \false@sw
  }{%
   \true@sw
  }%
 }{%
  \true@sw
 }{%
  \begin{thebibliography}{}\end{thebibliography}%
 }%
}%
\def\@doendnote#1#2{\bibitem{#1}#2}
\appdef\class@documenthook{%
\@ifxundefined\@evenhead{%
\@ifxundefined\@oddhead{%
\@ifxundefined\@oddfoot{%
\@ifxundefined\@evenfoot{%
\preprintsty@sw{%
\ps@preprint
}{%
\ps@article
}%
}{}%
}{}%
}{}%
}{}%
 \pagenumbering{arabic}%
 \normalsize
 \tabbingsep \labelsep
 \leftmargin\leftmargini
 \labelwidth\leftmargin\advance\labelwidth-\labelsep
 \let\@listi\@listI
 \@listi
 \def\curr@envir{document}%
 \mark@envir{\curr@envir}%
 \twocolumn@sw{%
  \open@column@mlt\tw@
  \set@colht
  \@floatplacement
  \@dblfloatplacement
  \sloppy
  \@booleanfalse\raggedcolumn@sw
  \let\set@listindent\set@listindent@
 }{%
  \open@column@one\@ne
  \set@colht
  \@floatplacement
  \@dblfloatplacement
  \sloppy
  \let\set@listindent\set@listindent@parshape
 }%
}%
\def\appendix{%
 \par
 \setcounter{section}\z@
 \setcounter{subsection}\z@
 \setcounter{subsubsection}\z@
 \def\thesubsection{\arabic{subsection}}%
 \def\thesubsubsection{\alph{subsubsection}}%
 \@addtoreset{equation}{section}%
 \def\theequation@prefix{\thesection}%
 \addtocontents{toc}{\protect\appendix}%
 \@ifstar{%
  \def\thesection{\unskip}%
  \def\theequation@prefix{A.}%
 }{%
  \def\thesection{\Alph{section}}%
 }%
}%
\def\title@column@default#1{%
 \minipagefootnote@init
 #1%
 \minipagefootnote@foot
}%
\def\close@column@default{%
 \newpage
}%
\def\galley@outdent{\rightmargin-\columnwidth\advance\rightmargin-\columnsep}%
\let\widetext@outdent\@empty
\newenvironment{widetext@galley}{%
  \list{}{%
    \topsep        \z@skip
    \listparindent \parindent
    \itemindent    \parindent
    \leftmargin    \z@
    \parsep        \z@\@plus\p@
    \widetext@outdent
    \relax
  }%
  \item\relax
}{
  \endlist
}%
\def\title@column@grid#1{%
 \minipagefootnote@init
  \onecolumngrid
  \begingroup
   \let\@footnotetext\frontmatter@footnotetext
   \let\set@footnotewidth\set@footnotewidth@ii
   #1%
  \endgroup
  \twocolumngrid
 \minipagefootnote@foot
}%
\def\close@column@grid{%
 \balancelastpage@sw{%
  \onecolumngrid
 }{}%
}%
\newenvironment{widetext@grid}{%
  \par\ignorespaces
  \setbox\widetext@top\vbox{%
   \vskip15\p@
   \hb@xt@\hsize{%
    \leaders\hrule\hfil
    \vrule\@height6\p@
   }%
   \vskip6\p@
  }%
  \setbox\widetext@bot\hb@xt@\hsize{%
    \vrule\@depth6\p@
    \leaders\hrule\hfil
  }%
  \onecolumngrid
  \dimen@\ht\widetext@top\advance\dimen@\dp\widetext@top
  \cleaders\box\widetext@top\vskip\dimen@
  \let\set@footnotewidth\set@footnotewidth@ii
}{%
  \par
  \setbox\widetext@bot\vbox{%
   \hb@xt@\hsize{\hfil\box\widetext@bot}%
   \vskip14\p@
  }%
  \dimen@\ht\widetext@bot\advance\dimen@\dp\widetext@bot
  \cleaders\box\widetext@bot\vskip\dimen@
  \twocolumngrid\global\@ignoretrue
  \@endpetrue
}%
\newbox\widetext@top
\newbox\widetext@bot
\appdef\class@documenthook{%
 \twocolumn@sw{%
  \let@environment{widetext}{widetext@grid}%
  \let\title@column\title@column@grid
  \let\close@column\close@column@grid
 }{%
  \let@environment{widetext}{widetext@galley}%
  \preprintsty@sw{%
  }{%
   \galley@sw{%
    \let\widetext@outdent\galley@outdent
   }{%
   }%
  }%
 }%
}%
\appdef\class@inithook{%
 \@ifxundefined\title@column{%
  \let\title@column\title@column@default
 }{}%
}%
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
\def\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
\def\notesname{Notes}
\def\partname{Part}
\def\tocname{Contents}
\def\lofname{List of Figures}
\def\lotname{List of Tables}
\def\refname{References}
\def\indexname{Index}
\def\figurename{FIG.}
\def\figuresname{Figures}%
\def\tablename{TABLE}
\def\tablesname{Tables}%
\def\abstractname{Abstract}
\def\appendixesname{Appendixes}%
\def\appendixname{Appendix}%
\def\acknowledgmentsname{Acknowledgments}
\def\journalname{??}
\def\copyrightname{??}
\def\andname{and}
\def\@pacs@name{PACS numbers: }%
\def\@keys@name{Keywords: }%
\def\ppname{pp}
\def\numbername{number}
\def\volumename{volume}
\def\Dated@name{Dated: }%
\def\Received@name{Received }%
\def\Revised@name{Revised }%
\def\Accepted@name{Accepted }%
\def\Published@name{Published }%
\def\address{\replace@command\address\affiliation}%
\def\altaddress{\replace@command\altaddress\altaffiliation}%
\newenvironment{references}{%
 \class@warn@end{The references environment is not supported; use thebibliography instead.}
 \gdef\references{\thebibliography{}}\references
}{%
 \endthebibliography
}%
\def\draft{%
 \class@warn@end{Command \string\draft\space is obsolete;^^JInvoke option draft instead.}%
 \@booleantrue\draft@sw
}%
\def\tighten{%
 \class@warn@end{Command \string\tighten\space is obsolete;^^JInvoke option tightenlines instead.}%
 \@booleantrue\tightenlines@sw
}%
\def\tableline{%
 \noalign{%
  \class@warn@end{Command \string\tableline\space is obsolete;^^JUse \string\colrule\space instead.}%
  \global\let\tableline\colrule
 }%
 \tableline
}%
\def\case{\replace@command\case\frac}%
\def\slantfrac{\replace@command\slantfrac\frac}%
\def\tablenote{\replace@command\tablenote\footnote}%
\def\tablenotemark{\replace@command\tablenotemark\footnotemark}%
\def\tablenotetext{\replace@command\tablenotetext\footnotetext}%
\DeclareRobustCommand\REV@text[1]{%
 \relax
 \ifmmode
  \mathchoice
   {\hbox{{\everymath{\displaystyle     }#1}}}%
   {\hbox{{\everymath{\textstyle        }#1}}}%
   {\hbox{{\everymath{\scriptstyle      }\let\f@size\sf@size\selectfont#1}}}%
   {\hbox{{\everymath{\scriptscriptstyle}\let\f@size\ssf@size\selectfont#1}}}%
  \glb@settings
 \else
  \mbox{#1}%
 \fi
}%
\DeclareRobustCommand\REV@bbox[1]{%
 \relax
 \ifmmode
  \mathchoice
   {\hbox{{\everymath{\displaystyle     }\boldmath$#1$}}}%
   {\hbox{{\everymath{\textstyle        }\boldmath$#1$}}}%
   {\hbox{{\everymath{\scriptstyle      }\boldmath$#1$}}}%
   {\hbox{{\everymath{\scriptscriptstyle}\boldmath$#1$}}}%
  \glb@settings
 \else
  \mbox{#1}%
 \fi
}%
\DeclareRobustCommand\REV@bm[1]{%
 \class@warn@end{To use \string\bm, please load the bm package!}%
 \global\let\bm\relax
}%
\def\FL{\obsolete@command\FL}%
\def\FR{\obsolete@command\FR}%
\def\narrowtext{\obsolete@command\narrowtext}%
\def\mediumtext{\obsolete@command\mediumtext}%
\newenvironment{quasitable}{%
 \let@environment{tabular}{longtable}%
}{%
}%
\appdef\class@inithook{%
 \@ifxundefined\text{\let\text\REV@text}{}%
 \@ifxundefined\bm{\let\bm\REV@bm}{}%
 \providecommand\bibinfo[2]{#2}%
 \providecommand\eprint[2][]{#2}%
 \providecommand\url[1]{#1}%
}%
\def\bbox#1{%
 \class@warn@end{\string\bbox\space is obsolete,^^Jload the bm package and use \string\bm\space instead.}%
 \global\let\bbox\relax
}%
\newenvironment{mathletters}{%
 \class@warn@end{Environment {mathletters} is obsolete;^^Jload the amsmath package and use {subequations}!}%
 \global\let\mathletters\@empty
}{%
}%
\def\eqnum#1{%
 \class@warn@end{\string\eqnum\space is obsolete, load the amsmath package and use \string\tag!}%
 \global\let\eqnum\@gobble
}%
\RequirePackage{revsymb}%
\appdef\class@inithook{\revsymb@inithook}%
\AtEndOfClass{%
 \InputIfFileExists{\jobname.rty}{}{}%
}%
\AtEndOfClass{%
 \@ifxundefined\secnumarabic@sw{\@booleanfalse\secnumarabic@sw}{}%
 \secnumarabic@sw{%
  \def\p@section        {}%
  \def\thesection       {\arabic{section}}%
  \def\p@subsection     {}%
  \def\thesubsection    {\thesection.\arabic{subsection}}%
  \def\p@subsubsection  {}%
  \def\thesubsubsection {\thesubsection.\arabic{subsubsection}}%
 }{}%
}%
\AtEndOfClass{%
 \appdef\class@documenthook{%
  \defaults@hook
 }%
}%
\@ifxundefined\defaults@hook{\let\defaults@hook\@empty}{}%
\endinput
%%
%% End of file `revtex4.cls'.
