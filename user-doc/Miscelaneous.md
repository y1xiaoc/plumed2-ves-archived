\page Miscelaneous Miscelaneous

- \subpage comments
- \subpage ContinuationLines
- \subpage VimSyntax
- \subpage includes
- \subpage load
- \subpage degub
- \subpage exchange-patterns
- \subpage mymodules
- \subpage misc

\page comments Comments

If you are an organised sort of person who likes to remember what the hell you were trying to do when you ran a 
particular simulation you might find it useful to put comments in your input file.  In PLUMED you can do this as 
comments can be added using a # sign.  On any given line everything after the # sign is ignored so 
erm... yes add lines of comments or trailing comments to your hearts content as shown below (using Shakespeare is optional):

\plumedfile
# This is the distance between two atoms:
DISTANCE ATOM=1,2 LABEL=d1
UPPER_WALLS ARG=d1 AT=3.0 KAPPA=3.0 LABEL=Snout # In this same interlude it doth befall.
# That I, one Snout by name, present a wall.
\endplumedfile
(see \ref DISTANCE and \ref UPPER_WALLS)

An alternative to including comments in this way is to use the command \subpage ENDPLUMED.  Everything in the PLUMED input after this
keyword will be ignored.

\page ContinuationLines Continuation lines

If your input lines get very long then editing them using vi and other such text editors becomes a massive pain in the arse.  
We at PLUMED are aware of this fact and thus have provided a way of doing line continuations so as to make your life that much 
easier - aren't we kind?  Well no not really, we have to use this code too.  Anyway, you can do continuations by using the "..." syntax
as this makes this: 

\plumedfile
DISTANCES ATOMS1=1,300 ATOMS2=1,400 ATOMS3=1,500 LABEL=dist
\endplumedfile
(see \ref DISTANCES)

equivalent to this:

\plumedfile
DISTANCES ...
  LABEL=dist
# we can also insert comments here
  ATOMS1=1,300
# multiple kewords per line are allowed
  ATOMS2=1,400 ATOMS3=1,500
#empty lines are also allowed

... DISTANCES
\endplumedfile

Notice that the closing `...` is followed by the word `DISTANCES`. This is optional, but might be
useful to find more easily which is the matching start of the statement. The following is equally correct
\plumedfile
DISTANCES ...
  LABEL=dist
# we can also insert comments here
  ATOMS1=1,300
# multiple kewords per line are allowed
  ATOMS2=1,400 ATOMS3=1,500
#empty lines are also allowed

...
\endplumedfile

Notice that PLUMED makes a check that the word following the closing `...` is actually identical to
the first word in the line with the first `...`. If not, it will throw an error.
Also notice that you might put more than one word in the first line. E.g.
\plumedfile
DISTANCES LABEL=dist ...
# we can also insert comments here
  ATOMS1=1,300
# multiple kewords per line are allowed
  ATOMS2=1,400 ATOMS3=1,500
#empty lines are also allowed
...
\endplumedfile
or, equivalently,
\plumedfile
dist: DISTANCES ...
# we can also insert comments here
  ATOMS1=1,300
# multiple kewords per line are allowed
  ATOMS2=1,400 ATOMS3=1,500
#empty lines are also allowed
...  
\endplumedfile


\page VimSyntax Using VIM syntax file

For the impatients:
- Add the following to your .vimrc file:
\verbatim
" This allows including the proper PLUMED syntax file:
:let &runtimepath.=','.$PLUMED_VIMPATH
" This makes autocompletion work in the expected way:
:set completeopt=longest,menuone
" This enables bindings of F2/F3/F4 to plumed specific commands:
:let plumed_shortcuts=1
\endverbatim
- When you open a PLUMED input file, you can enable syntax highlighting with:
\verbatim
:set ft=plumed
\endverbatim
This will also enable autocompletion. Use `<CTRL-X><CTRL-O>` to autocomplete a word.
- If you want to fold multiline statements, type
\verbatim
:setlocal foldmethod=syntax
\endverbatim
- While editing a plumed input file, you can use command `:PHelp` (or shortcut `<F2>`)
  to show in a split window a short help about the action defined in the line where the cursor is.
  Typing `:PHelp` again (or pushing `<F2>`) you will
  close that window. With `<CTRL-W><CTRL-W>` you go back and forth between the two windows.
- When you open a file starting with `#! FIELDS`, VIM will automatically understand it
  is a PLUMED outpt file (VIM filetype = plumedf) and will color fields and data columns with
  alternating colors. Typing `:PPlus` and `:PMinus` (or pushing `<F3>` and `<F4>`)
  you can move a highlighted column.

See below for more detailed instructions.

\par Configuration

When PLUMED is compiled, directories `help` and `syntax` will appear in `builddir/vim`.
They contain a VIM plugin that can be used to highlight proper PLUMED instructions in a PLUMED
input file and to quickly retrieve help.
There is also a file `builddir/vim/scripts.vim` that helps VIM in recognizing PLUMED output files.

\warning
Notice that these file do not appear if you are cross compiling.
In this case, you must copy the plugin files from another machine.

To make VIM aware of these files, you should copy them to your `$HOME/.vim` directory.
Later you can
enable plumed syntax with the command
\verbatim
:set ft=plumed
\endverbatim

If you work in an environment where several PLUMED versions are installed (e.g. using env modules),
we recommend the following procedure:
- Install PLUMED
- Add to your `.vimrc` file the following line:
\verbatim
:let &runtimepath.=','.$PLUMED_VIMPATH
\endverbatim

The modulefile provided with PLUMED should set the PLUMED_VIMPATH environemnt variable
to the proper path.
Thus, when working with a given PLUMED module loaded, you should be able to
enable to proper syntax by just typing
\verbatim
:set ft=plumed
\endverbatim
in VIM.
Notice that the variable PLUMED_VIMPATH is also set in the `sourceme.sh` script in the build directory.
This, if you modify your `.vimrc` file as suggested, you will be able to use the correct syntax both
when using an installed PLUMED and when running from a just compiled copy.

If you are tired of typing `:set ft=plumed`, you can use a modeline.
Add to your `.vimrc` file the following commands
\verbatim
:set modeline
:set modelines=5
\endverbatim
Then, at the beginning of your PLUMED input file, put the following comment:
\plumedfile
# vim:ft=plumed
d: DISTANCE ATOMS=1,2
RESTRAINT ARG=d AT=0.0 KAPPA=1.0
\endplumedfile
Now, every time you open this file, you will see it highlighted.

\par Syntax highlighting

The syntax file contains a definition of all possible PLUMED actions and keywords.
It is designed to allow for a quick validation of the PLUMED input file before running it.
As such, all the meaningful words in the input should be highlighted:
- Valid action names (such as `METAD`) and labels (such as `metad:` or `LABEL=metad`) will be
  highlighted in the brightest way (`Type` in VIM). Those are the most important words.
- Keyword and flag names (such as `ATOMS=` or `COMPONENTS` when part of the action \ref DISTANCE) will be highlighted with a different color
  (`Statement` in VIM).
- Values provided by users (such as the number of the atoms following `ATOMS=`) will be highlighted with a different color
  (`String` in VIM).
- Comments (see \ref comments) will be highlighted as comments (`Comment` in VIM).

If you see something that is not highlighted and appears in black, this is likely going to result in an error at runtime.
Think of this as a sort of preliminary spell-check.
For this checks to be effective, we recommend to use a syntax file generated with
exactly the same version of PLUMED that you are using.
In case you find that parts of an input file that is valid are not highlighted, then please
report it as a bug.
On the contrary, you cannot expect the VIM syntax file to recognize all possible errors
in a PLUMED input. Thus, a file for  which the highlighting looks correct might still contain errors.

\par Multi-line folding

Notice that syntax highlighting also allow VIM to properly fold multi-line actions.
Try to do the following:
- Open a PLUMED input file
- Enable PLUMED syntax
\verbatim
:set ft=plumed
\endverbatim
- Enable syntax-based folding
\verbatim
:setlocal foldmethod=syntax
\endverbatim

Now look at what happened to all the multi-line statements in PLUMED (i.e. those using
\ref ContinuationLines).  As you can see, they will be folded into single lines.
Folded lines can be expanded with `zo` and folded with `zc`. Look at VIM documentation to
learn more.
In case you want to use this feature, we suggest you to put both label
and action type on the first line of multi-line statements. E.g.
\plumedfile
m: METAD ...
  ARG=d
  HEIGHT=1.0
  SIGMA=0.5
  PACE=100
...
\endplumedfile
will be folded to
\verbatim
+--  6 lines: m: METAD ...------------------------------------------------------
\endverbatim
and
\plumedfile
METAD LABEL=m ...
  ARG=d
  HEIGHT=1.0
  SIGMA=0.5
  PACE=100
...
\endplumedfile
will be folded to
\verbatim
+--  6 lines: METAD LABEL=m ...-------------------------------------------------
\endverbatim
This will allow you to easily identify the folded lines by seeing the most important information,
that is the action type (`METAD`) and its label (`m`). This feature is convenient if
you want to browse files that contain a lot of actions defined on multiple lines.

\par Autocompletion

Another VIM feature that comes when you load PLUMED syntax is autocompletion of PLUMED
actions and keywords. Open your favorite PLUMED input file and set it to PLUMED syntax highlighting with
\verbatim
:set ft=plumed
\endverbatim
Now go into insert mode pressing `i` and type `DU` followed by `<CTRL+X><CTRL+O>`.
Here `<CTRL+X>` stands for autocompletion and `<CTRL+O>` for omnifunc autocompletion. You will see a short
menu listing the following actions
\verbatim
DUMPATOMS       
DUMPDERIVATIVES 
DUMPFORCES      
DUMPMASSCHARGE  
DUMPMULTICOLVAR 
DUMPPROJECTIONS 
\endverbatim
That is, all the actions starting with `DU`.
You can navigate it with up and down arrows so as to choose the
best match.

Notice that the default behavior of VIM is to use the first match by default.
In the first example (`DU<CTRL+X><CTRL+O`), it would be `DUMPATOMS`.
The following settings make it work as most of the people expect:
\verbatim
:set completeopt=longest,menuone
\endverbatim
With these settings, in the first example (`DU<CTRL+X><CTRL+O`) VIM will only complete up to the longest common part (`DUMP`).

As you can imagine,
if you use autocompletion after you have typed the word `DISTANCE` followed by a space you will see
a menu listing `LABEL=`, `COMPONENTS`, etc. Basically, all the keywords that are possibly used within a `DISTANCE` line
will be shown. This is very useful if you do not remember the exact name of the keywords associated with
a given action.

\par Quick help

You can also retrieve quick explanation of the input options for a specific action.
Try to do the following. Enable plumed syntax:
\verbatim
:set ft=plumed
\endverbatim
Then add the following line
\plumedfile
DISTANCE
\endplumedfile
Now, in normal mode, go with the cursor on the `DISTANCE` line and type
\verbatim
:PHelp
\endverbatim
A new split window should appear containing some documentation about the \ref DISTANCE collective variable.
You can go back and forth between the two windows with `<CTRL+W><CTRL+W>`, as usually in vim.
Notice that if you are in the help window and type `:PHelp` this window will be closed.

To make the navigation easier, you can add a shortcut in your .vimrc file. For example, adding:
\verbatim
: nmap <F2> : PHelp<CR>
\endverbatim
you should be able to open and close the manual hitting the F2 key.
This is done automatically in the PLUMED syntax file if you add `let plumed_shortcuts=1` to your
vimrc file.

\par Displaying output files

Most of the PLUMED output files look like this
\verbatim
#! FIELDS A B C
1 2 3
\endverbatim
This is useful since in the header you can see the name of the quantities that are printed in the
data lines. However, when you have an output file with many columns it might be a bit error prone to count them.
To simplify this, when PLUMED syntax for VIM is configured properly VIM should be able to:
- Detect that this file is a PLUMED output file with fields, automatically setting its type to `plumedf`. If not, just type
  `:set ft=plumedf`.
- Show this file with syntax highlighting to increase its readability.

Notice that the syntax file for the output files (`plumedf.vim`) is not the same one that is used for the PLUMED
input file (`plumed.vim`).

To make output files more readable, vim will show `FIELDS` and `SET` words in a different color,
and data columns with alternating colors (e.g. dark/light/dark/light).
The colors in the columns are consistent with those shown in the FIELD line.
In the example above, 1, 2, and 3 will be of the same color as A, B, and C respectively.
This should make it much easier to find which columns correspond to a given quantity.

It is also possible to highlight a specific field of the file. Typing
\verbatim
:5PCol
\endverbatim
you will highlight the fifth field. Notice that in the `FIELDS` line (the first line of the file)
the 7th word of the line will be highlighted, which is the one containing the name of the field.
This allows for easy matching of values shown
in the file and tags provided in the `FIELDS` line.
The highlighted column can be moved back and forth using `:PPlus` and `:PMinus`.
Adding a count to the command will move the highlighted column more. E.g. `:2PPlus` will move
the column to the right twice.

If you have a long output file, it might be convenient to split it with
`:split` so that one of the two windows will only show the header. The other window
can be used to navigate the file.


To make the navigation easier, you can add a shortcut in your .vimrc file. For example, adding:
\verbatim
: map <F3> :PMinus<CR>
: map <F4> :PPlus<CR>
\endverbatim
you should be able to move the highlight column using F3 and F4 buttons.
This is done automatically in the PLUMED syntax file if you add `let plumed_shortcuts=1` to your
vimrc file.

\page includes Including other files 

If, for some reason, you want to spread your PLUMED input over a number of files you can use \subpage INCLUDE as shown below:

\plumedfile
INCLUDE FILE=filename
\endplumedfile

So, for example, a single "plumed.dat" file:

\plumedfile
DISTANCE ATOMS=0,1 LABEL=dist
RESTRAINT ARG=dist
\endplumedfile

could be split up into two files as shown below:
 
\plumedfile
DISTANCE ATOMS=0,1 LABEL=dist
INCLUDE FILE=toBeIncluded.dat
\endplumedfile
plus a "toBeIncluded.dat" file
\plumedfile
RESTRAINT ARG=dist
\endplumedfile

However, when you do this it is important to recognise that \ref INCLUDE is a real directive that is only resolved
after all the \ref comments have been stripped and the \ref ContinuationLines have been unrolled.  This means it
is not possible to do things like:

\plumedfile
# this is wrong:
DISTANCE INCLUDE FILE=options.dat
RESTRAINT ARG=dist
\endplumedfile

\page load Loading shared libraries

You can introduce new functionality into PLUMED by placing it directly into the src directory and recompiling the 
PLUMED libraries.  Alternatively, if you want to keep your code independent from the rest of PLUMED (perhaps
so you can release it independely - we won't be offended), then you can create your own dynamic library.  To use this 
in conjuction with PLUMED you can then load it at runtime by using the \subpage LOAD keyword as shown below:

\plumedfile
LOAD FILE=library.so
\endplumedfile
 
N.B.  If your system uses a different suffix for dynamic libraries (e.g. macs use .dylib) then PLUMED will try to 
automatically adjust the suffix accordingly.

\page degub Debugging the code

The \subpage DEBUG action provides some functionality for debugging the code that may be useful if you are doing 
very intensive development of the code of if you are running on a computer with a strange architecture.

\page exchange-patterns Changing exchange patterns in replica exchange

Using the \subpage RANDOM_EXCHANGES keyword it is possible to make exchanges betweem randomly
chosen replicas. This is useful e.g. for bias exchange metadynamics \cite piana.

\page misc Frequently used tools

@DICTIONARY@
<TABLE ALIGN="center" FRAME="void" WIDTH="95%%" CELLPADDING="5%%">
<TR>
<TD WIDTH="5%"> 
\subpage Regex </TD><TD> </TD><TD> POSIX regular expressions can be used to select multiple actions when using ARG (i.e. \ref PRINT).
</TD>
</TR>
<TR>
<TD WIDTH="5%"> 
\subpage Files </TD><TD> </TD><TD> Dealing with Input/Outpt
</TD>
</TR>
</TABLE>

