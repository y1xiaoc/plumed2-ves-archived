/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionRegister.h"
#include "core/SetupMolInfo.h"

namespace PLMD {
namespace setup {

//+PLUMEDOC TOPOLOGY MOLINFO
/*
This command is used to provide information on the molecules that are present in your system.

The information on the molecules in your system can either be provided in the form of a pdb file
or as a set of lists of atoms that describe the various chains in your system. If a pdb file
is used plumed the MOLINFO command will endeavor to recognize the various chains and residues that
make up the molecules in your system using the chainIDs and resnumbers from the pdb file. You can
then use this information in later commands to specify atom lists in terms residues.  For example
using this command you can find the backbone atoms in your structure automatically.

\warning
Please be aware that the PDB parser in plumed is far from perfect. You should thus check the log file
and examine what plumed is actually doing whenenver you use the MOLINFO action.
Also make sure that the atoms are listed in the pdb with the correct order.
If you are using gromacs, the safest way is to use reference pdb file
generated with `gmx editconf -f topol.tpr -o reference.pdb`.

More information of the PDB parser implemented in PLUMED can be found \ref pdbreader "at this page".

Providing `MOLTYPE=protein`, `MOLTYPE=rna`, or `MOLTYPE=dna` will instruct plumed to look
for known residues from these three types of molecule. In other words, this is available for
historical reasons and to allow future extensions where alternative lists will be provided.
As of now, you can just ignore this keyoword.

Using MOLINFO with a protein's or nucleic acid's pdb extends the possibility of atoms selection using the @ special
symbol in the form

\verbatim
@"definition"-chainresiduenum
@"definition"-residuenum
\endverbatim

So for example

\verbatim
@psi-1 will select the atoms defining the psi torsion of residue 1
@psi-C1 will define the same torsion for residue 1 of chain C.
\endverbatim

In the following are listed the current available definitions:

For protein residues, the following groups are available:

\verbatim
@phi-#
@psi-#
@omega-#
@chi1-#
\endverbatim

that select the appropriate atoms that define each dihedral angle for residue #.

For DNA or RNA residues, the following groups are available:

\verbatim
# quadruplets for backbone dihedral angles
@alpha-#
@beta-#
@gamma-#
@delta-#
@epsilon-#
@zeta-#

# quadruplets for sugar dihedral angles
@v0-#
@v1-#
@v2-#
@v3-#
@v4-#

# quadruplet corresponding to the chi torsional angle
@chi-#

# backbone, sugar, and base heavy atoms
@back-#
@sugar-#
@base-#

# ordered triplets of atoms on the 6-membered ring of nucleobases
# namely:
#  C2/C4/C6 for pyrimidines
#  C2/C6/C4 for purines
@lcs-#
\endverbatim

Notice that `zeta` and `epsilon` groups should not be used on 3' end residue and `alpha` and `beta`
should not be used on 5' end residue.

Furthermore it is also possible to pick single atoms using the syntax
`@atom-chainresiduenum` or `@atom-residuenum`.

\warning If a residue-chain is repeated twice in the reference pdb only the first entry will be selected.

\bug At the moment the HA1 atoms in a GLY residues are treated as if they are the CB atoms. This may or
may not be true - GLY is problematic for secondary structure residues as it is achiral.

\bug If you use WHOLEMOLECULES RESIDUES=1-10 for a 18 amino acid protein
( 18 amino acids + 2 terminal groups = 20 residues ) the code will fail as it will not be able to
interpret terminal residue 1.

\par Examples

In the following example the MOLINFO command is used to provide the information on which atoms
are in the backbone of a protein to the ALPHARMSD CV.

\plumedfile
MOLINFO STRUCTURE=reference.pdb
ALPHARMSD RESIDUES=all TYPE=DRMSD LESS_THAN={RATIONAL R_0=0.08 NN=8 MM=12} LABEL=a
\endplumedfile

The following example prints the distance corresponding to the hydrogen bonds
in a GC Watson-Crick pair.

\plumedfile
MOLINFO STRUCTURE=reference.pdb
hb1: DISTANCE ATOMS=@N2-1,@O2-14
hb2: DISTANCE ATOMS=@N1-1,@N3-14
hb3: DISTANCE ATOMS=@O6-1,@N4-14
PRINT ARG=hb1,hb2,hb3
\endplumedfile

This example use MOLINFO to calculate torsions angles

\verbatim
MOLINFO MOLTYPE=protein STRUCTURE=myprotein.pdb
t1: TORSION ATOMS=@phi-3
t2: TORSION ATOMS=@psi-4
PRINT ARG=t1,t2 FILE=colvar STRIDE=10
\endverbatim

*/
//+ENDPLUMEDOC


/*
This action is defined in core/ as it is used by other actions.
Anyway, it is registered here, so that excluding this module from
compilation will exclude it from plumed.
*/

typedef PLMD::SetupMolInfo MolInfo;

PLUMED_REGISTER_ACTION(MolInfo,"MOLINFO")

}
}
