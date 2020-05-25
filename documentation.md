Welcome to the Eli Lilly Medchem Rules implementation.

BUILDING:

This software implements the rejection rules described by the 2012
publication in the Journal of Medicinal Chemistry by
R.F. Bruns and I.A. Watson (J. Med. Chem. 2012, 55, 9763-9772,
freely accessible under doi 10.1021/jm301008n).

Compiling should be as simple as entering 'make' in this directory.
Consider about 45 MB permanent memory to store the working program.
There is no separate installation, the scripts are designed to be run
from this location, although it would not be hard to change that. You
will need to figure out your LD_LIBRARY_PATH.

Once the executables are built, the rules can be invoked via the ruby
script Lilly_Medchem_Rules.rb in this directory.  That script will
locate the executables and data files it needs from the directory
hierarchy and run the rules.

EXECUTION:

The normal invocation will be of the form

ruby Lilly_Medchem_Rules.rb input.smi > okmedchem.smi

The molecules from input.smi will be scanned, and those that survive
will be written to stdout and captured in okmedchem.smi. Molecules
failing the rules will be written to one of four files, bad0.smi
through bad3.smi. Molecules in those files will be annotated with
the reason for their rejection.

In the preparation of the .smi files for this program, the user is
advised to refrain from an articifical aromatization programs like
openbabel may offer; this is quickly recognized by the use of lower
case characters e.g., c, n, o, p, s for atoms of carbon, nitrogen,
oxygen, phosphorous, and sulphur.

The output file will look like

ClC1=NC(Cl)=CC(=C1)C(=O)CN PBCHM19820997
O=C1NC(N)(C1)C1=CC=CC=C1 PBCHM19423780
S=C1N=C(N)NC(=N1)C1CCC1 PBCHM8009179 : D(30) thiocarbonyl_aromatic
BrC1=CC=C(S1)C1=CC(=NN1)N PBCHM16495471 : D(34) bromine
S1C=C(NC(=O)N)C(=C1)C(=O)OC PBCHM21568455 : D(60) thiophene_furan_n_acyl:ester
O=C(NC)C1=CN=C(NN)C=C1 PBCHM2826732 : D(75) N-N

The first token on each line is the smiles, followed by the molecule
name - whatever was in the input file. Some molecules will pass
unchanged through the rules, accruing no demerits - the first two
examples above. Other molecules that have passed, but have attracted
demerits will have the demerits shown in the form D(nn) above,
followed by the reason(s) for the demerits. If you don't care about
the demerits associated with passing molecules, invoke the script
with -noapdm and the demerit information will not be appended.

Note that the software is set up to ignore smiles it cannot interpret.
This might be a bad idea, potentially problematic structures should
generally be investigated. Check the file ok0.log after execution to 
see evidence of failed smiles interpretation.

Tool mc_summarise can be helpful in getting an overview of how the
rules might be impacting a particular set of molecules. An example
usage might be

bin/mc_summarise -T Reasons okmedchem.smi

QUERY FILES

The query files are on a modified Cerius-2 format. Today, there are
many better choices for file formats, and if the software were
being developed today, we would likely use JSON, Proto or similar.

SOFTWARE

This software has been extensively tested on a variety of *nix type
systems.  Today it is primarily used on RedHat systems, using either
gcc or the Intel Compiler.  In the past, the code has been
successfully deployed on cygwin, Windows Visual Studio, Sun, SGI, HP
and other systems. 

Work on the code commenced in early 1995.  At that time, commonly
accepted language features like STL, iterators, Boost and other
conveniences were either absent, or poorly standardised.  In fact the
mid 1990's, during which time most of the code was developed, was a
time of significant evolution of the language.  For that reason, there
are many constructs in the code that may seem strange today - user
written String and Vector classes for example.  But this was driven by
the seemingly overwhelming difficulties associated with trying to reconcile
different implementations, and changing standards of the time.  Today
many of these language features are stable both in time and across
implementations.  That was not always the case.  We are better off for
that evolution.

I have also learned a lot about OO design over the years, and many
of the ideas implemented here would benefit from refactoring to more
modern design sensibilities. There is however no motivation to do that.

Parsing of smarts and smiles is not elegant - it too would benefit
from a re-design.

The code is generally not thread safe - there is extensive use of file
scope static variables to indicate optional behaviours.  Parallel
processing for molecules is usually best done by separately processing
different chunks of molecules.  That said, we have successfully built
several multi-threaded applications, it all depends on which parts of
the code you use.

Speed was never a consideration for this code.  In a project driven
environment, the emphasis was always on getting functionality
available, with little thought to efficiency.  Besides, the 1990's was
the time when Intel would double the speed of computers every couple
of years.  For us, efficiency has nevertheless been quite adequate,
especially since most molecular processing tasks are pleasingly
parallel and well served by modern clusters and SMP systems.

Many class methods will start with

assert (ok())

where the ok() function is a function that checks the internal
integrity of the object.  If you compile with -DNDEBUG this will
compile out all these checks, and you will see a speed improvement of
maybe five percent.

A similar speed improvement is usually seen by compiling with the
Intel compiler rather than g++. There are flags in the Makefiles to
make this change.

While the Molecule object makes no claim to being a comprehensive 
Chemistry structure representation, it has provided the basis for
a great many useful chemically aware tools.

Design Philosophies

Molecules are represented in their Kekule form.  Atom objects do not
know they are part of a Molecule object.  Most processing of Atom
objects is via atom numbers in the parent molecule.  Rings consist of
ordered sets of atom numbers. Atom objects do not know anything about
Ring objects, although their corresponding atom numbers describe Ring
objects. Bond objects may know their ring membership.

The Molecule object is a lazy object.  It will not contain any derived
parameters unless requested.  For example, a ring determination will
not be performed unless some method is invoked that requires ring
perception.  Similarly for aromaticity, it will simply never be
computed unless it is requested. There is no concept of a molecule
being in a mutable or immutable state, it tries very hard
to respond properly to changes and queries. Almost certainly, some
efficiency is sacrificed because of this.

There are many limitations/bugs/blemishes.

Canonicalisation of cis-trans bonds does not work - the implementation
is wrong, but this is not a priority for us.

Aromaticity handling is a mess and should be re-done.  But the current
setup works pretty well, and redoing it would take a lot of work.

Some molecules are not canonicalised properly. This is a hard problem,
where seemingly all Molecular processing tools exhibit problematic
behaviour. Since most of the molecules that cause problems for this
implementation are not of pharmaceutical interest, we ignore these
difficulties. The tool tsmiles (cd Molecule && make tsmiles) is a tester 
that will automatically test canonicalisation and identify failures.
Suggest always using the -u option.

time Molecule/tsmiles -p 5 -a -t 0 -w 3 -m 2 -u -q test/example_molecules.smi

Takes around 1 minute to test these 36k molecules - no errors reported.

Getting consistent regular expression behaviour across implementations
proved impossible and I ultimately used the regexp source from grep. 
This is quite ugly and fragile, and should be fixed sometime.  The
whole regular expression area is weak, but not relied upon a great
deal.

As will be visible within the code, the code supports only a limited
number of file types, but for this distribution, we have included just
smiles and various mdl flavours.  The other file types are not used
extensively, and their handling is not always robust.  There are other
tools available for doing such conversions in a robust fashion.

Parsing of SDF files is very messy, but over the years, we have encountered
a vast number of variations in SDF files.

Toolkit usage:

The Molecule object is designed to be a functional tool for dealing
with commonly encountered, drug-like small molecules.  It is not
designed to handle all possible chemistry forms. 

There are a great many parameters that can be set which control how the
toolkit handles cases where reasonable people might need different
behaviour.  Among the most important is aromaticity.  The default
aromaticity definitions, named Pearlman (after Bob), are quite
restrictive in what is considered aromatic - only structures having
Kekule forms are considered aromatic. 

An alternate aromaticity definition, called Daylight (after the
company) takes a more liberal view of aromaticity.  It is recommended
that all new programmes begin with

 set_global_aromaticity_type(Daylight)

to set this as the default. By default, rings with 4 atoms (and 2 pi 
electrons) are not aromatic.  If this is desired, call 

 set_allow_two_electron_systems_to_be_aromatic(1)

which will enable this. There are other aromaticity definitions 
included that are variants on these. Aromaticity determination of
a Molecule will be governed by the currently applicable global
aromaticity definition.

Included in the source directory is a file, skeleton.cc, which gives
a starting place upon which Molecule applications can be built.
By default, it provides command line interfaces to a great many
optional settings. Linking will be similar to the other tools 
included.

Particular care is needed when a Molecule is altered. If you have
extracted a Ring or Bond via a pointer, and you change the molecule,
those pointers may no longer be valid. No warning, usually crashes.
Perform changes to the Molecule, then ask for Rings, Bonds, Atoms,
etc..



If you find bugs, we would be very interested in learning about those. 
Please contact ianiwatson@gmail.com  Or if you have a fix, then please
do send it along.

We hope this software proves useful, it has been very useful for us
and we are pleased to be able to share.

Longer term, the C++ software here will be replaced by the versions at
[LillyMol](https://github.com/EliLillyCo/LillyMol), which contains
more functionality. The query functionality will remain unchanged.


Please consult the LICENSE file for details of the license.
