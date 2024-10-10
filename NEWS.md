# Version 2, October 2024

This is a significant new release of the Lilly Medchem Rules. This brings
this repo into sync with how the rules had evolved within Lilly since
they were last synchronised.

The test suite has been re-implemented using Chembl structures.

There is a new option, `-nophosphorus` which is used to remove
all molecules contianing Phosphorus atoms. Most Phosphorus groups
were already matched by one or more rules, and several people
have mentioned that they usually discard Phosphorus anyway.

This release features some new rules and numerous fixes and tweaks to
existing rules, as well as some performance improvements.  One obvious
performance improvement that should have been done earlier is that all
queries are applied in order of prevalence within Chembl.  That way,
the queries that are most likely to remove a molecule are executed
first.  Once a molecule has been rejected, no further processing is
applied.

The default output format remains awkward, but most people using
the tool have worked out something that works in their context.

Your feedback continues to be both encouraged and most welcome, either via GitHub or to
ianiwatson@gmail.com.
