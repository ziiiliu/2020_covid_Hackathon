# Binding measurements of viral spike mutants

This dataset describes a saturation-mutagenesis experiment done on the receptor-binding domain
of the SARS-CoV-2 virus' spike protein. For every amino acid in the protein domain, mutant variants
were generated for every possible other amino acid (out of the 20 standard amino acids) 
that could exist at that position. Thus every possible single-position mutation is tested for its
ability to bind to the virus' known ACE2 receptor. The measurements for spike binding to ACE2-expressing
cells have already been normalized.

For each mutation, two measurements are made: how often that mutation is found when you analyze the 
very strongest binding mutants out of all of them ("Top"), and how often that mutation is found when
you analyze the weakest binding mutants out of all of them ("Bottom"). Although these metrics are correlated,
because of the complexitiies of measurement in an experiment like this they are not identical.

For each position in the protein, the amino acid that is naturally found in the spike's sequence (wild-type 
or "WT") is defined as having 0 enrichment in the top-binders and 0 enrichment in bottom-binders.

