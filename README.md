# panbiom
Estimation of a pan-microbiom.

For usage examples see the WIKI at https://github.com/timkahlke/panbiom/wiki


# Usage
panbiom.py [options] BIOM_FILE OUTPUT_FILE

# Options
-t, --treatments 
List of samples (treatments) that should be considered for analysis

-m, --abundance_minimum: Minimum abundance (between 0-1) an OTU must have to be considered present (default: 0.0)

-r, --replicate_threshold: If groups/replicates are defined in the treatments file (second column) the given abundance threshold has to be met in at least given number of replicates.

-p, --print_taxonomy: If set to False no taxonomy will be printed. Needed for biom files without taxonomy field.





