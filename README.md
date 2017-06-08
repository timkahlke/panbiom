# panbiom
Estimation of a pan-microbiom.

For usage examples see the WIKI at https://github.com/timkahlke/panbiom/wiki


# Usage
panbiom.py [options] BIOM_FILE OUTPUT_FILE

# Options
-t, --treatments 
List of samples (treatments) that should be considered for analysis

-m, --abundance_minimum: Minimum abundance (between 0-1) an OTU must have to be considered present (default: 0.0)

-p, --abundance_parameter: Defines whether the abundance threshold should be based on the complete data set (c) or on the defined treatments (t) (default: t)

-r, --replicate_threshold: If groups/replicates are defined in the treatments file (second column) the given abundance threshold has to be met in at least given number of replicates.





