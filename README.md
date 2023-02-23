# βNTI (Feature-level)
βNTI-feat is a metric that allows us to evaluate how a member of a community contributed to community structure. The scripts contained within this repository will help users calculate βNTI-feat, including null generation and null merging. The usage of these scripts will depend on your goal - specifically, βNTI-feat is equiped to handle data across various different scales. Whether looking at how one specific ASV contributes to community assembly in your dataset or how a molecular formula affects assemblage dynamics across treatment levels, βNTI-feat can help. The usage as follows:

- A file including "_create-nulls_" will generate the null values for a given comparison type, while "_merge-nulls_" will combine the created nulls and calculate βNTI-feat.

- "_new.R" performs βNTI-feat analyses on the whole dataset

- "_new_one_sample.R" performs βNTI-feat for pairwise comparisons against the specified sample

- "_new_within_group.R" performs βNTI-feat for a given treatment

- "_new_withing_subgroup.R" performs βNTI-feat for a level within a given treatment
