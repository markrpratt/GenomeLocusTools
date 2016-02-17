# GenomeLocusTools
a collection of matlab functions for performing genome locus set (GLS) algebra.  
coordinate intervals are stored internally in an Nx3 array of int32 where the first
number is the sequence name index (in segNames), the second is the start coordinate and the last is
the stop coordinate.  All coordinates are ones bases so this start coordinate is +1 the BED convention.
For the most part the tool set is reference unaware but I envision adding some support for this in the
future, particularly when it becomes helpful to know segment lengths and to standardize on segment names.

Mark Pratt
Feb 2016
