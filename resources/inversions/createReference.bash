awk 'BEGIN {print "ID", "bp", "str_allele1", "str_allele2", "CHR", "allele1", "allele2" } { print $4, $2, "N", "I", "chr"$1, "4064881669293999478", "4654963978389870087" }' resources/inversions/inversion_ranges.txt > resources/inversions/inversion_reference.txt
bgzip resources/inversions/inversion_reference.txt
