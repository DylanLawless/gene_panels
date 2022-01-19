#!/bin/bash

file=ge_panel_app_Primary_immunodeficiency_v2.1_20200224_coordinate_grch38p13

	# print
	# keep columns: chr start stop
	# remove headers
	# sort by chr then by position

cat $file.tsv | \
	cut -f 4,2,3 | \
	tail -n +2 | \
	sort -nk3,3 -nk1,1 \
	> $file.bed

