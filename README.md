# syntenyPlotteR_BETA
New drawing functions for syntenyPlotteR - IN ACTIVE DEVELOPMENT

WARNING! These functions are in active development and might contain bugs

## Microsynteny plots

Recommended for syntenic blocks <30Kbp 

The `draw.microsynteny` plot works similarly to `draw.linear` in the original syntenyPlotteR package (https://github.com/Farre-lab/syntenyPlotteR)

Example:

`draw.microsynteny("outputname","test_lengths.txt","test_alignment_1.txt",fileformat = "png",colours = c("red","blue"),w=13,h=5,opacity = .1,curve=.75,thickness=.5)`

## Linear plots 2.0

Additional parameter for `draw.linear` from the original syntenyPlotteR package (https://github.com/Farre-lab/syntenyPlotteR)

If working with small genomes, you can now alter the size of the gaps between the chromosomes to allow for more proportional sizes and therefore clearer plots

add the `insert_size = 6000000` (default) parameter to adjust insert size

Example: 

`draw.linear(output, sizefile, ..., fileformat = "png", colours = colours.default, w=13, h=5, opacity = .5,insert_size = 6000000)`

