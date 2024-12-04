# syntenyPlotteR_BETA
New drawing functions for syntenyPlotteR - IN ACTIVE DEVELOPMENT

WARNING! These functions are in active development and might contain bugs

## Microsynteny plots

Recommended for syntenic blocks <30Kbp 

The `draw.microsynteny` plot works similarly to `draw.linear` in the original syntenyPlotteR package (https://github.com/Farre-lab/syntenyPlotteR)

Example:

`draw.microsynteny("outputname","test_lengths.txt","test_alignment_1.txt",fileformat = "png",colours = c("red","blue"),w=13,h=5,opacity = .1,curve=.75,thickness=.5)`

## Linear plots 2.0

Optimisation for `draw.linear` from the original syntenyPlotteR package (https://github.com/Farre-lab/syntenyPlotteR)

Additional parameters: 

You can now alter the size of the gaps between the chromosomes to allow for more proportional sizes if working with larger or smaller genomes than generally expected 

add the `insert_size = 6000000` (default) parameter to adjust insert size

The sizes for the chromosome label and chromosome ID label can now be altered through the use of:
`chr.label.size = 2`
`sps.label.size = 2`

additionally the chromosome ID label angle can be altered along with the positioning (height) above the chromosome drawing:
`angle.chr.label = 45`
`chr.label.height = 0.2`

Finally, if you want to annotate the location of the centromere on the chromosome drawing, include a fourth column in the sizefile with the centromere position (see https://github.com/Farre-lab/syntenyPlotteR for more details on the sizefile), the colour can be adapted with: `centromere.colour = "red"`



Example: 

`draw.linear.2.0(output,sizefile,..., fileformat = "png", colours = colours.default, w=13, h=5, opacity = .5,insert.size = 6000000,chr.label.size = 4, sps.label.size = 7, angle.chr.label = 45, chr.label.height = 0.6)`


