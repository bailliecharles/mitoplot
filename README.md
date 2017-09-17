# mitoplot
plotting circular genomes with R

I recently needed to make visualisations of mitogenomes and struggled to find the right solution. There are many nice 
programs out there that do this job but they either needed .gbk file format, are unnecessarily complicated, or aesthetically 
limited (colors, too busy, too crammed). R has a number of cool packages for plotting data in circular format (which is what 
a mitogenome is, right?) but for genomic data these mainly seem to be used for things like visualising SNP locations, 
alignments, or linkage maps etc. None seemed flexible enough to me for plotting genes around a circle that could be any given 
length - probably (definitely!) me not seeing it! I bet the bioconductor folks have something awesome. 

Anyway, here's my hack using the {circlize} package. Its pretty neat and easy to learn with the brill documentation - took 
me an evening to put this together.  You might not want/ be able to include everything I did but if you're interested, I got 
coverage stats for each site using BBMap (its super fast). I also didn't use circlize to annotate anything (did mine in Inkscape) but you could easily add a legend, names, an arrow for translation direction etc. just see the documentation. My philosophy was 
to make something clear and easy to understand, not something with loads of detail because I find it just becomes cluttered and 
confusing. You will need/should, therefore, tweak the heights of the lines, the window and offsets etc. to suit the 
visualisation of your data. 

Essentially what I've done is very simple: plus and minus strands are plotted by track with intergenic regions not plotted as a 'gap', but as a filled segments the same color as the background. 
