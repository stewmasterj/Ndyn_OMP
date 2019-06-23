#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 5.0 patchlevel rc2    last modified 2014-08-28 
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2014
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	mailing list:     gnuplot-beta@lists.sourceforge.net
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
set output 'sene.png'
set term png font DejaVuSans 12
set yl 'Configuration Energy (eV)'
set xl 'X Strain'
p 'out' u ($1*0.1*1e-2/36.8):2 w l lw 3, 0.73025*x*x*36.8**3*0.74 w l
exit
