set yl 'Pressure (eV/Ang^3)'
set xl 'Volume (Ang^3)'
p 'S1/out' u 6:5 w l lw 3 t 'S1', \
'S2_low/out' u 6:5 w l lw 3 t 'S2 ', \
'S2/out' u 6:5 w l t 'S2 8x KE dens', \
'GP/out' u 6:5 w l t 'GP', 0 notitle
pause -1

