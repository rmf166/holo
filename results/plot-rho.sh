#!/bin/sh
rm *.pdf
sols="LD LC"
for sol in $sols
do
  for sw in 1 2 3 4
  do
    rm plot.p
    echo 'set autoscale' >> plot.p
    echo 'unset logscale; set logscale x' >> plot.p
    echo 'unset label' >> plot.p
    echo 'set xtic auto' >> plot.p
    echo 'set ytic auto' >> plot.p
    echo 'set title "'${sol}', s = '${sw}'"' >> plot.p
    echo 'set xlabel "{/Symbol t} (mfp)" enhanced' >> plot.p
    echo 'set ylabel "{/Symbol r}" enhanced' >> plot.p
    echo 'set yr [0:1]' >> plot.p
    echo 'set xr [0.01:20000]' >> plot.p
    echo 'plot    "numres-p1-s'${sw}'-'${sol}'.dat" using 1:2 notitle with points pointtype 1 ps 1 lc rgb "red", \' >> plot.p
    echo '        "numres-p1-s'${sw}'-'${sol}'.dat" using 1:3 notitle with points pointtype 1 ps 1 lc rgb "green", \' >> plot.p
    echo '        "numres-p1-s'${sw}'-'${sol}'.dat" using 1:4 notitle with points pointtype 1 ps 1 lc rgb "blue", \' >> plot.p
    echo '        "numres-p1-s'${sw}'-'${sol}'.dat" using 1:5 notitle with points pointtype 1 ps 1 lc rgb "violet", \' >> plot.p
    echo '        "result-p1-s'${sw}'-'${sol}'.dat" using 1:2 title "c=0.8"  with lines linetype 1 lc rgb "red", \' >> plot.p
    echo '        "result-p1-s'${sw}'-'${sol}'.dat" using 1:3 title "c=0.9"  with lines linetype 1 lc rgb "green", \' >> plot.p
    echo '        "result-p1-s'${sw}'-'${sol}'.dat" using 1:4 title "c=0.99" with lines linetype 1 lc rgb "blue", \' >> plot.p
    echo '        "result-p1-s'${sw}'-'${sol}'.dat" using 1:5 title "c=1.00" with lines linetype 1 lc rgb "violet"' >> plot.p
    echo 'set key left top' >> plot.p
    echo 'set terminal pdfcairo enhanced color dashed' >> plot.p
    echo 'set output "plot-'${sw}'-'${sol}'.pdf"' >> plot.p
    echo 'replot' >> plot.p
    echo 'set terminal x11' >> plot.p
    gnuplot plot.p
  done
done
pdfunite plot*LD.pdf plot*LC.pdf final.pdf
