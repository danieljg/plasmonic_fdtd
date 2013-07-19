set term png size 600,600
#THIS MUST MATCH WITH THE NUMBER OF SNAPSHOTS ON FDTD
snapshots=240

red = 1.0/20.0
yellow = 1.0/6.0
green = 1.0/3.0
blue = 2.0/3.0

set size ratio -1
#set xrange [0:1.620]
#set yrange [0:0.645]
set xlabel 'x (um)'
set ylabel 'z (um)'
set pm3d
unset surface
set view map

do for [idx=0:snapshots] {

 print 'plotting index '.idx.' of '.snapshots
 set palette model HSV defined ( 0 blue 0.7 0.3, 1 blue 0.7 1, 2 green 0.5 0.6, 3 green 0.7 0.3 )
outfile = sprintf('tmp/H%03.0f.png',idx)
 set output outfile
 set title 'Magnetic field'
 set cbrange [-1.:1.]
 unset logscale zcb
 splot 'tmp/H.dat' u ($1*1e4):($2*1e4):($3) index idx noti
 set output

 set palette model HSV defined ( 0 blue 0.7 0.3, 1.5 green 0.5 0.6, 2 yellow 0.5 0.7, 2.5 red 0.7 1, 3.5 red 0.8 0.9, 5 red 0.7 0.7 )
 outfile = sprintf('tmp/E%03.0f.png',idx)
 set title 'Electric field intensity'
 set output outfile
 set cbrange [0:1]
 splot 'tmp/E.dat' u ($1*1e4):($2*1e4):($3) index idx noti
 set output

}
#system('rm animation/*.eps')
#system('rm animation/*.pdf')
#system('mv *.png animation/')

