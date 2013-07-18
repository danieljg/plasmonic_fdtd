set term pngcairo size 800,480
#THIS MUST MATCH WITH THE NUMBER OF SNAPSHOTS ON FDTD
snapshots=240

green = 1.0/3.0
blue = 2.0/3.0
set palette model HSV defined ( 0 blue 0.7 0.3, 1 blue 0.7 1, 2 green 0.5 0.6, 3 green 0.7 0.3 )
set cbrange [-1.2:1.2]

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
outfile = sprintf('animation/H%03.0f.png',idx)
set output outfile
set title 'Magnetic field'
splot 'H.dat' u ($1*1e4):($2*1e4):($3) index idx noti
set output
outfile = sprintf('animation/Ex%03.0f.png',idx)
set title 'Electric field, x-component'
set output outfile
splot 'Ex.dat' u ($1*1e4):($2*1e4):($3) index idx noti
set output
outfile = sprintf('animation/Ez%03.0f.png',idx)
set title 'Electric field, z-component'
set output outfile
splot 'Ez.dat' u ($1*1e4):($2*1e4):($3) index idx noti
set output
}
#system('rm animation/*.eps')
#system('rm animation/*.pdf')
#system('mv *.png animation/')

