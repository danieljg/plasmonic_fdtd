#!/bin/bash
framespersecond=20
date
time ifort -o fdtd module.f90 fdtd.f90
time ./fdtd
#time gnuplot plot.p
#mencoder mf://animation/H*.png -mf fps=$framespersecond:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o fdtd1.avi > /dev/null
#mencoder mf://animation/Ex*.png -mf fps=$framespersecond:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o fdtd2.avi > /dev/null
#mencoder mf://animation/Ez*.png -mf fps=$framespersecond:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o fdtd3.avi > /dev/null

#rm *.dat
#rm animation/*.png
