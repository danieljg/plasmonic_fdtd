#!/bin/bash
framespersecond=4
date
time ifort -o fdtd module.f90 fdtd.f90 -r8 -parallel -ipo -fast
time ./fdtd
mencoder mf://tmp/H*.png -mf fps=$framespersecond:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o animation/H.avi > /dev/null
mencoder mf://tmp/E*.png -mf fps=$framespersecond:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o animation/E.avi > /dev/null

#rm tmp/*.dat
rm tmp/*.png
