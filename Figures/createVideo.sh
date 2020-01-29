#!/bin/bash
mencoder mf://$1/*.png -mf w=800:h=600:fps=32:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o $1$2.avi
cp daysFast.srt $1$2.srt
mencoder -sub $1$2.srt -fontconfig -font Arial -subfont-text-scale 4 -ovc lavc $1$2.avi -o $1$2sub.avi
