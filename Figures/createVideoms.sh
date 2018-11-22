#!/bin/bash
mencoder mf://$1/*.png -mf w=800:h=600:fps=20:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o $1$2.avi
cp ms.srt $1.srt
mencoder -sub $1.srt -fontconfig -font Arial -subfont-text-scale 4 -ovc lavc $1.avi -o $1sub.avi
