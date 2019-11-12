#!/usr/bin/expect -f
spawn ssh mcastro@129.20.25.215
expect "*assword:*"
send "L7\$vg9pf\r"
interact

