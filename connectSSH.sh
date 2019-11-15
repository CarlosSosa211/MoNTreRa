#!/usr/bin/expect -f
spawn ssh CarlosSosa@129.20.25.215
expect "*assword:*"
send "5OQnCeBe\r"
interact

