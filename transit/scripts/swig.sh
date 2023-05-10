#!/bin/bash
PYV=`python -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)";`

if [ $(echo "$PYV >= 3.0"|bc) -eq 1 ]
then
	swig -python -py3 -o python/transit_wrap.c src/transit.i
else
	swig -python -o python/transit_wrap.c src/transit.i
fi
