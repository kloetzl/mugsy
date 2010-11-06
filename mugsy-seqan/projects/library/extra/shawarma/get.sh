#!/bin/sh
if [ "`(which wget) 2> /dev/null | grep -v '^no'`" != "" ]; then
    wget $1 -O $2
elif [ "`(which curl) 2> /dev/null | grep -v '^no'`" != "" ]; then
    curl $1 -o $2
elif [ "`(which fetch) 2> /dev/null | grep -v '^no'`" != "" ]; then
    fetch $1 -o $2
else
    echo "No wget/curl/fetch found in $PATH -- aborting."
    exit 1
fi
