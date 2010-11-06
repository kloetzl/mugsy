#!/bin/sh
if [ "`(which wget) 2> /dev/null | grep -v '^no'`" != "" ]; then
    wget $1
elif [ "`(which curl) 2> /dev/null | grep -v '^no'`" != "" ]; then
    curl -O $1
elif [ "`(which fetch) 2> /dev/null | grep -v '^no'`" != "" ]; then
    fetch $1
else
    echo "No wget/curl/fetch found in $PATH -- aborting."
    exit 1
fi
