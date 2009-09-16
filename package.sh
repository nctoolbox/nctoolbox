#!/usr/bin/env bash
#
# Script to generate a neatly package zip archive for distributing
# ncdataset

# Brian Schlining - 20090701

TARGET=target/nctoolbox
mkdir -p $TARGET
cp -R data $TARGET/cdm
cp -R demos $TARGET/demos
cp -R java $TARGET/java
cp README $TARGET
cp setup_nctoolboxt.m $TARGET
rm -rf $(find $TARGET -name .svn)
ditto -c -k --keepParent -rsrc $TARGET nctoolbox-$(date "+%Y%m%d").zip
rm -rf target
