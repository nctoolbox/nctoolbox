#!/usr/bin/env bash
#
# Script to generate a neatly package zip archive for distributing
# ncdataset

# Brian Schlining - 20090701

TARGET=target/nctoolbox
if [ -d $TARGET ]; then
	rm -rf target
fi

mkdir -p $TARGET
cp -R cdm $TARGET/cdm
cp -R demos $TARGET/demos
cp -R java $TARGET/java
cp README $TARGET
cp setup_nctoolbox.m $TARGET
rm -rf $(find $TARGET -name .svn)
rm -rf $(find $TARGET -name .hg)
rm -rf $(find $TARGET -name *~)
rm -rf $(find $TARGET -name *DS_Store)

if [[ $OSTYPE == darwin* ]]; then
  ditto -c -k --keepParent --noextattr --norsrc $TARGET nctoolbox-$(date "+%Y%m%d").zip
elif [[ $OSTYPE == linux-gnu ]]; then
	myhome=${pwd}
	cd target
	ncname=nctoolbox-$(date "+%Y%m%d").zip
	zip -r $ncname nctoolbox
	mv $ncname ../$ncname
	cd $myhome
fi

rm -rf target

