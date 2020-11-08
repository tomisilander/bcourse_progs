#!/bin/bash

FLAGS="$@"

cd /root/bcourse_progs

objdir=/root/bcourse_progs/obj_files
bindir=/root/bcourse_progs/bin_files
mkdir -p $objdir $bindir

# build static library

for f in src/common/*.c ; do
	oname=`basename $f .c`
	gcc $FLAGS -Isrc/common -c $f -o $objdir/$oname.o
done
ar rvs $objdir/libcommon.a $objdir/*.o

# build programs
for f in src/*.c ; do
	bname=`basename $f .c`
	gcc $FLAGS -o $bindir/$bname -Isrc/common $f -lm -lJudy  $objdir/libcommon.a
done

rm -fr $objdir
