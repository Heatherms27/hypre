#!/bin/sh
#BHEADER**********************************************************************
# Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# This file is part of HYPRE.  See file COPYRIGHT for details.
#
# HYPRE is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.
#
# $Revision$
#EHEADER**********************************************************************

testname=`basename $0 .sh`

# Echo usage information
case $1 in
   -h|-help)
      cat <<EOF

   **** Only run this script on one of the tux machines. ****

   $0 [-h|-help] {src_dir}

   where: -h|-help   prints this usage information and exits
          {src_dir}  is the hypre source directory
          

   This script runs a number of tests suitable for the tux machines.

   Example usage: $0 ..

EOF
      exit
      ;;
esac

# Setup
test_dir=`pwd`
output_dir=`pwd`/$testname.dir
rm -fr $output_dir
mkdir -p $output_dir
src_dir=$1
shift

# Basic build and run tests
mo="-j test"
ro="-ams -ij -sstruct -struct"
eo=""

co=""
test.sh basictest.sh $src_dir -co: $co -mo: $mo
renametest.sh basictest $output_dir/basictest-default

co="--with-insure --enable-debug --with-print-errors"
test.sh basictest.sh $src_dir -co: $co -mo: $mo -ro: $ro -eo: $eo
renametest.sh basictest $output_dir/basictest--with-insure

co="--enable-global-partition --with-insure"
RO="-ams -ij -sstruct -struct -fac"
test.sh basictest.sh $src_dir -co: $co -mo: $mo -ro: $RO -eo: $eo
renametest.sh basictest $output_dir/basictest--enable-global-partition

co="--without-MPI"
test.sh basictest.sh $src_dir -co: $co -mo: $mo
renametest.sh basictest $output_dir/basictest--without-MPI

co="--with-strict-checking"
test.sh basictest.sh $src_dir -co: $co -mo: $mo
renametest.sh basictest $output_dir/basictest--with-strict-checking

co="--enable-shared"
test.sh basictest.sh $src_dir -co: $co -mo: $mo
renametest.sh basictest $output_dir/basictest--enable-shared

co="--enable-bigint --enable-debug"
test.sh basictest.sh $src_dir -co: $co -mo: $mo -ro: $ro -eo: -bigint
renametest.sh basictest $output_dir/basictest--enable-bigint

co="--enable-maxdim=4 --enable-debug"
test.sh basictest.sh $src_dir -co: $co -mo: $mo -eo: -maxdim
renametest.sh basictest $output_dir/basictest--enable-maxdim=4

co="--enable-complex --enable-maxdim=4 --enable-debug"
test.sh basictest.sh $src_dir -co: $co -mo: $mo -eo: -complex
# ignore complex compiler output for now
rm -fr basictest.dir/make.???
grep -v make.err basictest.err > basictest.err
renametest.sh basictest $output_dir/basictest--enable-complex

# Test babel build only if 'babel-runtime' directory is present
if [ -d $src_dir/babel-runtime ]; then
   co="--with-babel"
   MO="test"  # the -j option doesn't always work with the babel code
   test.sh basictest.sh $src_dir -co: $co -mo: $MO
   renametest.sh basictest $output_dir/basictest--with-babel
fi

# Test linking for different languages
link_opts="all++ all77"
for opt in $link_opts
do
   output_subdir=$output_dir/link$opt
   mkdir -p $output_subdir
   ./test.sh link.sh $src_dir $opt
   mv -f link.??? $output_subdir
done

# Test documentation build (only if 'docs_misc' directory is present)
if [ -d $src_dir/docs_misc ]; then
   ./test.sh docs.sh $src_dir
   mv -f docs.??? $output_dir
fi

# Check for 'int', 'double', and 'MPI_'
./test.sh check-int.sh $src_dir
mv -f check-int.??? $output_dir
./test.sh check-double.sh $src_dir
mv -f check-double.??? $output_dir
./test.sh check-mpi.sh $src_dir
mv -f check-mpi.??? $output_dir

# Echo to stderr all nonempty error files in $output_dir
for errfile in $( find $output_dir ! -size 0 -name "*.err" )
do
   echo $errfile >&2
done
