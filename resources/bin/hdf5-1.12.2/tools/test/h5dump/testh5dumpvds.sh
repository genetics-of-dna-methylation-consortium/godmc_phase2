#! /bin/sh
#
# Copyright by The HDF Group.
# Copyright by the Board of Trustees of the University of Illinois.
# All rights reserved.
#
# This file is part of HDF5.  The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the COPYING file, which can be found at the root of the source code
# distribution tree, or in https://www.hdfgroup.org/licenses.
# If you do not have access to either file, you may request a copy from
# help@hdfgroup.org.
#
# Tests for the h5dump tool with vds type files

srcdir=.

# Determine which filters are available
USE_FILTER_SZIP="no"
USE_FILTER_DEFLATE="yes"

TESTNAME=h5dump
EXIT_SUCCESS=0
EXIT_FAILURE=1

DUMPER=../../src/h5dump/h5dump                     # The tool name
DUMPER_BIN=`pwd`/$DUMPER          # The path of the tool binary

H5DIFF=../../src/h5diff/h5diff           # The h5diff tool name
H5DIFF_BIN=`pwd`/$H5DIFF          # The path of the h5diff  tool binary

H5IMPORT=../../src/h5import/h5import     # The h5import tool name
H5IMPORT_BIN=`pwd`/$H5IMPORT      # The path of the h5import  tool binary

RM='rm -rf'
CMP='cmp -s'
DIFF='diff -c'
CP='cp'
DIRNAME='dirname'
LS='ls'
AWK='awk'

nerrors=0
verbose=yes

# source dirs
SRC_TOOLS="$srcdir/../.."

SRC_TOOLS_TESTFILES="$SRC_TOOLS/testfiles"
# testfiles source dirs for tools
SRC_H5LS_TESTFILES="$SRC_TOOLS_TESTFILES"
SRC_H5DUMP_TESTFILES="$SRC_TOOLS_TESTFILES"
SRC_H5DUMP_ERRORFILES="$srcdir/errfiles"
SRC_H5DIFF_TESTFILES="$SRC_TOOLS/test/h5diff/testfiles"
SRC_H5COPY_TESTFILES="$SRC_TOOLS/test/h5copy/testfiles"
SRC_H5REPACK_TESTFILES="$SRC_TOOLS/test/h5repack/testfiles"
SRC_H5JAM_TESTFILES="$SRC_TOOLS/test/h5jam/testfiles"
SRC_H5STAT_TESTFILES="$SRC_TOOLS/test/h5stat/testfiles"
SRC_H5IMPORT_TESTFILES="$SRC_TOOLS/test/h5import/testfiles"

TEST_P_DIR=./testfiles
TESTDIR=./testfiles/vds
test -d $TEST_P_DIR || mkdir -p $TEST_P_DIR
test -d $TESTDIR || mkdir -p $TESTDIR

######################################################################
# test files
# --------------------------------------------------------------------
# All the test files copy from source directory to test directory
# NOTE: Keep this framework to add/remove test files.
#       Any test files from other tools can be used in this framework.
#       This list are also used for checking exist.
#       Comment '#' without space can be used.
# --------------------------------------------------------------------
LIST_HDF5_TEST_FILES="
$SRC_H5DUMP_TESTFILES/vds/1_a.h5
$SRC_H5DUMP_TESTFILES/vds/1_b.h5
$SRC_H5DUMP_TESTFILES/vds/1_c.h5
$SRC_H5DUMP_TESTFILES/vds/1_d.h5
$SRC_H5DUMP_TESTFILES/vds/1_e.h5
$SRC_H5DUMP_TESTFILES/vds/1_f.h5
$SRC_H5DUMP_TESTFILES/vds/1_vds.h5
$SRC_H5DUMP_TESTFILES/vds/2_a.h5
$SRC_H5DUMP_TESTFILES/vds/2_b.h5
$SRC_H5DUMP_TESTFILES/vds/2_c.h5
$SRC_H5DUMP_TESTFILES/vds/2_d.h5
$SRC_H5DUMP_TESTFILES/vds/2_e.h5
$SRC_H5DUMP_TESTFILES/vds/2_vds.h5
$SRC_H5DUMP_TESTFILES/vds/3_1_vds.h5
$SRC_H5DUMP_TESTFILES/vds/3_2_vds.h5
$SRC_H5DUMP_TESTFILES/vds/4_0.h5
$SRC_H5DUMP_TESTFILES/vds/4_1.h5
$SRC_H5DUMP_TESTFILES/vds/4_2.h5
$SRC_H5DUMP_TESTFILES/vds/4_vds.h5
$SRC_H5DUMP_TESTFILES/vds/5_a.h5
$SRC_H5DUMP_TESTFILES/vds/5_b.h5
$SRC_H5DUMP_TESTFILES/vds/5_c.h5
$SRC_H5DUMP_TESTFILES/vds/5_vds.h5
$SRC_H5DUMP_TESTFILES/vds/a.h5
$SRC_H5DUMP_TESTFILES/vds/b.h5
$SRC_H5DUMP_TESTFILES/vds/c.h5
$SRC_H5DUMP_TESTFILES/vds/d.h5
$SRC_H5DUMP_TESTFILES/vds/vds-percival-unlim-maxmin.h5
$SRC_H5DUMP_TESTFILES/vds/f-0.h5
$SRC_H5DUMP_TESTFILES/vds/f-3.h5
$SRC_H5DUMP_TESTFILES/vds/vds-eiger.h5
"

LIST_OTHER_TEST_FILES="
$SRC_H5DUMP_TESTFILES/vds/tvds-1.ddl
$SRC_H5DUMP_TESTFILES/vds/tvds-2.ddl
$SRC_H5DUMP_TESTFILES/vds/tvds-3_1.ddl
$SRC_H5DUMP_TESTFILES/vds/tvds-3_2.ddl
$SRC_H5DUMP_TESTFILES/vds/tvds-4.ddl
$SRC_H5DUMP_TESTFILES/vds/tvds-5.ddl
$SRC_H5DUMP_TESTFILES/vds/tvds_layout-1.ddl
$SRC_H5DUMP_TESTFILES/vds/tvds_layout-2.ddl
$SRC_H5DUMP_TESTFILES/vds/tvds_layout-3_1.ddl
$SRC_H5DUMP_TESTFILES/vds/tvds_layout-3_2.ddl
$SRC_H5DUMP_TESTFILES/vds/tvds_layout-4.ddl
$SRC_H5DUMP_TESTFILES/vds/tvds_layout-5.ddl
$SRC_H5DUMP_TESTFILES/vds/vds-first.ddl
$SRC_H5DUMP_TESTFILES/vds/vds-gap1.ddl
$SRC_H5DUMP_TESTFILES/vds/vds-gap2.ddl
$SRC_H5DUMP_TESTFILES/vds/vds_layout-eiger.ddl
$SRC_H5DUMP_TESTFILES/vds/vds_layout-maxmin.ddl
"

LIST_ERROR_TEST_FILES="
"

#
# copy test files and expected output files from source dirs to test dir
#
COPY_TESTFILES="$LIST_HDF5_TEST_FILES $LIST_OTHER_TEST_FILES $LIST_ERROR_TEST_FILES"

COPY_TESTFILES_TO_TESTDIR()
{
    # copy test files. Used -f to make sure get a new copy
    for tstfile in $COPY_TESTFILES
    do
        # ignore '#' comment
        echo $tstfile | tr -d ' ' | grep '^#' > /dev/null
        RET=$?
        if [ $RET -eq 1 ]; then
            # skip cp if srcdir is same as destdir
            # this occurs when build/test performed in source dir and
            # make cp fail
            SDIR=`$DIRNAME $tstfile`
            INODE_SDIR=`$LS -i -d $SDIR | $AWK -F' ' '{print $1}'`
            INODE_DDIR=`$LS -i -d $TESTDIR | $AWK -F' ' '{print $1}'`
            if [ "$INODE_SDIR" != "$INODE_DDIR" ]; then
              $CP -f $tstfile $TESTDIR
                if [ $? -ne 0 ]; then
                    echo "Error: FAILED to copy $tstfile ."

                    # Comment out this to CREATE expected file
                    exit $EXIT_FAILURE
                fi
            fi
        fi
    done
}

CLEAN_TESTFILES_AND_TESTDIR()
{
    # skip rm if srcdir is same as destdir
    # this occurs when build/test performed in source dir and
    # make cp fail
    SDIR=$SRC_H5DUMP_TESTFILES/vds
    INODE_SDIR=`$LS -i -d $SDIR | $AWK -F' ' '{print $1}'`
    INODE_DDIR=`$LS -i -d $TESTDIR | $AWK -F' ' '{print $1}'`
    if [ "$INODE_SDIR" != "$INODE_DDIR" ]; then
        $RM $TESTDIR
    fi
}

# Print a line-line message left justified in a field of 70 characters
# beginning with the word "Testing".
#
TESTING() {
   SPACES="                                                               "
   echo "Testing $* $SPACES" | cut -c1-70 | tr -d '\012'
}

# Source in the output filter function definitions.
. $srcdir/../../../bin/output_filter.sh

# Run a test and print PASS or *FAIL*.  If a test fails then increment
# the `nerrors' global variable and (if $verbose is set) display the
# difference between the actual output and the expected output. The
# expected output is given as the first argument to this function and
# the actual output file is calculated by replacing the `.ddl' with
# `.out'.  The actual output is not removed if $HDF5_NOCLEANUP has a
# non-zero value.
#
TOOLTEST() {
    expect="$TESTDIR/$1"
    actual="$TESTDIR/`basename $1 .ddl`.out"
    actual_err="$TESTDIR/`basename $1 .ddl`.err"
    actual_sav=${actual}-sav
    actual_err_sav=${actual_err}-sav
    shift

    # Run test.
    TESTING $DUMPER $@
    (
        cd $TESTDIR
        $RUNSERIAL $DUMPER_BIN "$@"
    ) >$actual 2>$actual_err

    # save actual and actual_err in case they are needed later.
    cp $actual $actual_sav
    STDOUT_FILTER $actual
    cp $actual_err $actual_err_sav
    STDERR_FILTER $actual_err
    cat $actual_err >> $actual

    if [ ! -f $expect ]; then
        # Create the expect file if it doesn't yet exist.
        echo " CREATED"
        cp $actual $expect
        echo "    Expected result (*.ddl) missing"
        nerrors="`expr $nerrors + 1`"
    elif $CMP $expect $actual; then
        echo " PASSED"
    else
        echo "*FAILED*"
        echo "    Expected result (*.ddl) differs from actual result (*.out)"
        nerrors="`expr $nerrors + 1`"
        test yes = "$verbose" && $DIFF $expect $actual |sed 's/^/    /'
    fi

    # Clean up output file
    if test -z "$HDF5_NOCLEANUP"; then
        rm -f $actual $actual_err $actual_sav $actual_err_sav $actual_ext
    fi
}


# same as TOOLTEST1 but compares generated file to expected output
#                   and compares the generated data file to the expected data file
# used for the binary tests that expect a full path in -o without -b
TOOLTEST2() {

    expectdata="$TESTDIR/$1"
    expect="$TESTDIR/`basename $1 .exp`.ddl"
    actualdata="$TESTDIR/`basename $1 .exp`.txt"
    actual="$TESTDIR/`basename $1 .exp`.out"
    actual_err="$TESTDIR/`basename $1 .exp`.err"
    shift

    # Run test.
    TESTING $DUMPER $@
    (
        cd $TESTDIR
        $RUNSERIAL $DUMPER_BIN "$@"
    ) >$actual 2>$actual_err
    cat $actual_err >> $actual

    if [ ! -f $expect ]; then
        # Create the expect file if it doesn't yet exist.
        echo " CREATED"
        cp $actual $expect
        echo "    Expected result (*.ddl) missing"
        nerrors="`expr $nerrors + 1`"
    elif $CMP $expect $actual; then
        if [ ! -f $expectdata ]; then
            # Create the expect data file if it doesn't yet exist.
            echo " CREATED"
            cp $actualdata $expectdata
            echo "    Expected data (*.exp) missing"
            nerrors="`expr $nerrors + 1`"
        elif $CMP $expectdata $actualdata; then
            echo " PASSED"
        else
            echo "*FAILED*"
            echo "    Expected datafile (*.exp) differs from actual datafile (*.txt)"
            nerrors="`expr $nerrors + 1`"
            test yes = "$verbose" && $DIFF $expectdata $actualdata |sed 's/^/    /'
        fi
    else
        echo "*FAILED*"
        echo "    Expected result (*.ddl) differs from actual result (*.out)"
        nerrors="`expr $nerrors + 1`"
        test yes = "$verbose" && $DIFF $expect $actual |sed 's/^/    /'
    fi

    # Clean up output file
    if test -z "$HDF5_NOCLEANUP"; then
        rm -f $actual $actualdata $actual_err
    fi

}

# same as TOOLTEST but filters error stack outp
# Extract file name, line number, version and thread IDs because they may be different
TOOLTEST3() {

    expect="$TESTDIR/$1"
    actual="$TESTDIR/`basename $1 .ddl`.out"
    actual_err="$TESTDIR/`basename $1 .ddl`.err"
    actual_ext="$TESTDIR/`basename $1 .ddl`.ext"
    actual_sav=${actual}-sav
    actual_err_sav=${actual_err}-sav
    shift

    # Run test.
    TESTING $DUMPER $@
    (
        cd $TESTDIR
        $RUNSERIAL $DUMPER_BIN "$@"
    ) >$actual 2>$actual_err

    # save actual and actual_err in case they are needed later.
    cp $actual $actual_sav
    STDOUT_FILTER $actual
    cp $actual_err $actual_err_sav
    STDERR_FILTER $actual_err

    # Extract file name, line number, version and thread IDs because they may be different
    sed -e 's/thread [0-9]*/thread (IDs)/' -e 's/: .*\.c /: (file name) /' \
        -e 's/line [0-9]*/line (number)/' \
        -e 's/v[1-9]*\.[0-9]*\./version (number)\./' \
        -e 's/[1-9]*\.[0-9]*\.[0-9]*[^)]*/version (number)/' \
        -e 's/H5Eget_auto[1-2]*/H5Eget_auto(1 or 2)/' \
        -e 's/H5Eset_auto[1-2]*/H5Eset_auto(1 or 2)/' \
        $actual_err > $actual_ext
    cat $actual_ext >> $actual

    if [ ! -f $expect ]; then
        # Create the expect file if it doesn't yet exist.
        echo " CREATED"
        cp $actual $expect
        echo "    Expected result (*.ddl) missing"
        nerrors="`expr $nerrors + 1`"
    elif $CMP $expect $actual; then
        echo " PASSED"
    else
        echo "*FAILED*"
        echo "    Expected result (*.ddl) differs from actual result (*.out)"
        nerrors="`expr $nerrors + 1`"
        test yes = "$verbose" && $DIFF $expect $actual |sed 's/^/    /'
    fi

    # Clean up output file
    if test -z "$HDF5_NOCLEANUP"; then
        rm -f $actual $actual_err $actual_sav $actual_err_sav
    fi

}

# same as TOOLTEST3 but filters error stack output and compares to an error file
# Extract file name, line number, version and thread IDs because they may be different
TOOLTEST4() {

    expect="$TESTDIR/$1"
    expect_err="$TESTDIR/`basename $1 .ddl`.err"
    actual="$TESTDIR/`basename $1 .ddl`.out"
    actual_err="$TESTDIR/`basename $1 .ddl`.oerr"
    actual_ext="$TESTDIR/`basename $1 .ddl`.ext"
    actual_sav=${actual}-sav
    actual_err_sav=${actual_err}-sav
    shift

    # Run test.
    TESTING $DUMPER $@
    (
        cd $TESTDIR
        $RUNSERIAL $DUMPER_BIN "$@"
    ) >$actual 2>$actual_err

    # save actual and actual_err in case they are needed later.
    cp $actual $actual_sav
    STDOUT_FILTER $actual
    cp $actual_err $actual_err_sav
    STDERR_FILTER $actual_err

    # Extract file name, line number, version and thread IDs because they may be different
    sed -e 's/thread [0-9]*/thread (IDs)/' -e 's/: .*\.c /: (file name) /' \
        -e 's/line [0-9]*/line (number)/' \
        -e 's/v[1-9]*\.[0-9]*\./version (number)\./' \
        -e 's/[1-9]*\.[0-9]*\.[0-9]*[^)]*/version (number)/' \
        -e 's/H5Eget_auto[1-2]*/H5Eget_auto(1 or 2)/' \
        -e 's/H5Eset_auto[1-2]*/H5Eset_auto(1 or 2)/' \
        $actual_err > $actual_ext
    #cat $actual_ext >> $actual

    if [ ! -f $expect ]; then
        # Create the expect file if it doesn't yet exist.
        echo " CREATED"
        cp $actual $expect
    elif $CMP $expect $actual; then
        if $CMP $expect_err $actual_ext; then
            echo " PASSED"
        else
            echo "*FAILED*"
            echo "    Expected result (*.err) differs from actual result (*.oerr)"
            nerrors="`expr $nerrors + 1`"
            test yes = "$verbose" && $DIFF $expect_err $actual_ext |sed 's/^/    /'
        fi
    else
        echo "*FAILED*"
        echo "    Expected result (*.ddl) differs from actual result (*.out)"
        nerrors="`expr $nerrors + 1`"
        test yes = "$verbose" && $DIFF $expect $actual |sed 's/^/    /'
    fi

    # Clean up output file
    if test -z "$HDF5_NOCLEANUP"; then
        rm -f $actual $actual_err $actual_sav $actual_err_sav
    fi
}

# Print a "SKIP" message
SKIP() {
    TESTING $DUMPER $@
    echo  " -SKIP-"
}

# Print a line-line message left justified in a field of 70 characters
#
PRINT_H5DIFF() {
    SPACES="                                                               "
    echo " Running h5diff $* $SPACES" | cut -c1-70 | tr -d '\012'
}


# Call the h5diff tool
#
DIFFTEST()
{
    PRINT_H5DIFF  $@
    (
        cd $TESTDIR
        $RUNSERIAL $H5DIFF_BIN "$@" -q
    )
    RET=$?

    if [ $RET != 0 ] ; then
         echo "*FAILED*"
         nerrors="`expr $nerrors + 1`"
    else
         echo " PASSED"
    fi
}

# Print a line-line message left justified in a field of 70 characters
# beginning with the word "Verifying".
#
PRINT_H5IMPORT() {
    SPACES="                                                               "
    echo " Running h5import $* $SPACES" | cut -c1-70 | tr -d '\012'
}

# Call the h5import tool
#
IMPORTTEST()
{
    # remove the output hdf5 file if it exists
    hdf5_file="$TESTDIR/$5"
    if [ -f $hdf5_file ]; then
        rm -f $hdf5_file
    fi

    PRINT_H5IMPORT  $@
    (
        cd $TESTDIR
        $RUNSERIAL $H5IMPORT_BIN "$@"
    )
    RET=$?

    if [ $RET != 0 ] ; then
         echo "*FAILED*"
         nerrors="`expr $nerrors + 1`"
    else
         echo " PASSED"
    fi
}


##############################################################################
##############################################################################
###        T H E   T E S T S                                            ###
##############################################################################
##############################################################################
# prepare for test
COPY_TESTFILES_TO_TESTDIR

####### test for dataset vds ######

# Data read
if test $USE_FILTER_DEFLATE = "yes" ; then
    TOOLTEST tvds-1.ddl --enable-error-stack 1_vds.h5
    TOOLTEST tvds-2.ddl --enable-error-stack 2_vds.h5
    TOOLTEST tvds-3_1.ddl --enable-error-stack 3_1_vds.h5
    TOOLTEST tvds-3_2.ddl --enable-error-stack 3_2_vds.h5
    TOOLTEST tvds-4.ddl --enable-error-stack 4_vds.h5
    TOOLTEST tvds-5.ddl --enable-error-stack 5_vds.h5
    TOOLTEST vds-first.ddl --vds-view-first-missing --enable-error-stack vds-percival-unlim-maxmin.h5
    TOOLTEST vds-gap1.ddl -d /VDS-Eiger --vds-gap-size=1 --enable-error-stack vds-eiger.h5
    TOOLTEST vds-gap2.ddl --vds-gap-size=2 --enable-error-stack vds-eiger.h5
fi

# Layout read
if test $USE_FILTER_DEFLATE = "yes" ; then
    TOOLTEST tvds_layout-1.ddl -p --enable-error-stack 1_vds.h5
    TOOLTEST tvds_layout-2.ddl -p --enable-error-stack 2_vds.h5
    TOOLTEST tvds_layout-3_1.ddl -p --enable-error-stack 3_1_vds.h5
    TOOLTEST tvds_layout-3_2.ddl -p --enable-error-stack 3_2_vds.h5
    TOOLTEST tvds_layout-4.ddl -p --enable-error-stack 4_vds.h5
    TOOLTEST tvds_layout-5.ddl -p --enable-error-stack 5_vds.h5
    TOOLTEST vds_layout-eiger.ddl -p --enable-error-stack vds-eiger.h5
    TOOLTEST vds_layout-maxmin.ddl -p --enable-error-stack vds-percival-unlim-maxmin.h5
fi

# Clean up temporary files/directories
CLEAN_TESTFILES_AND_TESTDIR

# Report test results and exit
if test $nerrors -eq 0 ; then
    echo "All $TESTNAME tests passed."
    exit $EXIT_SUCCESS
else
    echo "$TESTNAME tests failed with $nerrors errors."
    exit $EXIT_FAILURE
fi
