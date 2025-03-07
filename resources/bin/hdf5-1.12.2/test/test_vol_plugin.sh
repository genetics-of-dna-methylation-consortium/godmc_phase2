#! /bin/sh
#
# Copyright by The HDF Group.
# All rights reserved.
#
# This file is part of HDF5. The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the COPYING file, which can be found at the root of the source code
# distribution tree, or in https://www.hdfgroup.org/licenses.
# If you do not have access to either file, you may request a copy from
# help@hdfgroup.org.
#
# This shell script is for testing VOL connector plugins.
#
srcdir=.
TOP_BUILDDIR=..

EXIT_SUCCESS=0
EXIT_FAILURE=1

nerrors=0
verbose=yes
exit_code=$EXIT_SUCCESS

TEST_NAME=vol_plugin
TEST_BIN=`pwd`/$TEST_NAME
FROM_DIR=`pwd`/.libs
case $(uname) in
    CYGWIN* )
        NULL_VOL_PLUGIN="$FROM_DIR/cygnull_vol_connector*"
        ;;
    *)
        NULL_VOL_PLUGIN="$FROM_DIR/libnull_vol_connector*"
        ;;
esac
TEMP_PLUGIN_DIR=null_vol_plugin_dir
CP="cp -p"    # Use -p to preserve mode,ownership, timestamps
RM="rm -rf"

# Print a line-line message left justified in a field of 70 characters
# beginning with the word "Testing".
#
TESTING() {
    SPACES="                                                               "
    echo "Testing $* $SPACES" | cut -c1-70 | tr -d '\012'
}

# Main Body
# Create test directory if necessary.
test -d $TEMP_PLUGIN_DIR || mkdir -p $TEMP_PLUGIN_DIR
if [ $? != 0 ]; then
    echo "Failed to create VOL connector plugin test directory ($TEMP_PLUGIN_DIR)"
    exit $EXIT_FAILURE
fi

# Copy plugin for the tests.
$CP $NULL_VOL_PLUGIN $TEMP_PLUGIN_DIR
if [ $? != 0 ]; then
    echo "Failed to copy NULL VOL plugin ($NULL_VOL_PLUGIN) to test directory."
    exit $EXIT_FAILURE
fi

# setup plugin path
ENVCMD="env HDF5_PLUGIN_PATH=${TEMP_PLUGIN_DIR}"

# Run the test
$ENVCMD $TEST_BIN
if [ $? != 0 ]; then
    nerrors=`expr $nerrors + 1`
fi

# print results
if test $nerrors -ne 0 ; then
    echo "$nerrors errors encountered"
    exit_code=$EXIT_FAILURE
else
    echo "All VOL plugin tests passed."
    exit_code=$EXIT_SUCCESS
fi

# Clean up temporary files/directories and leave
$RM $TEMP_PLUGIN_DIR

exit $exit_code
