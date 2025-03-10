#! /bin/sh
#
# Copyright by The HDF Group.
# All rights reserved.
#
# This file is part of HDF5.  The full HDF5 copyright notice, including
# terms governing use, modification, and redistribution, is contained in
# the COPYING file, which can be found at the root of the source code
# distribution tree, or in https://www.hdfgroup.org/licenses.
# If you do not have access to either file, you may request a copy from
# help@hdfgroup.org.
#

top_builddir=../../..
top_srcdir=../../..
srcdir=.
IS_DARWIN="no"

TESTNAME=EX_Intro
EXIT_SUCCESS=0
EXIT_FAILURE=1

# Set up default variable values if not supplied by the user.
RM='rm -rf'
CMP='cmp'
DIFF='diff -c'
CP='cp'
DIRNAME='dirname'
LS='ls'
AWK='awk'

nerrors=0

# where the libs exist
HDFLIB_HOME="$top_srcdir/java/lib"
BLDDIR="."
BLDLIBDIR="$BLDDIR/testlibs"
HDFTEST_HOME="$top_srcdir/java/examples/intro"
JARFILE=jarhdf5-1.12.2.jar
TESTJARFILE=jarhdf5intro.jar
test -d $BLDLIBDIR || mkdir -p $BLDLIBDIR

######################################################################
# library files
# --------------------------------------------------------------------
# All the library files copy from source directory to test directory
# NOTE: Keep this framework to add/remove test files.
#       This list are also used for checking exist.
#       Comment '#' without space can be used.
# --------------------------------------------------------------------
LIST_LIBRARY_FILES="
$top_builddir/src/.libs/libhdf5.*
$top_builddir/java/src/jni/.libs/libhdf5_java.*
$top_builddir/java/src/$JARFILE
"
LIST_JAR_TESTFILES="
$HDFLIB_HOME/slf4j-api-1.7.33.jar
$HDFLIB_HOME/ext/slf4j-simple-1.7.33.jar
"
LIST_DATA_FILES="
$HDFTEST_HOME/../testfiles/examples.intro.H5_CreateDataset.txt
$HDFTEST_HOME/../testfiles/examples.intro.H5_CreateAttribute.txt
$HDFTEST_HOME/../testfiles/examples.intro.H5_CreateFile.txt
$HDFTEST_HOME/../testfiles/examples.intro.H5_CreateGroup.txt
$HDFTEST_HOME/../testfiles/examples.intro.H5_CreateGroupAbsoluteRelative.txt
$HDFTEST_HOME/../testfiles/examples.intro.H5_CreateGroupDataset.txt
$HDFTEST_HOME/../testfiles/examples.intro.H5_ReadWrite.txt
"

#
# copy files from source dirs to test dir
#
COPY_LIBFILES="$LIST_LIBRARY_FILES"
COPY_JARTESTFILES="$LIST_JAR_TESTFILES"

COPY_LIBFILES_TO_BLDLIBDIR()
{
    # copy test files. Used -f to make sure get a new copy
    for tstfile in $COPY_LIBFILES
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
            INODE_DDIR=`$LS -i -d $BLDLIBDIR | $AWK -F' ' '{print $1}'`
            if [ "$INODE_SDIR" != "$INODE_DDIR" ]; then
                $CP -f $tstfile $BLDLIBDIR
                if [ $? -ne 0 ]; then
                    echo "Error: FAILED to copy $tstfile ."

                    # Comment out this to CREATE expected file
                    exit $EXIT_FAILURE
                fi
            fi
        fi
    done
    if [ "$IS_DARWIN" = "yes" ]; then
       (cd $BLDLIBDIR; \
         install_name_tool -add_rpath @loader_path libhdf5_java.dylib; \
         exist_path=` otool -l libhdf5_java.dylib | grep libhdf5 | grep -v java | awk '{print $2}'`; \
         echo $exist_path; \
         install_name_tool -change $exist_path @rpath/libhdf5.dylib libhdf5_java.dylib)
    fi
    # copy jar files. Used -f to make sure get a new copy
    for tstfile in $COPY_JARTESTFILES
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
            INODE_DDIR=`$LS -i -d $BLDLIBDIR | $AWK -F' ' '{print $1}'`
            if [ "$INODE_SDIR" != "$INODE_DDIR" ]; then
                $CP -f $tstfile $BLDLIBDIR
                if [ $? -ne 0 ]; then
                    echo "Error: FAILED to copy $tstfile ."

                    # Comment out this to CREATE expected file
                    exit $EXIT_FAILURE
                fi
            fi
        fi
    done
}

CLEAN_LIBFILES_AND_BLDLIBDIR()
{
    # skip rm if srcdir is same as destdir
    # this occurs when build/test performed in source dir and
    # make cp fail
    SDIR=$HDFLIB_HOME
    INODE_SDIR=`$LS -i -d $SDIR | $AWK -F' ' '{print $1}'`
    INODE_DDIR=`$LS -i -d $BLDLIBDIR | $AWK -F' ' '{print $1}'`
    if [ "$INODE_SDIR" != "$INODE_DDIR" ]; then
        for tstfile in $COPY_JARTESTFILES
        do
            $RM $BLDLIBDIR/tstfile
        done
    fi
}

COPY_DATAFILES="$LIST_DATA_FILES"

COPY_DATAFILES_TO_BLDDIR()
{
    # copy test files. Used -f to make sure get a new copy
    for tstfile in $COPY_DATAFILES
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
            INODE_DDIR=`$LS -i -d $BLDDIR | $AWK -F' ' '{print $1}'`
            if [ "$INODE_SDIR" != "$INODE_DDIR" ]; then
                $CP -f $tstfile $BLDDIR
                if [ $? -ne 0 ]; then
                    echo "Error: FAILED to copy $tstfile ."

                    # Comment out this to CREATE expected file
                    exit $EXIT_FAILURE
                fi
            fi
        fi
    done
}

CLEAN_DATAFILES_AND_BLDDIR()
{
        $RM $BLDDIR/examples.intro.H5_*.txt
        $RM $BLDDIR/H5_*.out
        $RM $BLDDIR/H5_*.h5
}

# Print a line-line message left justified in a field of 70 characters
# beginning with the word "Testing".
#
TESTING() {
   SPACES="                                                               "
   echo "Testing $* $SPACES" | cut -c1-70 | tr -d '\012'
}

# where Java is installed (requires jdk1.7.x)
JAVAEXE=
JAVAEXEFLAGS=

###############################################################################
#            DO NOT MODIFY BELOW THIS LINE
###############################################################################

# prepare for test
COPY_LIBFILES_TO_BLDLIBDIR
COPY_DATAFILES_TO_BLDDIR

CPATH=".:"$BLDLIBDIR"/"$JARFILE":"$BLDLIBDIR"/slf4j-api-1.7.33.jar:"$BLDLIBDIR"/slf4j-simple-1.7.33.jar:"$TESTJARFILE""

TEST=/usr/bin/test
if [ ! -x /usr/bin/test ]
then
TEST=`which test`
fi

if $TEST -z "$CLASSPATH"; then
        CLASSPATH=""
fi
CLASSPATH=$CPATH":"$CLASSPATH
export CLASSPATH

if $TEST -n "$JAVAPATH" ; then
        PATH=$JAVAPATH":"$PATH
        export PATH
fi

if $TEST -e /bin/uname; then
   os_name=`/bin/uname -s`
elif $TEST -e /usr/bin/uname; then
   os_name=`/usr/bin/uname -s`
else
   os_name=unknown
fi

if $TEST -z "$LD_LIBRARY_PATH" ; then
        LD_LIBRARY_PATH=""
fi

case  $os_name in
    *)
    LD_LIBRARY_PATH=$BLDLIBDIR:$LD_LIBRARY_PATH
    ;;
esac

export LD_LIBRARY_PATH

echo "$JAVAEXE $JAVAEXEFLAGS -Xmx1024M -Dorg.slf4j.simpleLogger.defaultLog=trace -Djava.library.path=$BLDLIBDIR -cp $CLASSPATH examples.intro.H5_CreateDataset"
TESTING examples.intro.H5_CreateDataset
($RUNSERIAL $JAVAEXE $JAVAEXEFLAGS -Xmx1024M -Dorg.slf4j.simpleLogger.defaultLog=trace -Djava.library.path=$BLDLIBDIR -cp $CLASSPATH examples.intro.H5_CreateDataset > H5_CreateDataset.out)
if diff H5_CreateDataset.out examples.intro.H5_CreateDataset.txt > /dev/null; then
    echo "  PASSED      intro.H5_CreateDataset"
else
    echo "**FAILED**    intro.H5_CreateDataset"
    nerrors="`expr $nerrors + 1`"
fi

echo "$JAVAEXE $JAVAEXEFLAGS -Xmx1024M -Dorg.slf4j.simpleLogger.defaultLog=trace -Djava.library.path=$BLDLIBDIR -cp $CLASSPATH examples.intro.H5_CreateAttribute"
TESTING examples.intro.H5_CreateAttribute
($RUNSERIAL $JAVAEXE $JAVAEXEFLAGS -Xmx1024M -Dorg.slf4j.simpleLogger.defaultLog=trace -Djava.library.path=$BLDLIBDIR -cp $CLASSPATH examples.intro.H5_CreateAttribute > H5_CreateAttribute.out)
if diff H5_CreateAttribute.out examples.intro.H5_CreateAttribute.txt > /dev/null; then
    echo "  PASSED      intro.H5_CreateAttribute"
else
    echo "**FAILED**    intro.H5_CreateAttribute"
    nerrors="`expr $nerrors + 1`"
fi

echo "$JAVAEXE $JAVAEXEFLAGS -Xmx1024M -Dorg.slf4j.simpleLogger.defaultLog=trace -Djava.library.path=$BLDLIBDIR -cp $CLASSPATH examples.intro.H5_CreateFile"
TESTING examples.intro.H5_CreateFile
($RUNSERIAL $JAVAEXE $JAVAEXEFLAGS -Xmx1024M -Dorg.slf4j.simpleLogger.defaultLog=trace -Djava.library.path=$BLDLIBDIR -cp $CLASSPATH examples.intro.H5_CreateFile > H5_CreateFile.out)
if diff H5_CreateFile.out examples.intro.H5_CreateFile.txt > /dev/null; then
    echo "  PASSED      intro.H5_CreateFile"
else
    echo "**FAILED**    intro.H5_CreateFile"
    nerrors="`expr $nerrors + 1`"
fi

echo "$JAVAEXE $JAVAEXEFLAGS -Xmx1024M -Dorg.slf4j.simpleLogger.defaultLog=trace -Djava.library.path=$BLDLIBDIR -cp $CLASSPATH examples.intro.H5_CreateGroup"
TESTING examples.intro.H5_CreateGroup
($RUNSERIAL $JAVAEXE $JAVAEXEFLAGS -Xmx1024M -Dorg.slf4j.simpleLogger.defaultLog=trace -Djava.library.path=$BLDLIBDIR -cp $CLASSPATH examples.intro.H5_CreateGroup > H5_CreateGroup.out)
if diff H5_CreateGroup.out examples.intro.H5_CreateGroup.txt > /dev/null; then
    echo "  PASSED      intro.H5_CreateGroup"
else
    echo "**FAILED**    intro.H5_CreateGroup"
    nerrors="`expr $nerrors + 1`"
fi

echo "$JAVAEXE $JAVAEXEFLAGS -Xmx1024M -Dorg.slf4j.simpleLogger.defaultLog=trace -Djava.library.path=$BLDLIBDIR -cp $CLASSPATH examples.intro.H5_CreateGroupAbsoluteRelative"
TESTING examples.intro.H5_CreateGroupAbsoluteRelative
($RUNSERIAL $JAVAEXE $JAVAEXEFLAGS -Xmx1024M -Dorg.slf4j.simpleLogger.defaultLog=trace -Djava.library.path=$BLDLIBDIR -cp $CLASSPATH examples.intro.H5_CreateGroupAbsoluteRelative > H5_CreateGroupAbsoluteRelative.out)
if diff H5_CreateGroupAbsoluteRelative.out examples.intro.H5_CreateGroupAbsoluteRelative.txt > /dev/null; then
    echo "  PASSED      intro.H5_CreateGroupAbsoluteRelative"
else
    echo "**FAILED**    intro.H5_CreateGroupAbsoluteRelative"
    nerrors="`expr $nerrors + 1`"
fi

echo "$JAVAEXE $JAVAEXEFLAGS -Xmx1024M -Dorg.slf4j.simpleLogger.defaultLog=trace -Djava.library.path=$BLDLIBDIR -cp $CLASSPATH examples.intro.H5_CreateGroupDataset"
TESTING examples.intro.H5_CreateGroupDataset
($RUNSERIAL $JAVAEXE $JAVAEXEFLAGS -Xmx1024M -Dorg.slf4j.simpleLogger.defaultLog=trace -Djava.library.path=$BLDLIBDIR -cp $CLASSPATH examples.intro.H5_CreateGroupDataset > H5_CreateGroupDataset.out)
if diff H5_CreateGroupDataset.out examples.intro.H5_CreateGroupDataset.txt > /dev/null; then
    echo "  PASSED      intro.H5_CreateGroupDataset"
else
    echo "**FAILED**    intro.H5_CreateGroupDataset"
    nerrors="`expr $nerrors + 1`"
fi

echo "$JAVAEXE $JAVAEXEFLAGS -Xmx1024M -Dorg.slf4j.simpleLogger.defaultLog=trace -Djava.library.path=$BLDLIBDIR -cp $CLASSPATH examples.intro.H5_ReadWrite"
TESTING examples.intro.H5_ReadWrite
($RUNSERIAL $JAVAEXE $JAVAEXEFLAGS -Xmx1024M -Dorg.slf4j.simpleLogger.defaultLog=trace -Djava.library.path=$BLDLIBDIR -cp $CLASSPATH examples.intro.H5_ReadWrite > H5_ReadWrite.out)
if diff H5_ReadWrite.out examples.intro.H5_ReadWrite.txt > /dev/null; then
    echo "  PASSED      intro.H5_ReadWrite"
else
    echo "**FAILED**    intro.H5_ReadWrite"
    nerrors="`expr $nerrors + 1`"
fi

# Clean up temporary files/directories
CLEAN_LIBFILES_AND_BLDLIBDIR
CLEAN_DATAFILES_AND_BLDDIR

# Report test results and exit
if test $nerrors -eq 0 ; then
    echo "All $TESTNAME tests passed."
    exit $EXIT_SUCCESS
else
    echo "$TESTNAME tests failed with $nerrors errors."
    exit $EXIT_FAILURE
fi
