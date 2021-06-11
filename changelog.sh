#!/usr/bin/env bash
#
# Generates change log between two tags specified as arguments, e.g.:
#
#  ./changelog.sh 0.19.8 0.20.0
#
# will list changes between these two tag, or from a given tag and the
# present, e.g.
#
#  ./changelog.sh 0.20.8

if [ "$#" -eq 2 ]; then

    tag1=$1
    tag2=$2

    clog=docs/change_log_${tag1}_${tag2}.rst

    cat<<EOF > $clog
.. changelog created on `date`

.. include:: globals.rst

|hiper| pipeline change log from ${tag1} to ${tag2}
***************************************************

List of changes from git, newest first, with the commit keys linked to github:

EOF

    git log ${tag1}..${tag2} --pretty=format:'  * `%H <https://github.com/HiPERCAM/hipercam/commit/%H>`_ %s' >> $clog

    echo 'Written '$clog

elif [ "$#" -eq 1 ]; then

    tag=$1

    clog=docs/change_log_${tag}_to_now.rst

    cat<<EOF > $clog
.. changelog created on `date`

.. include:: globals.rst

|hiper| pipeline change log from ${tag2}
****************************************

List of changes from git, newest first, with the commit keys linked to github:

EOF

    git log ${tag}..HEAD --pretty=format:'  * `%H <https://github.com/HiPERCAM/hipercam/commit/%H>`_ %s' >> $clog

    echo 'Written '$clog

else

    echo "changelog takes either 1 or 2 arguments"

fi
