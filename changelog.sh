#!/usr/bin/env bash
#
# Generates change log between two tags specified as arguments
#
# Run e.g. as ./changelog.sh 0.19.8 0.20.0

tag1=$1
tag2=$2

clog=docs/change_log_${tag1}_${tag2}.rst

cat<<EOF > $clog
.. changelog created on `date`

.. include:: globals.rst

|hiper| pipeline change log from ${tag1} to ${tag2}
***************************************************

List of changes from git, oldest first:

EOF

git log ${tag1}..${tag2} --pretty=format:'  * `view commit &bull; <https://github.com/HiPERCAM/hipercam/commit/%H>`_ %s' --reverse >> $clog
