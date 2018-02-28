#!/bin/bash

# Run the TCS generated events through the experimental setup, acceptance tag.

echo $1

cd /home/vardan/g4-work/g4.10.03/
source setup_g4.sh

cd tcs/aux/tcsgen_data.test

# This part does not work.

#/usr/local/bin/root -b <<EOF
#{
#gSystem->Load("/home/vardan/g4-work/g4.10.03/tcs/aux/tcsgen_data.test/deepgen.C");
#      deepgen d;
##      d.Loop();
#}
#EOF

#/usr/local/bin/root -b -q run_deepgen.C'(''"'RootFiles/$1.root'"'')'

cd ../../tcs_setup-build
#./tcs run1.mac > out.run1

cd ../aux/acceptance_tag
/usr/local/bin/root -b -q ./tcs_ana.C'("'RootFiles.tracked/$1_tracked.root'"','"'RootFiles/$1.root'"','"'$1_tagged.root'")'
