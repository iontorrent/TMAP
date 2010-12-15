#!/bin/sh

rm -vr tmap-0*;

sh autogen.sh && ./configure && make dist;

TMAP=`ls -1 tmap-0.*.tar.gz | sed 's_.tar.gz__'`;
tar -zxvf ${TMAP}.tar.gz;

cd ${TMAP};

dh_make -e nils.homer@lifetech.com --createorig --single -c gpl2;

cd debian;

perl ../../scripts/debian/debian_mod_files.pl;

debuild -us -uc;

cd ../..;
