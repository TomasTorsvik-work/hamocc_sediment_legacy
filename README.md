# hamocc_sediment_legacy
This is a legacy sediment code for HAMOCC, converted from a mercurial repository (http://hg.gnu.org.ua/hgweb/burst-coupling). This source coude should be merged into the [BLOM](https://github.com/NorESMhub/BLOM) repository via a feature branch. __Please do not add further content directly to this repository.__

## Conversion of source code from mercurial to git
The procedure for converting the mercurial repository to git follows the outline in the [Migrating to Git](https://git-scm.com/book/en/v2/Git-and-Other-Systems-Migrating-to-Git) chapter from the [Pro Git book](https://git-scm.com/book/en/v2), written by Scott Chacon and Ben Straub and published by Apress.

The conversion was done by use of the [hg-fast-export](https://github.com/frej/fast-export.git) tool (version [44c50d0](https://github.com/frej/fast-export/tree/44c50d0fae5697ab9cdb05d0e484b97edda31042)).

__Essential steps:__

     $ cd /tmp
     $ git clone https://github.com/frej/fast-export.git
     $ hg clone http://hg.gnu.org.ua/hgweb/burst-coupling /tmp/hg-repo
     $ cd /tmp/hg-repo
     $ hg log | grep user: | sort | uniq | sed 's/user: *//' > ../authors
     $ git init /tmp/converted
     $ cd /tmp/converted
     $ /tmp/fast-export/hg-fast-export.sh -r /tmp/hg-repo -A /tmp/authors
     $ git remote add origin https://github.com/TomasTorsvik-work/hamocc_sediment_legacy.git
     $ git push origin --all

## LICENSE
The source code is licensed under the [GPLv3 license](https://github.com/TomasTorsvik-work/hamocc_sediment_legacy/blob/master/gpl.txt). Permision has been granted by Marco van Hulten to integrate this code into BLOM that is licensed under LGPLv3.
