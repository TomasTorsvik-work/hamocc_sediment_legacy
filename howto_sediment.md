How to use the sediment burst-coupling method in NorESM--OCL
=============================

The current version of this document is at
`oslo-s.uib.no:/vol/oslos/NAS29/mhu027/dokumentado/NorESM`.

This is a guide how to use the burst-coupling method.
It is tested on fram; adjust paths and other details to your environment.
We will assume that NorESM is set up properly.
If not done yet, set the shell variable `CASEROOT` to your case directory.
Such variable names are [CIME](http://esmci.github.io/cime/) convention.
If you use a Bourne shell compatible shell, do something like:

    export CASEROOT="${HOME}/NorESM/cases/NOINYOC_test"

At the moment, the most recent sediment code (*default* branch) can be found in
`fram:/cluster/projects/nn2980k/hulten/NorESM/sediment_code/` (or
`sediment_code_stable/` if you want the *stable* branch).
Check out the sediment code to your `${CASEROOT}/SourceMods/src.micom` or copy
over the files from the repository's working directory to `src.micom/`.
Another way could be:

    cd ~/src/
    module load Mercurial/4.3.3-intel-2017a-Python-2.7.13
    hg clone /cluster/projects/nn2980k/hulten/NorESM/sediment_code
    cd ${CASEROOT}/SourceMods/
    rmdir src.micom
    ln -s ~/src/sediment_code src.micom

Compile the model with the new SourceMods, and run for a couple of days
with daily output.
If that works, add the C preprocessor (cpp) directive `SED_OFFLINE` to include
the sediment burst-coupling code in the compilation.
I have this on line 87 of `${CASEROOT}/Buildconf/micom.buildexe.csh`:

    set cpp_ocn = "$cpp_ocn -DHAMOCC -DRESTART_BGC -DWLIN -DSED_OFFLINE"

During runtime, five additional namelist variables will be read.
I have this on lines 162--170 of `${CASEROOT}/Buildconf/micom.buildnml.csh`:

    &BGCNML
       ATM_CO2  = 278.000
      ,maxyear_sediment = 25
      ,maxyear_ocean    = 1
      ,nburst_last      = 0
      ,lsed_rclim  = .false.
      ,lsed_wclim  = .false.
      ,lsed_spinup = .false.
    /

Clean, compile and run for the same model duration as before.
The outputs should be identical, because---even though more code is compiled
(larger binary `ccsm.exe`)---the code is bypassed during runtime.
An explanation of the variables:

+ `ATM_CO2` is part of standard HAMOCC;
+ `maxyear_sediment` denotes the number of model years for the sediment
component (decoupled from the water column);
+ `maxyear_ocean` denotes the number of model years for the full MICOM/HAMOCC
model to run;
+ `nburst_last` should be set to burst iteration `nburst` from the end of the
previous simulation, if you your a restart from a burst coupling run (for
startup or the first burst coupling simulation use 0);
+ `lsed_rclim` is a logical/bool that tells us if we want to read a bottom-water
climatology from disk;
+ `lsed_wclim` tells us if a bottom-water climatology is written to disk during
full MICOM/HAMOCC;
+ `lsed_spinup` tells us if we want to run the sediment decoupled (spin up the
sediment).

Stand-alone sediment and full MICOM/HAMOCC are repeated until the last timestep
of MICOM/HAMOCC has been done (as defined by `STOP_N` and `STOP_OPTION` in
`${CASEROOT}/env_run.xml`).
We keep track of the iteration with the variable `nburst`.

\vskip -1cm
![Flow chart](SED_flow-crop.pdf)

In the flowchart we first check if the bottom water climatology should be read
("read forcing?"); this is decided by `lsed_rclim`.
If so, the climatology is read to be used as a forcing to spin up the sediment.
If not, we first need to create the climatology by running MICOM/HAMOCC for one
year and collect the relevant bottom-water variables.
Optionally, the climatology may be written to disk ("write forcing?"); this is
decided by `lsed_wclim`.
If `lsed_spinup .and. lcompleted_clim` ("spin up sediment?"), the sediment is
ran stand-alone, forced by the just calculated or read climatology.
When the sediment spin-up is done (after `maxyear_sediment`), it continues with
the full MICOM/HAMOCC.
This lasts for `maxyear_ocean` years; in the last year a new climatology is
built that is used as a forcing for the next stand-alone sediment simulation.
After each stand-alone sediment run, the date is restored to the value just
before the start of this run, implying that your whole job runs until
MICOM/HAMOCC has finished the last timestep of the model period defined by
`env_run.xml`.

![Burst coupling](burst-crop.pdf)

The schematic depth--time plot shows a simulation with `nburst = 4`
iterations.
It corresponds with `maxyear_sediment = 500` years of sediment and
`maxyear_ocean = 50000` years of MICOM/HAMOCC simulation.
Since `nburst` is a counter of which its maximum cannot be set explicitely, one
must set $4\times 500 = 2000$ years for the (coupled) model time in `env_run.xml`.

At the end of a job the sediment spin-up output is stored in the same directory
as the MICOM and HAMOCC output.
The stand-alone sediment output files are formatted as
`*.micom.hsedy.II.YYYYYY.nc` and `*.micom.hsedm.II.YYYYYY-MM.nc`, where II stands
for the iteration (`nburst`), YYYYYY for year and MM for month.
There can be similar filenames for decadal, centenial and millenial outputs.
The output frequency and variables that are to be stored for the
off-line sediment simulation are defined in a separate namelist group at
the bottom of `micom.buildnml.csh`, akin to the group for HAMOCC.
Only sediment variables are well defined, so it is not useful to store
water column variables for the off-line sediment.

## Caveats/bugs

At the time of writing there are several issues with this burst coupling.

- Output files for the off-line sediment sub-simulations don't have a
  correct time axis.
- Modelled (Julian) date is not consistent over off- and on-line
  sub-simulations; it is designed this way.
- One may expect `STOP_N/maxyear_ocean` iterations, but there will be,
  in fact, one more off-line sub-simulation; NorESM doesn't seem to
  finalise immediately after the last ocean timestep.
- ~~bug/FIXME: HAMOCC variables are not correctly passed between off-line
  and full!  (see note in calcurse 11 February 2019)~~
- Ocean--sediment fluxes stored to the bottom-water climatology are not
  identical to the same fields stored by HAMOCC.

