# LimitLFV13TeV

- set proper permission to all scripts by doing -> chmod u+x filename

- Create directories named root_out and txt_out.

- Run create_input_histos.cxx like this :
      root -l
      .x create_input_histos.cxx(8);
  This will create root files and txt files inside root_out and txt_out directories respectively.

- If not already present, create a directory named workdir
- Run organise_input_histos.sh like this :
      ./organise_input_histos.sh shape workdir

- cd workdir/ and check by submitting one condor job. If everything works fine then 
  submit all jobs ->  ./run_condor.sh

- When all jobs are done, combine all root files -> ./helper_get_observed_limit.sh

- Plot the final limit plot :
      root -l
      .x get_expected_limit.C("workdir")



Warning : Some paths, mass points and other things are hardcoded. Check all scripts before running. 