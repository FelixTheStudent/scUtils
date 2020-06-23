## Resubmission

This is a second resubmission. In this version I have:

* removed if(interactive()){} completely as suggested
* removed \dontrun and the example completely. This is going further
  than replacing it with \donttest as suggested, because the example is not relevant and
  not worth adding the microbenchmark package as a dependency. Thanks for the hint!


## Former resubmissions
This is a resubmission. In this version I have:

* removed inappropriate letter capitalization in DESCRIPTION file.
* replaced inappropriate quotes by undirected single quotes in DESCRIPTION file.
* removed typos and lack of clarity from function documentation in files
  plots.R and variance.R.
* removed \dontrun{} from all instances except one (where it's actually necessary)
* added Simon Anders as contributor to DESCRIPTION file.

I tested all changes again with devtools::check_rhub() and devtools::check_win_release().


## Test environments
* win-builder Windows, R-devel 
* R-hub builder, Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* R-hub builder, Ubuntu Linux 16.04 LTS, R-release, GCC
* R-hub builder, Fedora Linux, R-devel, clang, gfortran
* local ubuntu 18.04, R 4.0.0

## R CMD check results
There were no ERRORs and no WARNINGs.

There was 1 NOTE.

  * Maintainer: 'Felix Frauhammer <felixwertek@gmail.com>'
  New submission
  
  This is my first CRAN package.
  The email address is on my birth name, Wertek. I have used it
  for more than a decade and will continue using it perpetually.


## Downstream dependencies
There are currently no downstream dependencies for this package.
