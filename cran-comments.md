## Second resubmission

  1. **Win-builder** reported:
  
     ````markdown
     Flavor: r-devel-windows-ix86+x86_64
     Check: Overall checktime, Result: NOTE
     Overall checktime 11 min > 10 min
     ````
     **Uwe Ligges** asked: Pls reduce the overall check time a  bit.
    
     - **ACTION:** we reduced the check time to 538 seconds (~9 min.) according to 
                **Win-builder** <https://win-builder.r-project.org/>

Thank you for reviewing our submission!

## Fisrt resubmission

This is a first resubmission of a new package, responding to feedback from 
  **Gregor Seyer** on initial submission.

---

  1. You have examples for unexported functions.
  
     ````markdown
       `poly_cross_simulate()` in:
          sim_cross_one_informative_parent.Rd
       `poly_hmm_est()` in:
          sim_cross_two_informative_parents.Rd
     ````
     Please either omit these examples or export the functions.
    
     - **ACTION:** we ommited the examples.

---

  2. `\dontrun{}` should only be used if the example really cannot be executed 
    (e.g. because of missing additional software, missing API keys, ...) by 
    the user. That's why wrapping examples in `\dontrun{}` adds the comment 
    ('# Not run:') as a warning for the user.
    Does not seem necessary.
    Please unwrap the examples if they are executable in < 5 sec, or replace 
    `\dontrun{}` with `\donttest{}`.
    
      - **ACTION:** we provided smaller examples and removed all `\dontrun{}`. 
       For examples with execution time > 5 sec, we replace `\dontrun{}` with 
       `\donttest{}` . 

---

  3. You write information messages to the console that cannot be easily 
    suppressed. It is more R like to generate objects that can be used 
    to extract the information a user is interested in, and then print() 
    that object. Instead of `print()`/`cat()` rather use `message()`/`warning()`
    or `if(verbose)cat(..)` (or maybe `stop()`) if you really have to write 
    text to the console. (except for print, summary, interactive functions)
    
      - **ACTION:** we avoid using `cat()`. In some cases, we kept `cat()`, however 
        always preceded by `if(verbose)`. Now, messages to the console can be easily 
        suppressed.

---

   4. Please make sure that you do not change the user's options, par or 
     working directory. If you really have to do so within functions, please 
     ensure with an *immediate* call of `on.exit()` that the settings are reset 
     when the function is exited. e.g.:
     
      ````markdown
        oldpar <- par(no.readonly = TRUE)       # code line i
        on.exit(par(oldpar))                    # code line i + 1
        par(mfrow=c(2,2))                       # somewhere after
      ````
      If you're not familiar with the function, please check `?on.exit`. This 
      function makes it possible to restore options before exiting a function 
      even if the function breaks. Therefore it needs to be called immediately 
      after the option change within a function.
     
      - **ACTION:** we used the structure `on.exit(par(oldpar))` in all ocurrences 
      of `oldpar <- par(...)`.

---

  5. Please always add all authors, contributors and copyright holders in the 
    Authors@R field with the appropriate roles. From CRAN policies you agreed to:
    "The ownership of copyright and intellectual property rights of all components 
    of the package must be clear and unambiguous (including from the authors 
    specification in the DESCRIPTION file). Where code is copied (or derived) 
    from the work of others (including from R itself), care must be taken that 
    any copyright/license statements are preserved and authorship is not 
    misrepresented. Preferably, an ‘Authors@R’ would be used with ‘ctb’ roles 
    for the authors of such code. Alternatively, the ‘Author’ field should list 
    these authors as contributors. Where copyrights are held by an entity other 
    than the package authors, this should preferably be indicated via ‘cph’ 
    roles in the ‘Authors@R’ field, or using a ‘Copyright’ field (if necessary 
    referring to an inst/COPYRIGHTS file)."e.g.: Karl W Broman, Robert Gentleman, 
    Ross Ihaka, The R Foundation, The R Core Team"
    Please explain in the submission comments what you did about this issue.
    
     - **ACTION:** 
       - we removed Karl Broman's C functions from our package. Since these 
         functions were used only to allocate memory, we replaced them to the 
         operator "new" in C++. This modification can be seen in lines 190-197 
         and 439-436 in the file `est_map_hmm_given_prior.cpp`.
       - we included Robert Gentleman, Ross Ihaka, The R Foundation, and 
         The R Core Team as copyright holders in the ‘Authors@R’ field in the 
         DESCRIPTION file.

---

  6. The LICENSE file is only needed if you have additional restrictions to 
      the license which you have not? In that case omit the file and its 
      reference in the DESCRIPTION file.
      
     - **ACTION:** we omitted the file and its reference in the DESCRIPTION file

---

Other changes since first submission.

  - Include new test files 
  - Fix minor bugs 
  - Update documentation 
  - Update DESCRIPTION file 

Thank you for reviewing our submission!

## Test environments
* local R installation (macOS 10.15.6), R 4.1.0
* local R installation (Ubuntu 18.04), R 3.6.3
* Ubuntu 16.04 (on travis-ci), R 4.0.2
* Windows Server x64 (on appveyor), R 4.0.2
* Win-builder (3.6.3, 4.0.2, and devel)

## R CMD check results 

0 errors | 0 warnings | 2 note

* This is a new submission 

* installed size is 12.3Mb
  sub-directories of 1Mb or more:
    * R:      2.6Mb
    * data:   9.0Mb
    
## Downstream dependencies

 There are currently no downstream dependencies for this package
