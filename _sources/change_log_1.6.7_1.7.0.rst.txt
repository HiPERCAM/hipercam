.. changelog created on Fri  9 Oct 2026 12:43:50 BST

.. include:: globals.rst

|hiper| pipeline changes from v1.6.7 to v1.7.0
**************************************************

List of changes from git, newest first, with the commit keys linked to github:

  * `5d85591fc697f2743329e112aa7d756539f04cc1 <https://github.com/HiPERCAM/hipercam/commit/5d85591fc697f2743329e112aa7d756539f04cc1>`_ Merge pull request #162 from HiPERCAM/psfphot_v2
  * `b82d3c9a16f953dacc78ce5b227d7df511727a88 <https://github.com/HiPERCAM/hipercam/commit/b82d3c9a16f953dacc78ce5b227d7df511727a88>`_ scale bounds to binned pixels
  * `787c8d1d00f72123fd044f650bca1c604cc4b112 <https://github.com/HiPERCAM/hipercam/commit/787c8d1d00f72123fd044f650bca1c604cc4b112>`_ better bounds for beta parameter
  * `0657feddf9d03650dd672556ba15eafd7141b01a <https://github.com/HiPERCAM/hipercam/commit/0657feddf9d03650dd672556ba15eafd7141b01a>`_ initialise plot lims if any image plotted
  * `b1135cf3788e21621078bd1b18971751dbee39d1 <https://github.com/HiPERCAM/hipercam/commit/b1135cf3788e21621078bd1b18971751dbee39d1>`_ supply default if psfdevice missing from reduce file
  * `fddbe56dcd17a720692a6a828e6ca6374ea1197f <https://github.com/HiPERCAM/hipercam/commit/fddbe56dcd17a720692a6a828e6ca6374ea1197f>`_ fix bug in no reference aperture case
  * `26d986ab794f087be3d67f49b89e64c68273be92 <https://github.com/HiPERCAM/hipercam/commit/26d986ab794f087be3d67f49b89e64c68273be92>`_ remove store from arguments being passed to update_plots
  * `2f7cb69602b9ee2eb8faaa70227fc7d7be198e3b <https://github.com/HiPERCAM/hipercam/commit/2f7cb69602b9ee2eb8faaa70227fc7d7be198e3b>`_ better output from psfaper
  * `d527bc66dba0e87aa293ff55d5d796cb14a484fc <https://github.com/HiPERCAM/hipercam/commit/d527bc66dba0e87aa293ff55d5d796cb14a484fc>`_ clarifying docstring for create_psf_model
  * `cac230c42bf166d025af15e9dd4f248ed1f6d9fc <https://github.com/HiPERCAM/hipercam/commit/cac230c42bf166d025af15e9dd4f248ed1f6d9fc>`_ add shortcut for psfaper
  * `22110bdc10f74df4deb0d657c08439bdf6c70a4b <https://github.com/HiPERCAM/hipercam/commit/22110bdc10f74df4deb0d657c08439bdf6c70a4b>`_ Fix bug where we were fitting all apertures instead of references
  * `7156d1fab8f0a523bc07547fc5c38491bc388774 <https://github.com/HiPERCAM/hipercam/commit/7156d1fab8f0a523bc07547fc5c38491bc388774>`_ Update flag handling in PSF photometry
  * `f9a1eb25bcbdd25fb689fd8a076ddac2afdf2606 <https://github.com/HiPERCAM/hipercam/commit/f9a1eb25bcbdd25fb689fd8a076ddac2afdf2606>`_ better instructions for PSF apertures
  * `fa382cdebce1315cc330fb490b261d8f3bd5d101 <https://github.com/HiPERCAM/hipercam/commit/fa382cdebce1315cc330fb490b261d8f3bd5d101>`_ Merge pull request #178 from HiPERCAM/speedup
  * `8723214ca343434cbb49941dcea164154862458a <https://github.com/HiPERCAM/hipercam/commit/8723214ca343434cbb49941dcea164154862458a>`_ Add VSCode to ignore
  * `eeb2e200a774eb6b8b8ebb497c284d745c0b4994 <https://github.com/HiPERCAM/hipercam/commit/eeb2e200a774eb6b8b8ebb497c284d745c0b4994>`_ Vectorise to remove numba dependency
  * `86a82d2f57ae9193f8786732f5827d25f690c76b <https://github.com/HiPERCAM/hipercam/commit/86a82d2f57ae9193f8786732f5827d25f690c76b>`_ Remove non-C++ functions
  * `190509ae783df416b28e9f373b0b9a9f4c52be0b <https://github.com/HiPERCAM/hipercam/commit/190509ae783df416b28e9f373b0b9a9f4c52be0b>`_ Merge pull request #179 from HiPERCAM/setfringe_fix
  * `f07262d3c74ecea66ea184949d9c5a4585d9c705 <https://github.com/HiPERCAM/hipercam/commit/f07262d3c74ecea66ea184949d9c5a4585d9c705>`_ remove deprecation warning which hasn't been a thing since 2015
  * `b04a92a65fb0304a66fd98482f5b93090c8186e9 <https://github.com/HiPERCAM/hipercam/commit/b04a92a65fb0304a66fd98482f5b93090c8186e9>`_ proper handling of failed extraction
  * `556851b31edb7e73ff733982bedf5c3bebbd8559 <https://github.com/HiPERCAM/hipercam/commit/556851b31edb7e73ff733982bedf5c3bebbd8559>`_ bugfix in aperture display and suppression of FITS warnings
  * `5a96704dfc5d72671a209ead5cecb67e3ff6ea0f <https://github.com/HiPERCAM/hipercam/commit/5a96704dfc5d72671a209ead5cecb67e3ff6ea0f>`_ Merge branch 'master' into psfphot_v2
  * `b423f75b8799847664b30e727a4756f68eb0cdbd <https://github.com/HiPERCAM/hipercam/commit/b423f75b8799847664b30e727a4756f68eb0cdbd>`_ fix pbands prompts
  * `b7bf49beddd4e226f7822938a9b5819b2b9e6c4c <https://github.com/HiPERCAM/hipercam/commit/b7bf49beddd4e226f7822938a9b5819b2b9e6c4c>`_ Final C++ cleaning
  * `245c52be1cc5aaf00c795cb16c7adcab51210487 <https://github.com/HiPERCAM/hipercam/commit/245c52be1cc5aaf00c795cb16c7adcab51210487>`_ Add tol options
  * `d36db17790bf5b1bad34c2719eb2d30cee8b83c3 <https://github.com/HiPERCAM/hipercam/commit/d36db17790bf5b1bad34c2719eb2d30cee8b83c3>`_ Merge pull request #177 from HiPERCAM/hplot_bug
  * `0728579e1d92681e89d859ae0f6d68cb0c40364c <https://github.com/HiPERCAM/hipercam/commit/0728579e1d92681e89d859ae0f6d68cb0c40364c>`_ bypass callback if quit key pressed
  * `154f2ebb346d618f4cbcf2951e5118ecf19a2be1 <https://github.com/HiPERCAM/hipercam/commit/154f2ebb346d618f4cbcf2951e5118ecf19a2be1>`_ Add comparison test scripts for C++
  * `9a9cd5dbef80506fc3f98f4be50311ed578b3a98 <https://github.com/HiPERCAM/hipercam/commit/9a9cd5dbef80506fc3f98f4be50311ed578b3a98>`_ Move fully from Cython to Pybind11
  * `6a7a374ea0077a9f5f9417ad1af4dbc32df35dbb <https://github.com/HiPERCAM/hipercam/commit/6a7a374ea0077a9f5f9417ad1af4dbc32df35dbb>`_ More optimisations
  * `6e99f02397e41bcf8f28349e67758bafa50d7972 <https://github.com/HiPERCAM/hipercam/commit/6e99f02397e41bcf8f28349e67758bafa50d7972>`_ Add GIL release and OpenMP parallelisation
  * `215ce1d77ab95cf00e24dceae9b5edbdb3e7135a <https://github.com/HiPERCAM/hipercam/commit/215ce1d77ab95cf00e24dceae9b5edbdb3e7135a>`_ Optimise residual/jacobian in C++
  * `723335c1d061afbe7b52a9afea898ac7c6b28982 <https://github.com/HiPERCAM/hipercam/commit/723335c1d061afbe7b52a9afea898ac7c6b28982>`_ Some optimisations
  * `5c6bcb0d84194772962d925bd4dfa21de38c190b <https://github.com/HiPERCAM/hipercam/commit/5c6bcb0d84194772962d925bd4dfa21de38c190b>`_ Add C++ versions of fitting functions
  * `c5913dd01cd0897b6ee3255d9af63e53fbe5d2f2 <https://github.com/HiPERCAM/hipercam/commit/c5913dd01cd0897b6ee3255d9af63e53fbe5d2f2>`_ Optimisations in fit classes
  * `3200749fa6c714db982e1f98dec5a0ff97b6cff4 <https://github.com/HiPERCAM/hipercam/commit/3200749fa6c714db982e1f98dec5a0ff97b6cff4>`_ Bug fixing
  * `d36e86e2cc2b0b0c63f42b224dd47ba1e9c23063 <https://github.com/HiPERCAM/hipercam/commit/d36e86e2cc2b0b0c63f42b224dd47ba1e9c23063>`_ Merge pull request #174 from HiPERCAM/reduce_flags
  * `4657603b6f7b69a2bc9091eda5758736fbd4037c <https://github.com/HiPERCAM/hipercam/commit/4657603b6f7b69a2bc9091eda5758736fbd4037c>`_ plot non-linear points in reduce
  * `30055e84626e5179c6455f15475cc8b3a62e5aeb <https://github.com/HiPERCAM/hipercam/commit/30055e84626e5179c6455f15475cc8b3a62e5aeb>`_ Merge pull request #168 from HiPERCAM/plugins
  * `5027f47a3e0e56c131cb45c57c8463b6f07c8b6e <https://github.com/HiPERCAM/hipercam/commit/5027f47a3e0e56c131cb45c57c8463b6f07c8b6e>`_ turn on autoapi
  * `da4395f7edb742ca73af3420123564650bdcced7 <https://github.com/HiPERCAM/hipercam/commit/da4395f7edb742ca73af3420123564650bdcced7>`_ explain get_ccd_pars arguments
  * `f49188bb2b0519788033b762605d9f74a5aab06f <https://github.com/HiPERCAM/hipercam/commit/f49188bb2b0519788033b762605d9f74a5aab06f>`_ adopt mjd's suggestions
  * `7a5694343004fc7563dfbd4208e3aa0505b619ce <https://github.com/HiPERCAM/hipercam/commit/7a5694343004fc7563dfbd4208e3aa0505b619ce>`_ Update docs/plugins.rst
  * `473928c830476d15b3caad927bec28dcf71fce0c <https://github.com/HiPERCAM/hipercam/commit/473928c830476d15b3caad927bec28dcf71fce0c>`_ remove run from class
  * `4888f38d507eb7c52e6bea4f7e1ba8c1494699dc <https://github.com/HiPERCAM/hipercam/commit/4888f38d507eb7c52e6bea4f7e1ba8c1494699dc>`_ Merge pull request #171 from HiPERCAM/plog_fix
  * `f36b1e97c7b582709aa6f68fb73bf76716ad2a05 <https://github.com/HiPERCAM/hipercam/commit/f36b1e97c7b582709aa6f68fb73bf76716ad2a05>`_ rename param
  * `9005f9a630d1a388efc58fee81437c5c0c8e987d <https://github.com/HiPERCAM/hipercam/commit/9005f9a630d1a388efc58fee81437c5c0c8e987d>`_ pbands too
  * `a31b96c69bcb28f7e34fe6d1fdc7981ac1698fa7 <https://github.com/HiPERCAM/hipercam/commit/a31b96c69bcb28f7e34fe6d1fdc7981ac1698fa7>`_ add option to plot bad times in plog
  * `bcf80a97a1884f4f896e628a7bc4495871f9a41f <https://github.com/HiPERCAM/hipercam/commit/bcf80a97a1884f4f896e628a7bc4495871f9a41f>`_ change arguments to all relevant scripts
  * `6a0b9a9d9053b831517a0aa4f8935e8878c3fc13 <https://github.com/HiPERCAM/hipercam/commit/6a0b9a9d9053b831517a0aa4f8935e8878c3fc13>`_ first pass at plugins, new args for nrtplot
  * `2f6cbec7da1f9411c91d7880e723fa953932bcd8 <https://github.com/HiPERCAM/hipercam/commit/2f6cbec7da1f9411c91d7880e723fa953932bcd8>`_ fix bug found by @mkenne15 - psfplot must be defined
  * `363f3aa27250d23cca86146d4d27b4ba08b33484 <https://github.com/HiPERCAM/hipercam/commit/363f3aa27250d23cca86146d4d27b4ba08b33484>`_ dome docs
  * `fa789fc2d086e305c31f663209a7d15e6947af05 <https://github.com/HiPERCAM/hipercam/commit/fa789fc2d086e305c31f663209a7d15e6947af05>`_ better docs
  * `4ec7a165f4e3c541c01abe6182ab07300e5ae7a8 <https://github.com/HiPERCAM/hipercam/commit/4ec7a165f4e3c541c01abe6182ab07300e5ae7a8>`_ slight bugfix with pre-existing apertures
  * `d1171f890a7409f823ac7d72043aafce044fa54c <https://github.com/HiPERCAM/hipercam/commit/d1171f890a7409f823ac7d72043aafce044fa54c>`_ deploy psfaper script
  * `ec1782f84f8df97c36cc4063ffcd2f20cd8ea9bb <https://github.com/HiPERCAM/hipercam/commit/ec1782f84f8df97c36cc4063ffcd2f20cd8ea9bb>`_ two pass approach. fix PSF from reference stars
  * `80c981a20258cc130706e0716b72df5b450386c9 <https://github.com/HiPERCAM/hipercam/commit/80c981a20258cc130706e0716b72df5b450386c9>`_ gaussian PSF was not picking up settings from reference aperture fitting
  * `e1657eeb73a1b063d467fafa9ea3251b8b9b45bd <https://github.com/HiPERCAM/hipercam/commit/e1657eeb73a1b063d467fafa9ea3251b8b9b45bd>`_ missing import
  * `1218570fb8fdd04a33afb3d8a8dcf684dcfb51ed <https://github.com/HiPERCAM/hipercam/commit/1218570fb8fdd04a33afb3d8a8dcf684dcfb51ed>`_ use flags in PSF photom
  * `ce05b19224a18d42e46e666f63bc25899541739e <https://github.com/HiPERCAM/hipercam/commit/ce05b19224a18d42e46e666f63bc25899541739e>`_ add photutils as optional dependency
  * `2001abed97f9b614d940cbb9b981230be95bd612 <https://github.com/HiPERCAM/hipercam/commit/2001abed97f9b614d940cbb9b981230be95bd612>`_ improvements and bugfixes to genred for PSF photometry
  * `a9945c337708ef3c43389e46052497def2d803c5 <https://github.com/HiPERCAM/hipercam/commit/a9945c337708ef3c43389e46052497def2d803c5>`_ add plotting of residual PSF images
  * `e9859e8dc6125ccf9aa0a88c9da94305a9bf054f <https://github.com/HiPERCAM/hipercam/commit/e9859e8dc6125ccf9aa0a88c9da94305a9bf054f>`_ add gradient to moffat model
  * `919bfef7810cac314bf48f733f3a79a9bf71fbfc <https://github.com/HiPERCAM/hipercam/commit/919bfef7810cac314bf48f733f3a79a9bf71fbfc>`_ first pass at using new photutils models and routines