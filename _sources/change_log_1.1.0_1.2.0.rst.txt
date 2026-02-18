.. changelog created on Tue 29 Mar 12:30:37 BST 2022

.. include:: globals.rst

|hiper| pipeline changes from v1.1.0 to v1.2.0
**************************************************

List of changes from git, newest first, with the commit keys linked to github:

  * `edb84f5dacdc14fb553c8ec184097d930b49a72e <https://github.com/HiPERCAM/hipercam/commit/edb84f5dacdc14fb553c8ec184097d930b49a72e>`_ problem_runs.py, to find problematic runs
  * `126041afd410e6037d51e946ee90cacc39968b41 <https://github.com/HiPERCAM/hipercam/commit/126041afd410e6037d51e946ee90cacc39968b41>`_ avalanche.py, finds avalanche mode runs for ultraspec
  * `2f547e3b41a54f1faeda61b68150a4380b74264c <https://github.com/HiPERCAM/hipercam/commit/2f547e3b41a54f1faeda61b68150a4380b74264c>`_ documented new commands
  * `d0d824e66c99b11766fe1a63d3c206f4e3704ab9 <https://github.com/HiPERCAM/hipercam/commit/d0d824e66c99b11766fe1a63d3c206f4e3704ab9>`_ Added in pbands to standard setup files
  * `42393419dc42c44cff526421a613082cbd1c1d10 <https://github.com/HiPERCAM/hipercam/commit/42393419dc42c44cff526421a613082cbd1c1d10>`_ exploss, improved documentation
  * `051dc23ad7078515580cd52095f1707d9e21fd16 <https://github.com/HiPERCAM/hipercam/commit/051dc23ad7078515580cd52095f1707d9e21fd16>`_ pbands, new quick-look plot routine
  * `4ac93bb1daf0ca6d9367f1c37296e3f00866b326 <https://github.com/HiPERCAM/hipercam/commit/4ac93bb1daf0ca6d9367f1c37296e3f00866b326>`_ exploss, added some numerical info.
  * `1fba90acd8176bb782f495ecffb80bfa7f957baf <https://github.com/HiPERCAM/hipercam/commit/1fba90acd8176bb782f495ecffb80bfa7f957baf>`_ fringe-runs fix
  * `da8b1eb524bb6c448e1afc618d725a6304870b55 <https://github.com/HiPERCAM/hipercam/commit/da8b1eb524bb6c448e1afc618d725a6304870b55>`_ logsearch, bug fix, ignore others option
  * `5e2eafbfd8cdc60a2da88159dfce51298472ee84 <https://github.com/HiPERCAM/hipercam/commit/5e2eafbfd8cdc60a2da88159dfce51298472ee84>`_ caksearch.py, several upgrades
  * `00bec41db64beb12c382867a6e96e194f3873492 <https://github.com/HiPERCAM/hipercam/commit/00bec41db64beb12c382867a6e96e194f3873492>`_ exploss -- to help judge if readout noise is significant
  * `b1a7052bfff335c580ae9b461b95518eead36db5 <https://github.com/HiPERCAM/hipercam/commit/b1a7052bfff335c580ae9b461b95518eead36db5>`_ changes instrument column hlogger and logsearch
  * `b5c84998eb8830561b42d186c84875378743b6cb <https://github.com/HiPERCAM/hipercam/commit/b5c84998eb8830561b42d186c84875378743b6cb>`_ hlogger, added instrument name column
  * `fb2c616ac31d1799b8c617a70e2e61e18c19b0bb <https://github.com/HiPERCAM/hipercam/commit/fb2c616ac31d1799b8c617a70e2e61e18c19b0bb>`_ calsearch, new routine to find matching calibration runs
  * `befc6ccab161ad98a1734ddf8978fd6408584956 <https://github.com/HiPERCAM/hipercam/commit/befc6ccab161ad98a1734ddf8978fd6408584956>`_ logsearch, added option to skip "Others"
  * `9840e2e56b669db86287f3e0352d7d1323e39a3c <https://github.com/HiPERCAM/hipercam/commit/9840e2e56b669db86287f3e0352d7d1323e39a3c>`_ hlogger, updated to so with ultracam phase II names
  * `d3458c4f98cb1e0ff76585e29ced262dc534919f <https://github.com/HiPERCAM/hipercam/commit/d3458c4f98cb1e0ff76585e29ced262dc534919f>`_ calsearch, new routine to find calibration files
  * `07680306ea5a8e012fc584330431c8364ebfa307 <https://github.com/HiPERCAM/hipercam/commit/07680306ea5a8e012fc584330431c8364ebfa307>`_ logsearch, fixed bug in sql query
  * `ccb47e37823de53c335b08abef64d981e63ad2c1 <https://github.com/HiPERCAM/hipercam/commit/ccb47e37823de53c335b08abef64d981e63ad2c1>`_ hlogger, bug fixes for airmass and retry option above all
  * `1c52b6c73022a315ace77947d92e1d449f50c6f3 <https://github.com/HiPERCAM/hipercam/commit/1c52b6c73022a315ace77947d92e1d449f50c6f3>`_ utils.target_lookup, j or J now
  * `b8ab0685591811552e85def408e9dc3312715d8a <https://github.com/HiPERCAM/hipercam/commit/b8ab0685591811552e85def408e9dc3312715d8a>`_ hlog.Tseries.phase, sort now an option
  * `ccbb18f726053aea93b36863e6cf6359dd707d05 <https://github.com/HiPERCAM/hipercam/commit/ccbb18f726053aea93b36863e6cf6359dd707d05>`_ fixed problem with logsearch not recognising new downloads
  * `fac17c4a0a2242d449c5fae96343e356dc59067a <https://github.com/HiPERCAM/hipercam/commit/fac17c4a0a2242d449c5fae96343e356dc59067a>`_ hlogger, fixed check for number of entries in posdata
  * `981462dbdfc8afd3fe91c4301f21198f851fcfd1 <https://github.com/HiPERCAM/hipercam/commit/981462dbdfc8afd3fe91c4301f21198f851fcfd1>`_ logsearch, print statement added
  * `7bb17b9715c0bfedcd729e9cd90b6fd0ab92e2e9 <https://github.com/HiPERCAM/hipercam/commit/7bb17b9715c0bfedcd729e9cd90b6fd0ab92e2e9>`_ extinction.py, fixed problem with negative deviants
  * `46ce40f34ef61cf9eeaf55f7c53dbaa1f32c8fb9 <https://github.com/HiPERCAM/hipercam/commit/46ce40f34ef61cf9eeaf55f7c53dbaa1f32c8fb9>`_ logsearch, strange bug with password prompt
  * `3864112de52ea2c86ebe3b8055d0fb1355f9ca3b <https://github.com/HiPERCAM/hipercam/commit/3864112de52ea2c86ebe3b8055d0fb1355f9ca3b>`_ Reworked the useful/README a bit.
  * `6fc33f2f1c4020e8fbc0080a324f4c0e150e426d <https://github.com/HiPERCAM/hipercam/commit/6fc33f2f1c4020e8fbc0080a324f4c0e150e426d>`_ reorganised useful scripts directoy
  * `2d81bf8f70b28e8e0364fc56972c61ee24c11fe5 <https://github.com/HiPERCAM/hipercam/commit/2d81bf8f70b28e8e0364fc56972c61ee24c11fe5>`_ extinction.py, script for measuring extinction
  * `20d7466ddf03ab7e8480e47c1b6799d5b1c1ceeb <https://github.com/HiPERCAM/hipercam/commit/20d7466ddf03ab7e8480e47c1b6799d5b1c1ceeb>`_ hlog.py, removed irritating units from airmasses
  * `24026129e9c9b6412898d5637a7e85d2a15f7bde <https://github.com/HiPERCAM/hipercam/commit/24026129e9c9b6412898d5637a7e85d2a15f7bde>`_ hlog.py, secz to seczs bug fix
  * `baaa0e21db8cb7c0f31d33dae795d59c740fb992 <https://github.com/HiPERCAM/hipercam/commit/baaa0e21db8cb7c0f31d33dae795d59c740fb992>`_ hlog.py, added airmass computation
  * `f6b29811f5c3064c7107b9f744e1d9596cfac205 <https://github.com/HiPERCAM/hipercam/commit/f6b29811f5c3064c7107b9f744e1d9596cfac205>`_ hlogger, changed column name of Delta secz
  * `dee963dec8018902dc4ba24da89ded3edcc48f9a <https://github.com/HiPERCAM/hipercam/commit/dee963dec8018902dc4ba24da89ded3edcc48f9a>`_ hlogger.py, bugs fixed
  * `2ebc585b88ddbe990a91d870b4e6d7b5244c92a5 <https://github.com/HiPERCAM/hipercam/commit/2ebc585b88ddbe990a91d870b4e6d7b5244c92a5>`_ hlogger, added airmass info
  * `0e8e568103dd8c43a4e35c4709fa51dff8de4043 <https://github.com/HiPERCAM/hipercam/commit/0e8e568103dd8c43a4e35c4709fa51dff8de4043>`_ Added keyring to list of prerequisites.
  * `6beafa0d75f652fc0136d15e8fc6bbd5e214b331 <https://github.com/HiPERCAM/hipercam/commit/6beafa0d75f652fc0136d15e8fc6bbd5e214b331>`_ logsearch, automated the database download
  * `72c55e2d602cf73928e8cccd56247ea90a14fd21 <https://github.com/HiPERCAM/hipercam/commit/72c55e2d602cf73928e8cccd56247ea90a14fd21>`_ hpackage, generalised behaviour
  * `f8237920a0392728c4de40fb26c9c4cd2ff577ab <https://github.com/HiPERCAM/hipercam/commit/f8237920a0392728c4de40fb26c9c4cd2ff577ab>`_ hpackage, added option to pick up alternative reductions
  * `e6d3fa8fe552bb20dd1531579a7f5b1e1b8028b6 <https://github.com/HiPERCAM/hipercam/commit/e6d3fa8fe552bb20dd1531579a7f5b1e1b8028b6>`_ Tseries, add clear_bitmask
  * `e52af38020aa281d0fb3d983554c7030fbb4d169 <https://github.com/HiPERCAM/hipercam/commit/e52af38020aa281d0fb3d983554c7030fbb4d169>`_ Tseries.write, added option to skip masked data
  * `7876a627f87a23d8231dc45bb769133b831198cb <https://github.com/HiPERCAM/hipercam/commit/7876a627f87a23d8231dc45bb769133b831198cb>`_ Merge branch 'master' into trm-dev
  * `5bf9bee400e188d84ef9d7c5ec8bfdd8d4dc0acc <https://github.com/HiPERCAM/hipercam/commit/5bf9bee400e188d84ef9d7c5ec8bfdd8d4dc0acc>`_ hlogger, change to html column format
  * `a240e7536887fe671151ce14d350f7279c2ab8b0 <https://github.com/HiPERCAM/hipercam/commit/a240e7536887fe671151ce14d350f7279c2ab8b0>`_ logging doc change
  * `0912e281cc0c46fe0da885306f1ab32f1278e535 <https://github.com/HiPERCAM/hipercam/commit/0912e281cc0c46fe0da885306f1ab32f1278e535>`_ Merge pull request #92 from HiPERCAM/tseries_read_fix
  * `ee33c77d6dab7caa3626eff37399b941a2ba417a <https://github.com/HiPERCAM/hipercam/commit/ee33c77d6dab7caa3626eff37399b941a2ba417a>`_ fix typo in Tseries.read
  * `34f6d3519da2b42fd81aecf5a47362888c8d199a <https://github.com/HiPERCAM/hipercam/commit/34f6d3519da2b42fd81aecf5a47362888c8d199a>`_ joinup now joins up regardless of out of sync windows
  * `dab5608e407afb62fe0631a5dfece897d64b0447 <https://github.com/HiPERCAM/hipercam/commit/dab5608e407afb62fe0631a5dfece897d64b0447>`_ modified joinup to deal with asynchronous windows
  * `93f5bd57cf3c338aee988d9f1e9864f7c475f1a0 <https://github.com/HiPERCAM/hipercam/commit/93f5bd57cf3c338aee988d9f1e9864f7c475f1a0>`_ changed digest to harchive
  * `8d66fbf5ba18329ce23e44731d3baede52373bae <https://github.com/HiPERCAM/hipercam/commit/8d66fbf5ba18329ce23e44731d3baede52373bae>`_ put a try/except into hpackage to avoid joinup error
  * `2028dfa3e9a59fd787f8c2953a7a8ec8aef81357 <https://github.com/HiPERCAM/hipercam/commit/2028dfa3e9a59fd787f8c2953a7a8ec8aef81357>`_ generalised the logsearch a bit
  * `39a0deb6ae967e3afb4d03cb78179a680928fda3 <https://github.com/HiPERCAM/hipercam/commit/39a0deb6ae967e3afb4d03cb78179a680928fda3>`_ fixed problem with latest upgrade to logsearch
  * `991f3628cea8fbc5226dd8b80c3c0dd5fabbee03 <https://github.com/HiPERCAM/hipercam/commit/991f3628cea8fbc5226dd8b80c3c0dd5fabbee03>`_ added -m option to simply report failed names.
  * `a5a88de7515a891d8f0c3aa57a696e1701bd5362 <https://github.com/HiPERCAM/hipercam/commit/a5a88de7515a891d8f0c3aa57a696e1701bd5362>`_ added mjd to tdb correction to Tseries
  * `c078b3335e8437fe761d1f2332809e0804ed652b <https://github.com/HiPERCAM/hipercam/commit/c078b3335e8437fe761d1f2332809e0804ed652b>`_ added options to save results to a csv file to logsearch
  * `b1924bc0ac50a44fbb99d3f7b49bf5a0cea723e3 <https://github.com/HiPERCAM/hipercam/commit/b1924bc0ac50a44fbb99d3f7b49bf5a0cea723e3>`_ reduced minimum fwhm in genred
  * `e881e9bd98f65f9c1dcfacb3ae7fed85a104179b <https://github.com/HiPERCAM/hipercam/commit/e881e9bd98f65f9c1dcfacb3ae7fed85a104179b>`_ added plotting of the max counts in plog
  * `db1706632c88449e541a8f077f3454149d68802c <https://github.com/HiPERCAM/hipercam/commit/db1706632c88449e541a8f077f3454149d68802c>`_ cline now sets the default for the ignore value
  * `781022e2430bd0cda4a286b134b340afe8b04030 <https://github.com/HiPERCAM/hipercam/commit/781022e2430bd0cda4a286b134b340afe8b04030>`_ further adjustment to docs
  * `bc713d7b9f6d1e8ca3ee10e5a9a49e504a02c032 <https://github.com/HiPERCAM/hipercam/commit/bc713d7b9f6d1e8ca3ee10e5a9a49e504a02c032>`_ utils.py strips `STD_` from start of target names
  * `f6f2d350386d489cbd8c88d6cd04c589cd9bd97a <https://github.com/HiPERCAM/hipercam/commit/f6f2d350386d489cbd8c88d6cd04c589cd9bd97a>`_ doc updates to describe phase II process in the main
  * `bb76ebc1d62cdc41b28d0dcb168acec252e95255 <https://github.com/HiPERCAM/hipercam/commit/bb76ebc1d62cdc41b28d0dcb168acec252e95255>`_ further tweak following change to defect file
  * `721e56070fb4a03bf5dfb511a349828c5f6f4806 <https://github.com/HiPERCAM/hipercam/commit/721e56070fb4a03bf5dfb511a349828c5f6f4806>`_ new hipercam defects file following Vik's cleaning of CCD 2
  * `0c0c89d67383ee5a7d88bcf7bcd5328f1b03bf33 <https://github.com/HiPERCAM/hipercam/commit/0c0c89d67383ee5a7d88bcf7bcd5328f1b03bf33>`_ changelog for a single tag was not reporting its name
  * `5e6ab55b1121d449d7acbc0aed3b6cff9506d341 <https://github.com/HiPERCAM/hipercam/commit/5e6ab55b1121d449d7acbc0aed3b6cff9506d341>`_ many fixes for genred's template file option
  * `bab87543e02c3997e2dafb9cd7139ac1a4435c8d <https://github.com/HiPERCAM/hipercam/commit/bab87543e02c3997e2dafb9cd7139ac1a4435c8d>`_ expanded range of boolean inputs for reduce files
  * `a083e5f6da7a99e63615cf1a96f0a2d56b59ed89 <https://github.com/HiPERCAM/hipercam/commit/a083e5f6da7a99e63615cf1a96f0a2d56b59ed89>`_ fixed invalid keyword in template part of genred
  * `c8fc24b243ff9540f99d12bcddbb943143d3380b <https://github.com/HiPERCAM/hipercam/commit/c8fc24b243ff9540f99d12bcddbb943143d3380b>`_ got working version of logsearch
  * `80c9bfe8b4a0109aeee454f3e83de9fc0750349d <https://github.com/HiPERCAM/hipercam/commit/80c9bfe8b4a0109aeee454f3e83de9fc0750349d>`_ generalised the ignore option in cline
  * `437bbce5f6bd1ca901a075a6cf8058b679ba858e <https://github.com/HiPERCAM/hipercam/commit/437bbce5f6bd1ca901a075a6cf8058b679ba858e>`_ array index fix, redplt
  * `98e379ee0a3683b36756415c5d6bf6828556bd14 <https://github.com/HiPERCAM/hipercam/commit/98e379ee0a3683b36756415c5d6bf6828556bd14>`_ adjusted recent changes to redplt
  * `7b6b62da593acb9649baf25339603dd305df456f <https://github.com/HiPERCAM/hipercam/commit/7b6b62da593acb9649baf25339603dd305df456f>`_ fixed bugs to do with a single CCD
  * `eac66b1903acdb2a29d913ae8190a6e06f64f3d9 <https://github.com/HiPERCAM/hipercam/commit/eac66b1903acdb2a29d913ae8190a6e06f64f3d9>`_ horrible numpy/python issue in hlog.py
  * `f720d87afb5268bb91318065772151e80b2cb09a <https://github.com/HiPERCAM/hipercam/commit/f720d87afb5268bb91318065772151e80b2cb09a>`_ typo fixes to genred from James Wild
  * `1ffe1c1ee9528a9193cd6422b090fbe4a28c48de <https://github.com/HiPERCAM/hipercam/commit/1ffe1c1ee9528a9193cd6422b090fbe4a28c48de>`_ fixed some backwards compatibility issues
  * `8178c05c5752848e5ff9157cfcf03b82349a7bd9 <https://github.com/HiPERCAM/hipercam/commit/8178c05c5752848e5ff9157cfcf03b82349a7bd9>`_ placed limit on urllib3 version for security
  * `629764b779f1981b0857a8ff9f31b27f5f15db24 <https://github.com/HiPERCAM/hipercam/commit/629764b779f1981b0857a8ff9f31b27f5f15db24>`_ updated logsearch to work with database files
  * `eb703a1d056a3eec2a2a2a584b8f3d3d5b796300 <https://github.com/HiPERCAM/hipercam/commit/eb703a1d056a3eec2a2a2a584b8f3d3d5b796300>`_ changed how redplt makes figures
  * `6644a33292b9e02d29c6594894b10e070a5abf7a <https://github.com/HiPERCAM/hipercam/commit/6644a33292b9e02d29c6594894b10e070a5abf7a>`_ added "import warnings" to genred
  * `469888447d1d069a2ba67e7d506b7b9fe1e788d6 <https://github.com/HiPERCAM/hipercam/commit/469888447d1d069a2ba67e7d506b7b9fe1e788d6>`_ revised the way in which redplt makes plots
  * `5c00ae56b84a214f635cb9a111cd3de97c2e1f7f <https://github.com/HiPERCAM/hipercam/commit/5c00ae56b84a214f635cb9a111cd3de97c2e1f7f>`_ attempt to fix genred problem using template file
  * `4a7433a5abcdd4a878bf1bdd6fa0bf0b077466e1 <https://github.com/HiPERCAM/hipercam/commit/4a7433a5abcdd4a878bf1bdd6fa0bf0b077466e1>`_ fixed hlogger bug in accessing phase II data.
  * `d66e07d7ab00636815002b7f767f8999999c68c2 <https://github.com/HiPERCAM/hipercam/commit/d66e07d7ab00636815002b7f767f8999999c68c2>`_ added some info about plog in docs
  * `53a8ae0ec07199ff776ac71cdcf4c752c28f2f84 <https://github.com/HiPERCAM/hipercam/commit/53a8ae0ec07199ff776ac71cdcf4c752c28f2f84>`_ series of stupid bugs following withdrawal of Hlog.read
  * `51965e36cc507f860645027ee4b4608d5b5764d1 <https://github.com/HiPERCAM/hipercam/commit/51965e36cc507f860645027ee4b4608d5b5764d1>`_ docs update for ncal and setdefect
  * `b16a14257b2c33b0cc6ed3f7e92458abb0a83e9f <https://github.com/HiPERCAM/hipercam/commit/b16a14257b2c33b0cc6ed3f7e92458abb0a83e9f>`_ updated the way changelog works
  * `a9fe69118923283d63e5c216061ef4b2d4b9f795 <https://github.com/HiPERCAM/hipercam/commit/a9fe69118923283d63e5c216061ef4b2d4b9f795>`_ doc changes, now track api!
  * `fa4c9c34c24cf8e3c8521e927f3b11fb96139a33 <https://github.com/HiPERCAM/hipercam/commit/fa4c9c34c24cf8e3c8521e927f3b11fb96139a33>`_ documentation of ncal
  * `324a5737fc9b7771266d3d78af1dccf08d03b6a6 <https://github.com/HiPERCAM/hipercam/commit/324a5737fc9b7771266d3d78af1dccf08d03b6a6>`_ new routine ncal for noise calibration
  * `1ec0c2d1438bda92638746d5e8f68c3851e3cdc9 <https://github.com/HiPERCAM/hipercam/commit/1ec0c2d1438bda92638746d5e8f68c3851e3cdc9>`_ and more doc updates
  * `1ddd80fb4c91a932ac9e3e94525b40bc41549760 <https://github.com/HiPERCAM/hipercam/commit/1ddd80fb4c91a932ac9e3e94525b40bc41549760>`_ doc updates to reduction.rst
  * `67454220cf15f77190012960707b06f5c3bb9562 <https://github.com/HiPERCAM/hipercam/commit/67454220cf15f77190012960707b06f5c3bb9562>`_ doc updates
  * `a3cc8d4bd9d0c5cf53ae8d822171c203d00251f8 <https://github.com/HiPERCAM/hipercam/commit/a3cc8d4bd9d0c5cf53ae8d822171c203d00251f8>`_ a bit more documentation in setdefect
  * `699fe717ac727331d8606160528a578d44402168 <https://github.com/HiPERCAM/hipercam/commit/699fe717ac727331d8606160528a578d44402168>`_ hot pixels defects implemented properly
  * `237718b21e451014349ea452d4c9914528b6dcfc <https://github.com/HiPERCAM/hipercam/commit/237718b21e451014349ea452d4c9914528b6dcfc>`_ tweaked hist a bit
  * `470875a226a72877171a9394d5b44eaafd9b61f1 <https://github.com/HiPERCAM/hipercam/commit/470875a226a72877171a9394d5b44eaafd9b61f1>`_ fixed two bugs in genred for inst=other
  * `b9c50394cab8f865c6931fe92285c539158439d1 <https://github.com/HiPERCAM/hipercam/commit/b9c50394cab8f865c6931fe92285c539158439d1>`_ the new defect file
  * `359d6687f0a06837b3cf2e5da6a74ec8f22e0104 <https://github.com/HiPERCAM/hipercam/commit/359d6687f0a06837b3cf2e5da6a74ec8f22e0104>`_ updated hipercam defect file
  * `48be2fcd6b0f6c3ae12dae20ec5a93614f0b9b28 <https://github.com/HiPERCAM/hipercam/commit/48be2fcd6b0f6c3ae12dae20ec5a93614f0b9b28>`_ removed enforcement of plot limits from multiple routines
  * `4ddc1ab85e1a591a5322d1d143063bdd2f556298 <https://github.com/HiPERCAM/hipercam/commit/4ddc1ab85e1a591a5322d1d143063bdd2f556298>`_ micro re-format, jtrawl
  * `7be56e25dea4e561099a27dd0ed7800d3034b598 <https://github.com/HiPERCAM/hipercam/commit/7be56e25dea4e561099a27dd0ed7800d3034b598>`_ now looks for phase II position  data as a last resort
  * `79eaac7b639d67bfa932483cce3cede9b7ba7ab2 <https://github.com/HiPERCAM/hipercam/commit/79eaac7b639d67bfa932483cce3cede9b7ba7ab2>`_ removed full recompute option from hmeta
  * `bd93b605c348c33424a1d3f736ec6bf51f1be4c5 <https://github.com/HiPERCAM/hipercam/commit/bd93b605c348c33424a1d3f736ec6bf51f1be4c5>`_ fixed bugs in hmeta for hipercam
  * `cb7660e7b6799cde0b5896f71372243de563557c <https://github.com/HiPERCAM/hipercam/commit/cb7660e7b6799cde0b5896f71372243de563557c>`_ usual time stamp update
  * `2822c977bc20d8c492902c2e017c6692e87a0f86 <https://github.com/HiPERCAM/hipercam/commit/2822c977bc20d8c492902c2e017c6692e87a0f86>`_ Merge branch 'trm-dev'
  * `aec3ec2df7e5033c76ba91822bee17db6206333d <https://github.com/HiPERCAM/hipercam/commit/aec3ec2df7e5033c76ba91822bee17db6206333d>`_ another micro-tweak on genred
  * `4a9767db2eee9c9970fea689eb7f245cd4330bde <https://github.com/HiPERCAM/hipercam/commit/4a9767db2eee9c9970fea689eb7f245cd4330bde>`_ finer-grained control in genred
  * `8982ffb76fbea506132e0d4ec5b595a5f88c79c3 <https://github.com/HiPERCAM/hipercam/commit/8982ffb76fbea506132e0d4ec5b595a5f88c79c3>`_ small tweaks to default fit params in genred
  * `dae90f1971044da5eef3ef479d4b9d9e01c74688 <https://github.com/HiPERCAM/hipercam/commit/dae90f1971044da5eef3ef479d4b9d9e01c74688>`_ latest change log
  * `18657e2899421bbc1d3c79788409f7b183cc2cf2 <https://github.com/HiPERCAM/hipercam/commit/18657e2899421bbc1d3c79788409f7b183cc2cf2>`_ doc updates following latest release
