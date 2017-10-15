Setting up to use HiPERCAM
==========================

In order to use the HiPERCAM pipeline for normal reduction, you just need the
commands in your path. Try typing 'hplot<enter>' to see if you do. Apart from
this there are a few environment variables to know about.

Connecting to the server
========================

If using the pipeline in conjunction with the HiPERCAM instrument you need to
tell it where the server is which you do by setting an environment variable
called HIPERCAM_DEFAULT_URL. Typically I do this in my .cshrc file (C-shell,
use export with bash) with:

  setenv HIPERCAM_DEFAULT_URL 192.168.1.2:8007/

The value here is appropriate if you are connected to the HiPERCAM private
network. Note that there is no leading 'http://' unlike the equivalent in the
ULTRACAM pipeline because this address is used with both http and web sockets
and the appropriate signature is tacked on by whichever script is being used.
Also note the final '/' on the port number as in '8007/', which must be there.

Matplotlib-based plots
======================

It seems hard to get set the fontsize of matplotlib plots consistently across
devices. Therefore the following environment variable should allow the fonts
to be scaled in size to your preference:

  setenv HIPERCAM_MPL_FSCALE 0.5

which would make them all a factor of 2 smaller.
