# halfspace_WLC_simulations
Software for half-space (z>0) Worm-Like Chain simulations by Christian M. Kaiser and Dmitri Toptygin

This is a collection of one fortran function, two fortran subroutines, and two fortran programs, all of which were used to generate distance distributions and local concentrations for the article titled "Multivalent weak contacts shape chaperone-nascent protein interactions", authored by N. Rajasekaran, D. Toptygin, T.-W. Liao, V. J. Hilser, T. Ha, C. M. Kaiser, to be published in Proceedings of the National Academy of Sciences (2026). In that article half-space WLC simulations were intended to model the nascent polypeptide chain emerging from the ribosome exit tunnel. The ribosome surface was approximated as the Z=0 plane; the nascent chain was not allowed to enter the half-space Z≤0 that was occupied by the ribosome.

The numerical procedures are based on the theory described in the document ./doc/halfspace_WLC_simulations.pdf

The document ./doc/halfspace_WLC_simulations.pdf also describes how the programs can be used.
