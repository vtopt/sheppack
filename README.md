# ACM TOMS Algorithm 905: SHEPPACK: Modified Shepard Algorithm for Interpolation of Scattered Multivariate Data

SHEPPACK is a Fortran 95 package containing five versions of the modified Shepard algorithm: quadratic (Fortran 95 translations of Algorithms 660, 661, and 798), cubic (Fortran 95 translation of Algorithm 791), and linear variations of the original Shepard algorithm. An option to the linear Shepard code is a statistically robust fit, intended to be used when the data is known to contain outliers. SHEPPACK also includes a hybrid robust piecewise linear estimation algorithm RIPPLE (residual initiated polynomial-time piecewise linear estimation) intended for data from piecewise linear functions in arbitrary dimension m. The main goal of SHEPPACK is to provide users with a single consistent package containing most existing polynomial variations of Shepard’s algorithm. The algorithms target data of different dimensions. The linear Shepard algorithm, robust linear Shepard algorithm, and RIPPLE are the only algorithms in the package that are applicable to arbitrary dimensional data.

 - This code has been re-uploaded with the permission of Drs. William Thacker
   and Layne Watson.
   All comments and questions should be directed to them (see contact info at
   the bottom of this file).

## Organizational Details

The original source code, exactly as distributed by ACM TOMS, is included in
the ``src`` directory.
The ``src`` directory also contains its own ``README`` and build instructions.
Comments at the top of each subroutine document their proper usage.

Several minor modifications to the contents of ``src`` have been made:
 - The included ``Makefile`` has been slightly modified to run all tests
   when the ``make all`` command is run
 - All file extensions have been changed form ``.f95`` to ``.f90`` for
   compiler compatibility reasons.

## Reference and Contact

To cite this work, use:

```
    @article{alg905,
        author = {Thacker, William I. and Zhang, Jingwei and Watson, Layne T. and Birch, Jeffrey B. and Iyer, Manjula A. and Berry, Michael W.},
        title = {{Algorithm 905: SHEPPACK}: Modified {S}hepard Algorithm for Interpolation of Scattered Multivariate Data},
        year = {2010},
        volume = {37},
        number = {3},
        journal = {ACM Trans. Math. Softw.},
        articleno = {34},
        numpages = {20},
        doi = {10.1145/1824801.1824812}
    }
```

Inquiries should be directed to

William I. Thacker,  
Department of Computer Science, Winthrop University,  
Rock Hill, SC 29733;  
wthacker@winthrop.edu  

Layne T. Watson,  
Department of Computer Science, VPI&SU,  
Blacksburg, VA 24061-0106;  
(540) 231-7540;  
ltw@vt.edu  
