Starchart
=========

Starchart is a regression tree-based GPU hardware/software optimization tool developed by Wenhao Jia at Princeton University.

If you use Starchart, please cite the following paper.

> Starchart: Hardware and Software Optimization Using Recursive Partitioning Regression Trees  
> Wenhao Jia, Kelly A. Shaw, Margaret Martonosi  
> Proceedings of the 22nd International Conference on Parallel Architectures and Compilation Techniques (PACT '13)

Documentation
-------------

To use Starchart, you only need to download `starchart_1.0.tar.gz` and install it using the standard [R](http://www.r-project.org/) package installation method.

To modify Starchart source code for your own projects, it may be easier to start with `starchart.R` and make it into an R package if necessary.

For detailed help including a walkthrough tutorial, please refer to [Starchart's project webpage](http://www.princeton.edu/~wjia/starchart/).

Files
-----

* `starchart.R`: The original R source file that contains all Starchart code.
* `kmeans.csv`: The original CSV file that contains the kmeans example data set used in the tutorial.
* `starchart/`: Generated by R's package.skeleton function from starchart.R for producing the final binary package.
* `starchart_1.0.tar.gz`: The binary package that can be easily distributed and installed, generated from starchart/.

License
-------

Copyright 2014: Wenhao Jia (wjia@princeton.edu), Princeton University. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this list of conditions, and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions, and the following disclaimer in the documentation and/or other materials provided with the distribution.
* Neither the name of Princeton University nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY PRINCETON UNIVERSITY "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PRINCETON UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
