This repository contains c-source files for (dockerized) bcourse web service and a 
Docker container to compile them for it.

To build the docker image, say :
`docker build -t bc_progs .`

And to run run the docker image that compiles the programs, say 
`docker run --rm  -v $PWD:/root/bcourse_progs bc_progs`

You may also add gcc compilation flags like:
`docker run --rm  -v $PWD:/root/bcourse_progs bc_progs -Wall -O3`

After the compilation, there will be a a directory `bin_files` that contains
the binaries to be copied to bcourse/bin directory before running
the the [bcourse web_service](https://github.com/tomisilander/bcourse).
