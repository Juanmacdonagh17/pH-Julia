# This dockerfile is heavily influenced by the one created by Fons VanderPlas @ https://github.com/JuliaPluto/PlutoSliderServer.jl for Introduction to Computational Thinking @MIT
FROM julia:1.7-buster
RUN apt-get update -y && apt-get upgrade -y
RUN useradd -ms /bin/bash pluto
WORKDIR /home/pluto
USER pluto
COPY --chown=pluto . ${HOME}
RUN julia --project=${HOME}/env -e "import Pkg; Pkg.instantiate(); Pkg.precompile()"
CMD ["julia", "--project=/home/pluto/env","start.jl"]
