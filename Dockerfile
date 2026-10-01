# syntax=docker/dockerfile:1
#
# potpourri container: the conda environment from environment.yaml plus
# IPOPT (with MUMPS and ASL), CBC and SHOT installed system-wide.
# Usage is documented in docs/getting-started.md ("Option B - Docker").
#
# Layer order is deliberate: the slow solver builds come first and depend
# only on the pinned versions below, so editing the Python sources or
# environment.yaml never triggers a recompile.

# continuumio/miniconda3 is deprecated after Miniconda 26.7; anaconda/miniconda
# is its successor (Debian 13 "trixie", conda 26.7.1).
ARG MINICONDA_TAG=26.7.1
FROM anaconda/miniconda:${MINICONDA_TAG}

LABEL author="Steffen Kortmann"

# ---------------------------------------------------------------------------
# Pinned solver sources. Override at build time, e.g.
#   docker build --build-arg BUILD_JOBS=16 --build-arg SHOT_COMMIT=<sha> .
# ---------------------------------------------------------------------------
# IPOPT release (https://github.com/coin-or/Ipopt/releases); the same version
# is pinned for conda-forge in environment.yaml.
ARG IPOPT_VERSION=3.14.20
# coin-or-tools/ThirdParty-ASL and ThirdParty-Mumps build helpers. ASL gives
# IPOPT its AMPL (.nl) interface; MUMPS is the sparse linear solver IPOPT
# needs (the old recipe configured IPOPT without one).
ARG ASL_BUILD_VERSION=2.1.0
ARG MUMPS_BUILD_VERSION=3.0.14
# SHOT has no tagged release since 1.1.0 (2021); development happens on master,
# so a commit is pinned instead (master as of 2026-09-27).
ARG SHOT_COMMIT=c5cfcdb011ec33b96205a43312df3360fb840e2c
# Parallel build jobs for the C++/Fortran builds.
ARG BUILD_JOBS=8

ARG IPOPT_PREFIX=/opt/ipopt
ENV APP_HOME=/app
ENV DEBIAN_FRONTEND=noninteractive

# ---------------------------------------------------------------------------
# Toolchain and system libraries. CBC comes from Debian (2.10.12 in trixie):
# it is the Pyomo "cbc" solver and one of SHOT's MIP back-ends.
# ---------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        cmake \
        git \
        patch \
        pkg-config \
        wget \
        ca-certificates \
        liblapack-dev \
        libblas-dev \
        libmetis-dev \
        zlib1g-dev \
        libbz2-dev \
        coinor-cbc \
        coinor-libcbc-dev \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# IPOPT from source into ${IPOPT_PREFIX}: ASL, then MUMPS, then IPOPT itself,
# following https://coin-or.github.io/Ipopt/INSTALL.html. The install prefix
# is registered with the dynamic linker so that SHOT finds libipopt at run
# time without LD_LIBRARY_PATH.
# ---------------------------------------------------------------------------
RUN git clone --quiet --depth 1 --branch releases/${ASL_BUILD_VERSION} \
        https://github.com/coin-or-tools/ThirdParty-ASL /opt/src/ThirdParty-ASL \
    && cd /opt/src/ThirdParty-ASL \
    && ./get.ASL \
    && ./configure --prefix=${IPOPT_PREFIX} \
    && make -j${BUILD_JOBS} \
    && make install \
    && rm -rf /opt/src/ThirdParty-ASL

RUN git clone --quiet --depth 1 --branch releases/${MUMPS_BUILD_VERSION} \
        https://github.com/coin-or-tools/ThirdParty-Mumps /opt/src/ThirdParty-Mumps \
    && cd /opt/src/ThirdParty-Mumps \
    && ./get.Mumps \
    && ./configure --prefix=${IPOPT_PREFIX} \
    && make -j${BUILD_JOBS} \
    && make install \
    && rm -rf /opt/src/ThirdParty-Mumps

RUN mkdir -p /opt/src \
    && wget -q -O- https://github.com/coin-or/Ipopt/archive/refs/tags/releases/${IPOPT_VERSION}.tar.gz \
        | tar xz -C /opt/src \
    && cd /opt/src/Ipopt-releases-${IPOPT_VERSION} \
    && PKG_CONFIG_PATH=${IPOPT_PREFIX}/lib/pkgconfig ./configure --prefix=${IPOPT_PREFIX} \
    && make -j${BUILD_JOBS} \
    && make install \
    && echo "${IPOPT_PREFIX}/lib" > /etc/ld.so.conf.d/ipopt.conf \
    && ldconfig \
    && rm -rf /opt/src/Ipopt-releases-${IPOPT_VERSION}

# ---------------------------------------------------------------------------
# SHOT (https://github.com/coin-or/SHOT) with IPOPT, CBC and the bundled HiGHS;
# GAMS/CPLEX/Gurobi are switched off explicitly. The binary is called SHOT;
# Pyomo's generic ASL interface (``SolverFactory("shot")``) looks for a
# lowercase executable, hence the symlink.
# ---------------------------------------------------------------------------
RUN mkdir -p /opt/src/SHOT && cd /opt/src/SHOT \
    && git init --quiet \
    && git remote add origin https://github.com/coin-or/SHOT \
    && git fetch --quiet --depth 1 origin ${SHOT_COMMIT} \
    && git checkout --quiet FETCH_HEAD \
    && git submodule update --quiet --init --recursive --depth 1 \
    && cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DHAS_CBC=on -DCBC_DIR=/usr \
        -DHAS_IPOPT=on -DIPOPT_DIR=${IPOPT_PREFIX} \
        -DHAS_HIGHS=on \
        -DHAS_AMPL=on \
        -DHAS_GAMS=off -DHAS_CPLEX=off -DHAS_GUROBI=off \
    && cmake --build build -j${BUILD_JOBS} \
    && cmake --install build \
    # SHOT builds CppAD as a shared library inside its build tree and links
    # against it, but its install step does not copy it; without this the
    # installed binary fails to load libcppad_lib.so.
    && cp -a build/CppAD/lib*/libcppad_lib.so* /usr/local/lib/ \
    && ln -s /usr/local/bin/SHOT /usr/local/bin/shot \
    && ldconfig \
    && cd / && rm -rf /opt/src/SHOT \
    && ! ldd /usr/local/bin/SHOT | grep "not found"

# ---------------------------------------------------------------------------
# Conda environment. The Anaconda channels listed in environment.yaml require
# their Terms of Service to be accepted before conda will use them
# non-interactively.
# ---------------------------------------------------------------------------
# Where the base image installs conda: /opt/miniconda3 in anaconda/miniconda
# (the old continuumio/miniconda3 image used /opt/conda). The environment
# stage fails fast if this is wrong. Declared here, not at the top, because
# an ARG invalidates the build cache of every RUN that follows it.
ARG CONDA_ROOT=/opt/miniconda3
ARG CONDA_ENV=potpourri_env
WORKDIR ${APP_HOME}
COPY environment.yaml ${APP_HOME}/environment.yaml
RUN conda tos accept --override-channels \
        --channel https://repo.anaconda.com/pkgs/main \
        --channel https://repo.anaconda.com/pkgs/r \
    && conda env create --name ${CONDA_ENV} --file environment.yaml \
    && test -x ${CONDA_ROOT}/envs/${CONDA_ENV}/bin/python \
    && conda clean --all --yes \
    && conda init bash \
    && echo "conda activate ${CONDA_ENV}" >> /root/.bashrc

ENV PATH=${CONDA_ROOT}/envs/${CONDA_ENV}/bin:${IPOPT_PREFIX}/bin:$PATH

# ---------------------------------------------------------------------------
# The package itself, installed editable so that mounting a checkout over
# /app (see the docs) picks up local edits.
# ---------------------------------------------------------------------------
COPY . ${APP_HOME}
RUN pip install --no-deps --editable .

CMD ["bash"]
