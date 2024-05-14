# Use Miniconda as the base image
FROM continuumio/miniconda3

# Add metadata
LABEL author="Steffen Kortmann"

# Set the application directory in the container
ENV APP_HOME /app
WORKDIR $APP_HOME
COPY . $APP_HOME

# Update Conda in the base environment, create the specified environment
RUN conda update --name base --yes conda && \
    conda env create --name potpourri_env --file environment.yml

# Use conda run to execute commands within the Conda environment
RUN echo "source activate potpourri_env" > ~/.bashrc
ENV PATH /opt/conda/envs/potpourri_env/bin:$PATH

# Update and install system dependencies including CBC
RUN apt-get update && \
    apt-get install -y build-essential cmake git && \
    apt-get install -y coinor-cbc coinor-libcbc-dev

RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    g++ \
    gfortran \
    git \
    patch \
    wget \
    pkg-config \
    liblapack-dev \
    libmetis-dev

RUN cd /opt && \
    wget https://github.com/coin-or/Ipopt/archive/refs/tags/releases/3.14.16.tar.gz && \
    tar xzf 3.14.16.tar.gz && \
    mv Ipopt-releases-3.14.16 Ipopt && \
    cd Ipopt && \
    ./configure --prefix=/usr/local && \
    make && \
    make install

# Clone and build SHOT
RUN git clone https://www.github.com/coin-or/SHOT --recursive && \
    cd SHOT && \
    mkdir build && \
    cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release \
        -DHAS_CBC=on -DCBC_DIR=/opt/cbc \
        -DHAS_IPOPT=on -DIPOPT_DIR=/opt/Ipopt && \
    make && make install

# Set the default shell to run inside the 'potpourri_env' environment
SHELL ["conda", "run", "--name", "potpourri_env", "/bin/bash", "-c"]

# Add any additional commands or operations below
# Example: CMD ["python", "your_script.py"]
