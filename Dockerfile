FROM library/archlinux:latest
MAINTAINER Euxhen Hasanaj ehasanaj@cs.cmu.edu

RUN pacman -Syyu --noconfirm

RUN pacman -S --noconfirm\
        base-devel\
        git\
        vim
RUN pacman -Scc --noconfirm

RUN useradd -m nonroot
RUN echo "nonroot ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/nonroot
USER nonroot

# Install Miniconda
RUN mkdir /home/nonroot/downloads
RUN git clone https://aur.archlinux.org/miniconda3.git /home/nonroot/downloads/miniconda3
RUN cd /home/nonroot/downloads/miniconda3 && makepkg -scri --noconfirm
RUN echo "[ -f /opt/miniconda3/etc/profile.d/conda.sh ] && source /opt/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc
RUN source ~/.bashrc
RUN ls /opt

COPY env.yml /home/nonroot/downloads/

# Create conda environment
RUN /opt/miniconda3/bin/conda env create -f /home/nonroot/downloads/env.yml
ENV PATH /opt/miniconda3/bin:$PATH

# Cleanup
RUN conda clean -a && rm -rf /home/nonroot/downloads

# Init conda and close repo
RUN conda init bash
RUN mkdir /home/nonroot/cellar
ARG VER=unknown
RUN git clone https://github.com/ferrocactus/CellarV /home/nonroot/cellar

WORKDIR /home/nonroot/cellar
EXPOSE 8050
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "cellar", "python", "main.py"]
