FROM euxhen/bioenv
LABEL author="ehasanaj@cs.cmu.edu"

RUN mkdir /home/nonroot/cellar
ARG VER=unknown
RUN git clone https://github.com/euxhenh/cellar /home/nonroot/cellar

WORKDIR /home/nonroot/cellar
EXPOSE 8050
# ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "cellar", "gunicorn", "-t", "7200", "--threads", "8", "-b", "0.0.0.0:8050", "main:server"]
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "cellar", "python", "main.py"]