FROM python:2.7
ENV PYTHONUNBUFFERED 1

RUN mkdir /readuntil_tb
WORKDIR /readuntil_tb

## Install bwa
RUN wget https://github.com/lh3/bwa/archive/v0.7.13.tar.gz
RUN tar -xvf v0.7.13.tar.gz 
WORKDIR /readuntil_tb/bwa-0.7.13/
RUN make
RUN mv bwa /usr/bin/
WORKDIR /readuntil_tb/



ADD requirements.txt /readuntil_tb/
RUN apt-get update -y
RUN apt-get install -y emacs
RUN apt-get install -y netcdf-bin libhdf5-dev python-h5py python-numpy cython
RUN pip install --upgrade pip
RUN pip install -U -r requirements.txt
ADD . /readuntil_tb/
EXPOSE 8001
WORKDIR /readuntil_tb/nanonet/

RUN python setup.py install

WORKDIR /readuntil_tb/server/





CMD uwsgi --single-interpreter --enable-threads --processes 4 --socket :8001 --http-timeout 300 -w wsgi:app