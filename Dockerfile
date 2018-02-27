#Docker file for local building and serving only
FROM ubuntu:16.04
MAINTAINER James Scott-Brown <james@jamesscottbrown.com>

RUN apt-get update
RUN apt-get install wget -y


WORKDIR /root

# dReal
RUN wget https://github.com/dreal/dreal3/releases/download/v3.16.09.01/dReal-3.16.09.01-linux.tar.gz
RUN tar -zxvf dReal-3.16.09.01-linux.tar.gz

# MathSAT
RUN wget http://mathsat.fbk.eu/release/mathsat-5.5.1-linux-x86_64.tar.gz
RUN tar -zxvf mathsat-5.5.1-linux-x86_64.tar.gz

# iSAT binary
RUN wget http://www.avacs.org/fileadmin/tooladmin/open/isat-ode.tar.gz
RUN tar -zxvf isat-ode.tar.gz




## This site
RUN apt-get install -y python binutils g++ make sqlite3 python-pip git

RUN pip install --upgrade pip

ADD . /code
WORKDIR /code
RUN pip install -r requirements/dev.txt

ENV FLASK_APP=/code/autoapp.py
ENV FLASK_DEBUG=1
RUN flask db init
RUN flask db migrate
RUN flask db upgrade

CMD flask run --host=0.0.0.0


 


# root@f94d706413ea:~# dReal-3.16.09.01-linux/bin/dReal
# mathsat-5.5.1-linux-x86_64/bin/mathsat
