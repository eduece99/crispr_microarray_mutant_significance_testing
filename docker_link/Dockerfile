# syntax=docker/dockerfile:1

FROM python:3.8-slim-buster
WORKDIR /app
COPY requirements.txt requirements.txt


RUN apt-get update
RUN apt-get -y install gcc 
RUN apt-get -y install pkg-config
RUN apt-get -y install libcairo2-dev libdbus-1-3 libdbus-1-dev

RUN pip3 install -r requirements.txt


COPY . .

CMD [ "python3", "-m" , "flask", "run", "--host=0.0.0.0"]
