FROM debian:jessie
RUN apt-get update && apt-get install -y libjudy-dev gcc
CMD cd /root/bcourse_progs; bash build_all.sh
