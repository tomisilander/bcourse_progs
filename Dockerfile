FROM debian:jessie
RUN apt-get update && apt-get install -y libjudy-dev gcc
ENTRYPOINT ["bash", "/root/bcourse_progs/build_all.sh"]
