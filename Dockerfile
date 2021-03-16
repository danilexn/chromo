# Base image https://hub.docker.com/u/rocker/
FROM danilexn/tronos-base:latest

# copy necessary files
## app folder
COPY ./src/ /srv/shiny-server/tronos
COPY ./net/shiny-server.conf /etc/shiny/shiny-server.conf
# expose port
EXPOSE 3838
