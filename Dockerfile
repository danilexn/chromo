# Base image https://hub.docker.com/u/rocker/
FROM tronos-testing:latest

# copy necessary files
## app folder
COPY /Tronos ./tronos
# ## renv.lock file
COPY /Tronos/shiny-server.conf /etc/shiny-server/shiny-server.conf
COPY /Tronos/nginx.conf /etc/nginx/nginx.conf
# expose port
EXPOSE 3838