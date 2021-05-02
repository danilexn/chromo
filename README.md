# tronos-shiny

## Prerequisites

- Docker
- docker-compose

## Running ChroMo

Clone the repository in your local computer or server
<<<<<<< HEAD
<code>git clone https://github.com/danilexn/chromo</code>
=======

<code>git clone https://github.com/danilexn/tronos-shiny</code>
>>>>>>> b300e1cb9c641b97f96b5204623293f0e345c572

<code>cd chromo</code>

### ...with Docker

<code>docker build -t chromo .</code>

<code>docker run -p 127.0.0.1:8080:3838 chromo

Now, open a browser at http://127.0.0.1:8080, and the application should be ready

### ...with docker-compose
<<<<<<< HEAD
=======
<code>docker build -t tronos-shiny .</code>
>>>>>>> b300e1cb9c641b97f96b5204623293f0e345c572

<code>sudo docker-compose up -d</code>

Now, open a browser at http://127.0.0.1:80, and the application should be ready.
You can modify the nginx.conf file in the directory net to enable HTTPS.
Remember to mount the folder with certificates in docker-compose.yml, under the webserver section.
