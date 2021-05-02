# ChroMo

## Prerequisites

- Docker
- docker-compose
- ShinyProxy (optional)

## Running ChroMo

Clone the repository in your local computer or server
<code>git clone https://github.com/danilexn/chromo</code>

<code>cd chromo</code>

### ...with Docker

<code>docker build -t chromo .</code>

<code>docker run -p 127.0.0.1:8080:3838 chromo</code>

Now, open a browser at http://127.0.0.1:8080, and the application should be ready

### ...with docker-compose

<code>docker build -t chromo .</code>

<code>sudo docker-compose up -d</code>

Now, open a browser at http://127.0.0.1:80, and the application should be ready.
You can modify the nginx.conf file in the directory net to enable HTTPS.
Remember to mount the folder with certificates in docker-compose.yml, under the webserver section.
