# tronos-shiny
## Prerequisites
- Docker
- docker-compose

## Running tronos-shiny
Clone the repository in your local computer or server
<code>git clone https://github.com/danilexn/tronos-shiny</code>
<code>cd tronos-shiny</code>

### ...with Docker
<code>docker build -t tronos-shiny .</code>
<code>docker run -p 127.0.0.1:8080:3838 tronos-shiny</code>

Now, open a browser at http://127.0.0.1:8080, and the application should be ready

### ...with docker-compose
<code>sudo docker-compose up -d</code>

Now, open a browser at http://127.0.0.1:80, and the application should be ready.
You can modify the nginx.conf file in the directory net to enable HTTPS.
Remember to mount the folder with certificates in docker-compose.yml, under the webserver section.
