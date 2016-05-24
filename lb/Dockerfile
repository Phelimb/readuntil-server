FROM nginx
COPY readuntil.conf /etc/nginx/conf.d/default.conf.tmp

CMD sed  s/{{WEB_PORT_8001_TCP_ADDR}}/$WEB_PORT_8001_TCP_ADDR/g /etc/nginx/conf.d/default.conf.tmp > /etc/nginx/conf.d/default.conf.tmp2 && sed  s/{{WEB_PORT_8001_TCP_PORT}}/$WEB_PORT_8001_TCP_PORT/g /etc/nginx/conf.d/default.conf.tmp2 > /etc/nginx/conf.d/default.conf && nginx -g "daemon off;"

