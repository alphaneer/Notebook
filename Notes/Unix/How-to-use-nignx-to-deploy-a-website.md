---
title: 使用nginx作为网页服务器
tags: 网页服务器
notebook: *NIX基础
---
# 使用nginx作为网页服务器

## nginx是什么

Nginx是俄罗斯人编写的十分轻量级的HTTP服务器,Nginx，它的发音为“engine X”， 是一个高性能的HTTP和反向代理服务器，同时也是一个IMAP/POP3/SMTP 代理服务器．Nginx是由俄罗斯人 Igor Sysoev为俄罗斯访问量第二的 Rambler.ru站点开发.

Nginx以事件驱动的方式编写，所以有非常好的性能，同时也是一个非常高效的反向代理、负载平衡。其拥有匹配 Lighttpd的性能，同时还没有Lighttpd的内存泄漏问题，而且Lighttpd的mod_proxy也有一些问题并且很久没有更新。但是Nginx并不支持cgi方式运行，原因是可以减少因此带来的一些程序上的漏洞。所以必须使用**FastCGI**方式来执行PHP程序。

nginx做为HTTP服务器，有以下几项基本特性：

- 处理静态文件，索引文件以及自动索引；打开文件描述符缓冲
- 无缓存的反向代理加速，简单的负载均衡和容错
- FastCGI，简单的负载均衡和容错
- 模块化的结构。包括gzipping, byte ranges, chunked responses,以及 SSI-filter等filter。如果由FastCGI或其它代理服务器处理单页中存在的多个SSI，则这项处理可以并行运行，而不需要相互等待。
- Nginx专为性能优化而开发，性能是其最重要的考量,实现上非常注重效率。它支持内核Poll模型，能经受高负载的考验,有报告表明能支持高达 50,000个并发连接数。

Nginx具有很高的稳定性。其它HTTP服务器，当遇到访问的峰值，或者有人恶意发起慢速连接时，也很可能会导致服务器物理内存耗尽频繁交换，失去响应，只能重启服务器。例如当前apache一旦上到200个以上进程，web响应速度就明显非常缓慢了。而Nginx采取了分阶段资源分配技术，使得它的CPU与内存占用率非常低。nginx官方表示保持10,000个没有活动的连接，它只占2.5M内存，所以类似DOS这样的攻击对nginx来说基本上是毫无用处的。就稳定性而言,nginx比lighthttpd更胜一筹。

Nginx支持热部署。它的启动特别容易, 并且几乎可以做到7*24不间断运行，即使运行数个月也不需要重新启动。你还能够在不间断服务的情况下，对软件版本进行进行升级

## nignx安装

既可以简单的选择使用系统自带的管理工具apt/yum，也可以选择下载源代码进行编译。

nginx依赖于如下包：

- openssh
- zlib
- pcre

> nignx只是一个网页服务器，可以安装到任意地方

## nignx基本概念

**静态HTTP服务器**：Nginx是一个HTTP服务器，可以将服务器上的静态文件（如HTML、图片）通过HTTP协议展现给客户端。

**反向代理服务器**: 客户端本来可以直接通过HTTP协议访问某网站应用服务器，如果网站管理员在中间加上一个Nginx，客户端请求Nginx，Nginx请求应用服务器，然后将结果返回给客户端，此时Nginx就是反向代理服务器。

**负载均衡**: 当网站访问量非常大，也摊上事儿了。因为网站越来越慢，一台服务器已经不够用了。于是将相同的应用部署在多台服务器上，将大量用户的请求分配给多台机器处理。同时带来的好处是，其中一台服务器万一挂了，只要还有其他服务器正常运行，就不会影响用户使用。Nignx可以通过反响代理实现负载均衡。

![负载均衡](http://oex750gzt.bkt.clouddn.com/18-1-16/67715384.jpg)

**虚拟主机**：虽然不同域名被解析到相同IP，但是用户就能通过不同的域名打开不同的网站，互不影响，就像访问不同服务器一样，也就是虚拟主机。虚拟主机的**原理**是通过HTTP请求头中的Host是否匹配server_name来实现的，后续可以继续研究HTTP协议。

**FastCGI**: Nginx本身不支持PHP等语言，但是它可以通过FastCGI来将请求扔给某些语言或框架处理（例如PHP、Python、Perl）

## nignx常用命令

- 启动:`sudo nginx`
- 停止:`sudo nginx -s quit`
- 配置重载: `sudo nginx -s reload` 或 `sudo nginx reload`
- 指定配置文件: `sudo nginx -c /path/to/nginx.conf`
- 查看版本: `nginx -v` `nginx -V`
- 检查配置文件: `nginx -t`
- 帮助文件: `nginx -h`

## nignx配置文件

Nginx配置文件主要分成四部分：main（全局设置）、server（主机设置）、upstream（上游服务器设置，主要为反向代理、负载均衡相关配置）和 location（URL匹配特定位置后的设置），每部分包含若干个指令。main部分设置的指令将影响其它所有部分的设置；**server**部分的指令主要用于指定虚拟主机域名、IP和端口；upstream的指令用于设置一系列的后端服务器，设置反向代理及后端服务器的负载均衡；**location**部分用于匹配网页位置（比如，根目录“/”,“/images”,等等）。他们之间的关系式：server继承main，location继承server；upstream既不会继承指令也不会被继承。它有自己的特殊指令，不需要在其他地方的应用。

### 常用指令说明

#### location

http服务中，某些特定的URL对应的一系列配置项。

- `root /var/www/html`: 定义服务器的默认网站根目录位置。如果locationURL匹配的是子目录或文件，root没什么作用，一般放在server指令里面或/下。

- `index index.jsp index.html index.htm`:  定义路径下默认访问的文件名，一般跟着root放

- `proxy_pass http:/backend`: 请求转向backend定义的服务器列表，即反向代理，对应upstream负载均衡器。也可以`proxy_pass http://ip:port`。

关于location匹配规则的写法，可以说尤为关键且基础的，参考文章 nginx配置location总结及rewrite规则写法;

### 一个例子

下面的nginx.conf简单的实现nginx在前端做反向代理服务器的例子，处理js、png等静态文件，jsp等动态请求转发到其它服务器tomcat：

```bash
user  www www;
worker_processes  2;
error_log  logs/error.log;
#error_log  logs/error.log  notice;
#error_log  logs/error.log  info;

pid        logs/nginx.pid;
events {
    use epoll;
    worker_connections  2048;
}

http {
    include       mime.types;
    default_type  application/octet-stream;

    #log_format  main  '$remote_addr - $remote_user [$time_local] "$request" '
    #                  '$status $body_bytes_sent "$http_referer" '
    #                  '"$http_user_agent" "$http_x_forwarded_for"';

    #access_log  logs/access.log  main;

    sendfile        on;
    # tcp_nopush     on;

    keepalive_timeout  65;

  # gzip压缩功能设置
    gzip on;
    gzip_min_length 1k;
    gzip_buffers    4 16k;
    gzip_http_version 1.0;
    gzip_comp_level 6;
    gzip_types text/html text/plain text/css text/javascript application/json application/javascript application/x-javascript application/xml;
    gzip_vary on;

  # http_proxy 设置
    client_max_body_size   10m;
    client_body_buffer_size   128k;
    proxy_connect_timeout   75;
    proxy_send_timeout   75;
    proxy_read_timeout   75;
    proxy_buffer_size   4k;
    proxy_buffers   4 32k;
    proxy_busy_buffers_size   64k;
    proxy_temp_file_write_size  64k;
    proxy_temp_path   /usr/local/nginx/proxy_temp 1 2;

  # 设定负载均衡后台服务器列表
    upstream  backend  {
              #ip_hash;
              server   192.168.10.100:8080 max_fails=2 fail_timeout=30s ;
              server   192.168.10.101:8080 max_fails=2 fail_timeout=30s ;
    }

  # 很重要的虚拟主机配置
    server {
        listen       80;
        server_name  itoatest.example.com;
        root   /apps/oaapp;

        charset utf-8;
        access_log  logs/host.access.log  main;

        #对 / 所有做负载均衡+反向代理
        location / {
            root   /apps/oaapp;
            index  index.jsp index.html index.htm;

            proxy_pass        http://backend;
            proxy_redirect off;
            # 后端的Web服务器可以通过X-Forwarded-For获取用户真实IP
            proxy_set_header  Host  $host;
            proxy_set_header  X-Real-IP  $remote_addr;
            proxy_set_header  X-Forwarded-For  $proxy_add_x_forwarded_for;
            proxy_next_upstream error timeout invalid_header http_500 http_502 http_503 http_504;

        }

        #静态文件，nginx自己处理，不去backend请求tomcat
        location  ~* /download/ {
            root /apps/oa/fs;

        }
        location ~ .*\.(gif|jpg|jpeg|bmp|png|ico|txt|js|css)$
        {
            root /apps/oaapp;
            expires      7d;
        }
        location /nginx_status {
            stub_status on;
            access_log off;
            allow 192.168.10.0/24;
            deny all;
        }

        location ~ ^/(WEB-INF)/ {
            deny all;
        }
        #error_page  404              /404.html;

        # redirect server error pages to the static page /50x.html
        #
        error_page   500 502 503 504  /50x.html;
        location = /50x.html {
            root   html;
        }
    }

  ## 其它虚拟主机，server 指令开始
}
```

更多资料：

- <http://wiki.nginx.org/Pitfalls>
- <http://wiki.nginx.org/QuickStart>
- <http://wiki.nginx.org/Configuration>

## 参考资料

- [nginx服务器安装及配置文件详解](https://segmentfault.com/a/1190000002797601)
- [nginx配置location总结及rewrite规则写法](https://segmentfault.com/a/1190000002797606)