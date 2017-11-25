---
title: 使用远程服务器的几个技巧
tags: unix, SSH, 服务器
notebook: *NIX基础
---

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->

# 使用远程服务器的几个技巧

如果你有一台高性能的个人PC，那么绝大部分任务都能在本地解决。但是更为常见的情况是，你买了一台小型服务器或者有专门的服务器平台，然后在本地进行访问，那么这里有几个小技巧，你可以了解下。

## 保存常用SSH主机，避免重复输入

Linux下访问远程服务器的方法基本就是使用SSH(secure shell)。SSH其实是一种网络协议，用来计算机之间的加密登陆，保证传输过程中即便信息被截获了也无法解析出原始信息。SSH作为一种协议有多种使用方法，在Linux里面是OpenSSH，在Windows里面可以用putty或者Xshell。不过这里仅仅讨论Linux 里面的ssh。在Linux中访问远程主机的时候，大家肯定对如下指令不陌生:

```shell
ssh -p 22 xuzhougeng@10.10.87.36
```

其中`-p`指定端口号，如果远程服务器没有特殊说明，一般默认都是22，所以可以省去`-p 22`. 后面为用户@IP地址。第一次访问的时候会问你是否要将该主机的公钥加入信任名单中，当然是选`yes`了。

那么问题来了，能不能讨论不要输入"xuzhougegn@10.10.87.36"呢？方法当然是有的，你只需要创建`~/.ssh/config`文件，并添加主机信息

```shell
# 使用vi编辑器
# vi ~/.ssh/config
Host xzg
    HostName 10.10.87.36
    User xuzhougegn
    Port 22
```

然后就能以`ssh xzg`访问远程主机，而需要输入全部信息，又累还容易出错。

## 无需密码认证，快速登陆

上面的技巧使得你访问的时候不需要输入主机全称，但是依旧需要输入密码。如果密码比较长，那么人就容易出错，而且Linux输密码的时候啥都看不见，你都不知道自己输了多少个字符，如果你的电脑只有自己用，完全连密码输入这一步都可以省呀。

避免每次都要密码认证的方法就是使用SSH公钥。当你把个人电脑的SSH公钥存放到远程服务器的时候，远程服务器就完全信任了你，两个人之间就再也没有了隔阂。

首先是用`ssh-keygen`创建密钥（密码为空时，后续登陆时才能不需要输入密码）

```shell
$ ssh-keygen -b 2048
Generating public/private rsa key pair.
Enter file in which to save the key (/home/xzg/.ssh/id_rsa):
/home/xzg/.ssh/id_rsa already exists.
Enter passphrase (empty for no passphrase):
Enter same passphrase again:
Your identification has been saved in /home/xzg/.ssh/id_rsa.
Your public key has been saved in /home/xzg/.ssh/id_rsa.pub.
The key fingerprint is:
SHA256:/1dusQI7WRUfsf9pG3CjCstIzGHPdPyRCO+zVF/efVM xzg@DESKTOP-CNF0I9C
The key's randomart image is:
+---[RSA 2048]----+
|               ..|
|               o.|
|         .     .+|
|          + . ..o|
|       oS. = =.oE|
|      + =.o.o.*oX|
|       + +.+=o OO|
|      . o +=+.o.B|
|       . o oo..o |
+----[SHA256]-----+
```

然后是将自己的ssh公钥添加到远程服务器的`~/.ssh/authorized_keys`.这个方法比较多

- ssh-copy-id

```shell
$ ssh-copy-id xuzhougegn@10.10.87.36
/usr/bin/ssh-copy-id: INFO: Source of key(s) to be installed: "/home/xzg/.ssh/id_rsa.pub"
/usr/bin/ssh-copy-id: INFO: attempting to log in with the new key(s), to filter out any that are already installed
/usr/bin/ssh-copy-id: INFO: 1 key(s) remain to be installed -- if you are prompted now it is to install the new keys
xuzhougegn@10.10.87.36's password:

Number of key(s) added: 1

Now try logging into the machine, with:   "ssh 'xuzhougegn@10.10.87.36'"
and check to make sure that only the key(s) you wanted were added.
```

- 用`cat id_rsa.pub`显示，然后复制到远程服务器的`~/.ssh/authorized_keys`中
- `ssh xuzhougegn@10.10.87.36' mkdir -p .ssh && cat >> .ssh/authorized_keys' < ~/.ssh/id_rsa.pub`

## 任务挂起，安心关闭终端

远程操作时，一旦终端关闭，所有这个终端运行的进程都会收到`SIGHUP`信号，然后这些程序就会立即退出。如果你的命令需要运行好几个小时或者好几天，你肯定不愿意一直开着终端，因为网络问题功亏一篑。解决方法也是有的，而且还有好几种，这里就说说`nohup`和`screen`

`nohup`故名思意，就是不要hup，即能够捕捉到终端发出的SIGHUP信号并无视他，就不必担心自己的命令被终端关闭了。

```shell
# 实例
nohup bash snp_calling.sh > output.txt
```

如果你不知道自己的命令有多久，所以不知道要不要输入nohup, 其实用`screen`(或tmux)效果更好

```shell
# 开启一个screen
screen -S hisat2
# 使用ctrl +a ctrl +d 挂起
screen -list # 查看运行的sreen
There is a screen on:
        46953.hista2    (Detached)
1 Socket in /var/run/screen/S-xuzhougegn.
# 继续之前的screen
screen -r hisat2
```

## 总结

本文介绍三种连接远程服务器的小技巧

- 利用ssh config添加常用服务器
- 利用ssh公钥避免重复输入密码
- 利用screen/nohup 长时间运行程序