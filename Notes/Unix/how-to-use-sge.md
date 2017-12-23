---
title: 如何使用SGE提交任务
author: xuzhougeng
tag: UNIX
notebook: *NIX基础
---
# 如何使用SGE向运算节点提交任务

## SGE怎样工作

1. 接受用户投放的任务
1. 在任务运行以前，将任务放到一个存储区域
1. 发送任务到一个执行设备，并监控任务的运行
1. 运行结束写回结果并记录运行日志

## SGE命令的常见用法

1. 投递任务到指定队列all.q

```bash
# 方法一
qsub -cwd -l vf=*G -q all.q *.sh
# 方法二
qsub -cwd -S /bin/bash -l vf=*G -q all.q *.sh
```

>注： 方法一和方法二都可以投递任务到指定队列，但是方法一可能会输出警告信息“Warning: no access to tty (Bad file descriptor). Thus no job control in this shell.” 这是因为SGE默认使用的是tcsh，而\*.sh使用的是bash，所以应该在投递的时候指明命令解释器。若非要使用方法一的话，可以在脚本\*.sh的开头加上#$ -S /bin/bash。

参数讲解：

- `-cwd` 表示在当前路径下投递，sge的日志会输出到当前路径。
- `-l` vf=*G 任务的预估内存，内存估计的值应稍微大于真实的内存，内存预估偏小可能会导致节点跑挂。
- `-q` 指定要投递到的队列，如果不指定的话，SGE会在用户可使用的队列中选择一个满足要求的队列。

举例说明：

```bash
qsub -V -cwd -q wangjw -l vf=25g -j y -b y -sync y -N batch.sh
qsub -pe openmpi 4 -V -cwd -q wangjw -l vf=25g -j y -b y -sync y -N batch.sh
# -V: 将用户所有的当前变量输出成环境变量
# -j:
# -sync
```

1. 投递任务到指定节点

```bash
qsub -cwd -l vf=*G -l h=node1 *.sh
qsub -cwd -l vf=*G -l h=node1 -P project -q all.q *.sh
# -P 参数指明任务所属的项目
```

1. 查询任务
    qstat -f            查看所有任务
    qstat -j jobId   按任务id查看
    qstat -u user   按用户查看
    任务状态：
    qw    表示等待状态
    Eqw  投递任务出错
    r       表示任务正在运行
    dr     节点挂了之后，删除任务就会出现这个状态，只有节点重启之后，任务才会消失

1. 删除任务
    qdel -j 1111   删除任务号为1111的任务

1. 其他命令
    qrsh  与qsub相比，是交互式的投递任务，注意参数：
            -now yes|no   默认设置为yes 
                     若设置为yes，立即调度作业，如果没有可用资源，则拒绝作业，任务投递失败，任务状态为Eqw。
                     若设置为no，调度时如果没有可用资源，则将作业排入队列，等待调度。
            例子： qrsh -l vf=*G -q all.q -now no -w n *sh
    qacct  从集群日志中抽取任意账户信息
    qalter 更改已提交但正处于暂挂状态的作业的属性
    qconf 为集群和队列配置提供用户界面
    qhold 阻止已提交作业的执行
    qhost 显示SGE执行主机（即各个计算节点）的状态信息
    qlogin 启动telnet或类似的登录会话。

1. bash脚本与Linux环境变量
    为了防止脚本运行时找不到环境变量，在投递的bash脚本的前面最好加上以下两句话：
    #! /bin/bash
    #$ -S /bin/bash