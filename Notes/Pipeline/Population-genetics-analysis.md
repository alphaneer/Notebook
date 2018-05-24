# 群体演化分析

# 附录

EIG:

```bash
#OPENBLAS
wget https://github.com/xianyi/OpenBLAS/zipball/master 
mv master openblas.zip
unzip openblas
cd xianyi-OpenBLAS-5f998ef && make && make install PREFIX=$HOME/opt/biosoft/openblas
#GSL
wget -4 http://mirrors.ustc.edu.cn/gnu/gsl/gsl-latest.tar.gz
cd gsl-2.4
./configure --prefix=$HOME/opt/sysoft/gsl-2.4
make && make install
# EIG
git clone https://github.com/DReichLab/EIG.git
cd EIG/src
export CFLAGS="-I$HOME/opt/biosoft/openblas/include/ -I$HOME/opt/sysoft/gsl-2.4/include/"
export LDFLAGS="-L$HOME/opt/biosoft/openblas/lib -Wl,-R$HOME/opt/biosoft/openblas/lib -L$HOME/opt/sysoft/gsl-2.4/lib -Wl,-R$HOME/opt/sysoft/gsl-2.4/lib"
make
make install
```