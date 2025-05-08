#!/bin/bash
#																		                                                                            #
# Passo a Passo para compilação das bibliotecas e para compilação do modelo WRF com compilador gfortran.                                            #
#																		                                                                            #
# Este script esta preparado para instalar todas as bibliotecas necessárias para compilar o modelo WRF no diretório definido abaixo como DIR.       #
#																		                                                                            #	
# O script segue as "Variáveis de Controle" para cada etapa de instalação, assim se caso ocorrer algum erro durante o processamento de uma etapa    #
# basta o usuário alterar a flag "SIM" para "NAO" para o script seguir a partir da última configurada como "SIM".                                   #
#																		                                                                            #
# O script pode ser utilizado também como passo a passo para instalação manual, bastando o usuário seguir a lógica dos caminhos e dos comandos e    #
# executando em sua ordem no terminal do linux.                                                                                                     #
# O script encerra até o download das pastas com os arquivos do WPS (pré-processamento) e do WRF (processamento)                                    #                #    
# Ricardo Mollmann - criado em 05/2025                                                                                                              # 
#																		                                                                            #
# Referencial de instalação:                                                                                                                        #
# base WRF https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compilation_tutorial.php#STEP1                                                              #
# base BIBLIOTECAS https://forum.mmm.ucar.edu/threads/full-wrf-and-wps-installation-example-gnu.12385/                                              #
# download dos dados de contorno: https://www2.mmm.ucar.edu/wrf/src/wps_files/                                                                      #
# comando para baixa resolução: wget https://www2.mmm.ucar.edu/wrf/src/wps_files/geog_low_res_mandatory.tar.gz                                      #
# comando para alta resolução: wget https://www2.mmm.ucar.edu/wrf/src/wps_files/geog_high_res_mandatory.tar.gz                                      #
#####################################################################################################################################################
# 
set -e

# Caminho para as bibliotecas
mkdir $HOME/modelos
mkdir $HOME/modelos/wrf_bibliotecas
DIR=$HOME/modelos/wrf_bibliotecas
mkdir $DIR/grib2
mkdir $DIR/netcdf
INSTALL_PREFIX_GRIB2=$DIR/grib2
INSTALL_PREFIX_NETCDF=$DIR/netcdf
export CC=gcc
export CXX=g++
export FC=gfortran
export FCFLAGS="-m64 -fallow-argument-mismatch"
export F77=gfortran
export FFLAGS="-m64 -fallow-argument-mismatch"
export LDFLAGS="-L$NETCDF/lib -L$DIR/grib2/lib"
export CPPFLAGS="-I$NETCDF/include -I$DIR/grib2/include -fcommon"

# Variáveis de controle
SISTEMA="SIM"
ZLIB="SIM"     # Biblioteca de compressão básica utilizada por outras bibliotecas (como HDF5). Necessária para manipular arquivos compactados com eficiência.
HDF5="SIM"     # Formato de armazenamento hierárquico usado pelo NetCDF-C. Essencial para leitura e escrita de arquivos NetCDF-4.
NETCDF="SIM"   # Biblioteca principal para entrada e saída de dados no WRF/WPS. Armazena variáveis geofísicas com metadados em arquivos auto-descritivos.
NETCDFF="SIM"  # Interface Fortran da biblioteca NetCDF. Necessária para o WRF, que é escrito majoritariamente em Fortran.
MPICH="SIM"    # Implementação de MPI (Message Passing Interface). Permite execução paralela do WRF em múltiplos núcleos/processadores.
LIBPNG="SIM"   # Biblioteca para manipulação de imagens PNG. Usada principalmente durante a compilação do WPS para visualização de domínios.
JASPER="SIM"   # Biblioteca para compressão JPEG-2000. Necessária para lidar com arquivos GRIB2 no WPS, que utiliza essa compressão para campos meteorológicos.
WRF="SIM" 
WPS="SIM"

mkdir -p $INSTALL_PREFIX_GRIB2 $INSTALL_PREFIX_NETCDF
cd $DIR

if [[ $SISTEMA == "SIM" ]]; then
    echo "Instalando dependências do sistema..."
     sudo apt update
     sudo apt install -y build-essential gfortran m4 libcurl4-openssl-dev
fi

if [[ $ZLIB == "SIM" ]]; then
    echo "Instalando ZLIB"
    wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/zlib-1.2.11.tar.gz
    tar xzvf zlib-1.2.11.tar.gz
    cd zlib-1.2.11
    ./configure --prefix=$INSTALL_PREFIX_GRIB2
    make -j 4
    make install
    cd ..
    rm -rf zlib*
fi

if [[ $HDF5 == "SIM" ]]; then
    echo "Instalando HDF5"
    wget https://github.com/HDFGroup/hdf5/archive/hdf5-1_10_5.tar.gz
    tar -xzf hdf5-1_10_5.tar.gz
    cd hdf5-hdf5-1_10_5
    export CPPFLAGS="-I$INSTALL_PREFIX_GRIB2/include"
    export LDFLAGS="-L$INSTALL_PREFIX_GRIB2/lib"
    export LD_LIBRARY_PATH="$INSTALL_PREFIX_GRIB2/lib:$LD_LIBRARY_PATH"
    ./configure --prefix=$INSTALL_PREFIX_GRIB2 --with-zlib=$INSTALL_PREFIX_GRIB2 --enable-fortran --enable-cxx --enable-shared
    make -j 4
    make install
    cd ..
    rm -rf hdf5-*
fi

if [[ $NETCDF == "SIM" ]]; then
    echo "Instalando NetCDF-C"
    wget https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.2.tar.gz -O netcdf-c-4.9.2.tar.gz
    tar -xzf netcdf-c-4.9.2.tar.gz
    cd netcdf-c-4.9.2
    export CPPFLAGS="-I$INSTALL_PREFIX_GRIB2/include"
    export LDFLAGS="-L$INSTALL_PREFIX_GRIB2/lib"
    export LD_LIBRARY_PATH="$INSTALL_PREFIX_GRIB2/lib:$LD_LIBRARY_PATH"
    ./configure --prefix=$INSTALL_PREFIX_NETCDF --enable-netcdf-4 --enable-shared
    make -j 4
    make install
    cd ..
fi

if [[ $NETCDFF == "SIM" ]]; then
    echo "Instalando NetCDF-Fortran"
    wget https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.6.1.tar.gz -O netcdf-fortran-4.6.1.tar.gz
    tar -xzf netcdf-fortran-4.6.1.tar.gz
    cd netcdf-fortran-4.6.1
    export CPPFLAGS="-I$INSTALL_PREFIX_NETCDF/include"
    export LDFLAGS="-L$INSTALL_PREFIX_NETCDF/lib"
    export LD_LIBRARY_PATH="$INSTALL_PREFIX_NETCDF/lib:$LD_LIBRARY_PATH"
    ./configure --prefix=$INSTALL_PREFIX_NETCDF --enable-shared
    make -j 4
    make install
    cd ..
    rm -rf netcdf-*
fi

if [[ $MPICH == "SIM" ]]; then
    echo "Instalando MPICH"
    wget https://www.mpich.org/static/downloads/3.4.3/mpich-3.4.3.tar.gz
    tar xzvf mpich-3.4.3.tar.gz
    cd mpich-3.4.3
    ./configure --prefix=$DIR/mpich --with-device=ch3 FFLAGS="-fallow-argument-mismatch"
    make -j 4
    make install
    cd ..
    rm -rf mpich-3.4.3.tar.gz
fi

if [[ $LIBPNG == "SIM" ]]; then
    echo "Instalando LIBPNG"
    wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/libpng-1.2.50.tar.gz
    tar xzvf libpng-1.2.50.tar.gz
    cd libpng-1.2.50
    ./configure --prefix=$INSTALL_PREFIX_GRIB2
    make -j 4
    make install
    cd ..
    rm -rf libpng*
fi

if [[ $JASPER == "SIM" ]]; then
    echo "Instalando JASPER"
    wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/jasper-1.900.1.tar.gz
    tar xzvf jasper-1.900.1.tar.gz
    cd jasper-1.900.1
    ./configure --prefix=$INSTALL_PREFIX_GRIB2
    make -j 4
    make install
    cd ..
    rm -rf jasper*
    cd ..
fi

# Exportações finais
export NETCDF=$DIR/netcdf
export LD_LIBRARY_PATH=$NETCDF/lib:$DIR/grib2/lib:$DIR/mpich/lib:$LD_LIBRARY_PATH
export PATH=$NETCDF/bin:$DIR/mpich/bin:$PATH
export JASPERLIB=$DIR/grib2/lib
export JASPERINC=$DIR/grib2/include


if [[ $WRF == "SIM" ]]; then
    echo "Clonando e preparando WRF"
    git clone --recurse-submodule https://github.com/wrf-model/WRF.git
    echo "Entre na pasta WRF e execute ./configure (escolha opção 34 e depois 1), depois ./compile em_real"
fi

if [[ $WPS == "SIM" ]]; then
    echo "Clonando WPS"
    git clone https://github.com/wrf-model/WPS.git
    echo "Entre na pasta WPS, defina WRF_DIR corretamente, execute ./configure (opção 1) e ./compile"
    fi

