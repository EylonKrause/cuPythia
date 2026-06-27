#!/bin/bash
# Build the Rivet analysis toolchain (fastjet -> YODA -> Rivet) from source to ~/.local,
# NO sudo / NO pip, via --disable-pyext (avoids the missing Python.h). Matched to the
# already-installed HepMC3 3.03.01. Logs each phase; stops on first error (set -e).
set -e
PREFIX=/home/eylonk/.local
export PATH=$PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
cd /home/eylonk/rivetbuild
J=16
DL(){ rm -f "$2"; curl -L --retry 6 --retry-delay 5 --retry-all-errors "$1" -o "$2"; }   # fresh robust download (no resume)

echo "=== [1/3] fastjet 3.4.0 ==="
if [ -x $PREFIX/bin/fastjet-config ]; then
  echo "fastjet already installed: $($PREFIX/bin/fastjet-config --version) -- skipping"
else
  if [ ! -f fastjet-3.4.0.tar.gz ]; then DL http://fastjet.fr/repo/fastjet-3.4.0.tar.gz fastjet-3.4.0.tar.gz; fi
  rm -rf fastjet-3.4.0 && tar xf fastjet-3.4.0.tar.gz && cd fastjet-3.4.0
  ./configure --prefix=$PREFIX --enable-shared --disable-auto-ptr --enable-allcxxplugins > cfg.log 2>&1
  make -j$J > make.log 2>&1
  make install > install.log 2>&1
  cd ..
  echo "fastjet done: $($PREFIX/bin/fastjet-config --version 2>&1)"
fi

echo "=== [1a/3] fjcontrib 1.051 (Rivet 3.1.x requires the contrib headers + fragile shared lib) ==="
if [ -f $PREFIX/include/fastjet/contrib/SoftDrop.hh ]; then echo "fjcontrib already installed -- skipping"; else
  if [ ! -s fjcontrib-1.051.tar.gz ] || ! gzip -t fjcontrib-1.051.tar.gz 2>/dev/null; then DL "https://fastjet.hepforge.org/contrib/downloads/fjcontrib-1.051.tar.gz" fjcontrib-1.051.tar.gz; fi
  rm -rf fjcontrib-1.051 && tar xf fjcontrib-1.051.tar.gz && cd fjcontrib-1.051
  ./configure --fastjet-config=$PREFIX/bin/fastjet-config --prefix=$PREFIX > cfg.log 2>&1
  make -j$J > make.log 2>&1
  make install > install.log 2>&1
  make fragile-shared -j$J > fragile.log 2>&1
  make fragile-shared-install >> fragile.log 2>&1
  cd ..
  echo "fjcontrib done: $(ls $PREFIX/include/fastjet/contrib/SoftDrop.hh 2>&1)"
fi

echo "=== [1b/3] zlib (no system dev headers; needed by YODA + Rivet ref data) ==="
if [ -f $PREFIX/include/zlib.h ]; then echo "zlib already installed -- skipping"; else
  if [ ! -s zlib-1.3.1.tar.gz ] || ! gzip -t zlib-1.3.1.tar.gz 2>/dev/null; then
    DL "https://github.com/madler/zlib/releases/download/v1.3.1/zlib-1.3.1.tar.gz" zlib-1.3.1.tar.gz; fi
  rm -rf zlib-1.3.1 && tar xf zlib-1.3.1.tar.gz && cd zlib-1.3.1
  ./configure --prefix=$PREFIX > cfg.log 2>&1
  make -j$J > make.log 2>&1
  make install > install.log 2>&1
  cd ..
  echo "zlib done: $(ls $PREFIX/include/zlib.h 2>&1)"
fi

echo "=== [2/3] YODA 1.9.11 (--disable-pyext) ==="
if [ -x $PREFIX/bin/yoda-config ]; then echo "YODA already installed: $($PREFIX/bin/yoda-config --version) -- skipping"; else
  if [ ! -s YODA-1.9.11.tar.bz2 ]; then DL "https://yoda.hepforge.org/downloads/?f=YODA-1.9.11.tar.bz2" YODA-1.9.11.tar.bz2; fi
  rm -rf YODA-1.9.11 && tar xf YODA-1.9.11.tar.bz2 && cd YODA-1.9.11
  ./configure --prefix=$PREFIX --disable-pyext --with-zlib=$PREFIX > cfg.log 2>&1
  make -j$J > make.log 2>&1
  make install > install.log 2>&1
  cd ..
  echo "YODA done: $($PREFIX/bin/yoda-config --version 2>&1)"
fi

echo "=== [2b/3] HepMC3 3.3.1 rebuild WITH Search module (installs Relatives.h, needed by Rivet) ==="
if [ -f $PREFIX/include/HepMC3/Relatives.h ]; then echo "HepMC3 Search (Relatives.h) already present -- skipping"; else
  if [ ! -s HepMC3-3.3.1.tar.gz ] || ! gzip -t HepMC3-3.3.1.tar.gz 2>/dev/null; then DL "https://gitlab.cern.ch/hepmc/HepMC3/-/archive/3.3.1/HepMC3-3.3.1.tar.gz" HepMC3-3.3.1.tar.gz; fi
  rm -rf HepMC3-3.3.1 && tar xf HepMC3-3.3.1.tar.gz && cd HepMC3-3.3.1
  $PREFIX/bin/cmake -B bld -DCMAKE_INSTALL_PREFIX=$PREFIX -DHEPMC3_ENABLE_SEARCH=ON -DHEPMC3_ENABLE_PYTHON=OFF -DHEPMC3_ENABLE_ROOTIO=OFF -DHEPMC3_BUILD_EXAMPLES=OFF -DHEPMC3_ENABLE_TEST=OFF > cmake.log 2>&1
  $PREFIX/bin/cmake --build bld -j$J > cmakebuild.log 2>&1
  $PREFIX/bin/cmake --install bld > cmakeinstall.log 2>&1
  cd ..
  echo "HepMC3 Search done: $(ls $PREFIX/include/HepMC3/Relatives.h 2>&1)  $(ls $PREFIX/lib/libHepMC3search* 2>&1)"
fi

echo "=== [3/3] Rivet 3.1.11 (--disable-pyext) ==="
if [ ! -s Rivet-3.1.11.tar.gz ]; then DL "https://rivet.hepforge.org/downloads/?f=Rivet-3.1.11.tar.gz" Rivet-3.1.11.tar.gz; fi
rm -rf Rivet-3.1.11 && tar xf Rivet-3.1.11.tar.gz && cd Rivet-3.1.11
./configure --prefix=$PREFIX --with-yoda=$PREFIX --with-hepmc3=$PREFIX --with-fastjet=$PREFIX --with-zlib=$PREFIX --disable-pyext > cfg.log 2>&1
make -j$J > make.log 2>&1
make install > install.log 2>&1
cd ..
echo "Rivet done: $($PREFIX/bin/rivet-nopy --version 2>&1 || $PREFIX/bin/rivet --version 2>&1)"
echo "RIVET_BUILD_COMPLETE"
