# Instruction to install the Python package:
1. Download this Github repository by "git clone https://github.com/jinyan1214/ASREp.git" or download as zip file directly.
2. In your terminal/command prompt, navigate to the downloaded ASREp directory and run command "pip install ."
3. The package contains a c++ dynamic link library. The dynamic link library is pre-compiled for Windows and MacOS users and the pre-compiled library is located in ASREp/ASREp/ASREcpp/bin
4. If the pre-compiled c++ dynamic link library is not compatible with your computer or you would like to edit the source code of the c++ library, you can find compile commands below

# Commands to compile c++ dynamic link library
## On Windows using Visual Studio
```bash
cd ASREp/ASREp/ASREcpp/src
mkdir build
cd build
cmake .. -G "Visual Studio 17 2022" -A x64
cmake --build . --config Release --parallel
cmake --install .
```
## On macOS using Clang
Latest macOS has clang available from Apple's lightweight command line tools
```bash
cd ASREp/ASREp/ASREcpp/src
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=/path/to/your/clang -DCMAKE_CXX_COMPILER=/path/to/your/clang++ -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --config Release --parallel
cmake --install .
```

## (Legacy)On macOS using GCC
```bash
cd ASREp/ASREp/ASREcpp/src
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=/path/to/your/gcc/12.4.0/bin/gcc-12 -DCMAKE_CXX_COMPILER=/path/to/your/gcc/12.4.0/bin/g++-12 -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --config Release --parallel
cmake --install .
```

You can install on a Mac using homebrew

To install homebrew: please refer to: [https://brew.sh/] (https://brew.sh/)

To install gcc using homebrew: [https://formulae.brew.sh/formula/gcc] (https://formulae.brew.sh/formula/gcc)
