cd #FUJIN_ROOT
mkdir -p external_libraries
cd external_libraries
git clone https://github.com/gabime/spdlog.git
cd spdlog && mkdir build && cd build
cmake .. && make -j
