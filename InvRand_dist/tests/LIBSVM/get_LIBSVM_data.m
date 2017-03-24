cd tests/LIBSVM/data/

files = {'aloi', 'protein','gisette_scale', 'real-sim' };
%% Download and decompress the multiclass problems
system(['wget http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multiclass/' files{1}  '.bz2'])
system(['bzip2 -d '  files{1}  '.bz2'])
system(['wget http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multiclass/' files{2}  '.bz2'])
system(['bzip2 -d '  files{2}  '.bz2'])
%% Download and decompress the binary problems
system(['wget http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/' files{3}  '.bz2'])
system(['bzip2 -d '  files{3}  '.bz2'])
system(['wget http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/' files{4}  '.bz2'])
system(['bzip2 -d '  files{4}  '.bz2'])

cd ../../..;
