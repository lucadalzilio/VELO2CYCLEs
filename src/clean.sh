echo "=================================== "
echo ">>> clean folder"

rm -rf *lsf*
rm -rf *prn
rm -rf *1.gzip.h5
rm -rf *2.gzip.h5
rm -rf *3.gzip.h5
rm -rf *4.gzip.h5
rm -rf *5.gzip.h5
rm -rf *6.gzip.h5
rm -rf *7.gzip.h5
rm -rf *8.gzip.h5
rm -rf *9.gzip.h5
rm -rf stick*txt

echo ">>> move files"

mkdir h5_files
mv *h5 h5_files/

mkdir code/
mv *c code/
mv i2evps code/
mv in2evps code/
mv *sh code/

echo ">>> ok!"
echo "=================================== "
