# a script to make the Bins directory and 
# link bins from specified directory to it 
# if pass length threshold of 20000

mkdir Bins

for bin in `ls $1 | grep ".fa"`; do
    echo "checking: $1$bin"
    bases=$(grep -v ">" $1$bin | wc -c)
    if [ $bases -gt 19999 ]
    then
        ln -s $1$bin ./Bins/
    fi
done
