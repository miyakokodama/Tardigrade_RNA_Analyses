cd

for i in $(seq 1 18)
do
	wget ftp://20200812F20FTSEUET0118:RAMvqkE@5.57.48.133/F20FTSEUET0118_RAMvqkE/Clean/"$i"/"$i"_1.fq.gz
	sleep 1
	wget ftp://20200812F20FTSEUET0118:RAMvqkE@5.57.48.133/F20FTSEUET0118_RAMvqkE/Clean/"$i"/"$i"_2.fq.gz
	sleep 1
done
