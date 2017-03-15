PROJDIR='/users/mtaranov/dnaseEmbed_blood'
DATADIR=$PROJDIR/'data'
for chr in {1..22} X; do

echo 'chr'$chr
#separate chrs
zcat $DATADIR/E116.GM12878_Lymphoblastoid_Cells.ENCODE.Duke_Crawford.DNase-seq.merged.20bp.filt.50m_rep1-pr.IDR0.1.filt.narrowPeak.gz | awk '{if ($1==t) print $0}' t='chr'$chr > $DATADIR/E116.GM12878_'chr'${chr}

gzip $DATADIR/E116.GM12878_'chr'${chr}

#get sequences at dnase peaks
/users/mtaranov/local/anaconda2/bin/python $PROJDIR/get_seq_at_dnasePeaks.py --fasta_file /mnt/data/annotations/by_release/hg19.GRCh37/hg19.genome.fa --dnase_narrowPeak_file $DATADIR/E116.GM12878_'chr'${chr}.gz --window 2000 --output_name $DATADIR/seq_at_dnasePeaks_2kb_'chr'${chr}

#extract last layer embedding
/users/mtaranov/local/anaconda2/envs/keras_1/bin/python2.7 $PROJDIR/get_last_layer.py --json_model $DATADIR/json/record_1_model_GbsDh_modelJson.json --weights_file $DATADIR/weights/record_1_model_GbsDh_modelWeights.h5 --dnase_seq_file $DATADIR/seq_at_dnasePeaks_2kb_'chr'${chr}.npy --output_name $DATADIR/last_layer_embedding_'chr'${chr}

#compute t-sne
/users/mtaranov/local/anaconda2/envs/keras_1/bin/python2.7 $PROJDIR/tsne.py ${chr}

done
