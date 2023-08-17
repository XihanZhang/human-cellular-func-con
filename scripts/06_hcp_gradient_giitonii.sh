datadir=/gpfs/milgram/project/holmes/xz555/gradient_shift/data/Gradients_Margulies2016/fsaverage
outdir=${datadir}/nifti
cd ${datadir}

module load FreeSurfer

for ff in `ls *fsa6.func.gii`
do
	mri_convert -i ${ff} -o $outdir/${ff%\.*}.nii.gz

done
