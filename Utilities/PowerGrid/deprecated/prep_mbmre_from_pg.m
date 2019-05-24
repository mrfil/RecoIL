function mergePowerGridFileOutput(NSlices,NReps,NAvgs,NEchoes,NPhases)
% mergePowerGridFileOutput(NSlices,NReps,NAvgs,NEchoes,NPhases)
% Forms complex N-dimensional images from the output of PowerGrid. Does not
% handle multiband reshuffling if that is required.

for phs = 1:NPhases
for eco = 1:NEchoes
for avg = 1:NAvgs
for rep = 1:NReps
    for slc = 1:NSlices
            tmp1 = load_untouch_nii(sprintf('pcSENSE_Slice%i_Rep%i_Avg0_Echo0_Phase%i_mag.nii',ii-1,jj-1,kk-1),1,1,1,1);
            tmp2 = load_untouch_nii(sprintf('pcSENSE_Slice%i_Rep%i_Avg0_Echo0_Phase%i_phs.nii',ii-1,jj-1,kk-1),1,1,1,1);
            img(:,:,:,ii,jj,kk) = tmp1.img.*exp(1i*tmp2.img);
    end
end
end
end
end

save img.mat img

mkdir recon
!mv pcSENSE* recon/
